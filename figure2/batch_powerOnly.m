%% 批量处理脚本: 从文件夹批量加载 .set 文件，对全部epoch计算功率谱
%
% 功能说明:
%   - 批量处理指定文件夹下所有 .set 文件
%   - 自动适配不同导联数(19导/32导等)
%   - 逐通道计算功率谱，再对通道间结果取平均（兼容平均参考与非平均参考数据）
%   - 计算指标: 功率谱 (原始功率 + 归一化功率)
%   - 所有患者结果输出到同一个Excel文件
%
% 数据说明:
%   - 文件格式: EEGLAB .set (命名规则: 患者名字_pre_rejected_processed.set)
%   - 静息态数据，已做 1-40 Hz 带通滤波预处理
%   - 已分割为 2秒/epoch
%   - 对每位患者的全部epoch拼接后，逐通道计算功率谱再平均
%   - 频段: Delta(1-4), Theta(4-8), Alpha(8-13), Beta(13-30), Gamma(30-40)
%
% 注意:
%   - 使用逐通道计算功率谱再平均，而非先对信号取全导联均值再计算
%   - 原因: 平均参考预处理后全导联信号均值恒为0，直接对均值信号计算功率谱结果趋近于0
%
% 依赖:
%   EEGLAB (需要已添加到MATLAB路径中)
%   computeEEGPower.m
%
% 作者: Simba (批量处理版 - 仅功率谱)
% 日期: 2025

clear; clc; close all;

%% ==================== 用户配置区 ====================
% ---- 输入输出路径设置 ----
inputDir  = 'F:\EEG\LowVoltageResearch\PreprocessedData\lowThetaDelta_ICA';     % ★ 输入文件夹路径 (存放所有.set文件)
outputDir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\全频段功率谱';  % ★ 输出文件夹路径

% 检查输入目录
if ~exist(inputDir, 'dir')
    error('输入目录不存在: %s\n请修改 inputDir 变量', inputDir);
end

% 检查输出目录，不存在则创建
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('已创建输出目录: %s\n', outputDir);
end

% ---- 频谱参数 ----
cfg.fullBandRange = [1, 40];     % 归一化全频带范围 (Hz)
cfg.welchWin      = 2;           % Welch窗口长度 (秒)

cfg.freqBands = [1  4;           % Delta
                 4  8;           % Theta
                 8  13;          % Alpha
                 13 30;          % Beta
                 30 40];         % Gamma

cfg.bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
nBands = length(cfg.bandNames);

%% ==================== 扫描文件夹获取所有.set文件 ====================
fprintf('\n========================================\n');
fprintf('=== 扫描输入文件夹 ===\n');
fprintf('========================================\n');
fprintf('路径: %s\n', inputDir);

fileList = dir(fullfile(inputDir, '*.set'));
nFiles = length(fileList);

if nFiles == 0
    error('未在文件夹中找到任何 .set 文件');
end

fprintf('找到 %d 个 .set 文件:\n', nFiles);
for f = 1:nFiles
    fprintf('  [%d] %s\n', f, fileList(f).name);
end

%% ==================== 初始化汇总结果存储 ====================
% 定义表格列名
resultVarNames = {
    'Patient', 'nChannels', 'nEpochs', 'TotalDuration_sec', ...
    'FullBand_RawPower', ...
    'Delta_RawPower', 'Theta_RawPower', 'Alpha_RawPower', 'Beta_RawPower', 'Gamma_RawPower', ...
    'FullBand_NormPower', ...
    'Delta_NormPower', 'Theta_NormPower', 'Alpha_NormPower', 'Beta_NormPower', 'Gamma_NormPower'
};

% 初始化空表格
allPatientsResult = table('Size', [0, length(resultVarNames)], ...
    'VariableTypes', [{'string', 'double', 'double', 'double'}, repmat({'double'}, 1, length(resultVarNames)-4)], ...
    'VariableNames', resultVarNames);

% 记录处理状态
processStatus = cell(nFiles, 3);

%% ==================== 批量处理循环 ====================
fprintf('\n========================================\n');
fprintf('=== 开始批量处理 (%d 个文件) ===\n', nFiles);
fprintf('========================================\n');

for fileIdx = 1:nFiles
    fileName = fileList(fileIdx).name;
    filePath = fullfile(inputDir, fileName);

    fprintf('\n----------------------------------------\n');
    fprintf('[%d/%d] 处理: %s\n', fileIdx, nFiles, fileName);
    fprintf('----------------------------------------\n');

    try
        %% ===== 提取患者名 =====
        suffix = '_pre_rejected_processed.set';
        if endsWith(fileName, suffix, 'IgnoreCase', true)
            patientName = fileName(1:end - length(suffix));
        else
            parts = strsplit(fileName, '_');
            patientName = parts{1};
            warning('文件名不符合标准命名规则，提取患者名为: %s', patientName);
        end
        fprintf('患者: %s\n', patientName);

        %% ===== 加载数据 =====
        EEG = pop_loadset('filename', fileName, 'filepath', inputDir);

        % 提取基本信息
        fs        = EEG.srate;
        nChannels = EEG.nbchan;
        nEpochs   = EEG.trials;
        nPoints   = EEG.pnts;

        cfg.fs = fs;

        totalDuration = nEpochs * nPoints / fs;   % 总时长(秒)

        fprintf('采样率: %d Hz | 通道数: %d | Epoch数: %d\n', fs, nChannels, nEpochs);
        fprintf('总时长: %.1f 秒 (%.1f 分钟)\n', totalDuration, totalDuration / 60);

        %% ===== Nyquist检查 =====
        nyquistFreq = fs / 2;
        localFullBandRange = cfg.fullBandRange;
        if localFullBandRange(2) >= nyquistFreq
            localFullBandRange(2) = nyquistFreq - 1;
            warning('fullBandRange上限已调整为 %.1f Hz', localFullBandRange(2));
        end

        %% ===== 拼接所有epoch，逐通道计算功率谱再平均 =====
        % EEG.data 维度: [nChannels × nPoints × nEpochs]
        % reshape 后:   [nChannels × (nPoints * nEpochs)]
        allData = reshape(EEG.data, nChannels, []);
        signalLen = size(allData, 2);

        fprintf('导联数: %d | 拼接后信号长度: %d 样本\n', nChannels, signalLen);

        % 逐通道计算功率谱
        powerRaw_allChan      = zeros(nChannels, nBands);
        powerNorm_allChan     = zeros(nChannels, nBands);
        powerRaw_full_allChan = zeros(nChannels, 1);   % 全频段原始功率 (逐通道)

        for ch = 1:nChannels
            chSignal = double(allData(ch, :));   % 1 × (nPoints*nEpochs)
            [powerRaw_allChan(ch, :), powerNorm_allChan(ch, :)] = ...
                computeEEGPower(chSignal, fs, cfg.freqBands, localFullBandRange, cfg.welchWin);
            % 全频段: 将 localFullBandRange 作为单一频段传入，归一化分母与被积区间相同，结果=1
            [powerRaw_full_allChan(ch), ~] = ...
                computeEEGPower(chSignal, fs, localFullBandRange, localFullBandRange, cfg.welchWin);
        end

        % 通道间取平均
        powerRaw      = mean(powerRaw_allChan,      1);   % 1 × nBands
        powerNorm     = mean(powerNorm_allChan,     1);   % 1 × nBands
        powerRaw_full = mean(powerRaw_full_allChan, 1);   % 标量，全频段原始功率均值
        % 全频段归一化功率恒为1（归一化分母 = 全频段功率本身），记为1
        powerNorm_full = 1.0;

        fprintf('功率计算完成 (逐通道计算后平均，含全频段 %.0f-%.0f Hz)\n', ...
            localFullBandRange(1), localFullBandRange(2));

        %% ===== 存储结果 =====
        newRow = table(...
            string(patientName), nChannels, nEpochs, totalDuration, ...
            powerRaw_full, ...
            powerRaw(1), powerRaw(2), powerRaw(3), powerRaw(4), powerRaw(5), ...
            powerNorm_full, ...
            powerNorm(1), powerNorm(2), powerNorm(3), powerNorm(4), powerNorm(5), ...
            'VariableNames', resultVarNames);
        allPatientsResult = [allPatientsResult; newRow];

        %% ===== 记录成功状态 =====
        processStatus{fileIdx, 1} = fileName;
        processStatus{fileIdx, 2} = '成功';
        processStatus{fileIdx, 3} = sprintf('%d导, %d epoch, %.1f秒', nChannels, nEpochs, totalDuration);

    catch ME
        %% ===== 错误处理 =====
        warning('处理文件 %s 时出错: %s', fileName, ME.message);
        processStatus{fileIdx, 1} = fileName;
        processStatus{fileIdx, 2} = '失败';
        processStatus{fileIdx, 3} = ME.message;
    end
end

%% ==================== 导出Excel ====================
fprintf('\n========================================\n');
fprintf('=== 导出结果到Excel ===\n');
fprintf('========================================\n');

excelFileName = 'EEG_PowerSpectrum.xlsx';
excelFilePath = fullfile(outputDir, excelFileName);

% ----- Sheet 1: 处理状态 -----
statusTable = cell2table(processStatus, ...
    'VariableNames', {'FileName', 'Status', 'Note'});
writetable(statusTable, excelFilePath, 'Sheet', 'ProcessStatus');

% ----- Sheet 2: 分析参数 -----
paramData = {
    'AnalysisTime',      datestr(now, 'yyyy-mm-dd HH:MM:SS');
    'InputDir',          inputDir;
    'OutputDir',         outputDir;
    'TotalFiles',        nFiles;
    'SuccessCount',      sum(strcmp(processStatus(:,2), '成功'));
    'FullBandRange_Hz',  sprintf('%d-%d', cfg.fullBandRange(1), cfg.fullBandRange(2));
    'WelchWindow_sec',   cfg.welchWin;
    'ChannelAveraging',  'PerChannelPSD_ThenMean';
    'SegmentStrategy',   'AllEpochsConcatenated';
};
paramTable = cell2table(paramData, 'VariableNames', {'Parameter', 'Value'});
writetable(paramTable, excelFilePath, 'Sheet', 'Parameters');

% ----- Sheet 3: 频段定义 -----
freqDefData = cell(nBands + 1, 3);
freqDefData{1, 1} = 'FullBand';
freqDefData{1, 2} = cfg.fullBandRange(1);
freqDefData{1, 3} = cfg.fullBandRange(2);
for b = 1:nBands
    freqDefData{b + 1, 1} = cfg.bandNames{b};
    freqDefData{b + 1, 2} = cfg.freqBands(b, 1);
    freqDefData{b + 1, 3} = cfg.freqBands(b, 2);
end
freqDefTable = cell2table(freqDefData, 'VariableNames', {'BandName', 'LowFreq_Hz', 'HighFreq_Hz'});
writetable(freqDefTable, excelFilePath, 'Sheet', 'FreqBands');

% ----- Sheet 4: 功率谱结果 -----
if height(allPatientsResult) > 0
    writetable(allPatientsResult, excelFilePath, 'Sheet', 'PowerSpectrum');
end

fprintf('Excel已保存: %s\n', excelFilePath);

%% ==================== 最终汇总报告 ====================
fprintf('\n========================================\n');
fprintf('=== 批量处理完成 ===\n');
fprintf('========================================\n');

nSuccess = sum(strcmp(processStatus(:,2), '成功'));
nFailed  = sum(strcmp(processStatus(:,2), '失败'));

fprintf('总文件数:   %d\n', nFiles);
fprintf('成功处理:   %d\n', nSuccess);
fprintf('处理失败:   %d\n', nFailed);
fprintf('\n输出目录: %s\n', outputDir);
fprintf('Excel文件: %s\n', excelFileName);

fprintf('\n计算指标:\n');
fprintf('  1. 功率谱 (原始功率 + 归一化功率)\n');
fprintf('  2. 全频段功率 (FullBand %d-%d Hz, 原始功率 + 归一化功率)\n', ...
    cfg.fullBandRange(1), cfg.fullBandRange(2));

fprintf('\n频段定义:\n');
fprintf('  FullBand: %d-%d Hz\n', cfg.fullBandRange(1), cfg.fullBandRange(2));
for b = 1:nBands
    fprintf('  %s: %d-%d Hz\n', cfg.bandNames{b}, cfg.freqBands(b, 1), cfg.freqBands(b, 2));
end

fprintf('\n导联处理方式: 逐通道计算功率谱，再取通道间平均 (兼容平均参考与非平均参考)\n');
fprintf('分段策略:     全部epoch拼接，整体计算一次功率谱\n');
fprintf('========================================\n');

% 显示各文件处理状态
fprintf('\n各文件处理状态:\n');
fprintf('%-45s  %-6s  %s\n', '文件名', '状态', '备注');
fprintf('%s\n', repmat('-', 1, 80));
for f = 1:nFiles
    fprintf('%-45s  %-6s  %s\n', ...
            processStatus{f, 1}, processStatus{f, 2}, processStatus{f, 3});
end
fprintf('\n');