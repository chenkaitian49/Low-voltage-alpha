%% 批量处理脚本: 从文件夹批量加载 .set 文件，对全段数据计算后头部通道平均振幅
%
% 功能说明:
%   - 批量处理指定文件夹下所有 .set 文件
%   - 仅使用后头部指定通道 (P3/P4/O1/O2/T5/T6/Pz/PO3/PO4/PO7/PO8/Oz)
%   - 自动检测当前数据中存在哪些指定通道，只用实际存在的通道
%   - 对全部epoch一次性计算振幅，每个被试输出一行结果
%   - 所有患者结果输出到同一个Excel文件
%
% 数据说明:
%   - 文件格式: EEGLAB .set (命名规则: 患者名字_pre_rejected_processed.set)
%   - 静息态数据，已做 1-40 Hz 带通滤波预处理
%   - 已分割为 2秒/epoch
%   - 频段: Delta(1-4), Theta(4-8), Alpha(8-13), Beta(13-30), Gamma(30-40)
%
% 振幅指标说明:
%   - 对每个通道的每个epoch计算Hilbert包络均值
%   - 先做通道间平均，再做epoch间平均，得到该被试的振幅代表值
%   - 分别输出全频带振幅和各频段振幅
%
% 依赖:
%   EEGLAB (需要已添加到MATLAB路径中)
%   computeEEGAmplitude.m
%
% 作者: Simba
% 日期: 2025

clear; clc; close all;

%% ==================== 用户配置区 ====================
% ---- 输入输出路径设置 ----
inputDir  = 'F:\EEG\LowVoltageResearch\PreprocessedData\TestingData\RunICA';
outputDir = 'F:\EEG\LowVoltageResearch\OutputResult\外部验证数据\振幅';

% 检查输入目录
if ~exist(inputDir, 'dir')
    error('输入目录不存在: %s\n请修改 inputDir 变量', inputDir);
end

% 检查输出目录，不存在则创建
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('已创建输出目录: %s\n', outputDir);
end

% ---- 振幅计算参数 ----
cfg.freqBands = [1  4;    % Delta
                 4  8;    % Theta
                 8  13;   % Alpha
                 13 30;   % Beta
                 30 40];  % Gamma

cfg.bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
nBands = length(cfg.bandNames);

% ---- 后头部指定通道列表 ----
posterior_channels = {'P3', 'P4', 'O1', 'O2', 'T5', 'T6', 'Pz', ...
                      'PO3', 'PO4', 'PO7', 'PO8', 'Oz'};

%% ==================== 扫描文件夹 ====================
fprintf('\n========================================\n');
fprintf('=== 扫描输入文件夹 ===\n');
fprintf('========================================\n');
fprintf('路径: %s\n', inputDir);

fileList = dir(fullfile(inputDir, '*.set'));
nFiles   = length(fileList);

if nFiles == 0
    error('未在文件夹中找到任何 .set 文件');
end

fprintf('找到 %d 个 .set 文件:\n', nFiles);
for f = 1:nFiles
    fprintf('  [%d] %s\n', f, fileList(f).name);
end

%% ==================== 初始化结果存储 ====================
% 每个被试输出一行: Patient | nChannelsUsed | ChannelsUsed | nEpochs |
%                   Broadband + 5频段 振幅
resultVarNames = {
    'Patient', 'nChannelsUsed', 'ChannelsUsed', 'nEpochs', ...
    'Broadband_Amplitude', ...
    'Delta_Amplitude', 'Theta_Amplitude', 'Alpha_Amplitude', ...
    'Beta_Amplitude',  'Gamma_Amplitude'
};

allResults = table('Size', [0, length(resultVarNames)], ...
    'VariableTypes', [{'string','double','string','double'}, repmat({'double'}, 1, nBands+1)], ...
    'VariableNames', resultVarNames);

processStatus = cell(nFiles, 3);

%% ==================== 批量处理循环 ====================
fprintf('\n========================================\n');
fprintf('=== 开始批量处理 (%d 个文件) ===\n', nFiles);
fprintf('========================================\n');

for fileIdx = 1:nFiles
    fileName = fileList(fileIdx).name;

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
        EEG       = pop_loadset('filename', fileName, 'filepath', inputDir);
        fs        = EEG.srate;
        nEpochs   = EEG.trials;
        allLabels = {EEG.chanlocs.labels};

        cfg.fs = fs;

        fprintf('采样率: %d Hz | 通道数: %d | Epoch数: %d\n', fs, EEG.nbchan, nEpochs);

        %% ===== 后头部通道匹配 =====
        chanIdx    = [];
        usedLabels = {};
        for c = 1:length(posterior_channels)
            matchIdx = find(strcmpi(allLabels, posterior_channels{c}));
            if ~isempty(matchIdx)
                chanIdx(end+1)    = matchIdx(1);             %#ok<AGROW>
                usedLabels{end+1} = allLabels{matchIdx(1)};  %#ok<AGROW>
            end
        end

        if isempty(chanIdx)
            warning('未找到任何指定后头部通道，跳过: %s', fileName);
            processStatus{fileIdx, 1} = fileName;
            processStatus{fileIdx, 2} = '跳过';
            processStatus{fileIdx, 3} = '无匹配后头部通道';
            continue;
        end

        nChanUsed    = length(chanIdx);
        usedLabelStr = strjoin(usedLabels, ',');
        fprintf('后头部通道: %d / %d 个 → %s\n', nChanUsed, length(posterior_channels), usedLabelStr);

        %% ===== 全段振幅计算 =====
        % 逐epoch × 逐通道计算，先通道间平均，再epoch间平均
        ampPerEpoch = zeros(nEpochs, nBands + 1);

        for ep = 1:nEpochs
            ampPerChan = zeros(nChanUsed, nBands + 1);
            for ch = 1:nChanUsed
                chSignal = double(EEG.data(chanIdx(ch), :, ep));  % 1 × nPoints
                ampPerChan(ch, :) = computeEEGAmplitude(chSignal, fs, cfg);
            end
            ampPerEpoch(ep, :) = mean(ampPerChan, 1);  % 通道间平均
        end

        amplitude = mean(ampPerEpoch, 1);   % epoch间平均 → 1 × (nBands+1)

        fprintf('振幅计算完成 (共 %d epochs)\n', nEpochs);

        %% ===== 存储结果 =====
        newRow = table(...
            string(patientName), nChanUsed, string(usedLabelStr), nEpochs, ...
            amplitude(1), ...
            amplitude(2), amplitude(3), amplitude(4), amplitude(5), amplitude(6), ...
            'VariableNames', resultVarNames);
        allResults = [allResults; newRow]; %#ok<AGROW>

        processStatus{fileIdx, 1} = fileName;
        processStatus{fileIdx, 2} = '成功';
        processStatus{fileIdx, 3} = sprintf('%d后头部通道, %d epochs', nChanUsed, nEpochs);

    catch ME
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

excelFileName = 'Posterior_Amplitude.xlsx';
excelFilePath = fullfile(outputDir, excelFileName);

% Sheet 1: 处理状态
statusTable = cell2table(processStatus, 'VariableNames', {'FileName','Status','Note'});
writetable(statusTable, excelFilePath, 'Sheet', 'ProcessStatus');

% Sheet 2: 分析参数
paramData = {
    'AnalysisTime',    datestr(now, 'yyyy-mm-dd HH:MM:SS');
    'InputDir',        inputDir;
    'OutputDir',       outputDir;
    'TotalFiles',      nFiles;
    'SuccessCount',    sum(strcmp(processStatus(:,2), '成功'));
    'AmplitudeMethod', 'HilbertEnvelopeMean';
    'TargetChannels',  strjoin(posterior_channels, ', ');
    'Note',            '实际使用通道因数据集而异，见ChannelsUsed列';
};
writetable(cell2table(paramData, 'VariableNames', {'Parameter','Value'}), ...
    excelFilePath, 'Sheet', 'Parameters');

% Sheet 3: 频段定义
freqDefData = [cfg.bandNames', num2cell(cfg.freqBands)];
writetable(cell2table(freqDefData, 'VariableNames', {'BandName','LowFreq_Hz','HighFreq_Hz'}), ...
    excelFilePath, 'Sheet', 'FreqBands');

% Sheet 4: 主结果 (每被试一行)
if height(allResults) > 0
    writetable(allResults, excelFilePath, 'Sheet', 'Results');
end

fprintf('Excel已保存: %s\n', excelFilePath);

%% ==================== 汇总报告 ====================
fprintf('\n========================================\n');
fprintf('=== 批量处理完成 ===\n');
fprintf('========================================\n');

nSuccess = sum(strcmp(processStatus(:,2), '成功'));
nFailed  = sum(strcmp(processStatus(:,2), '失败'));
nSkipped = sum(strcmp(processStatus(:,2), '跳过'));

fprintf('总文件数: %d | 成功: %d | 失败: %d | 跳过: %d\n', nFiles, nSuccess, nFailed, nSkipped);
fprintf('输出: %s\n', excelFilePath);

fprintf('\n频段定义:\n');
for b = 1:nBands
    fprintf('  %s: %d-%d Hz\n', cfg.bandNames{b}, cfg.freqBands(b,1), cfg.freqBands(b,2));
end

fprintf('\n各文件处理状态:\n');
fprintf('%-45s  %-6s  %s\n', '文件名', '状态', '备注');
fprintf('%s\n', repmat('-',1,80));
for f = 1:nFiles
    fprintf('%-45s  %-6s  %s\n', processStatus{f,1}, processStatus{f,2}, processStatus{f,3});
end
fprintf('\n');