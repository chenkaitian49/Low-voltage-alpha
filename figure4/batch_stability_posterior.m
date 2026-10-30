%% 批量处理脚本: 后头部通道稳定性计算（整段信号）
%
% 功能说明:
%   - 批量处理指定文件夹下所有 .set 文件
%   - 只对后头部通道进行计算，取通道间均值
%   - 对每个患者的全部数据拼接为一段信号，计算一次稳定性
%   - 频段: DeltaTheta(1-8Hz), Alpha(8-13Hz)
%   - 所有患者结果输出到同一个Excel文件
%
% 数据说明:
%   - 文件格式: EEGLAB .set (命名规则: 患者名字_pre_rejected_processed.set)
%   - 静息态数据，已做预处理
%   - 已分割为 2秒/epoch
%
% 后头部通道:
%   P3, P4, O1, O2, T5, T6, Pz, PO3, PO4, PO7, PO8, Oz
%
% 依赖:
%   EEGLAB (需要已添加到MATLAB路径中)
%   computeEEGStability.m
%
% 作者: 改编自 example_usage_batch_fixed.m
% 日期: 2025

clear; clc; close all;

%% ==================== 用户配置区 ====================
% ---- 输入输出路径设置 ----
inputDir  = 'F:\EEG\LowVoltageResearch\PreprocessedData\TestingData\RunICA';
outputDir = 'F:\EEG\LowVoltageResearch\OutputResult\外部验证数据\主导节律稳定性原始数据（后头部1s窗宽）';

% 检查输入目录
if ~exist(inputDir, 'dir')
    error('输入目录不存在: %s\n请修改 inputDir 变量', inputDir);
end

% 检查输出目录，不存在则创建
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('已创建输出目录: %s\n', outputDir);
end

% ---- 后头部通道定义 ----
posterior_channels = {'P3', 'P4', 'O1', 'O2', 'T5', 'T6', 'Pz', ...
                      'PO3', 'PO4', 'PO7', 'PO8', 'Oz'};

% ---- 稳定性参数 ----
cfg.stabilityWin = 1;   % 稳定性计算滑动窗口 (秒)

cfg.freqBands = [1   8;   % DeltaTheta (Delta + Theta 合并)
                 8  13];  % Alpha

cfg.bandNames = {'DeltaTheta', 'Alpha'};
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
% 每个患者输出一行: 患者名、有效后头部通道数、总时长、两个频段稳定性
resultVarNames = {'Patient', 'nPostChan', 'TotalDuration_sec', ...
                  'DeltaTheta_Stability', 'Alpha_Stability'};

allResults = table('Size', [0, length(resultVarNames)], ...
    'VariableTypes', [{'string', 'double', 'double'}, repmat({'double'}, 1, 2)], ...
    'VariableNames', resultVarNames);

% 记录处理状态
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
        EEG = pop_loadset('filename', fileName, 'filepath', inputDir);

        fs        = EEG.srate;
        nEpochs   = EEG.trials;
        nPoints   = EEG.pnts;
        allLabels = {EEG.chanlocs.labels};
        totalDuration = nEpochs * nPoints / fs;  % 总时长 (秒)

        fprintf('采样率: %d Hz | 通道数: %d | Epoch数: %d\n', fs, EEG.nbchan, nEpochs);
        fprintf('总时长: %.1f 秒 (%.1f 分钟)\n', totalDuration, totalDuration / 60);

        %% ===== 确定后头部通道索引 =====
        postChanIdx = [];
        foundNames  = {};
        for c = 1:length(posterior_channels)
            idx = find(strcmpi(allLabels, posterior_channels{c}));
            if ~isempty(idx)
                postChanIdx(end+1) = idx(1); %#ok<AGROW>
                foundNames{end+1}  = allLabels{idx(1)}; %#ok<AGROW>
            else
                warning('通道 %s 未在数据中找到，已跳过', posterior_channels{c});
            end
        end

        if isempty(postChanIdx)
            error('未找到任何后头部通道，跳过此文件');
        end

        nPostChan = length(postChanIdx);
        fprintf('找到后头部通道 (%d/%d): %s\n', nPostChan, length(posterior_channels), ...
                strjoin(foundNames, ', '));

        %% ===== Nyquist检查 =====
        nyquistFreq = fs / 2;
        for b = 1:nBands
            if cfg.freqBands(b, 2) >= nyquistFreq
                error('频段 %s 上限 %.1f Hz 超过Nyquist频率 %.1f Hz，请检查采样率', ...
                      cfg.bandNames{b}, cfg.freqBands(b, 2), nyquistFreq);
            end
        end

        %% ===== 拼接全部数据为一段连续信号 =====
        % EEG.data: [nChan × nPoints × nEpochs]
        % 取后头部通道均值后展平为行向量: 1 × (nPoints × nEpochs)
        allData   = EEG.data(postChanIdx, :, :);    % [nPostChan × nPoints × nEpochs]
        chanAvg   = mean(allData, 1);               % [1 × nPoints × nEpochs]
        segSignal = double(reshape(chanAvg, 1, []));% [1 × (nPoints×nEpochs)]

        fprintf('信号长度: %d 个采样点 (%.1f 秒)\n', length(segSignal), length(segSignal)/fs);

        %% ===== 计算稳定性 =====
        % computeEEGStability 返回 [宽带, 频段1, ..., 频段N]
        % 用 (2:end) 去掉宽带，只保留频段部分
        stabAll   = computeEEGStability(segSignal, fs, cfg);
        stability = stabAll(2:end);   % [DeltaTheta_Stability, Alpha_Stability]

        fprintf('DeltaTheta 稳定性: %.4f\n', stability(1));
        fprintf('Alpha      稳定性: %.4f\n', stability(2));

        %% ===== 存储结果 =====
        newRow = table(...
            string(patientName), nPostChan, totalDuration, ...
            stability(1), stability(2), ...
            'VariableNames', resultVarNames);
        allResults = [allResults; newRow]; %#ok<AGROW>

        processStatus{fileIdx, 1} = fileName;
        processStatus{fileIdx, 2} = '成功';
        processStatus{fileIdx, 3} = sprintf('后头部%d导, %.1f秒', nPostChan, totalDuration);

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

excelFileName = 'Stability_Posterior_DeltaTheta_Alpha.xlsx';
excelFilePath = fullfile(outputDir, excelFileName);

% ----- Sheet 1: 处理状态 -----
statusTable = cell2table(processStatus, ...
    'VariableNames', {'FileName', 'Status', 'Note'});
writetable(statusTable, excelFilePath, 'Sheet', 'ProcessStatus');

% ----- Sheet 2: 分析参数 -----
paramData = {
    'AnalysisTime',        datestr(now, 'yyyy-mm-dd HH:MM:SS');
    'InputDir',            inputDir;
    'OutputDir',           outputDir;
    'TotalFiles',          nFiles;
    'SuccessCount',        sum(strcmp(processStatus(:,2), '成功'));
    'StabilityWindow_sec', cfg.stabilityWin;
    'ChannelScope',        strjoin(posterior_channels, ', ');
    'ChannelAveraging',    '后头部通道均值';
    'SignalScope',         '全部epoch拼接为一段连续信号';
};
paramTable = cell2table(paramData, 'VariableNames', {'Parameter', 'Value'});
writetable(paramTable, excelFilePath, 'Sheet', 'Parameters');

% ----- Sheet 3: 频段定义 -----
freqDefData = cell(nBands, 3);
for b = 1:nBands
    freqDefData{b, 1} = cfg.bandNames{b};
    freqDefData{b, 2} = cfg.freqBands(b, 1);
    freqDefData{b, 3} = cfg.freqBands(b, 2);
end
freqDefTable = cell2table(freqDefData, 'VariableNames', {'BandName', 'LowFreq_Hz', 'HighFreq_Hz'});
writetable(freqDefTable, excelFilePath, 'Sheet', 'FreqBands');

% ----- Sheet 4: 结果（每患者一行）-----
if height(allResults) > 0
    writetable(allResults, excelFilePath, 'Sheet', 'Results');
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

fprintf('\n频段定义:\n');
for b = 1:nBands
    fprintf('  %s: %d-%d Hz\n', cfg.bandNames{b}, cfg.freqBands(b, 1), cfg.freqBands(b, 2));
end
fprintf('\n后头部通道 (%d个): %s\n', length(posterior_channels), strjoin(posterior_channels, ', '));

fprintf('\n各文件处理状态:\n');
fprintf('%-45s  %-6s  %s\n', '文件名', '状态', '备注');
fprintf('%s\n', repmat('-', 1, 80));
for f = 1:nFiles
    fprintf('%-45s  %-6s  %s\n', ...
            processStatus{f, 1}, processStatus{f, 2}, processStatus{f, 3});
end
fprintf('\n');