function [powerRaw, powerNorm] = computeEEGPower(signal, fs, freqBands, fullBandRange, welchWin)
% computeEEGPower - 计算EEG/MUA信号在各频段的功率
%
% 功能描述:
%   使用Welch方法计算信号的功率谱密度(PSD)，并计算各指定频段的原始功率
%   和归一化功率(相对于全频带功率的百分比)。
%
% 语法:
%   [powerRaw, powerNorm] = computeEEGPower(signal, fs, freqBands)
%   [powerRaw, powerNorm] = computeEEGPower(signal, fs, freqBands, fullBandRange)
%   [powerRaw, powerNorm] = computeEEGPower(signal, fs, freqBands, fullBandRange, welchWin)
%
% 输入参数:
%   signal        - 1D数组，输入的EEG/MUA信号 (已预处理)
%   fs            - 标量，信号采样率 (Hz)，例如 256 Hz 或 500 Hz
%   freqBands     - Nx2矩阵，每行定义一个频段的[低频, 高频]边界 (Hz)
%                   例如: [0.5 4; 4 8; 8 14; 14 30; 30 50; 50 80]
%   fullBandRange - (可选) 1x2向量，归一化时使用的全频带范围 [低频, 高频] (Hz)
%                   默认值: [0.5, 80]，但不超过Nyquist频率
%                   EEG常用: [0.5, 45] 或 [0.5, 80]
%                   MUA常用: [0.5, 200] 或更高
%   welchWin      - (可选) 标量，Welch方法的窗口长度 (秒)
%                   默认值: [] (使用固定4段的方法，与原始代码兼容)
%                   推荐值: 2秒 (对应0.5Hz频率分辨率) 或 4秒 (0.25Hz分辨率)
%                   如果指定数值，则使用固定时间窗口方法
%
% 输出参数:
%   powerRaw  - 1xN数组，各频段的原始功率值 (单位: 信号单位^2/Hz)
%   powerNorm - 1xN数组，各频段的归一化功率 (相对全频带功率的百分比, 0-100)
%
% Welch参数说明:
%   - 当welchWin为空时: 使用固定段数方法 (4段, 50%重叠) [自定义的——窗口更长，频率分辨率更高]
%     窗口长度 L = M/2.5，频率分辨率随信号长度变化；同时通过手动设置nfft = fs * 2来保证显示分辨率
% P.S. 真正的物理频率分辨率由窗口长度L决定：Δf = fs/L。显示分辨率（零填充插值）。
%   - 当welchWin指定时: 使用固定时间窗口方法
%     频率分辨率 = 1/welchWin (如2秒窗口→0.5Hz分辨率)
%
% 依赖函数:
%   pwelch (Signal Processing Toolbox)
%   bandpower (Signal Processing Toolbox)
%
% 使用示例:
%   %% 示例1: 基本用法 (使用默认参数，兼容原代码)
%   freqBands = [0.5 4; 4 8; 8 14; 14 30; 30 50; 50 80];
%   [powerRaw, powerNorm] = computeEEGPower(EEG_signal, 256, freqBands);
%
%   %% 示例2: 指定全频带范围
%   fullBandRange = [0.5, 80];  
%   [powerRaw, powerNorm] = computeEEGPower(EEG_signal, 256, freqBands, fullBandRange);
%
%   %% 示例3: 使用固定2秒窗口
%   fullBandRange = [0.5, 80];
%   welchWin = 2;  % 2秒窗口，0.5Hz分辨率
%   [powerRaw, powerNorm] = computeEEGPower(EEG_signal, 256, freqBands, fullBandRange, welchWin);
%
%   %% 示例4: 只指定welchWin，使用默认fullBandRange
%   [powerRaw, powerNorm] = computeEEGPower(EEG_signal, 256, freqBands, [], 2);
%
% 作者: Simba
% 日期: 2025

%% 输入参数检查与默认值设置
if nargin < 3
    error('computeEEGPower:输入参数不足，需要 signal, fs, freqBands');
end

if ~isvector(signal)
    error('computeEEGPower:signal 必须是一维向量');
end

signal = signal(:)'; % 确保为行向量

if fs <= 0
    error('computeEEGPower:采样率 fs 必须为正数');
end

if size(freqBands, 2) ~= 2
    error('computeEEGPower:freqBands 必须是 Nx2 矩阵');
end

nBands = size(freqBands, 1);

% 默认全频带范围: [0.5, 80] Hz
if nargin < 4 || isempty(fullBandRange)
    % 检查默认值80Hz是否超过Nyquist频率
    if 80 >= fs/2
        error('computeEEGPower:采样率%.1fHz过低(Nyquist=%.1fHz)，请手动指定fullBandRange参数', fs, fs/2);
    end
    fullBandRange = [0.5, 80];
else
    % 用户指定了fullBandRange，检查并调整上限
    if fullBandRange(2) >= fs/2
        warning('computeEEGPower:fullBandRange上限 %.1f Hz 超过Nyquist频率 %.1f Hz，已自动调整', ...
                fullBandRange(2), fs/2);
        fullBandRange(2) = fs/2 - 0.5;
    end
end


% 默认welchWin为空，使用固定段数方法
if nargin < 5
    welchWin = [];
end

%% Welch参数设置
M = length(signal);

if isempty(welchWin)
    % 方法A: 固定段数 (原始代码逻辑，兼容性保留)
    segN = 4;                      % 4段
    L = fix(M / ((1 + segN) / 2)); % 窗口长度 = M/2.5
    noverlap = fix(0.5 * L);       % 50%重叠
    nfft = fs * 2;                 % FFT点数
    
    if L < 2
        warning('computeEEGPower:信号长度过短，结果可能不可靠');
        L = M;
        noverlap = 0;
    end
else
    % 方法B: 固定时间窗口
    L = round(welchWin * fs);      % 窗口长度 (样本数)
    noverlap = round(0.5 * L);     % 50%重叠
    nfft = max(L, 2^nextpow2(L));  % FFT点数，至少为窗口长度
    
    % 检查窗口长度是否合适
    if L > M
        warning('computeEEGPower:信号长度(%.2f秒)小于窗口长度(%.2f秒)，改用整段信号', ...
                M/fs, welchWin);
        L = M;
        noverlap = 0;
        nfft = max(L, 2^nextpow2(L));
    end
    
    % 检查最低频率是否可分辨
    freqRes = fs / L;  % 物理频率分辨率
    minFreqBand = min(freqBands(:, 1));
    if freqRes > minFreqBand
        warning('computeEEGPower:频率分辨率(%.2fHz)大于最低频段下限(%.2fHz)，建议增加welchWin', ...
                freqRes, minFreqBand);
    end
end

%% 计算功率谱密度
[pxx, f] = pwelch(signal, hamming(L), noverlap, nfft, fs);

%% 计算各频段功率
powerRaw = zeros(1, nBands);
for i = 1:nBands
    % 检查频段是否在有效范围内
    if freqBands(i, 1) < 0 || freqBands(i, 2) > fs/2
        warning('computeEEGPower:频段 [%.1f, %.1f] 超出Nyquist频率 %.1f Hz', ...
                freqBands(i, 1), freqBands(i, 2), fs/2);
    end
    powerRaw(i) = bandpower(pxx, f, freqBands(i, :), 'psd');
end

%% 计算归一化功率
powerFull = bandpower(pxx, f, fullBandRange, 'psd');

if powerFull > 0
    powerNorm = (powerRaw / powerFull) * 100;
else
    warning('computeEEGPower:全频带功率为零，归一化功率设为零');
    powerNorm = zeros(1, nBands);
end

end
