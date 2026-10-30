function amplitude = computeEEGAmplitude(signal, fs, cfg)
% computeEEGAmplitude - 计算EEG信号在各频段的振幅(Hilbert包络均值)
%
% 功能描述:
%   通过Hilbert变换计算信号的瞬时振幅包络，取包络均值作为振幅量化指标。
%   分别计算全频带(宽带)振幅和各指定频段的振幅。
%
% 原理说明:
%   1. Hilbert变换将实信号转换为解析信号: z(t) = x(t) + j*H{x(t)}
%   2. 瞬时振幅包络 = |z(t)| = sqrt(x(t)^2 + H{x(t)}^2)
%   3. 振幅指标 = mean(瞬时振幅包络)
%   4. 该指标反映信号的平均能量/强度水平
%
% 语法:
%   amplitude = computeEEGAmplitude(signal, fs, cfg)
%
% 输入参数:
%   signal - 1D数组，输入的EEG信号 (已预处理，如1-40Hz带通)
%   fs     - 标量，信号采样率 (Hz)
%   cfg    - 配置结构体，包含以下字段:
%            .freqBands - Nx2矩阵，频段定义 [低频, 高频] (Hz)
%                         例如: [1 4; 4 8; 8 13; 13 30; 30 40]
%
% 输出参数:
%   amplitude - 1x(N+1)数组，振幅指标 (单位与输入信号相同，通常为μV)
%               amplitude(1)     : 全频带(宽带)振幅 - 直接对原信号计算
%               amplitude(2:end) : 各频段振幅，顺序与freqBands一致
%               NaN表示无法计算(如滤波失败)
%
% 算法说明:
%   1. 全频带振幅: 直接对输入信号(已1-40Hz预处理)计算Hilbert包络均值
%   2. 各频段振幅: 先带通滤波到目标频段，再计算Hilbert包络均值
%   3. 使用4阶Butterworth带通滤波器 (butter + filtfilt零相位滤波)
%
% 依赖函数:
%   hilbert (Signal Processing Toolbox)
%   butter, filtfilt (Signal Processing Toolbox)
%
% 使用示例:
%   %% 示例1: 基本用法
%   cfg.freqBands = [1 4; 4 8; 8 13; 13 30; 30 40];
%   cfg.bandNames = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
%   
%   % 假设 EEG_signal 是预处理后的信号 (1-40Hz)，采样率 fs = 500 Hz
%   amplitude = computeEEGAmplitude(EEG_signal, 500, cfg);
%   
%   % 显示结果
%   fprintf('宽带振幅: %.4f μV\n', amplitude(1));
%   for i = 1:length(cfg.bandNames)
%       fprintf('%s 振幅: %.4f μV\n', cfg.bandNames{i}, amplitude(i+1));
%   end
%
%   %% 示例2: 批量处理多段数据
%   nSegments = 10;
%   allAmplitude = zeros(nSegments, 6);  % 宽带 + 5频段
%   for seg = 1:nSegments
%       allAmplitude(seg, :) = computeEEGAmplitude(segSignals{seg}, fs, cfg);
%   end
%
% 注意事项:
%   - 输入信号应已做预处理(如带通滤波)，全频带振幅基于此计算
%   - 各频段的频率范围应在信号的有效频带内
%   - 返回NaN表示该频段无法计算振幅
%
% 与功率的区别:
%   - 功率 (Power): 基于频谱分析，单位为 μV²/Hz，反映频域能量分布
%   - 振幅 (Amplitude): 基于时域包络，单位为 μV，反映信号强度
%   - 两者相关但不等价: Power ∝ Amplitude²
%
% 作者: Simba
% 日期: 2025
% 版本: 1.0

%% 输入参数检查
if nargin < 3
    error('computeEEGAmplitude:输入参数不足，需要 signal, fs, cfg');
end

if ~isvector(signal)
    error('computeEEGAmplitude:signal 必须是一维向量');
end

signal = signal(:)'; % 确保为行向量

if fs <= 0
    error('computeEEGAmplitude:采样率 fs 必须为正数');
end

% 检查cfg结构体必要字段
if ~isfield(cfg, 'freqBands')
    error('computeEEGAmplitude:cfg 缺少 freqBands 字段');
end

freqBands = cfg.freqBands;

if size(freqBands, 2) ~= 2
    error('computeEEGAmplitude:freqBands 必须是 Nx2 矩阵');
end

nBands = size(freqBands, 1);

%% 初始化输出
% 输出: [宽带振幅, 频段1振幅, ..., 频段N振幅]
amplitude = nan(1, nBands + 1);

%% 1. 计算全频带(宽带)振幅
% 直接对输入信号计算Hilbert包络 (信号已经是1-40Hz预处理后的)
try
    envelope = abs(hilbert(signal));
    amplitude(1) = mean(envelope);
catch ME
    warning('computeEEGAmplitude:宽带振幅计算失败: %s', ME.message);
end

%% 2. 计算各频段振幅
for i = 1:nBands
    lowFreq = freqBands(i, 1);
    highFreq = freqBands(i, 2);
    
    % 频率有效性检查
    if lowFreq <= 0 || highFreq >= fs/2
        warning('computeEEGAmplitude:频段 [%.1f, %.1f] Hz 无效，跳过', lowFreq, highFreq);
        continue;
    end
    
    if lowFreq >= highFreq
        warning('computeEEGAmplitude:频段下限 %.1f >= 上限 %.1f，跳过', lowFreq, highFreq);
        continue;
    end
    
    try
        % 计算归一化频率 (相对于Nyquist频率)
        Wn = [lowFreq, highFreq] / (fs / 2);
        
        % 4阶Butterworth带通滤波器
        filterOrder = 4;
        [b, a] = butter(filterOrder, Wn, 'bandpass');
        
        % 零相位滤波
        sigFilt = filtfilt(b, a, signal);
        
        % 计算Hilbert包络
        envelope = abs(hilbert(sigFilt));
        
        % 振幅 = 包络均值
        amplitude(i + 1) = mean(envelope);
        
    catch ME
        warning('computeEEGAmplitude:频段 %d [%.1f-%.1f Hz] 计算失败: %s', ...
                i, lowFreq, highFreq, ME.message);
    end
end

end
