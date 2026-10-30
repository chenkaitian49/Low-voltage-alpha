%% 四组患者EEG平均PSD对比图
%
% 功能说明:
%   - 从4个文件夹分别读取各组患者的 .set 文件
%   - 对每位患者: 先将所有电极在同一频率点上的PSD取均值 → 全脑平均PSD曲线
%   - 对每组: 收集所有患者的全脑平均PSD (转为dB) → 计算组均值曲线 + SEM阴影带
%   - 在同一坐标系内绘制4组的均值曲线 + SEM阴影带
%
% 修改说明 (v2):
%   - PSD统计在log₁₀(dB)空间进行，SEM上下对称
%   - 删除横向网格虚线
%   - 频段标签移至X轴上方，字体加大加粗，并加频段分界竖线
%
% 依赖:
%   EEGLAB (需要已添加到MATLAB路径中)
%   Signal Processing Toolbox (pwelch)
%
% 作者: Simba
% 日期: 2025

clear; clc; close all;

%% ==================== 用户配置区 ====================

% ---- 四组输入文件夹路径 (每组对应一个文件夹) ----
groupDirs = {
    'F:\EEG\LowVoltageResearch\PreprocessedData\highAlpha_ICA';
    'F:\EEG\LowVoltageResearch\PreprocessedData\highThetaDelta_ICA';
    'F:\EEG\LowVoltageResearch\PreprocessedData\lowAlpha_ICA1';
    'F:\EEG\LowVoltageResearch\PreprocessedData\lowThetaDelta_ICA';
};

% ---- 组名标签 (用于图例) ----
groupNames = {'HA', 'HTD', 'LA', 'LTD'};

% ---- 输出文件夹路径 ----
outputDir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\全频段功率谱\PSD示意图';

% ---- 分析参数 ----
freqRange = [1, 40];   % 显示的频率范围 (Hz)
welchWin  = 2;         % Welch窗口长度 (秒)，频率分辨率 = 1/welchWin Hz

% ---- 配色方案 (每组一种颜色) ----
groupColors = [
    182, 208, 147;   % #b6d093 - 第1组
    245, 207, 168;   % #f5cfa8 - 第2组
    148, 206, 229;   % #94cee5 - 第3组
    233, 162, 163;   % #e9a2a3 - 第4组
] / 255;

% SEM阴影透明度
shadingAlpha = 0.30;

%% ==================== 检查路径 ====================
nGroups = length(groupDirs);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('已创建输出目录: %s\n', outputDir);
end

for g = 1:nGroups
    if ~exist(groupDirs{g}, 'dir')
        error('第%d组输入目录不存在: %s\n请修改 groupDirs 变量', g, groupDirs{g});
    end
end

%% ==================== 逐组处理 ====================
groupPSD_matrix = cell(nGroups, 1);
freqVec = [];

for g = 1:nGroups
    fprintf('\n========================================\n');
    fprintf('=== 处理第 %d 组: %s ===\n', g, groupNames{g});
    fprintf('========================================\n');

    fileList = dir(fullfile(groupDirs{g}, '*.set'));
    nFiles   = length(fileList);

    if nFiles == 0
        error('第%d组文件夹中未找到任何 .set 文件: %s', g, groupDirs{g});
    end
    fprintf('找到 %d 个 .set 文件\n', nFiles);

    patientPSD_list = [];

    for fileIdx = 1:nFiles
        fileName = fileList(fileIdx).name;
        fprintf('  [%d/%d] %s ... ', fileIdx, nFiles, fileName);

        try
            % ---- 加载 .set 文件 ----
            EEG = pop_loadset('filename', fileName, 'filepath', groupDirs{g});

            fs        = EEG.srate;
            nChannels = EEG.nbchan;
            nEpochs   = EEG.trials;
            nPoints   = EEG.pnts;

            % ---- 拼接为连续信号 [channels × samples] ----
            if nEpochs > 1
                contData = reshape(EEG.data, nChannels, nPoints * nEpochs);
            else
                contData = EEG.data;
            end

            % ---- Welch参数 ----
            M        = size(contData, 2);
            L        = round(welchWin * fs);
            noverlap = round(0.5 * L);
            nfft     = max(L, 2^nextpow2(L));

            if L > M
                warning('患者 %s 信号过短，改用全段信号', fileName);
                L        = M;
                noverlap = 0;
                nfft     = max(L, 2^nextpow2(L));
            end

            % ---- 逐通道计算PSD，通道间取均值 ----
            [pxx_ch1, f_tmp] = pwelch(double(contData(1, :)), hamming(L), noverlap, nfft, fs);

            nFreqPoints   = length(f_tmp);
            chanPSD        = zeros(nChannels, nFreqPoints);
            chanPSD(1, :)  = pxx_ch1';

            for ch = 2:nChannels
                [pxx_ch, ~]    = pwelch(double(contData(ch, :)), hamming(L), noverlap, nfft, fs);
                chanPSD(ch, :) = pxx_ch';
            end

            % 全脑平均PSD (线性空间)
            brainAvgPSD = mean(chanPSD, 1);   % 1 × nFreqPoints

            % ---- 统一频率向量 ----
            if isempty(freqVec)
                freqVec = f_tmp';
            else
                if length(f_tmp) ~= length(freqVec) || max(abs(f_tmp' - freqVec)) > 1e-6
                    brainAvgPSD = interp1(f_tmp, brainAvgPSD, freqVec, 'linear', 'extrap');
                    warning('患者 %s 频率向量不一致，已插值对齐', fileName);
                end
            end

            patientPSD_list = [patientPSD_list; brainAvgPSD];   %#ok<AGROW>
            fprintf('完成 (%d导联, %d样本)\n', nChannels, M);

        catch ME
            fprintf('失败: %s\n', ME.message);
            warning('跳过文件 %s: %s', fileName, ME.message);
        end
    end

    if isempty(patientPSD_list)
        error('第%d组没有成功处理任何文件', g);
    end

    groupPSD_matrix{g} = patientPSD_list;   % [nPatients × nFreqPoints] 线性PSD
    fprintf('第%d组: 共%d位患者PSD已收集\n', g, size(patientPSD_list, 1));
end

%% ==================== dB空间计算均值与SEM ====================
fprintf('\n转换为dB并计算各组均值与SEM...\n');

freqMask = (freqVec >= freqRange(1)) & (freqVec <= freqRange(2));
fPlot    = freqVec(freqMask);

groupMean = cell(nGroups, 1);
groupSEM  = cell(nGroups, 1);

for g = 1:nGroups
    % 先转换为dB (10·log₁₀)，再计算统计量 → SEM在dB空间对称
    psdMat_dB    = 10 * log10(groupPSD_matrix{g}(:, freqMask));
    nPat         = size(psdMat_dB, 1);
    groupMean{g} = mean(psdMat_dB, 1);
    groupSEM{g}  = std(psdMat_dB, 0, 1) / sqrt(nPat);
    fprintf('第%d组 (%s): %d位患者\n', g, groupNames{g}, nPat);
end

%% ==================== 绘图 ====================
fprintf('\n绘制PSD对比图...\n');

% ---- 频段定义 (只标注 delta / theta / alpha) ----
% bandEdges 前4个点定义三个有标注的频段 (1-4, 4-8, 8-13)
% 最后保留40Hz边界仅供竖线和背景用
bandEdges  = [1, 4, 8, 13, 40];
bandLabels = {'\delta', '\theta', '\alpha'};   % 只标注前三段
nBands     = length(bandLabels);               % = 3

% ---- 创建图窗，留出X轴上方空间用于频段标签 ----
fig = figure('Units', 'centimeters', 'Position', [2, 2, 22, 15], ...
             'Color', 'white');

% 手动设置axes位置，在顶部留出空间给频段标签
ax = axes('Parent', fig, ...
          'Units', 'normalized', ...
          'Position', [0.10, 0.11, 0.86, 0.72]);   % [left, bottom, width, height]

hold(ax, 'on');

% ---- 先画频段背景色块（最底层）----
% delta(1-4Hz) + theta(4-8Hz): 暖黄阴影
% alpha(8-13Hz): 蓝紫阴影
% beta/gamma(13-40Hz): 纯白无阴影

% 先确定y范围（估算dB范围用于背景）
allMeans = cell2mat(groupMean');
yMin = min(allMeans(:)) - 5;
yMax = max(allMeans(:)) + 5;

% 三种背景色：delta/theta 共用中灰，alpha 用浅灰，其余白色
bgColorMap = {
    [0.82, 0.82, 0.82];   % delta  (1-4 Hz)  — 中灰
    [0.82, 0.82, 0.82];   % theta  (4-8 Hz)  — 中灰（同delta）
    [0.91, 0.91, 0.91];   % alpha  (8-13 Hz) — 浅灰
    [1.00, 1.00, 1.00];   % beta/gamma (13-40 Hz) — 白色（无阴影效果）
};
bgAlphaMap = [1.0, 1.0, 1.0, 1.0];   % 全部不透明，靠灰度深浅区分

% bandEdges 现在是 [1, 4, 8, 13, 40]，共4段
for b = 1:length(bandEdges)-1
    xL = bandEdges(b);
    xR = bandEdges(b + 1);
    fill(ax, [xL xR xR xL], [yMin yMin yMax yMax], ...
        bgColorMap{b}, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', bgAlphaMap(b), ...
        'HandleVisibility', 'off');
end

% ---- 绘制SEM阴影带 + 均值曲线 ----
lineHandles = gobjects(nGroups, 1);

for g = 1:nGroups
    mu  = groupMean{g};
    sem = groupSEM{g};
    col = groupColors(g, :);

    upper = mu + sem;
    lower = mu - sem;

    xFill = [fPlot, fliplr(fPlot)];
    yFill = [upper, fliplr(lower)];

    fill(ax, xFill, yFill, col, ...
        'FaceAlpha', shadingAlpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    % 均值曲线 (颜色比阴影深)
    lineHandles(g) = plot(ax, fPlot, mu, ...
        'Color', col * 0.72, ...
        'LineWidth', 2.2);
end

% ---- 频段分界竖线 (只在 4Hz、8Hz、13Hz 处画线) ----
for b = 2:nBands   % b=2→4Hz, b=3→8Hz, b=4→13Hz（nBands=3，所以到13Hz）
    xline(ax, bandEdges(b), '-', ...
        'Color', [0.78 0.78 0.78], ...
        'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end

% ---- 坐标轴范围与标签 ----
xlim(ax, [freqRange(1), freqRange(2)]);
ylim(ax, [yMin, yMax]);

xlabel(ax, 'Frequency (Hz)', 'FontSize', 13, 'FontName', 'Arial');
ylabel(ax, 'PSD (dB)', 'FontSize', 13, 'FontName', 'Arial');
% title(ax, 'Group-Averaged EEG Power Spectral Density', ...
%     'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

% ---- 关闭网格（不要横向虚线）----
set(ax, 'Box', 'off', ...
    'TickDir', 'out', ...
    'FontSize', 11, ...
    'FontName', 'Arial', ...
    'LineWidth', 1.0, ...
    'XGrid', 'off', ...
    'YGrid', 'off');   % ★ 完全关闭网格

% ---- 图例 ----
legend(ax, lineHandles, groupNames, ...
    'Location', 'northeast', ...
    'FontSize', 11, ...
    'FontName', 'Arial', ...
    'Box', 'off');

hold(ax, 'off');

% ---- 频段标签：绘制在X轴上方（axes外部，图窗坐标系）----
% 将每个频段中点的数据坐标转换为归一化图窗坐标
axPos  = ax.Position;   % [left, bottom, width, height] in normalized units
xRange = diff(freqRange);
yLabelNorm = axPos(2) + axPos(4) + 0.005;   % 紧贴axes顶部稍微上方

for b = 1:nBands
    xMid     = (bandEdges(b) + bandEdges(b+1)) / 2;
    % 将xMid从数据坐标转换为图窗归一化坐标
    xNorm = axPos(1) + axPos(3) * (xMid - freqRange(1)) / xRange;

    annotation('textbox', ...
        [xNorm - 0.04, yLabelNorm, 0.08, 0.06], ...
        'String', bandLabels{b}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 14, ...
        'FontName', 'Arial', ...
        'FontWeight', 'bold', ...
        'Color', [0.35 0.35 0.35], ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'Interpreter', 'tex');
end

%% ==================== 保存图片 ====================
outputFig_png = fullfile(outputDir, 'Group_PSD_Comparison.png');
outputFig_svg = fullfile(outputDir, 'Group_PSD_Comparison.svg');

exportgraphics(fig, outputFig_png, 'Resolution', 300, 'BackgroundColor', 'white');
print(fig, outputFig_svg, '-dsvg');

fprintf('\n图片已保存:\n  PNG: %s\n  SVG: %s\n', outputFig_png, outputFig_svg);
fprintf('完成!\n');