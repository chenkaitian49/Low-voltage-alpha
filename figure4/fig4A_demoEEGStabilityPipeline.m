function plotStabilityComparison_twoPatients()
% plotStabilityComparison_twoPatients
%   读取两个指定患者的 .set 文件，在同一个 figure 内生成三块子图:
%     左侧(跨两行): 全通道原始信号 (两患者均使用各自指定的时间窗);
%                   保留通道名; 右下角 L形标尺(开口朝左上): 横线向左+竖线向上
%     右上 (A): 全导联平均 波形 + 上包络 (两患者, 波幅放大);
%                   分离式标尺: 横向时间标尺位置不变，纵向波幅标尺移至绘图区左侧
%     右下 (B): movstd 对比; 图例直接显示全段稳定性 SI
%   无任何子图小标题, 无频段范围标注。每患者三图统一颜色。
%
% 依赖: EEGLAB (pop_loadset), Signal Processing Toolbox
% 日期: 2025

%% ========== 0. 用户配置区 (★ 改这里) ==========
patient1.filepath = 'F:\EEG\LowVoltageResearch\OutputResult\V2\稳定性示意图';
patient1.filename = 'ZhuYinShui_alphaband.set';
patient1.label    = 'High Alpha Stability';
patient1.band     = [8 13];
patient1.color    = [200 90 92]/255;         % 加深的红
patient1.SI       = 3.09;                    % ★ 全段稳定性 (直接写死)
patient1.timeRange = [5 10];                  % ★ 使用的数据时间窗 [起始秒 结束秒]
patient1.stabilityTimeRange = patient1.timeRange; % ★ 右下图自动跟随同一5秒时间窗

patient2.filepath = 'F:\EEG\LowVoltageResearch\OutputResult\V2\稳定性示意图';
patient2.filename = 'DuanBingJun_alphaband.set';
patient2.label    = 'Low Alpha Stability';
patient2.band     = [8 13];
patient2.color    = [25 95 145]/255;         % 更深的蓝
patient2.SI       = 2.15;                    % ★ 全段稳定性 (直接写死)
patient2.timeRange = [5 10];                  % ★ 使用的数据时间窗 [起始秒 结束秒]
patient2.stabilityTimeRange = patient2.timeRange; % ★ 右下图自动跟随同一5秒时间窗

% ---- 图片保存路径 ----
saveDir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\稳定性示意图\结果图';
savePrefix = 'StabilityCompare';

% ---- 分析参数 ----
stabilityWin = 1;

% ---- 比例尺 ----
rawScaleAmp_uV = 10;  % ★ 左侧波幅标尺（仅改变标尺，不改变波形倍率）
rawDisplayGain = 1/3; % ★ 固定左侧波形显示倍率，保持当前波形大小不变
avgScaleAmp_uV = 10;  % 右上平均波形波幅标尺
scaleTime_s     = 2;

% ---- 图A波幅放大系数 ----
ampGain = 3;

%% ========== 1. 准备输出目录 ==========
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
    fprintf('已创建输出目录: %s\n', saveDir);
end

%% ========== 2. 读取数据 ==========
[allData1, fs1, labels1] = loadPatient(patient1);
[allData2, fs2, ~]       = loadPatient(patient2);
if fs1 ~= fs2, warning('两患者采样率不同 (%g vs %g Hz)。', fs1, fs2); end
[data1, actualRange1] = selectTimeRange(allData1, fs1, patient1.timeRange, patient1.label);
[data2, actualRange2] = selectTimeRange(allData2, fs2, patient2.timeRange, patient2.label);
[stabilityData1, actualStabilityRange1] = selectTimeRange( ...
    allData1, fs1, patient1.stabilityTimeRange, patient1.label);
[stabilityData2, actualStabilityRange2] = selectTimeRange( ...
    allData2, fs2, patient2.stabilityTimeRange, patient2.label);
fs = fs1;
winSamples = round(fs * stabilityWin);

nCh = min(size(data1,1), size(data2,1));
if size(data1,1) ~= size(data2,1)
    warning('两患者通道数不同 (%d vs %d)，仅用前 %d 通道。', size(data1,1), size(data2,1), nCh);
end
chanLabels = labels1;

%% ========== 3. 全导联平均 + 截取 ==========
avg1 = mean(double(data1), 1);
avg2 = mean(double(data2), 1);
nSampAna = min(numel(avg1), numel(avg2));
avg1 = avg1(1:nSampAna);  avg2 = avg2(1:nSampAna);
t = (0:numel(avg1)-1) / fs;
Tana = t(end);

% 右下图使用与左侧/右上相同的5秒时间窗
avgStability1 = mean(double(stabilityData1), 1);
avgStability2 = mean(double(stabilityData2), 1);
nSampStability = min(numel(avgStability1), numel(avgStability2));
avgStability1 = avgStability1(1:nSampStability);
avgStability2 = avgStability2(1:nSampStability);
tStability = (0:nSampStability-1) / fs;
Tstability = tStability(end);

%% ========== 4. 流水线 (输入已按频段滤波, 仅求包络/movstd) ==========
[env1, ~, ~, ~] = pipeline(avg1, winSamples, t);
[env2, ~, ~, ~] = pipeline(avg2, winSamples, t);
[~, mov1, ~, tMov] = pipeline(avgStability1, winSamples, tStability);
[~, mov2, ~, ~]    = pipeline(avgStability2, winSamples, tStability);

%% ========== 5. 配色 ==========
c1 = patient1.color;  c2 = patient2.color;

%% ========== 6. figure ==========
fig = figure('Name','稳定性对比 (三合一)','Color','w','Position',[60 60 1300 850]);

% ===== 左列: 全通道原始信号 =====
axRaw = subplot(2,2,[1 3]); hold(axRaw,'on');
tRaw1 = (0:size(data1,2)-1)/fs1;
displaySplitTime = tRaw1(end);
tRaw2 = (0:size(data2,2)-1)/fs2 + displaySplitTime;
rawPlotDur = tRaw2(end);

seg1 = data1(1:nCh, :);  seg2 = data2(1:nCh, :);
chanGap = 6 * median(std([seg1, seg2], 0, 2));
if chanGap == 0 || isnan(chanGap), chanGap = 1; end

yTickPos = zeros(1, nCh);
h1 = []; h2 = [];
for ch = 1:nCh
    yoff = (nCh - ch) * chanGap;
    yTickPos(ch) = yoff;
    % 波形使用固定显示倍率；修改标尺数值不会改变波形大小
    sig1 = data1(ch,:) * rawDisplayGain;
    sig2 = data2(ch,:) * rawDisplayGain;
    hh1 = plot(axRaw, tRaw1, sig1 + yoff, 'Color', c1, 'LineWidth', 0.5);
    hh2 = plot(axRaw, tRaw2, sig2 + yoff, 'Color', c2, 'LineWidth', 0.5);
    if ch == 1, h1 = hh1; h2 = hh2; end
end
xline(axRaw, displaySplitTime, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

xlim(axRaw, [0 rawPlotDur]);
yLo = min(yTickPos)-2*chanGap;  yHi = max(yTickPos)+1.5*chanGap;
ylim(axRaw, [yLo yHi]);
axis(axRaw, 'off');

% 通道名
xLabelPos = -0.015 * rawPlotDur;
for ch = 1:nCh
    text(axRaw, xLabelPos, yTickPos(ch), chanLabels{ch}, ...
         'HorizontalAlignment','right','VerticalAlignment','middle', ...
         'FontSize',7,'Color',[0.2 0.2 0.2]);
end

% 左上角标题
text(axRaw, 0, yHi, 'Filtered signal', ...
     'HorizontalAlignment','left','VerticalAlignment','top', ...
     'FontWeight','bold','FontSize',9,'Color','k');

% ---- 分离式标尺 (横=时间, 竖=波幅), 各自独立位置 ----
% 时间标尺 (横线): 起点 (x,y), 向右延伸 scaleTime_s
hBarX = max(rawPlotDur-scaleTime_s-0.05*rawPlotDur, 0.05*rawPlotDur);
hBarY = yLo+0.4*chanGap;
drawScaleBar_H(axRaw, hBarX, hBarY, scaleTime_s, sprintf('%g s', scaleTime_s));
% 波幅标尺放在横向时间标尺的右上方，并留出间距，避免二者重叠
vBarX = min(hBarX + scaleTime_s + 0.65, rawPlotDur*0.98);
vBarY = hBarY + 0.55*chanGap;
drawScaleBar_V(axRaw, vBarX, vBarY, rawScaleAmp_uV*rawDisplayGain, ...
               sprintf('%g \\muV', rawScaleAmp_uV));

legend(axRaw, [h1 h2], {patient1.label, patient2.label}, ...
       'Location','northoutside','Orientation','horizontal');

% ===== 右上 (A): 全导联平均 波形 + 上包络 =====
axA = subplot(2,2,2); hold(axA,'on');
a1 = avg1*ampGain;  e1 = env1*ampGain;
a2 = avg2*ampGain;  e2 = env2*ampGain;
offset = 1.4 * (max(e1) + max(e2));

plot(axA, t, a1 + offset, 'Color', [c1 0.45], 'LineWidth', 0.4);
plot(axA, t, e1 + offset, 'Color', c1, 'LineWidth', 2.0);
plot(axA, t, a2, 'Color', [c2 0.45], 'LineWidth', 0.4);
plot(axA, t, e2, 'Color', c2, 'LineWidth', 2.0);

% 导联平均标签 (颜色已区分两组), 抬高以与包络线拉开间距
text(axA, 0.05, offset + max(e1)*1.55, 'Average across channels', ...
     'Color', 'k', 'FontWeight','bold','FontSize',9, 'VerticalAlignment','bottom');

xlim(axA,[0 Tana]);
yl = ylim(axA); ylim(axA, [yl(1)-0.12*range(yl), yl(2)+0.18*range(yl)]);
yl = ylim(axA);
axis(axA,'off');

% ---- 分离式标尺 (横=时间, 竖=波幅), 各自独立位置 ----
hBarX2 = Tana*0.05;  hBarY2 = yl(1)+0.05*range(yl);
drawScaleBar_H(axA, hBarX2, hBarY2, scaleTime_s, sprintf('%g s', scaleTime_s));
% 纵向波幅标尺与时间标尺分离，放在波形绘图区左侧（x < 0）。
% 保持横向时间标尺位置不动；drawScaleBar_V 中关闭对象裁剪，
% 使纵向标尺在 xlim 范围之外仍能正常显示和导出。
vBarX2 = -0.12*Tana;  vBarY2 = hBarY2;
drawScaleBar_V(axA, vBarX2, vBarY2, avgScaleAmp_uV*ampGain, sprintf('%g \\muV', avgScaleAmp_uV));

% ===== 右下 (B): movstd, 图例显示全段 SI =====
axB = subplot(2,2,4); hold(axB,'on');
hb1 = plot(axB, tMov, mov1, 'Color', c1, 'LineWidth', 1.8);
hb2 = plot(axB, tMov, mov2, 'Color', c2, 'LineWidth', 1.8);
xlabel(axB,'Time (s)'); ylabel(axB,'Moving standard deviation (μV)');
xlim(axB,[0 Tstability]);
legend(axB, [hb1 hb2], ...
    {sprintf('Stability = %.2f', patient1.SI), ...
     sprintf('Stability = %.2f', patient2.SI)}, ...
    'Location','best','FontSize',8);
grid(axB,'off'); box(axB,'off');

%% ========== 7. 保存 ==========
saveFigure(fig, saveDir, sprintf('%s_Combined', savePrefix));

%% ========== 8. 命令行输出 ==========
fprintf('\n===== 稳定性对比 (图例采用全段写死值) =====\n');
fprintf('%s : Stability(全段) = %.2f\n', patient1.label, patient1.SI);
fprintf('%s : Stability(全段) = %.2f\n', patient2.label, patient2.SI);
fprintf('%s : 使用时间窗 %.3f-%.3f s\n', patient1.label, actualRange1(1), actualRange1(2));
fprintf('%s : 使用时间窗 %.3f-%.3f s\n', patient2.label, actualRange2(1), actualRange2(2));
fprintf('%s : 右下图时间窗 %.3f-%.3f s\n', patient1.label, actualStabilityRange1(1), actualStabilityRange1(2));
fprintf('%s : 右下图时间窗 %.3f-%.3f s\n', patient2.label, actualStabilityRange2(1), actualStabilityRange2(2));
fprintf('图片已保存至: %s\n', saveDir);
fprintf('============================================\n');

end

% ====================================================================
function [data, actualRange] = selectTimeRange(data, fs, timeRange, patientLabel)
    % 按患者配置截取数据，并把截取后的第一个采样点作为绘图的 0 s。
    if numel(timeRange) ~= 2 || any(~isfinite(timeRange)) || ...
            timeRange(1) < 0 || timeRange(2) <= timeRange(1)
        error('%s 的 timeRange 必须是 [起始秒 结束秒]，且满足 0 <= 起始秒 < 结束秒。', patientLabel);
    end

    firstSample = round(timeRange(1) * fs) + 1;
    lastSample  = round(timeRange(2) * fs) + 1;
    if firstSample > size(data,2)
        error('%s 的起始时间 %.3f s 超出数据长度 %.3f s。', ...
              patientLabel, timeRange(1), (size(data,2)-1)/fs);
    end
    if lastSample > size(data,2)
        warning('%s 的结束时间 %.3f s 超出数据长度，已截到 %.3f s。', ...
                patientLabel, timeRange(2), (size(data,2)-1)/fs);
        lastSample = size(data,2);
    end

    data = data(:, firstSample:lastSample);
    actualRange = ([firstSample lastSample] - 1) / fs;
end

% ====================================================================
function drawScaleBar_H(ax, x, y, tLen, tLabel)
% 横向时间标尺: 从 (x,y) 向右延伸 tLen, 标签在下方
    plot(ax, [x x+tLen], [y y], 'k', 'LineWidth', 3);
    text(ax, x+tLen/2, y, tLabel, 'HorizontalAlignment','center', ...
         'VerticalAlignment','top','FontSize',9);
end

% ====================================================================
function drawScaleBar_V(ax, x, y, aLen, aLabel)
% 竖向波幅标尺: 从 (x,y) 向上延伸 aLen, 标签在左侧(竖排)
    plot(ax, [x x], [y y+aLen], 'k', 'LineWidth', 3, 'Clipping','off');
    text(ax, x-0.025*diff(xlim(ax)), y+aLen/2, aLabel, ...
         'HorizontalAlignment','right','VerticalAlignment','middle', ...
         'FontSize',9,'Rotation',90,'Clipping','off');
end

% ====================================================================
function [data, fs, labels] = loadPatient(p)
    if ~exist(fullfile(p.filepath, p.filename), 'file')
        error('找不到文件: %s', fullfile(p.filepath, p.filename));
    end
    EEG = pop_loadset('filename', p.filename, 'filepath', p.filepath);
    fs  = EEG.srate;
    labels = {EEG.chanlocs.labels};
    if isempty(labels)
        labels = arrayfun(@(x) sprintf('Ch%d',x), 1:EEG.nbchan, 'UniformOutput', false);
    end
    if ndims(EEG.data) == 3
        [nc, np, nt] = size(EEG.data);
        data = double(reshape(EEG.data, nc, np*nt));
    else
        data = double(EEG.data);
    end
    fprintf('已读取 %s : %d 通道, %g Hz, %d 点\n', p.filename, size(data,1), fs, size(data,2));
end

% ====================================================================
function [envelope, movstdVal, sdVal, tMov] = pipeline(signal, winSamples, t)
    % 输入信号已按频段滤波, 直接求希尔伯特包络
    envelope = abs(hilbert(signal));
    movstdVal = movstd(envelope, winSamples, 'Endpoints', 'discard');
    sdVal     = std(movstdVal);
    tMov = t( (1:numel(movstdVal)) + floor((winSamples-1)/2) );
end

% ====================================================================
function saveFigure(figH, saveDir, name)
    pngPath = fullfile(saveDir, [name '.png']);
    figPath = fullfile(saveDir, [name '.fig']);
    try
        exportgraphics(figH, pngPath, 'Resolution', 300);
    catch
        print(figH, pngPath, '-dpng', '-r300');
    end
    savefig(figH, figPath);
    fprintf('已保存: %s\n', pngPath);
end
