%% demo_EEG_Amplitude_real.m
% 使用两例患者真实EEG数据绘制宽带振幅示意图
%   图A — 宽带振幅: 两组信号上下排列，无坐标轴，scalebar风格
%          蓝色=原始信号, 绿色=上/下包络, 粉红虚线=振幅均值
%          上方=患者1, 下方=患者2
%
% 数据格式: EEGLAB .set 文件 (配套 .fdt 需在同一目录)
% 通  道: O1 (枕叶左侧)
% Epoch : 自动拼接全部 epoch → 连续信号
%
% 依赖文件 (需与本脚本放在同一目录):
%   computeEEGAmplitude.m
%   EEGLAB 工具箱 (已 addpath 或在 MATLAB 路径中)
%
% 作者: Simba

clear; clc; close all;

%% ===========================
% ★ 路径配置 (修改这里)
% ===========================

% 患者1 .set 文件完整路径 (显示在图A上方)
INPUT_SET_P1 = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅示意图\波幅示意图\ZhuYinShui_epoch18_21.set';   % ← 修改

% 患者2 .set 文件完整路径 (显示在图A下方)
INPUT_SET_P2 = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅示意图\波幅示意图\QiaoWenFei_epoch4_7.set';   % ← 修改

% 输出目录 (不存在时自动创建)
OUTPUT_DIR = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅示意图\波幅示意图';         % ← 修改

% 保存格式: 'png' | 'svg' | 'pdf' | 'tiff' (留空 '' 则不保存)
SAVE_FMT = 'png';

% 分辨率 DPI (仅 png/tiff 有效)
SAVE_DPI = 300;

%% ===========================
% 0. 配色方案
% ===========================
cGreen  = [182 208 147] / 255;   % #b6d093
cOrange = [245 207 168] / 255;   % #f5cfa8  (备用)
cBlue   = [148 206 229] / 255;   % #94cee5
cPink   = [233 162 163] / 255;   % #e9a2a3
cScale  = [0.20 0.20 0.20];

% 图A统一配色
cSig  = cBlue;
cEnv  = cGreen * 0.70;
cFill = cGreen;
cAmp  = cPink  * 0.80;

%% ===========================
% 1. 参数配置
% ===========================
cfg.freqBands    = [1  4; 4  8; 8  13; 13 30; 30 40];
cfg.bandNames    = {'Delta','Theta','Alpha','Beta','Gamma'};

%% ===========================
% 2. 载入数据 & 提取 O1 通道
% ===========================
[sig_p1, fs_p1, nEpochs_p1, label_p1] = loadO1Signal(INPUT_SET_P1);
[sig_p2, fs_p2, nEpochs_p2, label_p2] = loadO1Signal(INPUT_SET_P2);

% 采样率一致性检查
if fs_p1 ~= fs_p2
    error('两例患者采样率不一致 (%d Hz vs %d Hz)，请先重采样。', fs_p1, fs_p2);
end
fs = fs_p1;

% 时间轴
N_p1 = length(sig_p1);  dur_p1 = N_p1/fs;  t_p1 = (0:N_p1-1)/fs;
N_p2 = length(sig_p2);  dur_p2 = N_p2/fs;  t_p2 = (0:N_p2-1)/fs;

fprintf('\n患者1: %s | %d epochs | %.1f s\n', label_p1, nEpochs_p1, dur_p1);
fprintf('患者2: %s | %d epochs | %.1f s\n\n', label_p2, nEpochs_p2, dur_p2);

%% ===========================
% 3. 调用 computeEEGAmplitude
% ===========================
fprintf('计算宽带振幅...\n');
amp_p1 = computeEEGAmplitude(sig_p1, fs, cfg);
amp_p2 = computeEEGAmplitude(sig_p2, fs, cfg);

Amp_broad_p1 = amp_p1(1);
Amp_broad_p2 = amp_p2(1);

upperEnv_p1 =  abs(hilbert(sig_p1));
lowerEnv_p1 = -abs(hilbert(sig_p1));
upperEnv_p2 =  abs(hilbert(sig_p2));
lowerEnv_p2 = -abs(hilbert(sig_p2));

%% ===========================
% 4. 图A — 振幅示意 (scalebar风格)
%    上方: 患者1 / 下方: 患者2
% ===========================
showSec = 5;
idxA_p1 = 1 : min(round(showSec*fs), N_p1);
idxA_p2 = 1 : min(round(showSec*fs), N_p2);
t_A_p1  = t_p1(idxA_p1);
t_A_p2  = t_p2(idxA_p2);
showSec = min([showSec, t_A_p1(end), t_A_p2(end)]);   % 实际展示时长

yRange_p1 = max(upperEnv_p1(idxA_p1));
yRange_p2 = max(upperEnv_p2(idxA_p2));
gap       = yRange_p1 * 0.5;
offset_p1 = yRange_p2 + gap + yRange_p1;   % 患者1 (上方)
offset_p2 = 0;                              % 患者2 (下方)

figA = figure('Name','振幅计算示意图 (A)', 'Color','w', 'Position',[80 420 860 500]);
ax = axes('Position',[0.08 0.08 0.88 0.84]);
hold(ax,'on');
axis(ax,'off');

% ---- 患者1 (上) ----
fill(ax, [t_A_p1, fliplr(t_A_p1)], ...
     [upperEnv_p1(idxA_p1)+offset_p1, fliplr(lowerEnv_p1(idxA_p1)+offset_p1)], ...
     cFill, 'FaceAlpha',0.22, 'EdgeColor','none');
plot(ax, t_A_p1, upperEnv_p1(idxA_p1)+offset_p1, 'Color',cEnv, 'LineWidth',1.0);
plot(ax, t_A_p1, lowerEnv_p1(idxA_p1)+offset_p1, 'Color',cEnv, 'LineWidth',1.0);
plot(ax, t_A_p1, sig_p1(idxA_p1)+offset_p1,      'Color',cSig, 'LineWidth',2.0);
plot(ax, [t_A_p1(1) t_A_p1(end)], [ Amp_broad_p1+offset_p1,  Amp_broad_p1+offset_p1], ...
     '--', 'Color',cAmp, 'LineWidth',1.8);
plot(ax, [t_A_p1(1) t_A_p1(end)], [-Amp_broad_p1+offset_p1, -Amp_broad_p1+offset_p1], ...
     '--', 'Color',cAmp, 'LineWidth',1.8);
% 患者1标签 (信号上方，留足间距避免遮挡)
text(ax, t_A_p1(1), offset_p1 + yRange_p1*1.35, ...
     'Patient 1 (HA Group)', ...
     'FontSize',10, 'Color',[0.35 0.35 0.35], 'FontWeight','bold');

% ---- 患者2 (下) ----
fill(ax, [t_A_p2, fliplr(t_A_p2)], ...
     [upperEnv_p2(idxA_p2)+offset_p2, fliplr(lowerEnv_p2(idxA_p2)+offset_p2)], ...
     cFill, 'FaceAlpha',0.22, 'EdgeColor','none');
plot(ax, t_A_p2, upperEnv_p2(idxA_p2)+offset_p2, 'Color',cEnv, 'LineWidth',1.0);
plot(ax, t_A_p2, lowerEnv_p2(idxA_p2)+offset_p2, 'Color',cEnv, 'LineWidth',1.0);
plot(ax, t_A_p2, sig_p2(idxA_p2)+offset_p2,      'Color',cSig, 'LineWidth',2.0);
plot(ax, [t_A_p2(1) t_A_p2(end)], [ Amp_broad_p2+offset_p2,  Amp_broad_p2+offset_p2], ...
     '--', 'Color',cAmp, 'LineWidth',1.8);
plot(ax, [t_A_p2(1) t_A_p2(end)], [-Amp_broad_p2+offset_p2, -Amp_broad_p2+offset_p2], ...
     '--', 'Color',cAmp, 'LineWidth',1.8);
% 患者2标签 (信号上方，留足间距避免遮挡)
text(ax, t_A_p2(1), yRange_p2*1.35, ...
     'Patient 12 (LA Group)', ...
     'FontSize',10, 'Color',[0.35 0.35 0.35], 'FontWeight','bold');

% ---- Scalebar (右下角，紧贴 Patient 12 信号末端右侧) ----
% 时间标尺: 1 s；振幅标尺: 5 μV
scaleTime = 1;
scaleAmp  = 5;
lw_sb     = 2.2;

% 时间标尺：从信号末端向右偏移，垂直居中于 Patient 12 信号底部
sb_xend = showSec + 0.15;           % 紧贴信号右端
sb_x0   = sb_xend - scaleTime;     % 向左延伸 1 s
sb_y    = -yRange_p2 * 1;        % 下移，远离信号底部避免重叠

tick_t  = yRange_p2 * 0.06;
plot(ax, [sb_x0 sb_xend], [sb_y sb_y],               '-', 'Color',cScale, 'LineWidth',lw_sb);
plot(ax, [sb_x0  sb_x0],  [sb_y-tick_t sb_y+tick_t], '-', 'Color',cScale, 'LineWidth',lw_sb);
plot(ax, [sb_xend sb_xend],[sb_y-tick_t sb_y+tick_t], '-', 'Color',cScale, 'LineWidth',lw_sb);
text(ax, (sb_x0+sb_xend)/2, sb_y - yRange_p2*0.28, '1 s', ...
     'HorizontalAlignment','center', 'FontSize',11, 'FontWeight','bold', 'Color',cScale);

% 振幅标尺：紧贴时间标尺右侧
sb_xa = sb_xend + 0.15;
sb_y1 = sb_y + scaleAmp/2;
sb_y2 = sb_y - scaleAmp/2;
plot(ax, [sb_xa sb_xa],                [sb_y2 sb_y1], '-', 'Color',cScale, 'LineWidth',lw_sb);
plot(ax, [sb_xa-0.05 sb_xa+0.05], [sb_y1 sb_y1],      '-', 'Color',cScale, 'LineWidth',lw_sb);
plot(ax, [sb_xa-0.05 sb_xa+0.05], [sb_y2 sb_y2],      '-', 'Color',cScale, 'LineWidth',lw_sb);
text(ax, sb_xa+0.12, sb_y, sprintf('%d μV', scaleAmp), ...
     'VerticalAlignment','middle', 'FontSize',11, 'FontWeight','bold', 'Color',cScale);

% ---- 范围 & 图标签 ----
xlim(ax, [-0.5, sb_xa+0.8]);
ylim(ax, [sb_y2 - yRange_p2*0.5, offset_p1 + yRange_p1*1.55]);



%% ===========================
% 5. 保存图片
% ===========================
if ~isempty(SAVE_FMT)
    if ~exist(OUTPUT_DIR, 'dir')
        mkdir(OUTPUT_DIR);
        fprintf('已创建输出目录: %s\n', OUTPUT_DIR);
    end
    fileA = fullfile(OUTPUT_DIR, sprintf('O1_FigA_Amplitude.%s',  SAVE_FMT));
    exportgraphics(figA, fileA, 'Resolution', SAVE_DPI);
    fprintf('图A已保存: %s\n', fileA);
end

%% ===========================
% 6. 控制台打印
% ===========================
fprintf('\n================================================\n');
fprintf('  计算结果汇总  —  通道: O1\n');
fprintf('================================================\n');
fprintf('                    患者1          患者2\n');
fprintf('  文件            %-16s %-16s\n', label_p1,      label_p2);
fprintf('  Epoch 数        %6d         %6d\n',    nEpochs_p1,    nEpochs_p2);
fprintf('  总时长 (s)      %8.2f       %8.2f\n',  dur_p1,        dur_p2);
fprintf('  宽带振幅 (μV)   %8.4f       %8.4f\n',  Amp_broad_p1,  Amp_broad_p2);
fprintf('================================================\n\n');

%% ===========================
% 本地函数: 载入 .set 并提取 O1 通道 (拼接全部 epoch)
% ===========================
function [sig, fs, nEpochs, label] = loadO1Signal(setPath)
    [fPath, fName, fExt] = fileparts(setPath);
    if isempty(fPath), fPath = pwd; end
    label = fName;

    fprintf('载入: %s\n', setPath);
    EEG = pop_loadset('filename', [fName fExt], 'filepath', fPath);
    fs  = EEG.srate;

    % 定位 O1 通道
    chanLabels = {EEG.chanlocs.labels};
    o1Idx = find(strcmpi(chanLabels, 'O1'));
    if isempty(o1Idx)
        % 尝试模糊匹配 (兼容 'O1.' 'EEG O1' 等变体)
        o1Idx = find(~cellfun(@isempty, regexpi(chanLabels, '^O1')));
    end
    if isempty(o1Idx)
        error('文件 %s 中未找到 O1 通道。\n当前通道列表: %s', ...
              fName, strjoin(chanLabels, ', '));
    end
    o1Idx = o1Idx(1);
    fprintf('  → 通道: %s (第 %d 个)\n', chanLabels{o1Idx}, o1Idx);

    % 提取并拼接全部 epoch
    rawData = squeeze(EEG.data(o1Idx, :, :));   % [nPnts] 或 [nPnts × nEpochs]
    if ndims(EEG.data) == 3
        nEpochs = size(EEG.data, 3);
        sig = double(reshape(rawData, 1, []));   % 拼接为行向量
    else
        nEpochs = 1;
        sig = double(rawData(:)');
    end
end
