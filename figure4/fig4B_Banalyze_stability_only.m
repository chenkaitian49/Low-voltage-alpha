%% =========================================================
%  analyze_JointROC_3Vars_3D_sign_only_Var3_600dpi.m
%
%  功能：
%    1. 从 Excel 读取结局、基线诊断、主导频率及稳定性数据
%    2. 根据主导频率计算 Dominant Stability
%    3. 比较 UWS/MCS- 与 MCS+/EMCS 两组
%    4. 仅绘制并保存 boxplot_Var3_DominantStability.png
%
%  输出分辨率：600 dpi
%% =========================================================

clc; clear; close all;

%% =========================================================
%  用户配置区
%% =========================================================

INPUT_FILE = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
OUTPUT_DIR = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅+稳定性\稳定性';

% Excel 列名（须与表头完全一致）
COL_POST_DIAG  = 'postDiagnosis_VS_1MCS__2MCS__3EMCS_4_';
COL_PRE_DIAG   = 'preDiagnosis_VS_1MCS__2MCS__3EMCS_4_';
COL_FREQ       = 'Mean_Freq_Hz';
COL_ALPHA_STAB = 'Alpha_Stability';
COL_DELTA_STAB = 'DeltaTheta_Stability';

% 分组标准：post diagnosis >= 3 为 MCS+/EMCS
GRP2_THR = 3;

% 是否剔除基线已经为 MCS+/EMCS 的患者
EXCLUDE_PRE_MCSPLUS = true;
PRE_EXCL_THR = 3;

% 图形设置
VAR3_NAME = 'Dominant Stability';
OUTPUT_NAME = 'boxplot_Var3_DominantStability.png';
EXPORT_DPI = 600;

%% =========================================================
%  1. 检查路径并读取数据
%% =========================================================

if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
    fprintf('已创建输出目录：%s\n', OUTPUT_DIR);
end

if ~isfile(INPUT_FILE)
    error('找不到输入文件：%s\n请检查 INPUT_FILE 路径。', INPUT_FILE);
end

fprintf('正在读取文件：%s\n', INPUT_FILE);
try
    opts = detectImportOptions(INPUT_FILE);
    opts.VariableNamingRule = 'preserve';
    data = readtable(INPUT_FILE, opts);
catch
    data = readtable(INPUT_FILE);
end

% 检查必要列
required_cols = {COL_POST_DIAG, COL_PRE_DIAG, COL_FREQ, ...
                 COL_ALPHA_STAB, COL_DELTA_STAB};

for i = 1:numel(required_cols)
    if ~ismember(required_cols{i}, data.Properties.VariableNames)
        error('未找到列：%s。请核对 Excel 表头。', required_cols{i});
    end
end

%% =========================================================
%  2. 提取数据并计算 Dominant Stability
%% =========================================================

post_diag  = double(data.(COL_POST_DIAG));
pre_diag   = double(data.(COL_PRE_DIAG));
freq       = double(data.(COL_FREQ));
alpha_stab = double(data.(COL_ALPHA_STAB));
delta_stab = double(data.(COL_DELTA_STAB));

% Mean_Freq_Hz >= 8 Hz：采用 Alpha Stability
% Mean_Freq_Hz <  8 Hz：采用 Delta/Theta Stability
dominant_stability = NaN(size(freq));
is_alpha = freq >= 8;
dominant_stability(is_alpha)  = alpha_stab(is_alpha);
dominant_stability(~is_alpha) = delta_stab(~is_alpha);

% 保留结局、频率和对应稳定性均非缺失的患者
valid = ~isnan(post_diag) & ~isnan(freq) & ~isnan(dominant_stability);

% 剔除基线 MCS+/EMCS
if EXCLUDE_PRE_MCSPLUS
    exclude = ~isnan(pre_diag) & pre_diag >= PRE_EXCL_THR;
    fprintf('剔除基线 MCS+/EMCS：%d 例\n', sum(valid & exclude));
    valid = valid & ~exclude;
end

post_diag = post_diag(valid);
dominant_stability = dominant_stability(valid);

if isempty(dominant_stability)
    error('筛选后没有有效数据，无法绘图。');
end

% 0 = UWS/MCS-；1 = MCS+/EMCS
grp = double(post_diag >= GRP2_THR);

x_lo = dominant_stability(grp == 0);
x_hi = dominant_stability(grp == 1);

if isempty(x_lo) || isempty(x_hi)
    error('至少有一个结局组没有有效数据，无法进行两组比较。');
end

% Wilcoxon 秩和检验，用于显著性星号
p_value = ranksum(x_hi, x_lo);
fprintf('有效样本：UWS/MCS- = %d，MCS+/EMCS = %d；p = %.6f\n', ...
    numel(x_lo), numel(x_hi), p_value);

%% =========================================================
%  3. 仅绘制 Dominant Stability 箱线图
%% =========================================================

leg_lo = 'UWS/MCS-';
leg_hi = 'MCS+/EMCS';

% 配色
c_box_lo   = [148 206 229] / 255;   % 蓝色
c_box_hi   = [233 162 163] / 255;   % 红色
c_box_lo_e = c_box_lo * 0.72;
c_box_hi_e = c_box_hi * 0.72;

% 使用像素单位创建屏幕图窗，避免 Windows 高 DPI 缩放与
% exportgraphics(figure, ...) 的输出画布尺寸发生错位。
fig = figure( ...
    'Units', 'pixels', ...
    'Position', [100 100 640 720], ...
    'Color', 'w', ...
    'Renderer', 'opengl', ...
    'InvertHardcopy', 'off');

ax = axes('Parent', fig);
hold(ax, 'on');

% 箱线图骨架
boxplot(ax, [x_lo; x_hi], ...
    [zeros(numel(x_lo),1); ones(numel(x_hi),1)], ...
    'Labels', {leg_lo, leg_hi}, ...
    'Widths', 0.48, ...
    'Symbol', '', ...
    'Colors', 'k', ...
    'MedianStyle', 'line');

% 填充箱体颜色
box_handles = findobj(ax, 'Tag', 'Box');
if numel(box_handles) >= 2
    patch(ax, get(box_handles(1), 'XData'), get(box_handles(1), 'YData'), ...
        c_box_hi, 'FaceAlpha', 0.72, 'EdgeColor', c_box_hi_e, 'LineWidth', 1.4);
    patch(ax, get(box_handles(2), 'XData'), get(box_handles(2), 'YData'), ...
        c_box_lo, 'FaceAlpha', 0.72, 'EdgeColor', c_box_lo_e, 'LineWidth', 1.4);
end

% 加粗须线和中位线
whisker_h = findobj(ax, 'Tag', 'Whisker');
set(whisker_h, 'LineWidth', 1.1);

median_h = findobj(ax, 'Tag', 'Median');
set(median_h, 'Color', [0.20 0.20 0.20], 'LineWidth', 2.2);

% 抖动散点
jitter_w = 0.12;
rng(300, 'twister');
x_jit_lo = 1 + (rand(size(x_lo)) - 0.5) * 2 * jitter_w;
x_jit_hi = 2 + (rand(size(x_hi)) - 0.5) * 2 * jitter_w;

scatter(ax, x_jit_lo, x_lo, 28, c_box_lo, 'filled', ...
    'MarkerFaceAlpha', 0.70, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.7);

scatter(ax, x_jit_hi, x_hi, 28, c_box_hi, 'filled', ...
    'MarkerFaceAlpha', 0.70, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.7);

% 计算纵轴范围
all_values = [x_lo; x_hi];
y_top = max(all_values);
y_bot = min(all_values);
y_span = y_top - y_bot;

if y_span <= eps(max(abs(all_values)))
    y_span = max(abs(y_top) * 0.10, 1);
end

bar_y = y_top + y_span * 0.07;
sig_y = y_top + y_span * 0.10;

% 仅在 p < 0.05 时标注星号
if p_value < 0.001
    sig_str = '***';
elseif p_value < 0.01
    sig_str = '**';
elseif p_value < 0.05
    sig_str = '*';
else
    sig_str = '';
end

if ~isempty(sig_str)
    line(ax, [1 1 2 2], [bar_y sig_y sig_y bar_y], ...
        'Color', 'k', ...
        'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    text(ax, 1.5, sig_y + y_span * 0.02, sig_str, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    ylim_top = sig_y + y_span * 0.20;
else
    ylim_top = y_top + y_span * 0.12;
end

% 坐标轴格式
xlim(ax, [0.5 2.5]);
ylim(ax, [y_bot - y_span * 0.05, ylim_top]);
ylabel(ax, VAR3_NAME, 'FontSize', 10.5, 'FontWeight', 'bold');
set(ax, ...
    'XTick', [1 2], ...
    'XTickLabel', {leg_lo, leg_hi}, ...
    'FontSize', 9.5, ...
    'TickDir', 'out', ...
    'LineWidth', 1.0);

title(ax, '');
grid(ax, 'off');
box(ax, 'off');

% 以固定 4.0 × 4.5 英寸画布、600 dpi 导出 PNG。
% 最终图像尺寸约为 2400 × 2700 像素。
% 使用 print 而不是 exportgraphics(fig,...)，可避免部分 Windows/MATLAB
% 版本在高 DPI 显示缩放下出现“左上角小图 + 大面积黑底”的问题。
output_file = fullfile(OUTPUT_DIR, OUTPUT_NAME);

set(fig, ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 4.0 4.5], ...
    'PaperSize', [4.0 4.5]);
set(ax, 'Color', 'w');

drawnow;
print(fig, output_file, '-dpng', ...
    sprintf('-r%d', EXPORT_DPI), '-opengl');

fprintf('Dominant Stability 箱线图已保存：%s\n', output_file);
fprintf('导出分辨率：%d dpi（约 2400 × 2700 像素）\n', EXPORT_DPI);

close(fig);
