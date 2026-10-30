%% =========================================================
%  analyze_AmpConsistency.m
%
%  波幅目视/计算一致性分析：
%    1. Freq + Stability scatter with amplitude inconsistency markers (VhCl / VlCh)
%    2. Amp  + Stability scatter with VC consistency markers (VaCn / VnCa)
%
%  分组标准：
%    Post Diagnosis: VS/MCS- (1或2) vs MCS+/EMCS (3或4)
%
%  ★ 使用前仅需修改第一节"用户配置区"即可 ★
%
%  输出文件：
%    scatter_amp_FreqStab.png  — Freq + Stability scatter (amplitude inconsistency markers)
%    scatter_vc_AmpStab.png    — Amp  + Stability scatter (VC consistency markers)
%    results_summary_amp.txt   — amplitude consistency statistics
%% =========================================================

clc; clear; close all;

%% =========================================================
%  ★★★  用户配置区（仅需修改这里）★★★
%% =========================================================

INPUT_FILE = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
OUTPUT_DIR = 'F:\EEG\LowVoltageResearch\OutputResult\V2\目视和计算不一致比较';

% ── 列名（与 Excel 表头完全一致，区分大小写）──────────────
COL_POST_DIAG    = 'postDiagnosis_VS_1MCS__2MCS__3EMCS_4_';
COL_VAR1         = 'Mean_Freq_Hz';
COL_VAR2         = 'Broadband_Amplitude_Mean';
COL_ALPHA_STAB   = 'Alpha_Stability';
COL_DELTA_STAB   = 'DeltaTheta_Stability';
COL_VC_CONSIS    = 'VC_consistency_consis_0VaCn_1VnCa_2_';
COL_GROUP        = 'Group_HA_1HTD_2LA_3LTD_4_';
COL_NAME         = 'Name';                              % patient name/ID column
%  Dominant_Stability 根据 Mean_Freq_Hz 自动计算：
%    Mean_Freq_Hz >= 8 → Alpha_Stability，否则 → DeltaTheta_Stability

% ── 分组编码 ──────────────────────────────────────────────
GRP2_THR = 3;   % Post Diagnosis >= GRP2_THR → MCS+/EMCS (1)

% ── 图轴标注用变量名 ──────────────────────────────────────
VAR1_NAME = 'Dominant Freq (Hz)';
VAR2_NAME = 'Broadband Amp (\muV)';
VAR3_NAME = 'Dominant Stability';

% ── 波幅分界值 ────────────────────────────────────────────
% 留空时，根据目视高/低波幅人数比例自动确定阈值；
% 若需要固定阈值，也可直接填写数值（例如 AMP_THR = 3.9321）。
AMP_THR = [];   % [] = 与目视高/低波幅人数比例一致

%% =========================================================
%  以下代码无需修改
%% =========================================================

%% ── 0. 路径检查与创建 ────────────────────────────────────
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
    fprintf('Output directory created: %s\n', OUTPUT_DIR);
end
if ~isfile(INPUT_FILE)
    error('Input file not found: %s', INPUT_FILE);
end

%% ── 1. 读取数据 ──────────────────────────────────────────
fprintf('Reading: %s\n', INPUT_FILE);
try
    opts = detectImportOptions(INPUT_FILE);
    opts.VariableNamingRule = 'preserve';
    data = readtable(INPUT_FILE, opts);
catch
    data = readtable(INPUT_FILE);
end
fprintf('Detected %d columns, %d rows\n', width(data), height(data));

% ── 检查必要列 ────────────────────────────────────────────
required_cols = {COL_POST_DIAG, COL_VAR1, COL_VAR2, ...
                 COL_ALPHA_STAB, COL_DELTA_STAB, COL_VC_CONSIS, COL_GROUP, COL_NAME};
col_display   = {'COL_POST_DIAG','COL_VAR1','COL_VAR2', ...
                 'COL_ALPHA_STAB','COL_DELTA_STAB','COL_VC_CONSIS','COL_GROUP','COL_NAME'};
for ci = 1:length(required_cols)
    if ~ismember(required_cols{ci}, data.Properties.VariableNames)
        error('Column "%s" (%s) not found. Check column names.', ...
              required_cols{ci}, col_display{ci});
    end
end
fprintf('Column check passed.\n\n');

% ── 提取各列 ──────────────────────────────────────────────
post_diag  = double(data.(COL_POST_DIAG));
var1       = double(data.(COL_VAR1));
var2       = double(data.(COL_VAR2));
alpha_stab = double(data.(COL_ALPHA_STAB));
delta_stab = double(data.(COL_DELTA_STAB));
vc_consis  = double(data.(COL_VC_CONSIS));
group_vis  = double(data.(COL_GROUP));
pat_name   = data.(COL_NAME);

% ── 计算 Dominant_Stability ───────────────────────────────
var3 = zeros(size(var1));
is_alpha = var1 >= 8;
var3( is_alpha) = alpha_stab( is_alpha);
var3(~is_alpha) = delta_stab(~is_alpha);
fprintf('Dominant_Stability: %d alpha-dominant, %d delta/theta-dominant\n', ...
    sum(is_alpha), sum(~is_alpha));

% ── 保留有效行 ────────────────────────────────────────────
valid = ~isnan(var1) & ~isnan(var2) & ~isnan(var3) ...
      & ~isnan(alpha_stab) & ~isnan(delta_stab) ...
      & ~isnan(post_diag) & ~isnan(group_vis);
post_diag = post_diag(valid);
var1      = var1(valid);
var2      = var2(valid);
var3      = var3(valid);
vc_consis = vc_consis(valid);
group_vis = group_vis(valid);
pat_name  = pat_name(valid);
N = sum(valid);
fprintf('Valid samples: %d\n\n', N);

%% ── 2. 预后分组 ──────────────────────────────────────────
grp    = double(post_diag >= GRP2_THR);   % 1=MCS+/EMCS, 0=UWS/MCS-
leg_lo = 'UWS/MCS-';
leg_hi = 'MCS+/EMCS';

%% ── 2b. 按目视高/低波幅人数比例确定阈值 ─────────────────
% 目视高波幅：HA/HTD（编码 1/2）；目视低波幅：LA/LTD（编码 3/4）。
% 采用排序后目标分位点两侧数值的中点，使计算高波幅人数与目视高波幅人数一致。
vis_amp_high = ismember(group_vis, [1, 2]);
n_vis_amp_high = sum(vis_amp_high);

if isempty(AMP_THR)
    sorted_amp = sort(var2, 'descend');

    if n_vis_amp_high == 0
        % 没有目视高波幅样本：阈值放在最大值之上。
        AMP_THR = max(var2) + max(eps(max(var2)), eps);
    elseif n_vis_amp_high == N
        % 全部为目视高波幅：阈值取最小值，保证全部计为高波幅。
        AMP_THR = min(var2);
    else
        amp_high_edge = sorted_amp(n_vis_amp_high);
        amp_low_edge  = sorted_amp(n_vis_amp_high + 1);

        if amp_high_edge > amp_low_edge
            AMP_THR = (amp_high_edge + amp_low_edge) / 2;
        else
            % 并列值跨越分界时，单一数值阈值无法精确拆分并列样本。
            % 取并列值作为阈值，并在日志中报告实际人数差异。
            AMP_THR = amp_high_edge;
            warning(['Automatic amplitude threshold falls on tied values ' ...
                     '(%.6g); exact visual proportion cannot be separated ' ...
                     'with a single >= threshold.'], AMP_THR);
        end
    end
end

n_algo_amp_high = sum(var2 >= AMP_THR);
fprintf('Amplitude threshold: %.6f (visual high=%d/%d; calculated high=%d/%d)\n', ...
    AMP_THR, n_vis_amp_high, N, n_algo_amp_high, N);

%% ── 3. 颜色定义 ──────────────────────────────────────────
c_green = [182 208 147] / 255;   % VS/MCS-
c_red   = [233 162 163] / 255;   % MCS+/EMCS

%% ── 4. 日志文件 ──────────────────────────────────────────
log_file = fullfile(OUTPUT_DIR, 'results_summary_amp.txt');
fid = fopen(log_file, 'w');
fprintf(fid, 'Amplitude consistency analysis\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Valid samples N=%d\n', N);
fprintf(fid, 'Grouping: Post Diagnosis >= %d -> MCS+/EMCS (1), else VS/MCS- (0)\n', GRP2_THR);
fprintf(fid, 'Amplitude threshold: %.4f\n', AMP_THR);
fprintf(fid, 'Visual high/low count: %d/%d (%.2f%%/%.2f%%)\n', ...
    n_vis_amp_high, N - n_vis_amp_high, ...
    100 * n_vis_amp_high / N, 100 * (N - n_vis_amp_high) / N);
fprintf(fid, 'Calculated high/low count: %d/%d\n\n', ...
    n_algo_amp_high, N - n_algo_amp_high);

%% ── 5. 波幅一致性分类 ────────────────────────────────────
%
%  算法高/低：Broadband_Amplitude_Mean >= AMP_THR
%  目视高/低：Group HA(1)/HTD(2) = visual-high; LA(3)/LTD(4) = visual-low
%
%  VhCl (amp_incon=1): visual-high & calc-low
%  VlCh (amp_incon=2): visual-low  & calc-high
%  Consistent (amp_incon=0)

algo_amp_high = (var2 >= AMP_THR);

amp_incon = zeros(N, 1);
amp_incon( vis_amp_high & ~algo_amp_high) = 1;   % VhCl
amp_incon(~vis_amp_high &  algo_amp_high) = 2;   % VlCh

n_VhCl   = sum(amp_incon == 1);
n_VlCh   = sum(amp_incon == 2);
n_consis = sum(amp_incon == 0);

fprintf('Amplitude consistency:\n');
fprintf('  Consistent       : %d\n', n_consis);
fprintf('  VhCl (vis-H, calc-L): %d\n', n_VhCl);
fprintf('  VlCh (vis-L, calc-H): %d\n\n', n_VlCh);
fprintf(fid, 'Amplitude consistency:\n');
fprintf(fid, '  Consistent           : %d\n', n_consis);
fprintf(fid, '  VhCl (vis-H, calc-L): %d\n', n_VhCl);
fprintf(fid, '  VlCh (vis-L, calc-H): %d\n\n', n_VlCh);

%% ── 5b. 输出不一致患者名单 ───────────────────────────────

% VhCl 患者
names_VhCl = pat_name(amp_incon == 1);
names_VlCh = pat_name(amp_incon == 2);
names_VaCn = pat_name(vc_consis == 1);
names_VnCa = pat_name(vc_consis == 2);

% 对齐最大长度，空缺填 ""
max_len = max([numel(names_VhCl), numel(names_VlCh), ...
               numel(names_VaCn), numel(names_VnCa)]);
pad = @(x) [x; repmat({''}, max_len - numel(x), 1)];

out_tbl = table(pad(names_VhCl), pad(names_VlCh), ...
                pad(names_VaCn), pad(names_VnCa), ...
    'VariableNames', {'VhCl', 'VlCh', 'VaCn', 'VnCa'});

out_xlsx = fullfile(OUTPUT_DIR, 'inconsistent_patients.xlsx');
writetable(out_tbl, out_xlsx);
fprintf('Inconsistent patient list saved: inconsistent_patients.xlsx\n');
fprintf(fid, 'Inconsistent patient list saved: inconsistent_patients.xlsx\n');

%% ── 6. 绘图 ──────────────────────────────────────────────
%
%  图1: Freq (V1) + Dominant Stability (V3)  → scatter_amp_FreqStab.png
%       Overlay: amplitude inconsistency markers (VhCl / VlCh)
%       (corresponds to original scatter_amp_pair13)
%
%  图2: Amp  (V2) + Dominant Stability (V3)  → scatter_vc_AmpStab.png
%       Overlay: VC consistency markers (VaCn / VnCa)
%       (corresponds to original scatter_2d_pair23)

% ── 图1: Freq + Stability，波幅一致性标记 ─────────────────
fig1 = figure('Position',[100, 100, 580, 500], 'Color','w');
ax1  = axes('Parent', fig1);
hold(ax1, 'on');
draw_amp_scatter(ax1, var1, var3, grp, amp_incon, c_green, c_red, ...
    VAR1_NAME, VAR3_NAME, leg_lo, leg_hi, n_VhCl, n_VlCh);
exportgraphics(fig1, fullfile(OUTPUT_DIR, 'scatter_amp_FreqStab.png'), 'Resolution', 200);
fprintf('Saved: scatter_amp_FreqStab.png\n');

% ── 图2: Amp + Stability，VC一致性标记 ───────────────────
fig2 = figure('Position',[180, 160, 580, 500], 'Color','w');
ax2  = axes('Parent', fig2);
hold(ax2, 'on');
draw_vc_scatter(ax2, var2, var3, grp, vc_consis, c_green, c_red, ...
    VAR2_NAME, VAR3_NAME, leg_lo, leg_hi);
exportgraphics(fig2, fullfile(OUTPUT_DIR, 'scatter_vc_AmpStab.png'), 'Resolution', 200);
fprintf('Saved: scatter_vc_AmpStab.png\n');

%% ── 7. 关闭日志 ──────────────────────────────────────────
fclose(fid);
fprintf('\nResults saved to: %s\nDone.\n', log_file);

%% =========================================================
%  局部函数
%% =========================================================

function draw_amp_scatter(ax, xa, xb, grp, amp_incon, c_lo, c_hi, ...
        name_a, name_b, leg_lo, leg_hi, n_VhCl, n_VlCh)
%DRAW_AMP_SCATTER  2D scatter plot with amplitude inconsistency markers
%
%  Outcome groups (filled markers, background layer):
%    grp==0  VS/MCS-    green circle
%    grp==1  MCS+/EMCS  red triangle
%
%  Amplitude inconsistency (hollow markers, foreground layer):
%    amp_incon==1  VhCl  hollow circle,  dark orange border
%    amp_incon==2  VlCh  hollow diamond, dark purple border

    % ── Background: outcome groups ───────────────────────
    h_lo = scatter(ax, xa(grp==0), xb(grp==0), 45, c_lo, 'o', 'filled', ...
        'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
    h_hi = scatter(ax, xa(grp==1), xb(grp==1), 45, c_hi, '^', 'filled', ...
        'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);

    % ── Foreground: amplitude inconsistency markers ───────
    % VhCl (visual-high, calc-low): hollow circle, dark orange
    idx_VhCl = (amp_incon == 1);
    if any(idx_VhCl)
        h_VhCl = scatter(ax, xa(idx_VhCl), xb(idx_VhCl), 90, 'o', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.85 0.33 0.10], ...
            'LineWidth', 1.8);
    else
        h_VhCl = scatter(ax, NaN, NaN, 90, 'o', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.85 0.33 0.10], ...
            'LineWidth', 1.8);
    end

    % VlCh (visual-low, calc-high): hollow diamond, dark purple
    idx_VlCh = (amp_incon == 2);
    if any(idx_VlCh)
        h_VlCh = scatter(ax, xa(idx_VlCh), xb(idx_VlCh), 90, 'd', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.49 0.18 0.56], ...
            'LineWidth', 1.8);
    else
        h_VlCh = scatter(ax, NaN, NaN, 90, 'd', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.49 0.18 0.56], ...
            'LineWidth', 1.8);
    end

    % ── Axes labels ───────────────────────────────────────
    xlabel(ax, name_a, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(ax, name_b, 'FontSize', 10, 'FontWeight', 'bold');

    % ── Legend ────────────────────────────────────────────
    legend(ax, [h_lo, h_hi, h_VhCl, h_VlCh], ...
        {leg_lo, leg_hi, ...
         sprintf('VhCl (n=%d)', n_VhCl), ...
         sprintf('VlCh (n=%d)', n_VlCh)}, ...
        'Location', 'NorthEast', 'FontSize', 8);

    % ── Axes formatting ───────────────────────────────────
    grid(ax, 'off');
    box(ax, 'on');
    xlim(ax, [min(xa) - range(xa)*0.08, max(xa) + range(xa)*0.08]);
    ylim(ax, [min(xb) - range(xb)*0.08, max(xb) + range(xb)*0.08]);
end

function draw_vc_scatter(ax, xa, xb, grp, vc_consis, c_lo, c_hi, ...
        name_a, name_b, leg_lo, leg_hi)
%DRAW_VC_SCATTER  2D scatter plot with VC consistency markers
%
%  Outcome groups (filled markers, background layer):
%    grp==0  VS/MCS-    green circle
%    grp==1  MCS+/EMCS  red triangle
%
%  VC consistency (hollow markers, foreground layer):
%    vc_consis==1  VaCn  hollow circle,  dark blue border
%    vc_consis==2  VnCa  hollow diamond, dark red border

    % ── Background: outcome groups ───────────────────────
    h_lo = scatter(ax, xa(grp==0), xb(grp==0), 45, c_lo, 'o', 'filled', ...
        'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
    h_hi = scatter(ax, xa(grp==1), xb(grp==1), 45, c_hi, '^', 'filled', ...
        'MarkerFaceAlpha', 0.80, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);

    % ── Foreground: VC consistency markers ───────────────
    % VaCn (vc=1): hollow circle, dark blue
    idx_VaCn = (vc_consis == 1);
    if any(idx_VaCn)
        h_VaCn = scatter(ax, xa(idx_VaCn), xb(idx_VaCn), 80, 'o', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.08 0.17 0.55], ...
            'LineWidth', 1.8);
    else
        h_VaCn = scatter(ax, NaN, NaN, 80, 'o', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.08 0.17 0.55], ...
            'LineWidth', 1.8);
    end

    % VnCa (vc=2): hollow diamond, dark red
    idx_VnCa = (vc_consis == 2);
    if any(idx_VnCa)
        h_VnCa = scatter(ax, xa(idx_VnCa), xb(idx_VnCa), 80, 'd', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.60 0.05 0.05], ...
            'LineWidth', 1.8);
    else
        h_VnCa = scatter(ax, NaN, NaN, 80, 'd', ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.60 0.05 0.05], ...
            'LineWidth', 1.8);
    end

    % ── Axes labels ───────────────────────────────────────
    xlabel(ax, name_a, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(ax, name_b, 'FontSize', 10, 'FontWeight', 'bold');

    % ── Legend ────────────────────────────────────────────
    n_VaCn = sum(idx_VaCn);
    n_VnCa = sum(idx_VnCa);
    legend(ax, [h_lo, h_hi, h_VaCn, h_VnCa], ...
        {leg_lo, leg_hi, ...
         sprintf('VaCn (n=%d)', n_VaCn), ...
         sprintf('VnCa (n=%d)', n_VnCa)}, ...
        'Location', 'NorthEast', 'FontSize', 8);

    % ── Axes formatting ───────────────────────────────────
    grid(ax, 'off');
    box(ax, 'on');
    xlim(ax, [min(xa) - range(xa)*0.08, max(xa) + range(xa)*0.08]);
    ylim(ax, [min(xb) - range(xb)*0.08, max(xb) + range(xb)*0.08]);
end
