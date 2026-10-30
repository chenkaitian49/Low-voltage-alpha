%% Slope Plot: Diagnosis by Group — with Statistical Annotations
% Groups: HA=1, HTD=2, LA=3, LTD=4
% Colors: #b6d093 (HA), #f5cfa8 (HTD), #94cee5 (LA), #e9a2a3 (LTD)
%
% Three-step statistical framework (ordinal, non-parametric):
%   [1] KW on PRE      → baseline comparability
%   [2] Wilcoxon SR    → within-group intervention effect
%   [3] KW + Dunn's    → inter-group difference in CHANGE score
%
%   Dunn's p-values Bonferroni-corrected (6 pairs)
%   Stars: *p<0.05  **p<0.01  ***p<0.001  ns: not significant

clear; clc;

%% ===== USER-DEFINED PATHS =====
input_csv  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
output_dir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\量表评分+诊断（剔除无eeg）';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ===== Load Data =====
data = readtable(input_csv, 'PreserveVariableNames', true);

group_col  = double(data{:, 'Group(HA=1 HTD=2 LA=3 LTD=4)'});
pre_diag   = double(data{:, 'pre Diagnosis ( VS=1 MCS-=2 MCS+=3 EMCS=4)'});
post_diag  = double(data{:, 'post Diagnosis ( VS=1 MCS-=2 MCS+=3 EMCS=4)'});
delta_diag = post_diag - pre_diag;   % change score (Delta)

%% ===== Color & Label Definitions =====
colors = [
    182/255, 208/255, 147/255;   % HA   #b6d093
    245/255, 207/255, 168/255;   % HTD  #f5cfa8
    148/255, 206/255, 229/255;   % LA   #94cee5
    233/255, 162/255, 163/255;   % LTD  #e9a2a3
];
line_colors  = colors * 0.60;
group_labels = {'HA', 'HTD', 'LA', 'LTD'};
group_ids    = [1, 2, 3, 4];

x_pre  = [1, 3, 5, 7];
x_post = [2, 4, 6, 8];

all_xticks  = reshape([x_pre; x_post], 1, []);
tick_labels = {};
for g = 1:4
    tick_labels{end+1} = sprintf('Pre %s',  group_labels{g});
    tick_labels{end+1} = sprintf('Post %s', group_labels{g});
end

%% ================================================================
%% BLOCK A: Statistical Analysis + TXT export
%% ================================================================

txt_path = fullfile(output_dir, 'slope_diagnosis_stats.txt');
if exist(txt_path, 'file'), delete(txt_path); end
diary(txt_path);
diary on;

fprintf('Statistical Report: Diagnosis Level (Ordinal)\n');
fprintf('Generated: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Input:     %s\n', input_csv);
fprintf('============================================================\n\n');

%--- Extract per-group vectors ---%
pre_by_grp   = cell(4,1);
post_by_grp  = cell(4,1);
delta_by_grp = cell(4,1);
n_grp        = zeros(4,1);

fprintf('-- Descriptive Statistics --\n');
for g = 1:4
    idx = group_col == g;
    pre_by_grp{g}   = pre_diag(idx);
    post_by_grp{g}  = post_diag(idx);
    delta_by_grp{g} = delta_diag(idx);
    n_grp(g)        = sum(idx);
    fprintf('  Group %s: n=%d | Pre median=%.1f [IQR %.1f-%.1f] | Post median=%.1f [IQR %.1f-%.1f] | Delta median=%.1f\n', ...
        group_labels{g}, n_grp(g), ...
        median(pre_by_grp{g}),   quantile(pre_by_grp{g},  0.25), quantile(pre_by_grp{g},  0.75), ...
        median(post_by_grp{g}),  quantile(post_by_grp{g}, 0.25), quantile(post_by_grp{g}, 0.75), ...
        median(delta_by_grp{g}));
end

grp_vec = [ones(n_grp(1),1); 2*ones(n_grp(2),1); ...
           3*ones(n_grp(3),1); 4*ones(n_grp(4),1)];

%--------------------------------------------------------------------
% [1] KW on PRE — Baseline comparability
%--------------------------------------------------------------------
fprintf('\n-- [1] Kruskal-Wallis: Baseline (PRE) comparability --\n');
all_pre = vertcat(pre_by_grp{:});
[p_kw_pre, tbl_kw_pre, ~] = kruskalwallis(all_pre, grp_vec, 'off');
fprintf('  H = %.3f, df = %d, p = %.4f  %s\n', ...
    tbl_kw_pre{2,5}, tbl_kw_pre{2,3}, p_kw_pre, stars(p_kw_pre));
if p_kw_pre >= 0.05
    fprintf('  => Groups are comparable at baseline\n');
else
    fprintf('  => Baseline imbalance detected\n');
end

%--------------------------------------------------------------------
% [2] Wilcoxon Signed-Rank — Within-group Pre vs Post
%--------------------------------------------------------------------
fprintf('\n-- [2] Wilcoxon Signed-Rank: Within-group Pre vs Post --\n');
wilcox_p = NaN(4,1);
for g = 1:4
    d    = post_by_grp{g} - pre_by_grp{g};
    d_nz = d(d ~= 0);
    if length(d_nz) < 4
        fprintf('  %s: non-zero differences too few (n_nz=%d) — skipped\n', ...
            group_labels{g}, length(d_nz));
    else
        [p_w, ~] = signrank(pre_by_grp{g}, post_by_grp{g});
        wilcox_p(g) = p_w;
        fprintf('  %s: p = %.4f  %s\n', group_labels{g}, p_w, stars(p_w));
    end
end

%--------------------------------------------------------------------
% [3] KW + Dunn's on CHANGE score — Inter-group intervention effect
%--------------------------------------------------------------------
fprintf('\n-- [3] Kruskal-Wallis: Change score (Post - Pre) --\n');
all_delta = vertcat(delta_by_grp{:});
[p_kw_delta, tbl_kw_delta, ~] = kruskalwallis(all_delta, grp_vec, 'off');
fprintf('  H = %.3f, df = %d, p = %.4f  %s\n', ...
    tbl_kw_delta{2,5}, tbl_kw_delta{2,3}, p_kw_delta, stars(p_kw_delta));

dunn_delta_p = NaN(4,4);
if p_kw_delta < 0.05
    fprintf('  => Significant — Dunn''s post-hoc (Bonferroni):\n');
    dunn_delta_p = dunn_bonferroni(delta_by_grp, group_labels);
else
    fprintf('  => Not significant — Dunn''s post-hoc skipped\n');
end

fprintf('\n============================================================\n');
fprintf('Report saved to: %s\n', txt_path);
fprintf('============================================================\n');
diary off;
fprintf('Saved: slope_diagnosis_stats.txt\n');

%% ================================================================
%% BLOCK B: Figure
%% ================================================================
fig2 = figure('Name','Diagnosis','NumberTitle','off', ...
    'Color','white','Position',[50 50 1100 580]);

ax2 = axes('Parent', fig2);
hold(ax2, 'on');

diag_tick_labels = {'UWS', 'MCS-', 'MCS+', 'EMCS'};

% Horizontal guide lines at each diagnosis level
for lv = 1:4
    yline(ax2, lv, '--', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.8, 'Alpha', 0.9);
end

% Draw data per group
for g = 1:4
    idx = group_col == group_ids(g);
    pre = pre_diag(idx);
    post = post_diag(idx);
    n   = sum(idx);
    col = colors(g,:);
    dc  = line_colors(g,:);
    xp  = x_pre(g);
    xo  = x_post(g);

    % Transition lines (thickness & opacity proportional to patient count)
    pairs      = [pre, post];
    uniq_pairs = unique(pairs, 'rows');
    for u = 1:size(uniq_pairs, 1)
        cnt = sum(pre == uniq_pairs(u,1) & post == uniq_pairs(u,2));
        lw  = 1.0 + cnt * 0.75;
        alp = min(0.90, 0.25 + cnt * 0.12);
        plot(ax2, [xp, xo], [uniq_pairs(u,1), uniq_pairs(u,2)], '-', ...
            'Color', [dc, alp], 'LineWidth', lw);
    end

    % Jittered scatter (consistent lateral jitter per patient)
    rng(42 + g);
    jx = (rand(n,1) - 0.5) * 0.14;
    jy = (rand(n,1) - 0.5) * 0.14;
    scatter(ax2, xp + jx, pre  + jy, 40, col, 'filled', ...
        'MarkerFaceAlpha', 0.88, 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
    scatter(ax2, xo + jx, post + jy, 40, col, 'filled', ...
        'MarkerFaceAlpha', 0.88, 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);

    % Median line
    pm = median(pre);
    om = median(post);
    plot(ax2, [xp, xo], [pm, om], '-o', ...
        'Color', dc, 'LineWidth', 3.2, ...
        'MarkerFaceColor', dc, 'MarkerEdgeColor', 'w', 'MarkerSize', 10);
end

% Vertical group separators
for g = 1:3
    xsep = (x_post(g) + x_pre(g+1)) / 2;
    xline(ax2, xsep, '--', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.0, 'Alpha', 0.8);
end

%% ── Annotation Layer: Within-group Pre–Post comparisons only ───
% ── [2] Wilcoxon within-group brackets (dark grey, Pre <-> Post) ─
wilcox_y   = 4.30;
col_wilcox = [0.28 0.28 0.28];
for g = 1:4
    p = wilcox_p(g);
    if isnan(p), continue; end
    draw_bracket(ax2, x_pre(g), x_post(g), wilcox_y, ...
        stars(p), col_wilcox, 0.10, 12);
end

% Show only within-group significance brackets; between-group change-score
% comparisons remain in the statistical TXT report but are not plotted.
ylim(ax2, [0.6, 4.6]);

%% ── Axes Formatting (bold, no title, no legend) ─────────────────
xlim(ax2, [0.3, 8.2]);
set(ax2, ...
    'XTick',      all_xticks, ...
    'XTickLabel', tick_labels, ...
    'YTick',      [1 2 3 4], ...
    'YTickLabel', diag_tick_labels, ...
    'FontSize',   13, ...
    'FontName',   'Arial', ...
    'FontWeight', 'bold', ...
    'XColor',     [0.15 0.15 0.15], ...
    'YColor',     [0.15 0.15 0.15], ...
    'LineWidth',  1.6, ...
    'Box',        'off', ...
    'TickDir',    'out', ...
    'TickLength', [0.012 0.012]);
xtickangle(ax2, 35);
ylabel(ax2, 'Diagnosis Level', ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

hold(ax2, 'off');

exportgraphics(fig2, fullfile(output_dir, 'slope_diagnosis_stats.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');
fprintf('Saved: slope_diagnosis_stats.png\n');


%% ================================================================
%% LOCAL HELPER FUNCTIONS
%% ================================================================

function s = stars(p)
    if isnan(p) || isempty(p)
        s = '';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'ns';
    end
end

function s = p_fmt(p)
    if p < 0.001
        s = '<0.001***';
    elseif p < 0.01
        s = sprintf('=%.3f**', p);
    elseif p < 0.05
        s = sprintf('=%.3f*',  p);
    else
        s = sprintf('=%.3f (ns)', p);
    end
end

function draw_bracket(ax, x1, x2, y, label_str, col, tick_h, fsize)
    if nargin < 7, tick_h = 0.10; end
    if nargin < 8, fsize  = 12;   end
    lw = 1.4;
    plot(ax, [x1, x1], [y - tick_h, y], '-', 'Color', col, 'LineWidth', lw);
    plot(ax, [x2, x2], [y - tick_h, y], '-', 'Color', col, 'LineWidth', lw);
    plot(ax, [x1, x2], [y,          y], '-', 'Color', col, 'LineWidth', lw);
    text(ax, (x1 + x2)/2, y + 0.06, label_str, ...
        'HorizontalAlignment', 'center', 'FontSize', fsize, ...
        'FontName', 'Arial', 'FontWeight', 'bold', 'Color', col);
end

function p_mat = dunn_bonferroni(data_by_grp, group_labels)
% Dunn's test with Bonferroni correction (manual, with tie correction)
% Returns symmetric 4x4 matrix of Bonferroni-corrected p-values
%
% Formula (Dunn 1964):
%   SE_ij  = sqrt[ tie_factor * N(N+1)/12 * (1/n_i + 1/n_j) ]
%   z_ij   = (Rbar_i - Rbar_j) / SE_ij
%   p_bonf = min(1, p_raw * C),  C = G*(G-1)/2

    G     = length(data_by_grp);
    allv  = vertcat(data_by_grp{:});
    N     = length(allv);
    r_all = tiedrank(allv);

    n_grp = cellfun(@length, data_by_grp);
    cum   = [0; cumsum(n_grp)];
    R_bar = zeros(G, 1);
    for g = 1:G
        R_bar(g) = mean(r_all(cum(g)+1 : cum(g)+n_grp(g)));
    end

    [~, ~, ic] = unique(allv);
    tie_cnt    = accumarray(ic, 1);
    T_val      = sum(tie_cnt.^3 - tie_cnt);
    tie_factor = 1 - T_val / (N^3 - N);

    n_pairs = G*(G-1)/2;
    p_mat   = NaN(G, G);
    pairs   = nchoosek(1:G, 2);

    for pc = 1:size(pairs, 1)
        i  = pairs(pc,1);  j = pairs(pc,2);
        ni = n_grp(i);    nj = n_grp(j);

        SE     = sqrt(tie_factor * N*(N+1)/12 * (1/ni + 1/nj));
        z      = (R_bar(i) - R_bar(j)) / SE;
        p_raw  = 2 * (1 - normcdf(abs(z)));
        p_bonf = min(1, p_raw * n_pairs);

        p_mat(i,j) = p_bonf;
        p_mat(j,i) = p_bonf;

        fprintf('    %s vs %s: z=%+.3f, p_raw=%.4f, p_bonf=%.4f  %s\n', ...
            group_labels{i}, group_labels{j}, z, p_raw, p_bonf, stars(p_bonf));
    end
end