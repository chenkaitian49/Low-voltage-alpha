%% Post-Outcome Comparison: Overall Pearson Chi-Square + Pairwise Fisher + Bonferroni + Stacked Bar
%  Outcome: post Diagnosis categorised as
%           VS/MCS-   (post Diagnosis = 1 or 2)
%           MCS+/EMCS (post Diagnosis = 3 or 4)
%
%  Groups : HA=1, HTD=2, LA=3, LTD=4
%  Colors : #b6d093 (HA), #f5cfa8 (HTD), #94cee5 (LA), #e9a2a3 (LTD)
%
%  Statistics:
%   [1] Overall 4-group Pearson chi-square (2×4 contingency table)
%   [2] Pairwise two-sided Fisher's exact tests × 6 + Bonferroni correction
%   [3] Stacked bar figure matching the reference style

clear; clc;

%% ===== USER-DEFINED PATHS =====
input_csv  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
output_dir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\量表评分+诊断（剔除无eeg）';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ===== Load Data =====
data = readtable(input_csv, 'PreserveVariableNames', true);

group_col = double(data{:, 'Group(HA=1 HTD=2 LA=3 LTD=4)'});
post_diag = double(data{:, 'post Diagnosis ( VS=1 MCS-=2 MCS+=3 EMCS=4)'});

group_labels = {'HA', 'HTD', 'LA', 'LTD'};   % 直接对应 Excel: HA=1 HTD=2 LA=3 LTD=4
group_short  = {'HA', 'HTD', 'LA', 'LTD'};
group_ids    = [1, 2, 3, 4];
n_groups     = 4;

%% ===== Data validation =====
% Only patients with valid group coding (1-4) and valid follow-up diagnosis (1-4)
% are included. This prevents missing diagnosis values from being incorrectly
% classified as unfavorable outcomes.
valid_idx = ismember(group_col, group_ids) & ismember(post_diag, 1:4) & ...
            ~isnan(group_col) & ~isnan(post_diag);

if any(~valid_idx)
    fprintf('\nWarning: %d rows excluded because Group or post Diagnosis was missing/invalid.\n', ...
        sum(~valid_idx));
end

group_col = group_col(valid_idx);
post_diag = post_diag(valid_idx);

%% ===== Binary Outcome: 0 = VS/MCS-  |  1 = MCS+/EMCS =====
outcome_bin = double(post_diag >= 3);   % 1 if MCS+ or EMCS

%% ===== Count per group =====
n_total   = zeros(n_groups, 1);
n_good    = zeros(n_groups, 1);   % MCS+/EMCS
n_poor    = zeros(n_groups, 1);   % VS/MCS-

for g = 1:n_groups
    idx         = group_col == group_ids(g);
    n_total(g)  = sum(idx);
    n_good(g)   = sum(outcome_bin(idx) == 1);
    n_poor(g)   = sum(outcome_bin(idx) == 0);
end

fprintf('\n=== Outcome Counts ===\n');
fprintf('%-18s  %5s  %10s  %10s\n', 'Group','n','VS/MCS-','MCS+/EMCS');
fprintf('%s\n', repmat('-',1,48));
for g = 1:n_groups
    fprintf('%-18s  %5d  %10d  %10d\n', ...
        group_labels{g}, n_total(g), n_poor(g), n_good(g));
end

%% ===== [1] Overall Chi-Square (2 × 4) =====
contingency = [n_poor'; n_good'];   % 2 × 4 matrix
[chi2_overall, p_overall] = chi2test_2xK(contingency);

fprintf('\n=== [1] Overall Chi-Square (2×4) ===\n');
fprintf('χ² = %.3f, df = %d, p = %.4f  %s\n', ...
    chi2_overall, (2-1)*(n_groups-1), p_overall, sig_star(p_overall));

%% ===== [2] Pairwise Fisher's Exact Test + Bonferroni =====
pairs     = nchoosek(1:n_groups, 2);
n_pairs   = size(pairs, 1);          % 6 pairwise comparisons
alpha_bon = 0.05 / n_pairs;          % equivalent raw-p threshold = 0.0083

fprintf('\n=== [2] Pairwise two-sided Fisher exact tests + Bonferroni ===\n');
fprintf('Bonferroni: 6 comparisons; adjusted significance level = %.4f\n', alpha_bon);
fprintf('%-10s  %-10s  %12s  %12s  %6s\n', ...
    'GroupA','GroupB','p_raw','p_adj','Sig');
fprintf('%s\n', repmat('-',1,58));

p_raw_mat = NaN(n_groups, n_groups);
p_adj_mat = NaN(n_groups, n_groups);

for pc = 1:n_pairs
    ga = pairs(pc,1);
    gb = pairs(pc,2);

    % 2×2 contingency table: rows = outcome, columns = phenotype groups
    ct = [n_poor(ga), n_poor(gb); ...
          n_good(ga), n_good(gb)];

    % Two-sided Fisher's exact test is preferred for post-hoc 2×2 tables
    % because some cells/expected counts are small (e.g. LTD favorable = 0).
    [~, p_p] = fishertest(ct, 'Tail', 'both');

    % Bonferroni correction across all six pairwise comparisons
    p_adj_p = min(p_p * n_pairs, 1.0);

    p_raw_mat(ga,gb) = p_p;
    p_raw_mat(gb,ga) = p_p;
    p_adj_mat(ga,gb) = p_adj_p;
    p_adj_mat(gb,ga) = p_adj_p;

    fprintf('%-10s  %-10s  %12.6g  %12.6g  %6s\n', ...
        group_short{ga}, group_short{gb}, p_p, p_adj_p, ...
        sig_star(p_adj_p));
end

%% ===== [3] Figure: Stacked Bar + Significance Brackets =====

% Color definitions
% VS/MCS-   → each group's light color (same as group color)
% MCS+/EMCS → each group's dark (salmon-like) version
% Reference colors per group (matching the codes):
% HA  #b6d093  |  HTD #f5cfa8  |  LA #94cee5  |  LTD #e9a2a3
% For the stacked bars we use a unified 2-color scheme
% (blue = VS/MCS-, salmon = MCS+/EMCS) matching the reference figure:
col_poor = [148/255, 188/255, 216/255];   % blue  (VS/MCS-)
col_good = [233/255, 148/255, 150/255];   % salmon (MCS+/EMCS)

fig = figure('Name','Post Outcome by Group','NumberTitle','off', ...
    'Color','white','Position',[80 80 780 560]);

ax = axes('Parent', fig);
hold(ax, 'on');

x_pos = 1:n_groups;
bar_w = 0.52;

% Draw stacked bars using matrix input (standard MATLAB idiom)
% Each row = one group; columns = [VS/MCS-,  MCS+/EMCS]
bar_data = [n_poor, n_good];           % n_groups × 2
hb = bar(ax, x_pos, bar_data, bar_w, 'stacked', ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.2);
hb(1).FaceColor = col_poor;            % VS/MCS-   (bottom, blue)
hb(2).FaceColor = col_good;            % MCS+/EMCS (top, salmon)

%% --- Percentage labels above each bar ---
for g = 1:n_groups
    pct = 100 * n_good(g) / n_total(g);
    text(ax, x_pos(g), n_total(g) + 0.6, ...
        sprintf('%.1f%%', pct), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontSize',   12, ...
        'FontName',   'Arial', ...
        'FontWeight', 'bold', ...
        'Color',      [0.20 0.20 0.20]);
end

%% --- Significance Brackets ---
% Only draw brackets for pairs with p_adj < 0.05
% Stagger heights to avoid overlap

y_top   = max(n_total) + 5;   % +5 为百分比标签 + 括号留出足够空间
brk_gap = 3.2;                 % 括号层间距
tick_h  = 0.7;                 % 括号竖线高度

% Determine bracket levels (sort by span to reduce crossing)
sig_pairs = zeros(0, 2);   % pre-allocate with correct column count
for pc = 1:n_pairs
    ga = pairs(pc,1);  gb = pairs(pc,2);
    if p_adj_mat(ga,gb) < 0.05
        sig_pairs(end+1,:) = [ga, gb]; %#ok<AGROW>
    end
end

% Sort by span (shortest first → lower brackets)
if ~isempty(sig_pairs)
    spans = sig_pairs(:,2) - sig_pairs(:,1);
    [~, si] = sort(spans);
    sig_pairs = sig_pairs(si, :);

    for lv = 1:size(sig_pairs, 1)
        ga  = sig_pairs(lv,1);
        gb  = sig_pairs(lv,2);
        by  = y_top + (lv - 1) * brk_gap;
        slabel = sig_star(p_adj_mat(ga,gb));

        % Bracket lines
        plot(ax, [x_pos(ga), x_pos(ga)], [by - tick_h, by], '-', ...
            'Color', [0.20 0.20 0.20], 'LineWidth', 1.5);
        plot(ax, [x_pos(gb), x_pos(gb)], [by - tick_h, by], '-', ...
            'Color', [0.20 0.20 0.20], 'LineWidth', 1.5);
        plot(ax, [x_pos(ga), x_pos(gb)], [by, by], '-', ...
            'Color', [0.20 0.20 0.20], 'LineWidth', 1.5);
        text(ax, (x_pos(ga)+x_pos(gb))/2, by + 0.25, slabel, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 15, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'Color', [0.15 0.15 0.15]);
    end
end

%% --- Legend (use bar handles directly) ---
lgd = legend(ax, [hb(1), hb(2)], {'UWS/MCS-', 'MCS+/EMCS'}, ...
    'Location', 'northeast', 'FontSize', 12, 'FontName', 'Arial', ...
    'Box', 'off', 'FontWeight', 'bold');
% 图例上移一点，避免与上方显著性标记线过近
lgd.Units = 'normalized';
lgd.Position(2) = min(lgd.Position(2) + 0.03, 1 - lgd.Position(4));

%% --- Axes formatting ---
xlim(ax, [0.4, n_groups + 0.6]);

% ylim下限设为-0.1，防止柱子底部边框与X轴重叠产生细线
y_bottom = -0.1;
if ~isempty(sig_pairs)
    ylim(ax, [y_bottom, y_top + size(sig_pairs,1) * brk_gap + 2]);
else
    ylim(ax, [y_bottom, y_top + 3]);
end

% Y轴刻度从0开始，不显示负数刻度
y_tick_max = ceil(max(n_total) / 5) * 5;
ax.YTick = 0 : 5 : y_tick_max;

set(ax, ...
    'XTick',       x_pos, ...
    'XTickLabel',  group_labels, ...
    'FontSize',    13, ...
    'FontName',    'Arial', ...
    'FontWeight',  'bold', ...
    'XColor',      [0.15 0.15 0.15], ...
    'YColor',      [0.15 0.15 0.15], ...
    'LineWidth',   1.6, ...
    'Box',         'off', ...
    'TickDir',     'out', ...
    'TickLength',  [0.012 0.012]);

% xtickangle(ax, 30);  % 短标签无需旋转
ylabel(ax, 'Number of patients', ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

hold(ax, 'off');

%% --- Save figure ---
out_fig = fullfile(output_dir, 'outcome_fisher_bonferroni_stackedbar.png');
exportgraphics(fig, out_fig, 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('\nSaved figure: %s\n', out_fig);

%% ===== TXT Report =====
txt_path = fullfile(output_dir, 'outcome_fisher_bonferroni_stats.txt');
fid = fopen(txt_path, 'w');

fprintf(fid, 'Follow-up Outcome — Statistical Report\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Outcome categories:\n');
fprintf(fid, '  VS/MCS-   : post Diagnosis = 1 (VS) or 2 (MCS-)\n');
fprintf(fid, '  MCS+/EMCS : post Diagnosis = 3 (MCS+) or 4 (EMCS)\n\n');

fprintf(fid, '1. Counts\n');
fprintf(fid, '%-18s  %5s  %10s  %10s  %10s\n', ...
    'Group','n','VS/MCS-','MCS+/EMCS','%MCS+/EMCS');
fprintf(fid, '%s\n', repmat('-',1,58));
for g = 1:n_groups
    fprintf(fid, '%-18s  %5d  %10d  %10d  %9.1f%%\n', ...
        group_labels{g}, n_total(g), n_poor(g), n_good(g), ...
        100*n_good(g)/n_total(g));
end

fprintf(fid, '\n2. Overall Chi-Square (2×4)\n');
fprintf(fid, '   χ² = %.3f, df = %d, p = %.4f  %s\n\n', ...
    chi2_overall, (2-1)*(n_groups-1), p_overall, sig_star(p_overall));

fprintf(fid, '3. Pairwise two-sided Fisher exact tests + Bonferroni\n');
fprintf(fid, '   Number of pairwise comparisons = %d; Bonferroni raw-p threshold = %.4f\n', ...
    n_pairs, alpha_bon);
fprintf(fid, '%-10s  %-10s  %12s  %12s  %6s\n', ...
    'GroupA','GroupB','p_raw','p_adj','Sig');
fprintf(fid, '%s\n', repmat('-',1,58));
for pc = 1:n_pairs
    ga = pairs(pc,1);  gb = pairs(pc,2);
    fprintf(fid, '%-10s  %-10s  %12.6g  %12.6g  %6s\n', ...
        group_short{ga}, group_short{gb}, ...
        p_raw_mat(ga,gb), p_adj_mat(ga,gb), ...
        sig_star(p_adj_mat(ga,gb)));
end
fprintf(fid, '\nNote: pairwise tests are two-sided Fisher exact tests.\n');
fprintf(fid, '      Bonferroni p_adj = min(p_raw × %d, 1.0).\n', n_pairs);
fprintf(fid, '      Significance is determined from the adjusted p value: * p<0.05, ** p<0.01, *** p<0.001.\n');

fclose(fid);
fprintf('Saved statistics: %s\n', txt_path);

%% ================================================================
%% LOCAL HELPER FUNCTIONS
%% ================================================================

function [chi2, p] = chi2test_2xK(ct)
% Chi-square test for any 2×K contingency table
    col_sum  = sum(ct, 1);
    row_sum  = sum(ct, 2);
    N        = sum(col_sum);
    expected = (row_sum * col_sum) / N;   % outer product → 2×K matrix
    % Guard against zero expected cells (structural zeros → skip)
    mask = expected > 0;
    chi2 = sum(sum((ct(mask) - expected(mask)).^2 ./ expected(mask)));
    df   = (size(ct,1)-1) * (size(ct,2)-1);
    p    = 1 - chi2cdf(chi2, df);
end

function s = sig_star(p)
    if isnan(p) || p >= 0.05
        s = 'ns';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    else
        s = '*';
    end
end
