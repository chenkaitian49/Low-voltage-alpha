%% Slope Plot: CRS-R Score by Group — with Mixed-Effects Two-Way ANOVA
%  Statistical model : Two-way ANOVA (Group × Time) with mixed-effects / REML
%  Post-hoc          : Simple Effects (within each group, Pre vs Post)
%                      + Bonferroni correction for 4 comparisons (one per group)
%  Effect sizes      : Cohen's dz (paired simple effects)
%                      + approximate partial eta-squared (F-based) for ANOVA terms
%
%  Groups: HA=1, HTD=2, LA=3, LTD=4
%  Colors: #b6d093 (HA), #f5cfa8 (HTD), #94cee5 (LA), #e9a2a3 (LTD)

clear; clc;

%% ===== USER-DEFINED PATHS =====
input_csv  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
output_dir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\量表评分+诊断（剔除无eeg）';

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ===== Load Data =====
data = readtable(input_csv, 'PreserveVariableNames', true);

group_col = double(data{:, 'Group(HA=1 HTD=2 LA=3 LTD=4)'});
pre_crs   = double(data{:, 'pre CRS-R'});
post_crs  = double(data{:, 'post CRS-R'});

n_total   = numel(pre_crs);

%% ===== Build Long-Format Table for Mixed-Effects ANOVA =====
% Each subject contributes 2 rows (Pre and Post).
% Subject ID is used as the random effect to model repeated measures.

subject_id = repmat((1:n_total)', 2, 1);           % Subject IDs (repeated)
score      = [pre_crs; post_crs];                  % Outcome variable
time_var   = [zeros(n_total,1); ones(n_total,1)];  % 0=Pre, 1=Post
group_var  = [group_col; group_col];               % Group labels (1-4)

% Convert to categorical for fitlme
time_cat  = categorical(time_var,  [0 1],     {'Pre','Post'});
group_cat = categorical(group_var, [1 2 3 4], {'HA','HTD','LA','LTD'});

tbl = table(subject_id, score, time_cat, group_cat, ...
    'VariableNames', {'Subject','CRS_R','Time','Group'});

%% ===== Mixed-Effects Two-Way ANOVA (REML) =====
% Fixed effects : Group, Time, Group×Time interaction
% Random effect : Intercept per subject  (1|Subject)
% Estimation    : REML (default in fitlme)

fprintf('\n=== Mixed-Effects Two-Way ANOVA (REML) ===\n');
lme = fitlme(tbl, 'CRS_R ~ Group * Time + (1|Subject)', ...
    'FitMethod', 'REML');
disp(lme);

% ANOVA table (F-tests for fixed effects)
anova_tbl = anova(lme, 'DFMethod', 'Satterthwaite');
disp(anova_tbl);

%% ===== Approximate Partial Eta-Squared for Fixed Effects =====
% A mixed-effects model fitted by REML has no classical sum-of-squares
% partition, so a conventional partial eta-squared cannot be computed
% directly. We use the standard F-based approximation:
%       eta_p^2 ≈ (F * df1) / (F * df1 + df2)
% NOTE: df2 are Satterthwaite approximations under the mixed model, so
%       these values should be reported as estimates only.

F_vals   = anova_tbl.FStat;
df1_vals = anova_tbl.DF1;
df2_vals = anova_tbl.DF2;
eta_p2   = (F_vals .* df1_vals) ./ (F_vals .* df1_vals + df2_vals);

% Term names: anova(lme) returns them as a categorical/cell column;
% convert to cellstr for robust indexing across MATLAB versions.
term_names = cellstr(string(anova_tbl.Term));
n_terms    = numel(eta_p2);

fprintf('\n=== Approximate Partial Eta-Squared (F-based estimate) ===\n');
for i = 1:n_terms
    fprintf('%-14s  eta_p^2 = %.4f\n', term_names{i}, eta_p2(i));
end

%% ===== Simple Effects: Pre vs Post within each group =====
% Strategy: fit the same LME model; then use coefTest or manually
% extract paired t-statistics with Bonferroni correction.
%
% We use a direct paired t-test per group (equivalent to the simple
% effect under the mixed model) and apply Bonferroni × 4 correction.
% Effect size: Cohen's dz = t / sqrt(n) = mean(diff) / std(diff).

group_labels = {'HA', 'HTD', 'LA', 'LTD'};
group_ids    = [1, 2, 3, 4];
n_groups     = 4;
alpha_bonf   = 0.05 / n_groups;   % Bonferroni-corrected α = 0.0125

fprintf('\n=== Simple Effects: Pre vs Post by Group (Bonferroni corrected) ===\n');
fprintf('Bonferroni-corrected α = %.4f  (0.05 / %d)\n\n', alpha_bonf, n_groups);
fprintf('%-6s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %6s\n', ...
    'Group','n','Pre_M','Post_M','t','dz','p_raw','p_adj','','Sig');
fprintf('%s\n', repmat('-',1,90));

p_raw    = zeros(n_groups, 1);
t_stat   = zeros(n_groups, 1);
df_stat  = zeros(n_groups, 1);
dz_stat  = zeros(n_groups, 1);    % Cohen's dz
sig_bonf = false(n_groups, 1);

for g = 1:n_groups
    idx  = group_col == group_ids(g);
    pre  = pre_crs(idx);
    post = post_crs(idx);
    n    = sum(idx);

    [~, p, ~, stats] = ttest(pre, post);   % paired t-test
    t_stat(g)  = stats.tstat;
    df_stat(g) = stats.df;
    p_raw(g)   = p;
    dz_stat(g) = stats.tstat / sqrt(n);    % Cohen's dz
    p_adj      = min(p * n_groups, 1.0);   % Bonferroni adjusted p
    sig_bonf(g) = p_adj < 0.05;

    fprintf('%-6s  %5d  %8.2f  %8.2f  %8.3f  %8.3f  %8.4f  %8.4f  %8s  %s\n', ...
        group_labels{g}, n, ...
        mean(pre), mean(post), ...
        t_stat(g), dz_stat(g), p_raw(g), p_adj, '', ...
        ternary(sig_bonf(g), '*', 'ns'));
end

%% ===== Color & Layout Definitions =====
colors = [
    182/255, 208/255, 147/255;   % HA   #b6d093
    245/255, 207/255, 168/255;   % HTD  #f5cfa8
    148/255, 206/255, 229/255;   % LA   #94cee5
    233/255, 162/255, 163/255;   % LTD  #e9a2a3
];
line_colors  = colors * 0.60;

x_pre  = [1, 3, 5, 7];
x_post = [2, 4, 6, 8];

all_xticks  = reshape([x_pre; x_post], 1, []);
tick_labels = {};
for g = 1:4
    tick_labels{end+1} = sprintf('Pre %s',  group_labels{g});
    tick_labels{end+1} = sprintf('Post %s', group_labels{g});
end

%% ================================================================
%% FIGURE: CRS-R Slope Plot with Significance Brackets
%% ================================================================
fig1 = figure('Name','CRS-R with Stats','NumberTitle','off', ...
    'Color','white','Position',[50 50 1100 600]);

ax1 = axes('Parent', fig1);
hold(ax1, 'on');

y_min = max(0, min([pre_crs; post_crs]) - 2);
y_max = max([pre_crs; post_crs]) + 2;

%% --- Draw data (same as original) ---
for g = 1:4
    idx = group_col == group_ids(g);
    pre = pre_crs(idx);
    post = post_crs(idx);
    n   = sum(idx);
    col = colors(g,:);
    dc  = line_colors(g,:);
    xp  = x_pre(g);
    xo  = x_post(g);

    rng(42 + g);
    jx = (rand(n,1) - 0.5) * 0.18;

    % Connecting lines
    for i = 1:n
        plot(ax1, [xp + jx(i), xo + jx(i)], [pre(i), post(i)], '-', ...
            'Color', [dc, 0.32], 'LineWidth', 0.9);
    end

    % Scatter dots
    scatter(ax1, xp + jx, pre,  40, col, 'filled', ...
        'MarkerFaceAlpha', 0.88, 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);
    scatter(ax1, xo + jx, post, 40, col, 'filled', ...
        'MarkerFaceAlpha', 0.88, 'MarkerEdgeColor', 'w', 'LineWidth', 0.6);

    % Mean ± SD lines
    pm = mean(pre);  ps = std(pre);
    om = mean(post); os = std(post);
    errorbar(ax1, xp, pm, ps, 'Color', dc, 'LineWidth', 1.8, ...
        'CapSize', 7, 'LineStyle', 'none');
    errorbar(ax1, xo, om, os, 'Color', dc, 'LineWidth', 1.8, ...
        'CapSize', 7, 'LineStyle', 'none');
    plot(ax1, [xp, xo], [pm, om], '-o', ...
        'Color', dc, 'LineWidth', 3.2, ...
        'MarkerFaceColor', dc, 'MarkerEdgeColor', 'w', 'MarkerSize', 10);
end

%% --- Vertical separators ---
for g = 1:3
    xsep = (x_post(g) + x_pre(g+1)) / 2;
    xline(ax1, xsep, '--', 'Color', [0.78 0.78 0.78], 'LineWidth', 1.0, 'Alpha', 0.8);
end

%% --- Significance Brackets ---
% 所有组括号统一高度：取四组最高原始数据点的全局最大值 + gap
bracket_tickh = 0.4;
bracket_gap   = 0.8;   % 最高数据点到括号底端的间距

% 计算全局最高点，所有组括号对齐到同一高度
global_top = max([pre_crs; post_crs]);
y_brk_all  = global_top + bracket_gap;

for g = 1:4
    xp = x_pre(g);
    xo = x_post(g);

    p_adj_g = min(p_raw(g) * n_groups, 1.0);
    if     p_adj_g < 0.001;  sig_label = '***';
    elseif p_adj_g < 0.01;   sig_label = '**';
    elseif p_adj_g < 0.05;   sig_label = '*';
    else;                     sig_label = 'ns';
    end

    if ~strcmp(sig_label, 'ns')
        plot(ax1, [xp, xp], [y_brk_all - bracket_tickh, y_brk_all], '-', ...
            'Color', [0.25 0.25 0.25], 'LineWidth', 1.4);
        plot(ax1, [xp, xo], [y_brk_all, y_brk_all], '-', ...
            'Color', [0.25 0.25 0.25], 'LineWidth', 1.4);
        plot(ax1, [xo, xo], [y_brk_all, y_brk_all - bracket_tickh], '-', ...
            'Color', [0.25 0.25 0.25], 'LineWidth', 1.4);

        text(ax1, (xp+xo)/2, y_brk_all + 0.15, sig_label, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'Color', [0.15 0.15 0.15], 'Clipping', 'off');
    end
end

%% --- Axes formatting ---
xlim(ax1, [0.3, 8.2]);
set(ax1, ...
    'XTick',      all_xticks, ...
    'XTickLabel', tick_labels, ...
    'FontSize',   13, ...
    'FontName',   'Arial', ...
    'FontWeight', 'bold', ...
    'XColor',     [0.15 0.15 0.15], ...
    'YColor',     [0.15 0.15 0.15], ...
    'LineWidth',  1.6, ...
    'Box',        'off', ...
    'TickDir',    'out', ...
    'TickLength', [0.012 0.012]);

% Y 轴下限固定，上限自动（为显著性标记留空间）
ylim(ax1, [y_min, Inf]);

% 刻度只标到 23（CRS-R 满分），23 以上不显示刻度
y_ticks_base = y_min : 5 : 20;
y_ticks_show = unique([y_ticks_base, 23]);
ax1.YTick      = y_ticks_show;
ax1.YTickLabel = arrayfun(@num2str, y_ticks_show, 'UniformOutput', false);
xtickangle(ax1, 35);
ylabel(ax1, 'CRS-R Total Score', ...
    'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');
% title(ax1, 'CRS-R Score: Pre vs. Post by Group', ...
%     'FontSize', 15, 'FontName', 'Arial', 'FontWeight', 'bold', ...
%     'Color', [0.15 0.15 0.15]);

legend(ax1, 'off');

hold(ax1, 'off');

%% --- Save figure ---
out_fig = fullfile(output_dir, 'slope_crsr_stats.png');
exportgraphics(fig1, out_fig, 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('\nSaved figure: %s\n', out_fig);

%% ================================================================
%% TXT: Full Statistics Report
%% ================================================================
txt_path = fullfile(output_dir, 'crsr_stats_full.txt');
fid = fopen(txt_path, 'w');

fprintf(fid, 'CRS-R Total Score — Full Statistical Report\n');
fprintf(fid, '============================================\n\n');

%% Descriptive stats
fprintf(fid, '1. Descriptive Statistics\n');
fprintf(fid, '%-6s  %5s  %-18s  %-18s  %-8s\n', ...
    'Group','n','Pre (Mean±SD)','Post (Mean±SD)','Change');
fprintf(fid, '%s\n', repmat('-',1,62));
for g = 1:4
    idx  = group_col == group_ids(g);
    pre  = pre_crs(idx);
    post = post_crs(idx);
    n    = sum(idx);
    fprintf(fid, '%-6s  %5d  %-18s  %-18s  %+.2f\n', ...
        group_labels{g}, n, ...
        sprintf('%.2f±%.2f', mean(pre), std(pre)), ...
        sprintf('%.2f±%.2f', mean(post), std(post)), ...
        mean(post)-mean(pre));
end

%% ANOVA table
fprintf(fid, '\n\n2. Two-Way ANOVA (Mixed-Effects, REML)\n');
fprintf(fid, '   Model: CRS_R ~ Group * Time + (1|Subject)\n\n');
fprintf(fid, '%s\n', evalc('disp(anova_tbl)'));

%% Approximate partial eta-squared
fprintf(fid, '   Approximate partial eta-squared (F-based estimate):\n');
fprintf(fid, '   eta_p^2 = (F*df1)/(F*df1 + df2); df2 via Satterthwaite, estimate only.\n');
for i = 1:n_terms
    fprintf(fid, '   %-14s  eta_p^2 = %.4f\n', term_names{i}, eta_p2(i));
end

%% Simple effects
fprintf(fid, '\n3. Simple Effects: Pre vs Post within each Group\n');
fprintf(fid, '   (Bonferroni correction for %d comparisons; adjusted α = %.4f)\n\n', ...
    n_groups, alpha_bonf);
fprintf(fid, '%-6s  %5s  %8s  %8s  %8s  %8s  %8s  %6s\n', ...
    'Group','n','t','df','dz','p_raw','p_adj','Sig');
fprintf(fid, '%s\n', repmat('-',1,76));
for g = 1:4
    idx  = group_col == group_ids(g);
    n    = sum(idx);
    p_adj_g = min(p_raw(g) * n_groups, 1.0);
    if     p_adj_g < 0.001;  sl = '***';
    elseif p_adj_g < 0.01;   sl = '**';
    elseif p_adj_g < 0.05;   sl = '*';
    else;                     sl = 'ns';
    end
    fprintf(fid, '%-6s  %5d  %8.3f  %8.1f  %8.3f  %8.4f  %8.4f  %4s\n', ...
        group_labels{g}, n, t_stat(g), df_stat(g), dz_stat(g), ...
        p_raw(g), p_adj_g, sl);
end

fprintf(fid, '\nNote: Paired t-test used for simple effects.\n');
fprintf(fid, '      Bonferroni-adjusted p = min(p_raw × %d, 1.0)\n', n_groups);
fprintf(fid, '      Effect size: Cohen''s dz = t / sqrt(n) for paired comparisons.\n');
fprintf(fid, '      Significance: * p<0.05  ** p<0.01  *** p<0.001  ns: not significant\n');

fclose(fid);
fprintf('Saved statistics: %s\n', txt_path);

%% ===== Helper function =====
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end