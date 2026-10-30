%% =========================================================
%  基线资料组间比较 (4组: HA=1, HTD=2, LA=3, LTD=4)
%  连续变量: Age, Duration  → Kruskal-Wallis 检验 + 事后两两比较
%  分类变量: Sex, Etiology  → 卡方检验 (Chi-square)
%  输出: 命令行汇总表
% ==========================================================

clear; clc;

%% ---------- 0. 路径配置 ----------
input_csv  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
output_dir = 'F:\EEG\LowVoltageResearch\OutputResult\V2\量表评分+诊断（剔除无eeg）';

% 若输出文件夹不存在则自动创建
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 将命令行输出同步写入 log 文件
log_file = fullfile(output_dir, 'baseline_comparison_log.txt');
diary(log_file);
diary on;

%% ---------- 1. 读取数据 ----------
filename = input_csv;
T = readtable(filename, 'PreserveVariableNames', true);

% 列名映射
grp_col  = 'Group(HA=1 HTD=2 LA=3 LTD=4)';
sex_col  = 'Sex ( female=0 male=1)';
age_col  = 'Age';
etio_col = 'Etiology (trauma=1 stroke=2 ischemia=3)';
dur_col  = 'Duration';
crs_col  = 'pre CRS-R';

Group    = double(T{:, grp_col});
Sex      = double(T{:, sex_col});
Age      = double(T{:, age_col});
Etiology = double(T{:, etio_col});
Duration = double(T{:, dur_col});
preCRS   = double(T{:, crs_col});

groups     = [1 2 3 4];
groupNames = {'HA (1)', 'HTD (2)', 'LA (3)', 'LTD (4)'};
nGroups    = numel(groups);

%% ---------- 2. 辅助函数 ----------
% 中位数 [Q1, Q3]
medIQR = @(x) [median(x), quantile(x, 0.25), quantile(x, 0.75)];

%% ---------- 3. 打印汇总表头 ----------
fprintf('\n============================================================\n');
fprintf('         基线资料组间比较 — 四组临床数据\n');
fprintf('============================================================\n\n');

%% ========== 4. 连续变量: Age ==========
fprintf('【Age】  中位数 [Q1, Q3]\n');
fprintf('%-12s  N    Median   Q1     Q3\n', 'Group');
fprintf('%s\n', repmat('-',1,45));

age_data = cell(nGroups, 1);
for i = 1:nGroups
    idx = Group == groups(i);
    d   = Age(idx);
    age_data{i} = d;
    m   = medIQR(d);
    fprintf('%-12s  %-4d  %-7.1f  %-6.1f  %-6.1f\n', ...
        groupNames{i}, sum(idx), m(1), m(2), m(3));
end

% Kruskal-Wallis
[p_age, tbl_age, stats_age] = kruskalwallis( ...
    [age_data{1}; age_data{2}; age_data{3}; age_data{4}], ...
    [ones(size(age_data{1})); 2*ones(size(age_data{2})); ...
     3*ones(size(age_data{3})); 4*ones(size(age_data{4}))], 'off');
fprintf('\nKruskal-Wallis: H = %.3f,  df = %d,  p = %.4f\n', ...
    tbl_age{2,5}, nGroups-1, p_age);

if p_age < 0.05
    fprintf('→ 差异显著 (p < 0.05)，进行事后两两比较 (Dunn / multcompare):\n');
    c_age = multcompare(stats_age, 'Display', 'off');
    pairs = {groupNames{1},groupNames{2}; groupNames{1},groupNames{3}; ...
             groupNames{1},groupNames{4}; groupNames{2},groupNames{3}; ...
             groupNames{2},groupNames{4}; groupNames{3},groupNames{4}};
    fprintf('   %-14s vs %-14s  p-adj\n', 'Group A', 'Group B');
    for r = 1:size(c_age,1)
        fprintf('   %-14s vs %-14s  %.4f %s\n', ...
            pairs{r,1}, pairs{r,2}, c_age(r,6), ...
            pStar(c_age(r,6)));
    end
else
    fprintf('→ 差异不显著 (p ≥ 0.05)\n');
end

%% ========== 5. 连续变量: Duration ==========
fprintf('\n------------------------------------------------------------\n');
fprintf('【Duration (months)】  中位数 [Q1, Q3]\n');
fprintf('%-12s  N    Median   Q1     Q3\n', 'Group');
fprintf('%s\n', repmat('-',1,45));

dur_data = cell(nGroups, 1);
for i = 1:nGroups
    idx = Group == groups(i);
    d   = Duration(idx);
    dur_data{i} = d;
    m   = medIQR(d);
    fprintf('%-12s  %-4d  %-7.1f  %-6.1f  %-6.1f\n', ...
        groupNames{i}, sum(idx), m(1), m(2), m(3));
end

[p_dur, tbl_dur, stats_dur] = kruskalwallis( ...
    [dur_data{1}; dur_data{2}; dur_data{3}; dur_data{4}], ...
    [ones(size(dur_data{1})); 2*ones(size(dur_data{2})); ...
     3*ones(size(dur_data{3})); 4*ones(size(dur_data{4}))], 'off');
fprintf('\nKruskal-Wallis: H = %.3f,  df = %d,  p = %.4f\n', ...
    tbl_dur{2,5}, nGroups-1, p_dur);

if p_dur < 0.05
    fprintf('→ 差异显著 (p < 0.05)，进行事后两两比较:\n');
    c_dur = multcompare(stats_dur, 'Display', 'off');
    pairs = {groupNames{1},groupNames{2}; groupNames{1},groupNames{3}; ...
             groupNames{1},groupNames{4}; groupNames{2},groupNames{3}; ...
             groupNames{2},groupNames{4}; groupNames{3},groupNames{4}};
    fprintf('   %-14s vs %-14s  p-adj\n', 'Group A', 'Group B');
    for r = 1:size(c_dur,1)
        fprintf('   %-14s vs %-14s  %.4f %s\n', ...
            pairs{r,1}, pairs{r,2}, c_dur(r,6), ...
            pStar(c_dur(r,6)));
    end
else
    fprintf('→ 差异不显著 (p ≥ 0.05)\n');
end

%% ========== 5b. 连续变量: pre CRS-R (One-way ANOVA) ==========
fprintf('\n------------------------------------------------------------\n');
fprintf('【pre CRS-R】  均值 ± 标准差\n');
fprintf('%-12s  N    Mean     SD\n', 'Group');
fprintf('%s\n', repmat('-',1,40));

crs_data = cell(nGroups, 1);
crs_grp  = cell(nGroups, 1);
for i = 1:nGroups
    idx = Group == groups(i);
    d   = preCRS(idx);
    crs_data{i} = d;
    crs_grp{i}  = i * ones(size(d));
    fprintf('%-12s  %-4d  %-7.2f  %-6.2f\n', ...
        groupNames{i}, sum(idx), mean(d), std(d));
end

crs_all = [crs_data{1}; crs_data{2}; crs_data{3}; crs_data{4}];
grp_all = [crs_grp{1};  crs_grp{2};  crs_grp{3};  crs_grp{4}];

% 方差齐性检验 (Levene, 基于绝对离差的单因素ANOVA等价实现)
grpMeans = cellfun(@mean, crs_data);
absdev   = abs(crs_all - grpMeans(grp_all));
p_levene = anova1(absdev, grp_all, 'off');
fprintf('\nLevene 方差齐性检验: p = %.4f  → %s\n', p_levene, ...
    ternary(p_levene >= 0.05, '方差齐 (满足ANOVA前提)', '方差不齐 (结果需谨慎)'));

% One-way ANOVA
[p_crs, tbl_crs, stats_crs] = anova1(crs_all, grp_all, 'off');
df_b = tbl_crs{2,3};   % 组间自由度
df_w = tbl_crs{3,3};   % 组内自由度
F_crs = tbl_crs{2,5};  % F 统计量
ss_b = tbl_crs{2,2};   % 组间平方和
ss_t = tbl_crs{4,2};   % 总平方和
eta2 = ss_b / ss_t;    % η²
fprintf('\nOne-way ANOVA: F(%d,%d) = %.3f,  p = %.4f,  η² = %.3f\n', ...
    df_b, df_w, F_crs, p_crs, eta2);

if p_crs < 0.05
    fprintf('→ 差异显著 (p < 0.05)，进行事后两两比较 (Tukey HSD):\n');
    c_crs = multcompare(stats_crs, 'Display', 'off');   % 默认 Tukey-Kramer
    pairs = {groupNames{1},groupNames{2}; groupNames{1},groupNames{3}; ...
             groupNames{1},groupNames{4}; groupNames{2},groupNames{3}; ...
             groupNames{2},groupNames{4}; groupNames{3},groupNames{4}};
    fprintf('   %-14s vs %-14s  p-adj\n', 'Group A', 'Group B');
    for r = 1:size(c_crs,1)
        fprintf('   %-14s vs %-14s  %.4f %s\n', ...
            pairs{r,1}, pairs{r,2}, c_crs(r,6), ...
            pStar(c_crs(r,6)));
    end
else
    fprintf('→ 差异不显著 (p ≥ 0.05)\n');
end

%% ========== 6. 分类变量: Sex ==========
fprintf('\n------------------------------------------------------------\n');
fprintf('【Sex】  n (%%)\n');
fprintf('%-12s  N    Female n(%%)    Male n(%%)\n', 'Group');
fprintf('%s\n', repmat('-',1,50));

sex_tbl = zeros(nGroups, 2);   % col1=female(0), col2=male(1)
for i = 1:nGroups
    idx = Group == groups(i);
    sex_tbl(i,1) = sum(Sex(idx) == 0);
    sex_tbl(i,2) = sum(Sex(idx) == 1);
    n = sum(idx);
    fprintf('%-12s  %-4d  %-14s  %s\n', groupNames{i}, n, ...
        sprintf('%d (%.1f%%)', sex_tbl(i,1), 100*sex_tbl(i,1)/n), ...
        sprintf('%d (%.1f%%)', sex_tbl(i,2), 100*sex_tbl(i,2)/n));
end

[~, p_sex, stats_sex] = chi2cont(sex_tbl);
fprintf('\nChi-square: χ² = %.3f,  df = %d,  p = %.4f  %s\n', ...
    stats_sex.chi2stat, stats_sex.df, p_sex, pStar(p_sex));
if p_sex < 0.05
    fprintf('→ 差异显著 (p < 0.05)\n');
else
    fprintf('→ 差异不显著 (p ≥ 0.05)\n');
end

%% ========== 7. 分类变量: Etiology ==========
fprintf('\n------------------------------------------------------------\n');
fprintf('【Etiology】  n (%%)\n');
fprintf('%-12s  N    Trauma n(%%)    Stroke n(%%)    Ischemia n(%%)\n', 'Group');
fprintf('%s\n', repmat('-',1,65));

etio_tbl = zeros(nGroups, 3);  % col1=trauma(1), col2=stroke(2), col3=ischemia(3)
for i = 1:nGroups
    idx = Group == groups(i);
    etio_tbl(i,1) = sum(Etiology(idx) == 1);
    etio_tbl(i,2) = sum(Etiology(idx) == 2);
    etio_tbl(i,3) = sum(Etiology(idx) == 3);
    n = sum(idx);
    fprintf('%-12s  %-4d  %-16s  %-15s  %s\n', groupNames{i}, n, ...
        sprintf('%d (%.1f%%)', etio_tbl(i,1), 100*etio_tbl(i,1)/n), ...
        sprintf('%d (%.1f%%)', etio_tbl(i,2), 100*etio_tbl(i,2)/n), ...
        sprintf('%d (%.1f%%)', etio_tbl(i,3), 100*etio_tbl(i,3)/n));
end

[~, p_etio, stats_etio] = chi2cont(etio_tbl);
fprintf('\nChi-square: χ² = %.3f,  df = %d,  p = %.4f  %s\n', ...
    stats_etio.chi2stat, stats_etio.df, p_etio, pStar(p_etio));
if p_etio < 0.05
    fprintf('→ 差异显著 (p < 0.05)\n');
else
    fprintf('→ 差异不显著 (p ≥ 0.05)\n');
end

%% ---------- 汇总 ----------
fprintf('\n============================================================\n');
fprintf('汇总 (p 值)\n');
fprintf('%-20s  p\n', '变量');
fprintf('%s\n', repmat('-',1,32));
fprintf('%-20s  %.4f  %s\n', 'Age',      p_age,  pStar(p_age));
fprintf('%-20s  %.4f  %s\n', 'Duration', p_dur,  pStar(p_dur));
fprintf('%-20s  %.4f  %s\n', 'pre CRS-R', p_crs, pStar(p_crs));
fprintf('%-20s  %.4f  %s\n', 'Sex',      p_sex,  pStar(p_sex));
fprintf('%-20s  %.4f  %s\n', 'Etiology', p_etio, pStar(p_etio));
fprintf('============================================================\n');
fprintf('显著性标注: * p<0.05  ** p<0.01  *** p<0.001\n\n');
fprintf('结果已保存至: %s\n', log_file);
diary off;

%% =========================================================
%  局部辅助函数
% ==========================================================

function [tbl, p, stats] = chi2cont(observed)
% 对任意 r×c 列联表执行卡方检验
    [r, c] = size(observed);
    N      = sum(observed(:));
    rowSum = sum(observed, 2);
    colSum = sum(observed, 1);
    expected = rowSum * colSum / N;
    chi2stat = sum((observed(:) - expected(:)).^2 ./ expected(:));
    df       = (r-1)*(c-1);
    p        = 1 - chi2cdf(chi2stat, df);
    tbl      = table(chi2stat, df, p);
    stats.chi2stat = chi2stat;
    stats.df       = df;
    stats.p        = p;
end

function s = pStar(p)
    if p < 0.001,     s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'ns';
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end