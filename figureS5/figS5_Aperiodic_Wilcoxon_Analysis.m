%% 独立样本 Wilcoxon 秩和检验（Mann-Whitney U）+ 皮尔逊相关分析
% 分析变量: Mean_Aperiodic_Exponent, Mean_Aperiodic_Offset
% 比较方式1（波幅）: HA+HTD vs LA+LTD（高波幅 vs 低波幅）
% 比较方式2（频率）: HA+LA vs HTD+LTD（Alpha Group (HA+LA) vs Theta/Delta Group (HTD+LTD)）
% 相关分析: 所有患者合并为1组，比较 Exponent 与 Offset 之间的相关性
%
% 说明：本版本将原独立样本 t 检验替换为 Wilcoxon 秩和检验（ranksum，
%       即 Mann-Whitney U 检验）。非参数方法不依赖方差齐性假设，
%       故不再进行 Levene 方差齐性检验。描述统计改用中位数 (Q1, Q3)。

clear; clc; close all;

%% 0. 路径设置
inputFile  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
outputPath = 'F:\EEG\LowVoltageResearch\OutputResult\非周期成分';

if ~isfile(inputFile)
    error('输入文件不存在: %s\n请检查 inputFile 路径。', inputFile);
end
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
    fprintf('已创建输出文件夹: %s\n', outputPath);
end

%% 1. 读取数据
fprintf('正在读取文件: %s\n', inputFile);
try
    opts = detectImportOptions(inputFile);
    opts.VariableNamingRule = 'preserve';
    data = readtable(inputFile, opts);
catch
    data = readtable(inputFile);
end

fprintf('检测到的列名:\n');
disp(data.Properties.VariableNames');

% 自动定位 Group 列
groupColName = '';
colNames = data.Properties.VariableNames;
for i = 1:length(colNames)
    if contains(colNames{i}, 'Group', 'IgnoreCase', true)
        groupColName = colNames{i};
        break;
    end
end
if isempty(groupColName)
    error('未找到包含 "Group" 的列名，请检查Excel文件。');
end
fprintf('使用分组列: "%s"\n', groupColName);

groupNum   = data.(groupColName);
groupOrder = {'HA', 'HTD', 'LA', 'LTD'};
groupCodes = [1, 2, 3, 4];

group = cell(height(data), 1);
for i = 1:height(data)
    idx = find(groupCodes == groupNum(i), 1);
    if ~isempty(idx)
        group{i} = groupOrder{idx};
    else
        group{i} = 'Unknown';
    end
end

% 分析变量
varNames  = {'Mean_Aperiodic_Exponent', 'Mean_Aperiodic_Offset'};
varLabels = {'Mean Aperiodic Exponent',  'Mean Aperiodic Offset'};
varData   = {data.Mean_Aperiodic_Exponent, data.Mean_Aperiodic_Offset};

fprintf('\n各组样本数:\n');
for g = 1:4
    fprintf('  %s: %d 例\n', groupOrder{g}, sum(strcmp(group, groupOrder{g})));
end

%% 2. 比较方式定义
comparisons(1).name        = 'Amplitude';
comparisons(1).group1Tags  = {'HA',  'HTD'};
comparisons(1).group2Tags  = {'LA',  'LTD'};
comparisons(1).group1Name  = 'High Amplitude (HA+HTD)';
comparisons(1).group2Name  = 'Low Amplitude (LA+LTD)';
comparisons(1).xTickLabels = {'High Amplitude', 'Low Amplitude'};

comparisons(2).name        = 'Rhythm';
comparisons(2).group1Tags  = {'HA',  'LA'};
comparisons(2).group2Tags  = {'HTD', 'LTD'};
comparisons(2).group1Name  = 'Alpha Group (HA+LA)';
comparisons(2).group2Name  = 'Theta/Delta Group (HTD+LTD)';
comparisons(2).xTickLabels = {'Alpha Group', 'Theta/Delta Group'};

%% 3. 颜色定义（统一使用指定色系）
% 指定色系: #b6d093, #f5cfa8, #94cee5, #e9a2a3
c1 = [182/255, 208/255, 147/255];  % #b6d093 绿色
c2 = [245/255, 207/255, 168/255];  % #f5cfa8 橙色
c3 = [148/255, 206/255, 229/255];  % #94cee5 蓝色
c4 = [233/255, 162/255, 163/255];  % #e9a2a3 红色

% 比较方式1（波幅）: 蓝色 vs 橙色
colorsAmp  = [c3; c2];
% 比较方式2（频率）: 绿色 vs 红色
colorsFreq = [c1; c4];
colorMap   = {colorsAmp, colorsFreq};

%% 4. 结果容器（Wilcoxon 秩和检验）
ranksumResults = table('Size', [0, 11], ...
    'VariableTypes', {'string','string','double','double','double','double', ...
                      'double','double','double','string','string'}, ...
    'VariableNames', {'Variable','Comparison', ...
                      'N_Group1','N_Group2','Median_Group1','Median_Group2', ...
                      'RankSum_W','z_value','p_value','Significance','Test_Type'});

descStats = table();

%% 5. 主循环（Wilcoxon 秩和检验）
for v = 1:length(varNames)
    varName  = varNames{v};
    varLabel = varLabels{v};
    y = varData{v};
    if ~isnumeric(y); y = str2double(y); end

    for c = 1:length(comparisons)
        comp   = comparisons(c);
        colors = colorMap{c};
        edgeC  = colors * 0.65;

        % 提取两组（去NaN）
        idx1 = ismember(group, comp.group1Tags) & ~isnan(y);
        idx2 = ismember(group, comp.group2Tags) & ~isnan(y);
        y1   = y(idx1);
        y2   = y(idx2);
        n1   = length(y1);
        n2   = length(y2);

        % ----------------------------------------------------------
        % Wilcoxon 秩和检验（Mann-Whitney U test）
        % ----------------------------------------------------------
        [p_val, ~, stats_w] = ranksum(y1, y2);

        % 秩和统计量 W（第一组秩和）与正态近似 z 值
        W_stat = stats_w.ranksum;
        if isfield(stats_w, 'zval')
            z_val = stats_w.zval;          % 大样本时 ranksum 自动给出 z
        else
            z_val = NaN;                   % 小样本（精确法）无 z 值
        end

        if     p_val < 0.001; sig = '***';
        elseif p_val < 0.01;  sig = '**';
        elseif p_val < 0.05;  sig = '*';
        else;                  sig = 'ns';
        end

        fprintf('\n[%s | %s]\n  Wilcoxon rank-sum: W = %.1f, z = %.4f, p = %.4f  %s\n', ...
            varName, comp.name, W_stat, z_val, p_val, sig);

        % 写入结果表
        newRow = {varName, comp.name, n1, n2, ...
                  round(median(y1),4), round(median(y2),4), ...
                  round(W_stat,1), round(z_val,4), p_val, sig, 'Wilcoxon rank-sum'};
        ranksumResults = [ranksumResults; newRow];

        % 描述统计（中位数 + 四分位距）
        for gi = 1:2
            if gi == 1; gData = y1; gName = comp.group1Name;
            else;        gData = y2; gName = comp.group2Name; end
            newDesc = table({varName}, {comp.name}, {gName}, length(gData), ...
                round(median(gData),4), round(prctile(gData,25),4), round(prctile(gData,75),4), ...
                'VariableNames', {'Variable','Comparison','Group','N','Median','Q1','Q3'});
            descStats = [descStats; newDesc];
        end

        % 存入 plotData
        plotData(v, c).y1       = y1;
        plotData(v, c).y2       = y2;
        plotData(v, c).n1       = n1;
        plotData(v, c).n2       = n2;
        plotData(v, c).p_val    = p_val;
        plotData(v, c).sig      = sig;
        plotData(v, c).colors   = colors;
        plotData(v, c).edgeC    = edgeC;
        plotData(v, c).xLabels  = comp.xTickLabels;
        plotData(v, c).varLabel = varLabel;
        plotData(v, c).varName  = varName;
        plotData(v, c).compName = comp.name;
    end
end

%% 5b. 绘图（箱线图）
for v = 1:length(varNames)
    fig = figure('Position', [100, 100, 860, 520], 'Color', 'w');

    for c = 1:length(comparisons)
        pd     = plotData(v, c);
        colors = pd.colors;
        edgeC  = pd.edgeC;
        y1     = pd.y1;
        y2     = pd.y2;
        n1     = pd.n1;
        n2     = pd.n2;
        p_val  = pd.p_val;
        sig    = pd.sig;

        ax = subplot(1, 2, c);

        allY     = [y1; y2];
        groupVec = [ones(n1,1); 2*ones(n2,1)];

        bp = boxplot(allY, groupVec, ...
            'Labels', pd.xLabels, ...
            'Colors', [0.3 0.3 0.3], ...
            'Widths', 0.50, ...
            'symbol', '', ...
            'MedianStyle', 'line', ...
            'OutlierSize', 4);
        set(bp, 'LineWidth', 0.8);
        hold on;

        % 箱体填色
        boxes = findobj(ax, 'Tag', 'Box');
        for j = 1:length(boxes)
            ci = length(boxes) - j + 1;
            if ci >= 1 && ci <= size(colors,1)
                patch(get(boxes(j),'XData'), get(boxes(j),'YData'), ...
                    colors(ci,:), 'FaceAlpha', 0.75, ...
                    'EdgeColor', edgeC(ci,:), 'LineWidth', 1.0);
            end
        end

        % 散点抖动
        for gi = 1:2
            if gi == 1; gData = y1; else; gData = y2; end
            rng(gi + v*10 + c*100);
            xj = gi + 0.18*(rand(length(gData),1) - 0.5);
            scatter(xj, gData, 22, colors(gi,:), 'filled', ...
                'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
        end

        % 显著性标注
        yL  = ylim;
        yR  = yL(2) - yL(1);
        yBr = yL(2) + yR * 0.05;
        bH  = yR * 0.018;

        if p_val < 0.05
            plot([1,1,2,2], [yBr, yBr+bH, yBr+bH, yBr], ...
                'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
            text(1.5, yBr + bH*2.2, sig, ...
                'HorizontalAlignment','center','FontSize',12, ...
                'FontWeight','bold','Color',[0.15 0.15 0.15]);
            ylim([yL(1), yBr + yR*0.15]);
        end

        ylabel(pd.varLabel, 'FontSize',12,'FontWeight','bold', ...
            'Interpreter','none','Color',[0.2 0.2 0.2]);
        xlabel(''); title('');
        set(ax, 'XTickLabelRotation',0,'FontSize',11,'FontName','Arial', ...
            'FontWeight','bold','LineWidth',0.8,'TickDir','out', ...
            'TickLength',[0.015 0.015]);
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontWeight = 'bold';
        grid off; box off;
        set(ax,'Color','w');
        hold off;
    end

    set(fig,'Color','w');
    figName = fullfile(outputPath, varNames{v});
    print(fig, figName, '-dpng', '-r300');
    saveas(fig, [figName, '.fig']);
    fprintf('图片已保存: %s.png\n', figName);
end

%% ================================================================
%% 6. 皮尔逊相关分析 —— 全样本（所有患者合并为1组）
%%    Exponent vs Offset
%% ================================================================
fprintf('\n====== 皮尔逊相关分析（全样本）======\n');

% 提取两个变量，去除任一变量为NaN的行
exp_all = data.Mean_Aperiodic_Exponent;
off_all = data.Mean_Aperiodic_Offset;
if ~isnumeric(exp_all); exp_all = str2double(exp_all); end
if ~isnumeric(off_all); off_all = str2double(off_all); end

validIdx   = ~isnan(exp_all) & ~isnan(off_all);
exp_valid  = exp_all(validIdx);
off_valid  = off_all(validIdx);
n_corr     = sum(validIdx);

% 皮尔逊相关系数及显著性检验
[r_val, p_corr] = corr(exp_valid, off_valid, 'Type', 'Pearson');

% 95% 置信区间（Fisher z 变换）
z_r    = 0.5 * log((1 + r_val) / (1 - r_val));
se_z   = 1 / sqrt(n_corr - 3);
z_low  = z_r - 1.96 * se_z;
z_high = z_r + 1.96 * se_z;
r_low  = (exp(2*z_low)  - 1) / (exp(2*z_low)  + 1);
r_high = (exp(2*z_high) - 1) / (exp(2*z_high) + 1);

if     p_corr < 0.001; sig_corr = '***';
elseif p_corr < 0.01;  sig_corr = '**';
elseif p_corr < 0.05;  sig_corr = '*';
else;                   sig_corr = 'ns';
end

fprintf('  N = %d\n', n_corr);
fprintf('  Pearson r = %.4f\n', r_val);
fprintf('  95%% CI: [%.4f, %.4f]\n', r_low, r_high);
fprintf('  p = %.4f  %s\n', p_corr, sig_corr);

% 线性回归（用于绘图拟合线）
p_fit   = polyfit(exp_valid, off_valid, 1);
x_range = linspace(min(exp_valid), max(exp_valid), 200);
y_fit   = polyval(p_fit, x_range);

% 按原始分组着色（4组分别用4种颜色）
groupColors = [c3; c1; c2; c4];  % HA=蓝, HTD=绿, LA=橙, LTD=红
% groupCodes = [1, 2, 3, 4] 对应 HA, HTD, LA, LTD

%% 6b. 皮尔逊相关散点图
fig_corr = figure('Position', [200, 200, 620, 560], 'Color', 'w');
ax_corr  = axes('Parent', fig_corr);
hold(ax_corr, 'on');

% 分组绘制散点，赋予不同颜色
legendHandles = gobjects(4,1);
for g = 1:4
    gTag    = groupOrder{g};
    gIdx    = strcmp(group(validIdx), gTag);
    gColor  = groupColors(g,:);
    gEdge   = gColor * 0.65;
    hs = scatter(ax_corr, exp_valid(gIdx), off_valid(gIdx), 45, ...
        'filled', ...
        'MarkerFaceColor', gColor, ...
        'MarkerEdgeColor', gEdge, ...
        'MarkerFaceAlpha', 0.85, ...
        'LineWidth', 0.5);
    legendHandles(g) = hs;
end

% 拟合线（深灰色实线）
plot(ax_corr, x_range, y_fit, '-', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 1.8);

% 95% 置信带（基于预测方差）
n_fit  = n_corr;
x_mean = mean(exp_valid);
s_res  = sqrt(sum((off_valid - polyval(p_fit, exp_valid)).^2) / (n_fit - 2));
se_fit = s_res * sqrt(1/n_fit + (x_range - x_mean).^2 / sum((exp_valid - x_mean).^2));
t_crit = tinv(0.975, n_fit - 2);
ci_up  = y_fit + t_crit * se_fit;
ci_lo  = y_fit - t_crit * se_fit;

fill(ax_corr, [x_range, fliplr(x_range)], [ci_up, fliplr(ci_lo)], ...
    [0.6 0.6 0.6], 'FaceAlpha', 0.18, 'EdgeColor', 'none');

% 图例
legend(ax_corr, legendHandles, groupOrder, ...
    'Location', 'southeast', ...
    'FontSize', 10, ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'Box', 'off');

xlabel(ax_corr, 'Exponent', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'FontName', 'Arial', 'Color', [0.2 0.2 0.2]);
ylabel(ax_corr, 'Offset', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'FontName', 'Arial', 'Color', [0.2 0.2 0.2]);

set(ax_corr, 'FontSize', 11, 'FontName', 'Arial', 'FontWeight', 'bold', ...
    'LineWidth', 0.8, 'TickDir', 'out', ...
    'TickLength', [0.015 0.015], 'Color', 'w');
grid(ax_corr, 'off');
box(ax_corr, 'off');
hold(ax_corr, 'off');

set(fig_corr, 'Color', 'w');
figCorrName = fullfile(outputPath, 'Pearson_Exponent_vs_Offset');
print(fig_corr, figCorrName, '-dpng', '-r300');
saveas(fig_corr, [figCorrName, '.fig']);
fprintf('相关分析图已保存: %s.png\n', figCorrName);

%% 7. 保存 Excel（含皮尔逊相关结果 Sheet）
excelFile = fullfile(outputPath, 'Wilcoxon_Aperiodic_Results.xlsx');
writetable(ranksumResults, excelFile, 'Sheet', 'Wilcoxon');
writetable(descStats,      excelFile, 'Sheet', 'Descriptive');

% 皮尔逊相关结果表
% p_value 用科学计数法字符串保存，避免极小值被 round 截为 0
p_corr_str = sprintf('%.4e', p_corr);
corrResult = table({'Mean_Aperiodic_Exponent'}, {'Mean_Aperiodic_Offset'}, ...
    n_corr, round(r_val,4), {p_corr_str}, {sig_corr}, ...
    round(r_low,4), round(r_high,4), ...
    'VariableNames', {'Variable_X','Variable_Y','N','Pearson_r','p_value', ...
                      'Significance','CI_lower_95','CI_upper_95'});
writetable(corrResult, excelFile, 'Sheet', 'Pearson_Correlation');
fprintf('\nExcel 已保存: %s\n', excelFile);

%% 8. 结果摘要
fprintf('\n====== Results Summary ======\n');
disp(ranksumResults);
fprintf('\n--- 皮尔逊相关 ---\n');
disp(corrResult);
fprintf('\nWilcoxon rank-sum test (Mann-Whitney U):\n');
fprintf('  非参数方法，不依赖方差齐性假设；描述统计采用中位数 (Q1, Q3)。\n');
fprintf('Significance: * p<0.05  ** p<0.01  *** p<0.001  ns=n.s.\n');
fprintf('\nAnalysis completed!\n');