%% Broadband Amplitude 组间分析
% 功能：IQR离群点剔除、两组独立 t 检验（HA vs LA, HTD vs LTD）、美化可视化
% 输入：合并结果v4_filled.xlsx（单Sheet，Group列: HA=1, HTD=2, LA=3, LTD=4）
% 输出：统计结果Excel + 双子图箱线图PNG/FIG

clear; clc; close all;

%% ==================== 路径设置 ====================
inputPath  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
outputPath = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅\波幅4组统计结果+绘图';

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

%% ==================== 参数设置 ====================
IQR_multiplier = 3;
alpha = 0.05;

groupNames = {'High Alpha', 'High Theta/Delta', 'Low Alpha', 'Low Theta/Delta'};
groupAbbr  = {'HA', 'HTD', 'LA', 'LTD'};
nGroups    = 4;

%% ==================== 1. 读取数据 ====================
T = readtable(inputPath, 'Sheet', 1);

groupData_raw  = cell(1, nGroups);
groupNames_raw = cell(1, nGroups);   % 每组对应的 Name 列表
for g = 1:nGroups
    idx = T{:, 1} == g;
    groupData_raw{g}  = T.Broadband_Amplitude_Mean(idx);
    groupNames_raw{g} = T.Name(idx);
end

% 需要高亮标记的受试者名称
highlightNames = {'ZhuYinShui', 'QiaoWenFei'};

fprintf('================================================================================\n');
fprintf('                    Broadband Amplitude 组间分析报告\n');
fprintf('================================================================================\n\n');
fprintf('【数据概览】\n');
totalN = 0;
for g = 1:nGroups
    fprintf('  %s: %d 例\n', groupNames{g}, length(groupData_raw{g}));
    totalN = totalN + length(groupData_raw{g});
end
fprintf('  总计: %d 例\n', totalN);

%% ==================== 2. IQR 离群点检测（全局边界）====================
fprintf('\n【离群点检测】方法: IQR (倍数=%.1f)\n', IQR_multiplier);

allMetricData = vertcat(groupData_raw{:});
Q1 = prctile(allMetricData, 25);
Q3 = prctile(allMetricData, 75);
IQR_val = Q3 - Q1;
lowerBound = Q1 - IQR_multiplier * IQR_val;
upperBound = Q3 + IQR_multiplier * IQR_val;
fprintf('  全局边界: [%.4f, %.4f]\n', lowerBound, upperBound);

outlierCount = zeros(1, nGroups);
groupData_clean  = cell(1, nGroups);
groupNames_clean = cell(1, nGroups);
outlierDetails = {};
detailIdx = 1;

for g = 1:nGroups
    data = groupData_raw{g};
    isOut = (data < lowerBound) | (data > upperBound);
    outlierCount(g) = sum(isOut);
    groupData_clean{g}  = data(~isOut);
    groupNames_clean{g} = groupNames_raw{g}(~isOut);
    if any(isOut)
        outlierIdx = find(isOut);
        for i = 1:length(outlierIdx)
            outlierDetails{detailIdx, 1} = groupNames{g};
            outlierDetails{detailIdx, 2} = 'Broadband Amplitude';
            outlierDetails{detailIdx, 3} = outlierIdx(i);
            outlierDetails{detailIdx, 4} = data(outlierIdx(i));
            outlierDetails{detailIdx, 5} = lowerBound;
            outlierDetails{detailIdx, 6} = upperBound;
            detailIdx = detailIdx + 1;
        end
    end
    fprintf('  %s: 剔除 %d 例，保留 %d 例\n', groupAbbr{g}, outlierCount(g), length(groupData_clean{g}));
end

% 描述统计
groupMean = zeros(1, nGroups);
groupStd  = zeros(1, nGroups);
groupN    = zeros(1, nGroups);
for g = 1:nGroups
    groupMean(g) = mean(groupData_clean{g}, 'omitnan');
    groupStd(g)  = std(groupData_clean{g},  'omitnan');
    groupN(g)    = length(groupData_clean{g});
end

%% ==================== 3. 两组独立检验：HA vs LA, HTD vs LTD ====================
% 比较对: {HA(1) vs LA(3)}, {HTD(2) vs LTD(4)}
comparePairs = {[1, 3], [2, 4]};
pairLabels   = {'HA vs LA', 'HTD vs LTD'};

fprintf('\n【组间比较分析结果】（两独立样本 Wilcoxon 秩和检验）\n\n');

pairResults = {};
for k = 1:2
    g1 = comparePairs{k}(1);
    g2 = comparePairs{k}(2);
    d1 = groupData_clean{g1};
    d2 = groupData_clean{g2};

    [p_val, ~, stats_out] = ranksum(d1, d2);

    % Cohen's d（效应量）
    pooled_std = sqrt(((length(d1)-1)*std(d1)^2 + (length(d2)-1)*std(d2)^2) / ...
                      (length(d1) + length(d2) - 2));
    cohens_d = (mean(d1) - mean(d2)) / pooled_std;

    if p_val < 0.001; sigStr = '***';
    elseif p_val < 0.01; sigStr = '**';
    elseif p_val < 0.05; sigStr = '*';
    else; sigStr = 'n.s.'; end

    fprintf('  %s\n', pairLabels{k});
    fprintf('    %s: Mean=%.4f, SD=%.4f, n=%d\n', groupAbbr{g1}, groupMean(g1), groupStd(g1), groupN(g1));
    fprintf('    %s: Mean=%.4f, SD=%.4f, n=%d\n', groupAbbr{g2}, groupMean(g2), groupStd(g2), groupN(g2));
    fprintf('    Wilcoxon p = %.4f, Cohen''s d = %.3f  %s\n\n', p_val, cohens_d, sigStr);

    pairResults{k, 1}  = groupAbbr{g1};
    pairResults{k, 2}  = groupAbbr{g2};
    pairResults{k, 3}  = groupN(g1);
    pairResults{k, 4}  = groupMean(g1);
    pairResults{k, 5}  = groupStd(g1);
    pairResults{k, 6}  = groupN(g2);
    pairResults{k, 7}  = groupMean(g2);
    pairResults{k, 8}  = groupStd(g2);
    pairResults{k, 9}  = p_val;
    pairResults{k, 10} = cohens_d;
    pairResults{k, 11} = sigStr;
end

%% ==================== 4. 双子图可视化 ====================
% 颜色: HA=#b6d093(绿), HTD=#f5cfa8(橙), LA=#94cee5(蓝), LTD=#e9a2a3(粉)
colors = [182/255, 208/255, 147/255;   % HA
          245/255, 207/255, 168/255;   % HTD
          148/255, 206/255, 229/255;   % LA
          233/255, 162/255, 163/255];  % LTD

% 子图布局: 左 = HA/LA, 右 = HTD/LTD
subplotPairs  = {[1, 3], [2, 4]};   % group indices
subplotAbbr   = {{'HA', 'LA'}, {'HTD', 'LTD'}};
subplotTitles = {'Alpha Groups', 'Theta/Delta Groups'};

figure('Position', [150, 80, 760, 520], 'Color', 'w');

axHandles = gobjects(1, 2);

for k = 1:2
    g1 = subplotPairs{k}(1);
    g2 = subplotPairs{k}(2);

    d1 = groupData_clean{g1};
    d2 = groupData_clean{g2};

    plotData   = [d1; d2];
    plotLabels = [ones(length(d1),1); 2*ones(length(d2),1)];

    axHandles(k) = subplot(1, 2, k);
    ax = axHandles(k);

    bp = boxplot(plotData, plotLabels, ...
        'Colors', [0.3 0.3 0.3], 'Symbol', '', 'Widths', 0.52);
    set(bp, 'LineWidth', 0.9);
    set(findobj(ax, 'Tag', 'Median'), 'LineWidth', 1.0, 'Color', [0.15 0.15 0.15]);

    hold on;

    % 填充箱体颜色（boxplot 绘图顺序为倒序）
    h = findobj(ax, 'Tag', 'Box');
    groupOrder = [g2, g1];   % boxplot 从右到左绘制
    for j = 1:length(h)
        gi = groupOrder(j);
        patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(gi, :), ...
            'FaceAlpha', 0.65, 'EdgeColor', colors(gi,:) * 0.72, 'LineWidth', 1.0);
    end

    % 散点抖动
    localGroups = [g1, g2];
    for li = 1:2
        gi = localGroups(li);
        data_g  = groupData_clean{gi};
        names_g = groupNames_clean{gi};
        rng(gi*100 + k*10 + li);   % 固定随机种子保证抖动可复现
        jitter  = 0.13 * (rand(size(data_g)) - 0.5);

        % 判断哪些点需要高亮
        isHL = ismember(names_g, highlightNames);

        % 普通散点（白色边框）
        if any(~isHL)
            scatter(li + jitter(~isHL), data_g(~isHL), 40, colors(gi,:), 'filled', ...
                'MarkerFaceAlpha', 0.72, 'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
        end

        % 高亮散点（灰色边框标记）
        if any(isHL)
            scatter(li + jitter(isHL), data_g(isHL), 32, colors(gi,:), 'filled', ...
                'MarkerFaceAlpha', 0.72, 'MarkerEdgeColor', [0.6, 0.6, 0.6], 'LineWidth', 1.2);
        end
    end

    % 显著性标注
    p_val  = cell2mat(pairResults(k, 9));
    sigStr = pairResults{k, 11};
    if ~strcmp(sigStr, 'n.s.')
        ylims  = ylim;
        yRange = ylims(2) - ylims(1);
        yBar   = ylims(2) + 0.04 * yRange;
        yTick  = yBar - 0.06 * yRange;
        plot([1,1], [yTick, yBar], 'k-', 'LineWidth', 1.0);
        plot([1,2], [yBar, yBar],  'k-', 'LineWidth', 1.0);
        plot([2,2], [yTick, yBar], 'k-', 'LineWidth', 1.0);
        text(1.5, yBar + 0.025*yRange, sigStr, ...
            'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        ylim([ylims(1), yBar + 0.12*yRange]);
    end

    xlim([0.4, 2.6]);
    set(ax, 'XTick', [1 2], 'XTickLabel', subplotAbbr{k}, ...
        'FontSize', 11, 'FontWeight', 'bold', 'LineWidth', 0.8, ...
        'Box', 'off', 'TickDir', 'out', 'TickLength', [0.018 0.018], ...
        'XColor', [0.18 0.18 0.18], 'YColor', [0.18 0.18 0.18]);

    if k == 1
        ylabel('Broadband Amplitude (µV)', 'FontSize', 12, 'FontWeight', 'bold');
    end

end

% 统一两子图 Y 轴范围
allY = [axHandles(1).YLim, axHandles(2).YLim];
yMin = min(allY); yMax = max(allY);
set(axHandles(1), 'YLim', [yMin, yMax]);
set(axHandles(2), 'YLim', [yMin, yMax]);
set(axHandles(2), 'YTickLabel', []);   % 右图不重复 Y 刻度标签

% 调整子图间距
set(axHandles(1), 'Position', [0.12, 0.13, 0.36, 0.76]);
set(axHandles(2), 'Position', [0.56, 0.13, 0.36, 0.76]);

pngFile = fullfile(outputPath, 'Broadband_Amplitude_Comparison.png');
figFile = fullfile(outputPath, 'Broadband_Amplitude_Comparison.fig');
print(gcf, pngFile, '-dpng', '-r200');
saveas(gcf, figFile);
fprintf('\n图已保存: %s\n', pngFile);

%% ==================== 5. 保存统计结果 Excel ====================
outputFile = fullfile(outputPath, 'Broadband_Amplitude_Stats_Results.xlsx');
if exist(outputFile, 'file'); delete(outputFile); end

% Sheet 1: 描述统计
descTable = table(groupNames', groupN', groupMean', groupStd', ...
    'VariableNames', {'Group', 'N', 'Mean', 'SD'});
writetable(descTable, outputFile, 'Sheet', 'Descriptive_Stats');

% Sheet 2: 组间比较结果
compTable = cell2table(pairResults, ...
    'VariableNames', {'Group1', 'Group2', ...
                      'N_G1', 'Mean_G1', 'SD_G1', ...
                      'N_G2', 'Mean_G2', 'SD_G2', ...
                      'p_value', 'Cohens_d', 'Significance'});
writetable(compTable, outputFile, 'Sheet', 'Pairwise_Comparison');

% Sheet 3: 离群点详情
if ~isempty(outlierDetails)
    outlierTable = cell2table(outlierDetails, ...
        'VariableNames', {'Group', 'Metric', 'RowIndex', 'Value', 'Lower_Bound', 'Upper_Bound'});
    writetable(outlierTable, outputFile, 'Sheet', 'Outlier_Details');
end

fprintf('统计结果已保存: %s\n', outputFile);
fprintf('\n================================================================================\n');
fprintf('                              分析完成!\n');
fprintf('================================================================================\n');