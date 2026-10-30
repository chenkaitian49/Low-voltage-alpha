%% Wilcoxon秩和检验 —— Mean_Freq_Hz
% 分析变量: Mean_Freq_Hz
% 分组变量: Group (HA=1, HTD=2, LA=3, LTD=4)
% 图1左坐标系: HA vs HTD（仅在这两组间做Wilcoxon检验）
% 图1右坐标系: LA vs LTD（仅在这两组间做Wilcoxon检验）

clear; clc; close all;

%% 0. 设置输入和输出路径
% ====== 在此修改输入文件路径 ======
inputFile = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
% ====== 在此修改输出文件夹路径 ======
outputPath = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅\频率4组统计结果+绘图';
% =====================================

% ====== 重点标记患者设置 ======
highlightRows  = [18, 99];
highlightNames = {'SunEnQiang', 'GaoXingLong'};
% ==============================

if ~isfile(inputFile)
    error('输入文件不存在: %s\n请检查 inputFile 路径是否正确。', inputFile);
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

% --- 读取分组变量 ---
groupColName = '';
colNames = data.Properties.VariableNames;
for i = 1:length(colNames)
    if contains(colNames{i}, 'Group', 'IgnoreCase', true)
        groupColName = colNames{i};
        break;
    end
end
if isempty(groupColName)
    error('未找到包含 "Group" 的列名，请检查Excel文件的列名。');
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

% --- 读取分析变量 ---
Mean_Freq_Hz = data.Mean_Freq_Hz;
varNames  = {'Mean_Freq_Hz'};
varLabels = {'Dominant Freq (Hz)'};
varData   = {Mean_Freq_Hz};

fprintf('\n各组样本数:\n');
for g = 1:4
    fprintf('  %s (编码%d): %d 例\n', groupOrder{g}, groupCodes(g), sum(strcmp(group, groupOrder{g})));
end

%% 2. 结果表格初始化
% 每个子图对应一次独立Wilcoxon检验，共2次（HA-HTD 和 LA-LTD）
wilcoxonResults = cell(0, 6);

%% 3. 颜色与样式定义
% 颜色顺序对应 HA / HTD / LA / LTD
colors = [182/255 208/255 147/255;   % HA  - 绿
          245/255 207/255 168/255;   % HTD - 橙
          148/255 206/255 229/255;   % LA  - 蓝
          233/255 162/255 163/255];  % LTD - 红
edgeColors = [colors(:,1)*0.65, colors(:,2)*0.65, colors(:,3)*0.65];

highlightSize      = 32;
highlightEdgeColor = [0.6, 0.6, 0.6];
highlightLineWidth = 1.2;

% 两个子图的分组配置
% panelGroups{p} = {组名列表}，colorIdx{p} = 对应colors行号
panelGroups   = {{'HA','HTD'}, {'LA','LTD'}};
panelColorIdx = {[1,2],        [3,4]};
panelTitles   = {'HA vs HTD',  'LA vs LTD'};

%% 4. 主循环：每个分析变量生成一张双子图
for v = 1:length(varNames)
    varName  = varNames{v};
    varLabel = varLabels{v};
    y = varData{v};
    if ~isnumeric(y)
        y = str2double(y);
    end

    validIdx    = ~isnan(y);
    y_valid     = y(validIdx);
    group_valid = group(validIdx);
    originalRows = find(validIdx);

    % ── 生成一张宽图，包含左右两个子图 ──────────────────────────────────
    fig = figure('Position', [50, 50, 760, 480], 'Color', 'w');

    for p = 1:2   % p=1: HA/HTD，p=2: LA/LTD
        pGroups   = panelGroups{p};      % 本子图的组名，如 {'HA','HTD'}
        pColIdx   = panelColorIdx{p};    % 对应colors的行索引

        % 筛选本子图数据
        panelMask = ismember(group_valid, pGroups);
        y_panel   = y_valid(panelMask);
        g_panel   = group_valid(panelMask);
        origRows_panel = originalRows(panelMask);

        if isempty(y_panel)
            warning('子图%d 变量 %s 没有有效数据', p, varName);
            continue;
        end

        % ── Wilcoxon秩和检验（仅本子图两组之间）──────────────────────────
        data1 = y_panel(strcmp(g_panel, pGroups{1}));
        data2 = y_panel(strcmp(g_panel, pGroups{2}));

        [p_wilcox, ~] = ranksum(data1, data2);
        medDiff = median(data1) - median(data2);

        if p_wilcox < 0.001
            sig_str = '***';
        elseif p_wilcox < 0.01
            sig_str = '**';
        elseif p_wilcox < 0.05
            sig_str = '*';
        else
            sig_str = 'ns';
        end

        % 记录到结果表
        rowIdx = size(wilcoxonResults, 1) + 1;
        wilcoxonResults{rowIdx, 1} = varName;
        wilcoxonResults{rowIdx, 2} = pGroups{1};
        wilcoxonResults{rowIdx, 3} = pGroups{2};
        wilcoxonResults{rowIdx, 4} = round(medDiff, 4);
        wilcoxonResults{rowIdx, 5} = round(p_wilcox, 6);
        wilcoxonResults{rowIdx, 6} = sig_str;

        % ── 子图绘制 ──────────────────────────────────────────────────────
        ax = subplot(1, 2, p, 'Parent', fig);

        % 将组名映射为子图内局部编号 1/2
        groupNum_plot = zeros(size(g_panel));
        for g = 1:2
            groupNum_plot(strcmp(g_panel, pGroups{g})) = g;
        end

        bp = boxplot(y_panel, groupNum_plot, ...
            'Labels', pGroups, ...
            'Colors', [0.3 0.3 0.3], ...
            'Widths', 0.55, ...
            'symbol', '', ...
            'MedianStyle', 'line', ...
            'OutlierSize', 4);
        set(bp, 'LineWidth', 0.8);
        hold on;

        % 箱体填色
        boxes = findobj(ax, 'Tag', 'Box');
        for j = 1:length(boxes)
            localColorIdx = pColIdx(length(boxes) - j + 1);
            patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), ...
                colors(localColorIdx, :), ...
                'FaceAlpha', 0.75, ...
                'EdgeColor', edgeColors(localColorIdx, :), ...
                'LineWidth', 1.0);
        end

        % 散点叠加
        for g = 1:2
            gIdxInPanel = find(strcmp(g_panel, pGroups{g}));
            gData = y_panel(gIdxInPanel);

            if ~isempty(gData)
                rng(pColIdx(g));   % 固定随机种子（用全局颜色索引确保可复现）
                x_jitter = g + 0.18 * (rand(length(gData), 1) - 0.5);

                % 判断重点患者
                isHighlight = false(length(gIdxInPanel), 1);
                for h = 1:length(highlightRows)
                    posInOriginal = find(origRows_panel == highlightRows(h), 1);
                    if ~isempty(posInOriginal)
                        posInGroup = find(gIdxInPanel == posInOriginal, 1);
                        if ~isempty(posInGroup)
                            isHighlight(posInGroup) = true;
                        end
                    end
                end

                % 普通点
                normalIdx = ~isHighlight;
                if any(normalIdx)
                    scatter(x_jitter(normalIdx), gData(normalIdx), 22, ...
                        colors(pColIdx(g), :), 'filled', ...
                        'MarkerFaceAlpha', 0.85, ...
                        'MarkerEdgeColor', 'w', ...
                        'LineWidth', 0.4);
                end

                % 重点患者点
                hlIdx = find(isHighlight);
                for h = 1:length(hlIdx)
                    scatter(x_jitter(hlIdx(h)), gData(hlIdx(h)), highlightSize, ...
                        colors(pColIdx(g), :), 'filled', ...
                        'MarkerFaceAlpha', 1.0, ...
                        'MarkerEdgeColor', highlightEdgeColor, ...
                        'LineWidth', highlightLineWidth);
                end
            end
        end

        % ── 显著性标注横线 ─────────────────────────────────────────────
        yLimits = ylim;
        yRange  = yLimits(2) - yLimits(1);

        yMax = 13;   % Y轴上限固定为13 Hz
        yMin = yLimits(1);

        if p_wilcox < 0.05
            yOffset   = yMax - (yMax - yMin) * 0.06;
            barHeight = (yMax - yMin) * 0.015;
            plot([1, 1, 2, 2], ...
                 [yOffset, yOffset+barHeight, yOffset+barHeight, yOffset], ...
                 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
            text(1.5, yOffset + barHeight*1.8, sig_str, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 11, 'FontWeight', 'bold', ...
                'Color', [0.15 0.15 0.15]);
        end
        ylim([yMin, yMax]);

        % ── 坐标轴样式 ─────────────────────────────────────────────────
        title('', 'Parent', ax);
        if p == 1
            ylabel(varLabel, 'FontSize', 12, 'FontWeight', 'bold', ...
                'Interpreter', 'none', 'Color', [0.2 0.2 0.2], 'Parent', ax);
        else
            ylabel('', 'Parent', ax);
        end
        xlabel('', 'Parent', ax);

        set(ax, 'XTickLabelRotation', 0, ...
            'FontSize', 11, ...
            'FontName', 'Arial', ...
            'FontWeight', 'bold', ...
            'LineWidth', 0.8, ...
            'TickDir', 'out', ...
            'TickLength', [0.015 0.015]);
        ax.YAxis.FontWeight = 'bold';
        ax.XAxis.FontWeight = 'bold';
        grid off;
        box off;
        set(ax, 'Color', 'w');
        hold off;
    end

    set(fig, 'Color', 'w');

    % 保存图片
    outName = fullfile(outputPath, varNames{v});
    print(fig, outName, '-dpng', '-r300');
    saveas(fig, [outName, '.fig']);
    fprintf('图片已保存: %s.png\n', outName);
end

%% 5. 保存结果到Excel
if ~isempty(wilcoxonResults)
    wilcoxonTable = cell2table(wilcoxonResults, 'VariableNames', ...
        {'Variable', 'Group1', 'Group2', 'MedianDiff', 'p_Wilcoxon', 'Significance'});
else
    wilcoxonTable = table();
    warning('Wilcoxon检验结果为空');
end

% 描述统计（中位数 + 四分位距）
descStats = table();
for v = 1:length(varNames)
    y = varData{v};
    if ~isnumeric(y); y = str2double(y); end
    for g = 1:4
        gData = y(strcmp(group, groupOrder{g}));
        gData = gData(~isnan(gData));
        if ~isempty(gData)
            q1 = prctile(gData, 25);
            q3 = prctile(gData, 75);
            newRow = table({varNames{v}}, {groupOrder{g}}, length(gData), ...
                round(median(gData),4), round(q1,4), round(q3,4), round(q3-q1,4), ...
                'VariableNames', {'Variable', 'Group', 'N', 'Median', 'Q1', 'Q3', 'IQR'});
            descStats = [descStats; newRow];
        end
    end
end

excelFile = fullfile(outputPath, 'Wilcoxon_Freq_Results.xlsx');
if ~isempty(wilcoxonTable) && height(wilcoxonTable) > 0
    writetable(wilcoxonTable, excelFile, 'Sheet', 'Wilcoxon');
end
writetable(descStats, excelFile, 'Sheet', 'Descriptive');
fprintf('Excel结果已保存至: %s\n', excelFile);

%% 6. 显示结果
fprintf('\n========== Wilcoxon Results (Mean_Freq_Hz) ==========\n');
if ~isempty(wilcoxonTable) && height(wilcoxonTable) > 0
    disp(wilcoxonTable);
end
fprintf('\nSignificance: * p<0.05, ** p<0.01, *** p<0.001\n');
fprintf('\n========== Output Files ==========\n');
fprintf('1. %s.png\n', fullfile(outputPath, 'Mean_Freq_Hz'));
fprintf('2. %s.fig\n', fullfile(outputPath, 'Mean_Freq_Hz'));
fprintf('3. %s\n', excelFile);
fprintf('\nAnalysis completed!\n');