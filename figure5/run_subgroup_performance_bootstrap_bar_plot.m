function run_subgroup_performance_bootstrap_bar_plot
%RUN_SUBGROUP_PERFORMANCE_BOOTSTRAP_BAR_PLOT
% 启动时显示实际执行文件，便于排除同名旧文件或 MATLAB 路径缓存。
fprintf('Running MATLAB file: %s%s', mfilename('fullpath'), newline);
% 使用冻结的联合模型，在天坛开发集与丰台独立测试集中，
% 按初始诊断、病因、年龄和病程进行亚组分析，并绘制敏感度、特异度、
% 准确率三联点线图（点为 bootstrap 均值，水平线为 bootstrap SD）。
%
% 重要说明：
% 1) 模型参数和分类阈值均直接读取冻结的 .mat 文件，不重新拟合模型，
%    也不在测试集中重新选择阈值。
% 2) 年龄分界值默认取“开发集纳入分析患者”的年龄中位数，并将同一分界值
%    应用于测试集，以保证两个队列的亚组定义一致。
% 3) 病程分界值取开发集病程中位数并四舍五入为整数，同一分界值应用于测试集。
% 4) 初始诊断仅比较 UWS 与 MCS-；冻结模型原本即排除基线 MCS+/EMCS。
%
% MATLAB 建议版本：R2020a 或更高（低版本也包含 print 导出回退方案）。

%% ======================== 1. 输入与输出路径 =========================
% 如需使用绝对路径，可将上面四行替换为类似：
DEV_XLSX  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
TEST_XLSX = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\测试集结果.xlsx';
MODEL_MAT = 'F:\EEG\LowVoltageResearch\OutputResult\V3\3模型对比\模型\joint_model_frozen.mat';
OUTPUT_DIR = 'F:\EEG\LowVoltageResearch\OutputResult\V3\3模型对比\分亚组比较\boostrap_Mean+SD';

OUTPUT_PNG  = fullfile(OUTPUT_DIR, 'subgroup_performance_bootstrap_mean_sd.png');
OUTPUT_TIFF = fullfile(OUTPUT_DIR, 'subgroup_performance_bootstrap_mean_sd.tif');
OUTPUT_PDF  = fullfile(OUTPUT_DIR, 'subgroup_performance_bootstrap_mean_sd.pdf');
OUTPUT_FIG  = fullfile(OUTPUT_DIR, 'subgroup_performance_bootstrap_mean_sd.fig');
OUTPUT_XLSX = fullfile(OUTPUT_DIR, 'subgroup_performance_bootstrap_mean_sd.xlsx');

%% ========================== 2. 可调参数 =============================
AGE_CUTOFF_MODE = 'development_median';  % 可选：'development_median' 或数值
DURATION_CUTOFF_MODE = 'development_median'; % 可选：'development_median' 或数值
REFERENCE_VALUE = NaN;                  % 设为 NaN 可取消竖向参考线
CI_LEVEL = 0.95;
BOOTSTRAP_N = 2000;                    % 建议至少 2000 次
BOOTSTRAP_SEED = 20260723;             % 固定随机种子，保证结果可重复

% 绘图样式
COLOR_DEV  = [148 206 229]/255;   % #94cee5 蓝  
COLOR_TEST = [233/255, 162/255, 163/255];  % #e9a2a3  红 
COLOR_IMPROVED = [182 208 147]/255;       % #b6d093 有改善：绿色
COLOR_NOT_IMPROVED = [205 205 205]/255;   % 无改善：灰色
DEV_COHORT_LABEL = 'M_joint_Tiantan';
TEST_COHORT_LABEL = 'M_joint_Fengtai (external validation)';
FONT_NAME = 'Arial';
FONT_SIZE = 10;
LINE_WIDTH = 1.25;
Y_OFFSET = 0.20;          % 增大两个队列结果的上下间距
MARKER_SIZE = 7;
COMPOSITION_BAR_WIDTH = 0.50; % 加高患者结局构成的堆叠条
CAP_HALF_HEIGHT = 0.065;

%% =========================== 3. 检查文件 ============================
inputFiles = {DEV_XLSX, TEST_XLSX, MODEL_MAT};
for i = 1:numel(inputFiles)
    if ~isfile(inputFiles{i})
        error('找不到输入文件：%s\n请在代码第 1 部分修改输入路径。', inputFiles{i});
    end
end
if ~isfolder(OUTPUT_DIR)
    mkdir(OUTPUT_DIR);
end

%% =========================== 4. 读取模型 ============================
S = load(MODEL_MAT);
if isfield(S, 'model')
    model = S.model;
else
    names = fieldnames(S);
    validNames = names(~startsWith(names, '__'));
    if numel(validNames) ~= 1
        error('无法自动确定冻结模型变量。请确认 .mat 文件中包含名为 model 的结构体。');
    end
    model = S.(validNames{1});
end

requiredModelFields = {'intercept','beta','mu','sg','cutoff','cfg', ...
    'PRE_EXCL_THR','GRP_THR','freqGate'};
for i = 1:numel(requiredModelFields)
    if ~isfield(model, requiredModelFields{i})
        error('冻结模型缺少字段：model.%s', requiredModelFields{i});
    end
end

%% =========================== 5. 读取数据 ============================
% 兼容新旧 MATLAB：
% 新版优先保留 Excel 原始列名；旧版不支持 VariableNamingRule 时自动回退。
TdevRaw  = readTableVersionRobust(DEV_XLSX);
TtestRaw = readTableVersionRobust(TEST_XLSX);

Ddev  = applyFrozenJointModel(TdevRaw,  model, DEV_COHORT_LABEL);
Dtest = applyFrozenJointModel(TtestRaw, model, TEST_COHORT_LABEL);

% 年龄分界值：只由开发集确定，然后固定应用于测试集
if isnumeric(AGE_CUTOFF_MODE)
    ageCutoff = double(AGE_CUTOFF_MODE);
elseif strcmpi(AGE_CUTOFF_MODE, 'development_median')
    ageCutoff = median(Ddev.age);
else
    error('AGE_CUTOFF_MODE 必须为 ''development_median'' 或一个数值。');
end

% 病程分界值：取开发集病程中位数并四舍五入为整数，然后固定应用于测试集。
if isnumeric(DURATION_CUTOFF_MODE)
    durationCutoff = round(double(DURATION_CUTOFF_MODE));
    durationCutoffSource = 'User-specified numeric value, rounded to integer';
elseif strcmpi(DURATION_CUTOFF_MODE, 'development_median')
    durationCutoff = round(median(Ddev.duration));
    durationCutoffSource = 'Rounded development-set median';
else
    error('DURATION_CUTOFF_MODE 必须为 ''development_median'' 或一个数值。');
end

fprintf('\n========== Subgroup analysis ==========%s', newline);
fprintf('Development set: raw N = %d; eligible complete cases = %d%s', ...
    height(TdevRaw), Ddev.n, newline);
fprintf('External test set: raw N = %d; eligible complete cases = %d%s', ...
    height(TtestRaw), Dtest.n, newline);
fprintf('Frozen cutoff = %.6f%s', double(model.cutoff), newline);
fprintf('Age cutoff (development-set median) = %s years%s', ...
    formatNumber(ageCutoff), newline);
fprintf('Duration cutoff (rounded development-set median) = %s months%s\n', ...
    formatNumber(durationCutoff), newline);

% 验证开发集整体结果是否与冻结模型记录一致（如字段存在）
overallDev = calculateMetrics(Ddev.y, Ddev.pred, CI_LEVEL);
overallTest = calculateMetrics(Dtest.y, Dtest.pred, CI_LEVEL);
rng(BOOTSTRAP_SEED, 'twister');
overallBootDev = bootstrapMetrics(Ddev.y, Ddev.pred, BOOTSTRAP_N);
overallBootTest = bootstrapMetrics(Dtest.y, Dtest.pred, BOOTSTRAP_N);
if isfield(model, 'sens_dev_at_cutoff') && isfield(model, 'spec_dev_at_cutoff') && ...
        isfield(model, 'acc_dev_at_cutoff')
    recorded = [double(model.sens_dev_at_cutoff), ...
                double(model.spec_dev_at_cutoff), ...
                double(model.acc_dev_at_cutoff)];
    recalculated = [overallDev.sensitivity, overallDev.specificity, overallDev.accuracy];
    if any(abs(recorded - recalculated) > 1e-10)
        warning(['使用当前 Excel 数据重新计算的开发集整体指标与冻结模型中记录值不一致。' ...
            '请检查输入文件、变量名和数据版本。']);
    else
        fprintf('Frozen-model check passed: sensitivity %.3f, specificity %.3f, accuracy %.3f.%s\n', ...
            recalculated(1), recalculated(2), recalculated(3), newline);
    end
end

%% ======================= 6. 定义临床亚组 ===========================
% 每一行：Section, Label, Variable, Operator, Cutoff/Code
ageLowLabel  = sprintf('< %s years',  formatNumber(ageCutoff));
ageHighLabel = sprintf('≥ %s years', formatNumber(ageCutoff));
durLowLabel  = sprintf('< %s months',  formatNumber(durationCutoff));
durHighLabel = sprintf('≥ %s months', formatNumber(durationCutoff));

defs = {
    'Initial diagnosis', 'UWS',                'preDiag',  'eq', 1;
    'Initial diagnosis', 'MCS-',               'preDiag',  'eq', 2;
    'Etiology',          'Traumatic',          'etiology', 'eq', 1;
    'Etiology',          'Stroke',             'etiology', 'eq', 2;
    'Etiology',          'Anoxia',              'etiology', 'eq', 3;
    'Age',               ageLowLabel,          'age',      'lt', ageCutoff;
    'Age',               ageHighLabel,         'age',      'ge', ageCutoff;
    'Disease duration',  durLowLabel,          'duration', 'lt', durationCutoff;
    'Disease duration',  durHighLabel,         'duration', 'ge', durationCutoff;
    };

nGroups = size(defs, 1);
Mdev  = repmat(emptyMetricStruct(), nGroups, 1);
Mtest = repmat(emptyMetricStruct(), nGroups, 1);
Bdev  = repmat(emptyBootstrapMetricStruct(), nGroups, 1);
Btest = repmat(emptyBootstrapMetricStruct(), nGroups, 1);

for i = 1:nGroups
    maskDev  = buildMask(Ddev,  defs{i,3}, defs{i,4}, defs{i,5});
    maskTest = buildMask(Dtest, defs{i,3}, defs{i,4}, defs{i,5});
    yDevGroup = Ddev.y(maskDev);
    predDevGroup = Ddev.pred(maskDev);
    yTestGroup = Dtest.y(maskTest);
    predTestGroup = Dtest.pred(maskTest);
    Mdev(i)  = calculateMetrics(yDevGroup,  predDevGroup,  CI_LEVEL);
    Mtest(i) = calculateMetrics(yTestGroup, predTestGroup, CI_LEVEL);
    Bdev(i)  = bootstrapMetrics(yDevGroup,  predDevGroup,  BOOTSTRAP_N);
    Btest(i) = bootstrapMetrics(yTestGroup, predTestGroup, BOOTSTRAP_N);
end

%% ======================= 7. 导出数值结果 ============================
resultTable = buildResultTable(defs, Mdev, Mtest, ...
    DEV_COHORT_LABEL, TEST_COHORT_LABEL);
bootstrapTable = buildBootstrapResultTable(defs, ...
    [overallDev; Mdev], [overallTest; Mtest], ...
    [overallBootDev; Bdev], [overallBootTest; Btest], ...
    DEV_COHORT_LABEL, TEST_COHORT_LABEL);
cohortTable = table( ...
    {DEV_COHORT_LABEL; TEST_COHORT_LABEL}, ...
    [height(TdevRaw); height(TtestRaw)], ...
    [Ddev.nAfterBaselineExclusion; Dtest.nAfterBaselineExclusion], ...
    [Ddev.n; Dtest.n], ...
    [sum(Ddev.y); sum(Dtest.y)], ...
    [sum(~Ddev.y); sum(~Dtest.y)], ...
    [overallDev.sensitivity; overallTest.sensitivity], ...
    [overallDev.specificity; overallTest.specificity], ...
    [overallDev.accuracy; overallTest.accuracy], ...
    'VariableNames', {'Cohort','Raw_N','N_after_baseline_exclusion','Complete_case_N', ...
    'Outcome_positive_N','Outcome_negative_N','Sensitivity','Specificity','Accuracy'});

settingsTable = table( ...
    string({'Frozen probability cutoff'; 'Age cutoff (years)'; 'Duration cutoff (months)'; ...
            'Confidence interval'; 'Age cutoff source'; 'Duration cutoff source'; 'Bootstrap iterations'; ...
            'Bootstrap seed'; 'Bootstrap method'}), ...
    string({formatNumber(double(model.cutoff)); formatNumber(ageCutoff); ...
            formatNumber(durationCutoff); sprintf('%.0f%% Wilson', CI_LEVEL*100); ...
            'Development-set median'; durationCutoffSource; formatNumber(BOOTSTRAP_N); ...
            formatNumber(BOOTSTRAP_SEED); ...
            'Patient-level, outcome-stratified; frozen model and cutoff'}), ...
    'VariableNames', {'Setting','Value'});

if isfile(OUTPUT_XLSX)
    delete(OUTPUT_XLSX);
end
writetable(resultTable,  OUTPUT_XLSX, 'Sheet', 'Subgroup_metrics');
writetable(bootstrapTable, OUTPUT_XLSX, 'Sheet', 'Bootstrap_mean_SD');
writetable(cohortTable,  OUTPUT_XLSX, 'Sheet', 'Cohort_summary');
writetable(settingsTable, OUTPUT_XLSX, 'Sheet', 'Settings');

%% ============== 8. 绘制 bootstrap 均值 ± SD 点线图 ===============
% 点表示 bootstrap 均值，水平线及两端端帽表示 ±1 bootstrap SD。
rowY = [15, 12, 11, 8.5, 7.5, 6.5, 4, 3, 0.5, -0.5];
plotLabels = [{''}; defs(:,2)];
BdevPlot = [overallBootDev; Bdev];
BtestPlot = [overallBootTest; Btest];
MdevPlot = [overallDev; Mdev];
MtestPlot = [overallTest; Mtest];
nPlotGroups = numel(BdevPlot);

sectionNames = {'Diagnosis','Etiology','Age','Duration'};
sectionY = [13.0, 9.5, 5.0, 1.5];
yLimits = [-1.25, 15.75];

metricMeanFields = {'sensitivityMean','specificityMean','accuracyMean'};
metricSDFields = {'sensitivitySD','specificitySD','accuracySD'};
metricLabels = {'Sensitivity','Specificity','Accuracy'};

% 前 3 个子图显示模型性能；第 4 个子图显示各组有/无改善的人数构成。
% 加宽画布，为每个性能子图右侧的 Mean 数值列留出空间。
fig = figure('Color','w', 'Units','centimeters', 'Position',[1, 2, 48, 21]);
tl = tiledlayout(fig, 1, 4, 'TileSpacing','none', 'Padding','compact');
% 顶部预留独立图例带；底部和左侧预留刻度及亚组标签空间。
tl.Units = 'normalized';
tl.Position = [0.05, 0.05, 0.945, 0.88];

for p = 1:3
    ax = nexttile(tl, p);
    hold(ax, 'on');
    ax.FontName = FONT_NAME;
    ax.FontSize = FONT_SIZE;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.Layer = 'top';
    ax.YLim = yLimits;
    % 横轴从 0 开始，使纵轴与水平柱之间没有空隙。
    ax.XLim = [0, 1.22];
    ax.XTick = 0:0.2:1.0;
    drawLimitedBottomAxis(ax, 0:0.2:1.0, false, yLimits, ...
        FONT_NAME, FONT_SIZE);

    [yTickAscending, yTickOrder] = sort(double(rowY(:)), 'ascend');
    yTickLabels = plotLabels(yTickOrder);
    ax.TickLabelInterpreter = 'none';
    if p == 1
        % 亚组名称只在第一个性能面板显示，避免 3 个面板重复。
        set(ax, 'YTick', yTickAscending, 'YTickLabel', yTickLabels);
        text(ax, -0.065, rowY(1), 'Overall', ...
            'FontName',FONT_NAME, 'FontSize',FONT_SIZE, 'FontWeight','bold', ...
            'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
            'Interpreter','none', 'Clipping','off');
    else
        set(ax, 'YTick', [], 'YTickLabel', {});
    end

    if isfinite(REFERENCE_VALUE)
        line(ax, [REFERENCE_VALUE REFERENCE_VALUE], yLimits, ...
            'Color',[0.45 0.70 0.88], 'LineStyle','--', 'LineWidth',1.0, ...
            'HandleVisibility','off');
    end
    devMean = arrayfun(@(x) x.(metricMeanFields{p}), BdevPlot);
    devSD = arrayfun(@(x) x.(metricSDFields{p}), BdevPlot);
    testMean = arrayfun(@(x) x.(metricMeanFields{p}), BtestPlot);
    testSD = arrayfun(@(x) x.(metricSDFields{p}), BtestPlot);

    % 点位为 bootstrap 均值；水平线及两端端帽表示 mean ± 1 bootstrap SD。
    for i = 1:nPlotGroups
        drawBootstrapSD(ax, devMean(i), devSD(i), rowY(i)+Y_OFFSET, ...
            COLOR_DEV, LINE_WIDTH, CAP_HALF_HEIGHT);
        drawBootstrapSD(ax, testMean(i), testSD(i), rowY(i)-Y_OFFSET, ...
            COLOR_TEST, LINE_WIDTH, CAP_HALF_HEIGHT);
        if ~isfinite(devMean(i))
            text(ax, 0.50, rowY(i)+Y_OFFSET, 'NA', 'Color',COLOR_DEV, ...
                'FontWeight','bold', 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontSize',FONT_SIZE-1);
        end
        if ~isfinite(testMean(i))
            text(ax, 0.50, rowY(i)-Y_OFFSET, 'NA', 'Color',COLOR_TEST, ...
                'FontWeight','bold', 'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', 'FontSize',FONT_SIZE-1);
        end

        % 固定在每个子图右侧仅显示 bootstrap Mean；SD 仍由误差线表示。
        text(ax, 1.055, rowY(i)+Y_OFFSET, formatMean(devMean(i)), ...
            'Color',[0 0 0], 'FontName',FONT_NAME, 'FontSize',FONT_SIZE-1, ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
            'Interpreter','none', 'Clipping','on');
        text(ax, 1.055, rowY(i)-Y_OFFSET, formatMean(testMean(i)), ...
            'Color',[0 0 0], 'FontName',FONT_NAME, 'FontSize',FONT_SIZE-1, ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
            'Interpreter','none', 'Clipping','on');
    end

    % 均值标记最后绘制，确保圆点/三角形位于 SD 线之上。
    hDev = plot(ax, devMean, rowY+Y_OFFSET, 'o', ...
        'LineStyle','none', 'MarkerSize',MARKER_SIZE, ...
        'MarkerFaceColor',COLOR_DEV, 'MarkerEdgeColor',COLOR_DEV);
    hTest = plot(ax, testMean, rowY-Y_OFFSET, '^', ...
        'LineStyle','none', 'MarkerSize',MARKER_SIZE, ...
        'MarkerFaceColor',COLOR_TEST, 'MarkerEdgeColor',COLOR_TEST);


    text(ax, 1.055, 15.55, 'Mean', ...
        'FontName',FONT_NAME, 'FontSize',FONT_SIZE-1, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
        'Interpreter','none', 'Clipping','on');

    if p == 1
        % 第一个面板紧邻画布左边界，标题稍向右移以避免 Diagnosis 被裁切。
        sectionLabelX = -0.020;
        for s = 1:numel(sectionNames)
            text(ax, sectionLabelX, sectionY(s), sectionNames{s}, ...
                'FontName',FONT_NAME, 'FontSize',FONT_SIZE+0.5, 'FontWeight','bold', ...
                'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
                'Interpreter','none', 'Clipping','off');
        end
    end

    % 指标小标题放在上方；数值刻度仍保留在横轴下方。
    title(ax, metricLabels{p}, ...
        'FontName',FONT_NAME, 'FontSize',FONT_SIZE+1, ...
        'FontWeight','bold', 'Interpreter','none');

end

%% ========== 8.1 各组总人数及有/无改善人数占比（最右侧） ==========
axComp = nexttile(tl, 4);
hold(axComp, 'on');
axComp.FontName = FONT_NAME;
axComp.FontSize = FONT_SIZE;
axComp.LineWidth = 0.8;
axComp.TickDir = 'out';
axComp.Box = 'off';
axComp.Layer = 'top';
axComp.YLim = yLimits;
% 横轴从 0 开始，使纵轴与堆叠条之间没有空隙。
axComp.XLim = [0, 1.46];
axComp.XTick = 0:0.2:1.0;
axComp.XTickLabel = arrayfun(@(x) sprintf('%d%%', round(100*x)), ...
    axComp.XTick, 'UniformOutput', false);
drawLimitedBottomAxis(axComp, 0:0.2:1.0, true, yLimits, ...
    FONT_NAME, FONT_SIZE);

% 与左侧性能图逐行对齐，因此不再重复显示相同的 Y 轴文字。
set(axComp, 'YTick', [], 'YTickLabel', {});
axComp.TickLabelInterpreter = 'none';

text(axComp, 1.025, 15.55, 'Number', ...
    'FontName',FONT_NAME, 'FontSize',FONT_SIZE-1, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
    'Interpreter','none', 'Clipping','on');

devN = arrayfun(@(x) x.n, MdevPlot);
testN = arrayfun(@(x) x.n, MtestPlot);
devImprovedN = arrayfun(@(x) x.positiveN, MdevPlot);
testImprovedN = arrayfun(@(x) x.positiveN, MtestPlot);

devImprovedProp = safeProportion(devImprovedN, devN);
testImprovedProp = safeProportion(testImprovedN, testN);
devNotImprovedProp = safeProportion(devN-devImprovedN, devN);
testNotImprovedProp = safeProportion(testN-testImprovedN, testN);

% 将开发集/测试集交错排列，使其与左侧 3 个性能子图中的上下顺序一致。
compositionY = reshape([rowY+Y_OFFSET; rowY-Y_OFFSET], 1, []);
improvedProp = reshape([devImprovedProp(:)'; testImprovedProp(:)'], 1, []);
notImprovedProp = reshape([devNotImprovedProp(:)'; testNotImprovedProp(:)'], 1, []);
totalN = reshape([devN(:)'; testN(:)'], 1, []);

hComposition = barh(axComp, compositionY, ...
    [improvedProp(:), notImprovedProp(:)], 'stacked', ...
    'BarWidth',COMPOSITION_BAR_WIDTH, ...
    'EdgeColor','none');
hComposition(1).FaceColor = COLOR_IMPROVED;
hComposition(2).FaceColor = COLOR_NOT_IMPROVED;

for i = 1:numel(compositionY)
    if totalN(i) <= 0
        text(axComp, 0.50, compositionY(i), 'NA', ...
            'Color',[0.35 0.35 0.35], 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'FontSize',FONT_SIZE-1);
        text(axComp, 1.025, compositionY(i), '0 (0/0)', ...
            'FontName',FONT_NAME, 'FontSize',FONT_SIZE-1, ...
            'HorizontalAlignment','left', 'VerticalAlignment','middle');
        continue;
    end

    improvedCount = round(improvedProp(i) * totalN(i));
    notImprovedCount = totalN(i) - improvedCount;
    text(axComp, 1.025, compositionY(i), ...
        sprintf('%d (%d/%d)', totalN(i), improvedCount, notImprovedCount), ...
        'FontName',FONT_NAME, 'FontSize',FONT_SIZE-1, ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
        'Interpreter','none', 'Clipping','on');
end

title(axComp, 'Patient proportion', ...
    'FontName',FONT_NAME, 'FontSize',FONT_SIZE+1, ...
    'FontWeight','bold', 'Interpreter','none');

% 使用覆盖整张 figure 的透明坐标层创建统一图例，使其按四联图整体居中。
legendAxes = axes(fig, 'Position',[0 0 1 1], 'Visible','off', ...
    'Color','none', 'HitTest','off', 'HandleVisibility','off');
hold(legendAxes, 'on');
hDevLegend = plot(legendAxes, NaN, NaN, 'o', ...
    'LineStyle','none', 'MarkerSize',MARKER_SIZE, ...
    'MarkerFaceColor',COLOR_DEV, 'MarkerEdgeColor',COLOR_DEV);
hTestLegend = plot(legendAxes, NaN, NaN, '^', ...
    'LineStyle','none', 'MarkerSize',MARKER_SIZE, ...
    'MarkerFaceColor',COLOR_TEST, 'MarkerEdgeColor',COLOR_TEST);
hImprovedLegend = plot(legendAxes, NaN, NaN, 's', ...
    'LineStyle','none', 'MarkerSize',MARKER_SIZE+1, ...
    'MarkerFaceColor',COLOR_IMPROVED, 'MarkerEdgeColor',COLOR_IMPROVED);
hNotImprovedLegend = plot(legendAxes, NaN, NaN, 's', ...
    'LineStyle','none', 'MarkerSize',MARKER_SIZE+1, ...
    'MarkerFaceColor',COLOR_NOT_IMPROVED, 'MarkerEdgeColor',COLOR_NOT_IMPROVED);
lgd = legend(legendAxes, ...
    [hDevLegend hTestLegend hImprovedLegend hNotImprovedLegend], ...
    {DEV_COHORT_LABEL, TEST_COHORT_LABEL, 'Favorable outcome', 'Unfavorable outcome'}, ...
    'Location','north', 'Orientation','horizontal', ...
    'Box','off', 'Interpreter','none', 'FontName',FONT_NAME, ...
    'FontSize',FONT_SIZE);
try
    lgd.ItemTokenSize = [18, 10];
catch
    % 较旧 MATLAB 版本可能不支持 ItemTokenSize。
end

drawnow;
% 解除自动 Location 布局后，以整张 figure 的中心（x=0.5）固定图例。
lgd.Units = 'normalized';
lgd.Location = 'none';
legendPosition = lgd.Position;
legendPosition(1) = 0.5 - legendPosition(3)/2;
% 将图例顶部固定在 figure 的 99.5% 高度，避免与各面板小标题重叠。
legendPosition(2) = 0.995 - legendPosition(4);
lgd.Position = legendPosition;

%% =========================== 9. 导出图片 ============================
savefig(fig, OUTPUT_FIG);

% 优先使用 exportgraphics；旧版 MATLAB 自动回退至 print
try
    exportgraphics(fig, OUTPUT_PNG,  'Resolution',600);
    exportgraphics(fig, OUTPUT_TIFF, 'Resolution',600);
    exportgraphics(fig, OUTPUT_PDF,  'ContentType','vector');
catch ME
    warning('exportgraphics 失败，改用 print 导出：%s', ME.message);
    set(fig, 'PaperPositionMode','auto');
    print(fig, OUTPUT_PNG,  '-dpng',  '-r600');
    print(fig, OUTPUT_TIFF, '-dtiff', '-r600');
    print(fig, OUTPUT_PDF,  '-dpdf',  '-painters');
end

fprintf('结果已输出至：%s%s', OUTPUT_DIR, newline);
fprintf('Figure: %s%s', OUTPUT_PNG, newline);
fprintf('Metrics: %s%s\n', OUTPUT_XLSX, newline);

end

%% ===================================================================
%                           本地函数
% ====================================================================
function T = readTableVersionRobust(inputFile)
%READTABLEVERSIONROBUST 兼容不同 MATLAB 版本读取 Excel。
%
% 优先保留 Excel 原始表头。对于不支持 VariableNamingRule 的旧版本，
% 自动回退到 PreserveVariableNames 或默认 readtable。旧版自动修改的
% 原始表头通常会保存在 T.Properties.VariableDescriptions 中，后续
% resolveTableVariableName 会据此恢复映射。

% 关闭旧版本自动修改变量名时的提示；若该 warning ID 不存在也不会报错。
try
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
catch
end

% 方案 1：较新版本，detectImportOptions 支持 VariableNamingRule。
try
    opts = detectImportOptions(inputFile, 'VariableNamingRule', 'preserve');
    T = readtable(inputFile, opts);
    return;
catch
end

% 方案 2：部分版本的 readtable 可直接接受 VariableNamingRule。
try
    T = readtable(inputFile, 'VariableNamingRule', 'preserve');
    return;
catch
end

% 方案 3：较旧版本可能使用 PreserveVariableNames。
try
    T = readtable(inputFile, 'PreserveVariableNames', true);
    return;
catch
end

% 方案 4：detectImportOptions 存在，但不支持上述命名参数。
try
    opts = detectImportOptions(inputFile);
    try
        opts.PreserveVariableNames = true;
    catch
        % 当前版本没有 PreserveVariableNames 属性，继续使用默认 opts。
    end
    T = readtable(inputFile, opts);
    return;
catch
end

% 最后回退：使用该 MATLAB 版本的默认读取方式。
T = readtable(inputFile);
end

function actualName = resolveTableVariableName(T, requestedName)
%RESOLVETABLEVARIABLENAME 将模型中记录的原始列名映射到当前表变量名。

actualName = '';
requestedName = char(requestedName);
variableNames = T.Properties.VariableNames;

% 1) 实际变量名完全一致。
idx = find(strcmp(variableNames, requestedName), 1);
if ~isempty(idx)
    actualName = variableNames{idx};
    return;
end

% 2) 忽略大小写匹配。
idx = find(strcmpi(variableNames, requestedName), 1);
if ~isempty(idx)
    actualName = variableNames{idx};
    return;
end

% 3) 旧版 readtable 常将原始 Excel 表头保存在 VariableDescriptions。
try
    descriptions = T.Properties.VariableDescriptions;
catch
    descriptions = {};
end
if iscell(descriptions) && numel(descriptions) == numel(variableNames)
    idx = find(strcmp(descriptions, requestedName), 1);
    if isempty(idx)
        idx = find(strcmpi(descriptions, requestedName), 1);
    end
    if ~isempty(idx)
        actualName = variableNames{idx};
        return;
    end
end

% 4) 按 MATLAB 的合法变量名转换规则匹配。
try
    validRequestedName = matlab.lang.makeValidName(requestedName);
    idx = find(strcmp(variableNames, validRequestedName), 1);
    if isempty(idx)
        idx = find(strcmpi(variableNames, validRequestedName), 1);
    end
    if ~isempty(idx)
        actualName = variableNames{idx};
        return;
    end
catch
end

% 5) 最后进行宽松规范化匹配（去除空格、下划线及常见标点）。
normalizedRequested = normalizeVariableName(requestedName);
normalizedNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    normalizedNames{i} = normalizeVariableName(variableNames{i});
end
idx = find(strcmp(normalizedNames, normalizedRequested));
if numel(idx) == 1
    actualName = variableNames{idx};
end
end

function s = normalizeVariableName(s)
%NORMALIZEVARIABLENAME 用于新旧 MATLAB 列名之间的宽松匹配。
s = char(s);
try
    s = matlab.lang.makeValidName(s);
catch
end
s = lower(s);
s = regexprep(s, '[\s_\-\(\)\[\]\{\}/\\\.\,\:\;]+', '');
end

function D = applyFrozenJointModel(Traw, model, cohortName)
% 根据冻结参数生成概率与二分类预测，完全不重新拟合。

cfg = model.cfg;

% cfg 中保存的是建模时使用的列名。旧版 MATLAB 读取 Excel 时可能自动
% 修改不符合变量命名规则的表头，因此这里同时依据实际变量名、
% VariableDescriptions 中保存的原始表头及 makeValidName 结果进行匹配。
requestedColumns = { ...
    char(cfg.COL_AGE), ...
    char(cfg.COL_DURATION), ...
    char(cfg.COL_PRE_CRS), ...
    char(cfg.COL_PRE_DIAG), ...
    char(cfg.COL_POST_DIAG), ...
    char(cfg.COL_ETIOLOGY), ...
    char(cfg.COL_VAR1), ...
    char(cfg.COL_VAR2), ...
    char(cfg.COL_ALPHA_STAB), ...
    char(cfg.COL_DELTA_STAB)};

resolvedColumns = cell(size(requestedColumns));
missingColumns = {};
for c = 1:numel(requestedColumns)
    resolvedColumns{c} = resolveTableVariableName(Traw, requestedColumns{c});
    if isempty(resolvedColumns{c})
        missingColumns{end+1} = requestedColumns{c}; %#ok<AGROW>
    end
end
if ~isempty(missingColumns)
    error(['%s 缺少必要变量：%s\n当前 Excel 被 MATLAB 读取后的变量名为：%s'], ...
        cohortName, strjoin(missingColumns, ', '), ...
        strjoin(Traw.Properties.VariableNames, ', '));
end

colAge       = resolvedColumns{1};
colDuration  = resolvedColumns{2};
colPreCRS    = resolvedColumns{3};
colPreDiag   = resolvedColumns{4};
colPostDiag  = resolvedColumns{5};
colEtiology  = resolvedColumns{6};
colFreq      = resolvedColumns{7};
colBroadAmp  = resolvedColumns{8};
colAlphaStab = resolvedColumns{9};
colDeltaStab = resolvedColumns{10};

age       = numericColumn(Traw, colAge);
duration  = numericColumn(Traw, colDuration);
preCRS    = numericColumn(Traw, colPreCRS);
preDiag   = numericColumn(Traw, colPreDiag);
postDiag  = numericColumn(Traw, colPostDiag);
etiology  = numericColumn(Traw, colEtiology);
freq      = numericColumn(Traw, colFreq);
broadAmp  = numericColumn(Traw, colBroadAmp);
alphaStab = numericColumn(Traw, colAlphaStab);
deltaStab = numericColumn(Traw, colDeltaStab);

baselineEligible = isfinite(preDiag) & preDiag < double(model.PRE_EXCL_THR);
complete = isfinite(age) & isfinite(duration) & isfinite(preCRS) & ...
    isfinite(preDiag) & isfinite(postDiag) & isfinite(etiology) & ...
    isfinite(freq) & isfinite(broadAmp) & isfinite(alphaStab) & isfinite(deltaStab);
keep = baselineEligible & complete;

nAfterBaselineExclusion = sum(baselineEligible);
if any(baselineEligible & ~complete)
    warning('%s：%d 例满足基线诊断条件，但因模型变量缺失被排除。', ...
        cohortName, sum(baselineEligible & ~complete));
end

age       = age(keep);
duration  = duration(keep);
preCRS    = preCRS(keep);
preDiag   = preDiag(keep);
postDiag  = postDiag(keep);
etiology  = etiology(keep);
freq      = freq(keep);
broadAmp  = broadAmp(keep);
alphaStab = alphaStab(keep);
deltaStab = deltaStab(keep);

% 主导频率稳定性：频率达到 freqGate 时使用 Alpha_Stability，
% 否则使用 DeltaTheta_Stability；此逻辑与冻结模型构建一致。
domStability = deltaStab;
useAlpha = freq >= double(model.freqGate);
domStability(useAlpha) = alphaStab(useAlpha);

% 病因参照组为 trauma=1；模型包含 stroke 和 ischemia 两个哑变量。
if isfield(model, 'etioRef') && double(model.etioRef) ~= 1
    error('当前代码按 etioRef=1（trauma）实现，但冻结模型 etioRef=%g。', double(model.etioRef));
end

X = [age, duration, preCRS, preDiag, ...
     double(etiology == 2), double(etiology == 3), ...
     freq, broadAmp, domStability];

mu = double(model.mu(:))';
sg = double(model.sg(:))';
beta = double(model.beta(:));
if size(X,2) ~= numel(mu) || numel(mu) ~= numel(sg) || numel(beta) ~= numel(mu)
    error('冻结模型参数维度与当前构造的 9 个特征不一致。');
end
if any(~isfinite(sg) | sg == 0)
    error('冻结模型标准差 model.sg 含 0 或非有限值。');
end

% 兼容较旧 MATLAB 版本，使用 bsxfun 完成标准化
Z = bsxfun(@rdivide, bsxfun(@minus, X, mu), sg);
linearPredictor = double(model.intercept) + Z * beta;
probability = 1 ./ (1 + exp(-linearPredictor));
prediction = probability >= double(model.cutoff);
outcome = postDiag >= double(model.GRP_THR);

D = struct();
D.cohort = cohortName;
D.n = numel(outcome);
D.nAfterBaselineExclusion = nAfterBaselineExclusion;
D.age = age;
D.duration = duration;
D.preDiag = preDiag;
D.etiology = etiology;
D.y = logical(outcome);
D.p = probability;
D.pred = logical(prediction);
end

function x = numericColumn(T, variableName)
% 将 numeric/logical/cell/string/categorical 列稳健转换为 double。
v = T.(variableName);
if isnumeric(v) || islogical(v)
    x = double(v);
elseif iscell(v)
    x = str2double(string(v));
elseif isstring(v) || ischar(v) || iscategorical(v)
    x = str2double(string(v));
else
    error('变量 %s 无法转换为数值型。', variableName);
end
x = x(:);
end

function mask = buildMask(D, variableName, operatorName, value)
x = D.(variableName);
switch lower(operatorName)
    case 'eq'
        mask = x == value;
    case 'lt'
        mask = x < value;
    case 'le'
        mask = x <= value;
    case 'gt'
        mask = x > value;
    case 'ge'
        mask = x >= value;
    otherwise
        error('不支持的亚组运算符：%s', operatorName);
end
mask = logical(mask(:));
end

function M = calculateMetrics(y, pred, ciLevel)
y = logical(y(:));
pred = logical(pred(:));
if numel(y) ~= numel(pred)
    error('y 与 pred 长度不一致。');
end

M = emptyMetricStruct();
M.n = numel(y);
M.positiveN = sum(y);
M.negativeN = sum(~y);
M.tp = sum(pred & y);
M.fn = sum(~pred & y);
M.tn = sum(~pred & ~y);
M.fp = sum(pred & ~y);

[M.sensitivity, M.sensitivityL, M.sensitivityU] = ...
    wilsonInterval(M.tp, M.positiveN, ciLevel);
[M.specificity, M.specificityL, M.specificityU] = ...
    wilsonInterval(M.tn, M.negativeN, ciLevel);
[M.accuracy, M.accuracyL, M.accuracyU] = ...
    wilsonInterval(M.tp + M.tn, M.n, ciLevel);
end

function M = emptyMetricStruct()
M = struct( ...
    'n',0,'positiveN',0,'negativeN',0,'tp',0,'fn',0,'tn',0,'fp',0, ...
    'sensitivity',NaN,'sensitivityL',NaN,'sensitivityU',NaN, ...
    'specificity',NaN,'specificityL',NaN,'specificityU',NaN, ...
    'accuracy',NaN,'accuracyL',NaN,'accuracyU',NaN);
end

function B = bootstrapMetrics(y, pred, nBoot)
%BOOTSTRAPMETRICS 患者层面、按真实结局分层的非参数 bootstrap。
% 模型预测和分类阈值保持冻结；每次只重新抽取患者。
% 本函数假定 Excel 中每一行对应一名独立患者；若同一患者有多行记录，
% 应改为按患者 ID 进行整簇 bootstrap。
y = logical(y(:));
pred = logical(pred(:));
if numel(y) ~= numel(pred)
    error('y 与 pred 长度不一致。');
end
if ~isscalar(nBoot) || nBoot < 2 || nBoot ~= round(nBoot)
    error('BOOTSTRAP_N 必须是大于等于 2 的整数。');
end

B = emptyBootstrapMetricStruct();
B.nBoot = nBoot;
B.n = numel(y);
B.positiveN = sum(y);
B.negativeN = sum(~y);

sensValues = NaN(nBoot,1);
specValues = NaN(nBoot,1);
accValues = NaN(nBoot,1);
predPositive = pred(y);
predNegative = pred(~y);

for b = 1:nBoot
    correctN = 0;
    if B.positiveN > 0
        sampledPositive = predPositive(randi(B.positiveN, B.positiveN, 1));
        sensValues(b) = mean(sampledPositive);
        correctN = correctN + sum(sampledPositive);
    end
    if B.negativeN > 0
        sampledNegative = predNegative(randi(B.negativeN, B.negativeN, 1));
        specValues(b) = mean(~sampledNegative);
        correctN = correctN + sum(~sampledNegative);
    end
    if B.n > 0
        accValues(b) = correctN / B.n;
    end
end

[B.sensitivityMean, B.sensitivitySD] = finiteMeanSD(sensValues);
[B.specificityMean, B.specificitySD] = finiteMeanSD(specValues);
[B.accuracyMean, B.accuracySD] = finiteMeanSD(accValues);
end

function B = emptyBootstrapMetricStruct()
B = struct( ...
    'nBoot',0,'n',0,'positiveN',0,'negativeN',0, ...
    'sensitivityMean',NaN,'sensitivitySD',NaN, ...
    'specificityMean',NaN,'specificitySD',NaN, ...
    'accuracyMean',NaN,'accuracySD',NaN);
end

function [valueMean, valueSD] = finiteMeanSD(values)
values = values(isfinite(values));
if isempty(values)
    valueMean = NaN;
    valueSD = NaN;
else
    valueMean = mean(values);
    valueSD = std(values, 0);
end
end

function [p, lower, upper] = wilsonInterval(k, n, ciLevel)
% Wilson score interval；95% CI 时 z=1.959963984540054。
if n <= 0
    p = NaN;
    lower = NaN;
    upper = NaN;
    return;
end
p = k / n;
alpha = 1 - ciLevel;
% 避免依赖 Statistics and Machine Learning Toolbox：
% 使用 erfcinv 计算标准正态分位数。
z = -sqrt(2) * erfcinv(2 * (1 - alpha/2));
denominator = 1 + z^2/n;
center = (p + z^2/(2*n)) / denominator;
halfWidth = z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / denominator;
lower = max(0, center - halfWidth);
upper = min(1, center + halfWidth);
end

function drawLimitedBottomAxis(ax, tickValues, usePercentLabels, yLimits, ...
    fontName, fontSize)
%DRAWLIMITEDBOTTOMAXIS 重画只到 1 为止的底部横轴及刻度。
% 主坐标轴仍保留额外的右侧空间用于 Mean/N 数值列，但横轴线不会延伸过去。
ax.XColor = 'none';
axisY = yLimits(1);
tickLength = 0.10;
labelY = axisY - tickLength - 0.035;

line(ax, [tickValues(1) tickValues(end)], [axisY axisY], ...
    'Color',[0 0 0], 'LineWidth',0.8, 'HandleVisibility','off', ...
    'Clipping','off');
for i = 1:numel(tickValues)
    x = tickValues(i);
    line(ax, [x x], [axisY axisY-tickLength], ...
        'Color',[0 0 0], 'LineWidth',0.8, 'HandleVisibility','off', ...
        'Clipping','off');
    if usePercentLabels
        tickLabel = sprintf('%d%%', round(100*x));
    elseif abs(x-round(x)) < 1e-10
        tickLabel = sprintf('%.0f', x);
    else
        tickLabel = sprintf('%.1f', x);
    end
    text(ax, x, labelY, tickLabel, ...
        'Color',[0 0 0], 'FontName',fontName, 'FontSize',fontSize, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'Interpreter','none', 'Clipping','off');
end
end

function drawBootstrapSD(ax, valueMean, valueSD, y, colorValue, ...
    lineWidth, capHalfHeight)
%DRAWBOOTSTRAPSD 在水平柱末端绘制 mean ± 1 bootstrap SD。
if ~isfinite(valueMean) || ~isfinite(valueSD)
    return;
end
lowerValue = valueMean - valueSD;
upperValue = valueMean + valueSD;
line(ax, [lowerValue upperValue], [y y], 'Color',colorValue, ...
    'LineWidth',lineWidth, 'HandleVisibility','off');
line(ax, [lowerValue lowerValue], [y-capHalfHeight y+capHalfHeight], ...
    'Color',colorValue, 'LineWidth',lineWidth, 'HandleVisibility','off');
line(ax, [upperValue upperValue], [y-capHalfHeight y+capHalfHeight], ...
    'Color',colorValue, 'LineWidth',lineWidth, 'HandleVisibility','off');
end

function txt = formatMean(valueMean)
%FORMATMEAN 将 bootstrap 均值格式化为右侧数值标签；SD 由误差线表示。
if ~isfinite(valueMean)
    txt = 'NA';
else
    txt = sprintf('%.3f', valueMean);
end
end

function p = safeProportion(numerator, denominator)
%SAFEPROPORTION 避免空亚组除以 0；空亚组返回 0，随后单独标记为 NA。
p = zeros(size(denominator));
valid = denominator > 0;
p(valid) = numerator(valid) ./ denominator(valid);
end

function T = buildResultTable(defs, Mdev, Mtest, devCohortLabel, testCohortLabel)
nGroups = size(defs,1);
rows = cell(nGroups*2, 22);
r = 0;
for i = 1:nGroups
    r = r + 1;
    rows(r,:) = metricRow(defs(i,:), devCohortLabel, Mdev(i));
    r = r + 1;
    rows(r,:) = metricRow(defs(i,:), testCohortLabel, Mtest(i));
end

varNames = {'Section','Subgroup','Cohort','N','Outcome_positive_N','Outcome_negative_N', ...
    'TP','FN','TN','FP', ...
    'Sensitivity','Sensitivity_LCL','Sensitivity_UCL', ...
    'Specificity','Specificity_LCL','Specificity_UCL', ...
    'Accuracy','Accuracy_LCL','Accuracy_UCL', ...
    'Subgroup_variable','Operator','Cutoff_or_code'};
T = cell2table(rows, 'VariableNames', varNames);
end

function T = buildBootstrapResultTable(defs, MdevPlot, MtestPlot, ...
    BdevPlot, BtestPlot, devCohortLabel, testCohortLabel)
%BUILDBOOTSTRAPRESULTTABLE 导出原始点估计、bootstrap 均值与 bootstrap SD。
sectionLabels = [{'Overall'}; defs(:,1)];
subgroupLabels = [{'Overall'}; defs(:,2)];
nPlotGroups = numel(MdevPlot);
rows = cell(nPlotGroups*2, 17);
r = 0;
for i = 1:nPlotGroups
    r = r + 1;
    rows(r,:) = bootstrapMetricRow(sectionLabels{i}, subgroupLabels{i}, ...
        devCohortLabel, MdevPlot(i), BdevPlot(i));
    r = r + 1;
    rows(r,:) = bootstrapMetricRow(sectionLabels{i}, subgroupLabels{i}, ...
        testCohortLabel, MtestPlot(i), BtestPlot(i));
end

varNames = {'Section','Subgroup','Cohort','N','Outcome_positive_N', ...
    'Outcome_negative_N','Original_sensitivity','Bootstrap_sensitivity_mean', ...
    'Bootstrap_sensitivity_SD','Original_specificity', ...
    'Bootstrap_specificity_mean','Bootstrap_specificity_SD', ...
    'Original_accuracy','Bootstrap_accuracy_mean','Bootstrap_accuracy_SD', ...
    'Bootstrap_iterations','Bootstrap_method'};
T = cell2table(rows, 'VariableNames', varNames);
end

function row = bootstrapMetricRow(sectionName, subgroupName, cohortName, M, B)
row = {string(sectionName), string(subgroupName), string(cohortName), ...
    M.n, M.positiveN, M.negativeN, ...
    M.sensitivity, B.sensitivityMean, B.sensitivitySD, ...
    M.specificity, B.specificityMean, B.specificitySD, ...
    M.accuracy, B.accuracyMean, B.accuracySD, B.nBoot, ...
    string('Patient-level outcome-stratified; frozen model and cutoff')};
end

function row = metricRow(defRow, cohortName, M)
row = {string(defRow{1}), string(defRow{2}), string(cohortName), ...
    M.n, M.positiveN, M.negativeN, M.tp, M.fn, M.tn, M.fp, ...
    M.sensitivity, M.sensitivityL, M.sensitivityU, ...
    M.specificity, M.specificityL, M.specificityU, ...
    M.accuracy, M.accuracyL, M.accuracyU, ...
    string(defRow{3}), string(defRow{4}), double(defRow{5})};
end

function txt = formatNumber(x)
if abs(x - round(x)) < 1e-10
    txt = sprintf('%.0f', x);
else
    txt = sprintf('%.1f', x);
end
end
