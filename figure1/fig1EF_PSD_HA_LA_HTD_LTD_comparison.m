%% =========================================================
%  EEG_HA_LA_HTD_LTD_comparison.m
%
%  Figure 1 — HA (Group=1) vs LA (Group=3):
%             Alpha_RawPower (absolute band power)
%  Figure 2 — HTD (Group=2) vs LTD (Group=4):
%             ThetaDelta_RawPower (absolute band power)
%  Relative-power metrics are retained for statistics but are not plotted.
%
%  Statistics: Wilcoxon rank-sum test
%  Data verified: HA=29, HTD=31, LA=22, LTD=19
%% =========================================================

clear; clc; close all;

%% =========================================================
%  Path configuration  ← modify these two lines as needed
% ---------------------------------------------------------
INPUT_FILE   = 'F:\EEG\LowVoltageResearch\OutputResult\V2\全频段功率谱\临床资料_合并功率谱.xlsx';
OUTPUT_DIR   = 'F:\EEG\LowVoltageResearch\OutputResult\V2\全频段功率谱\高低波幅PSD比较';
% ---------------------------------------------------------
OUTPUT_STATS = 'EEG_HA_LA_HTD_LTD_stats.xlsx';
%% =========================================================

% Create output directory if needed
if ~exist(OUTPUT_DIR, 'dir')
    mkdir(OUTPUT_DIR);
    fprintf('Output directory created: %s\n', OUTPUT_DIR);
end

%% --- 1. Read data (full raw including header row) ---
if ~exist(INPUT_FILE, 'file')
    error('Input file not found: %s\nPlease check the path.', INPUT_FILE);
end

fprintf('Reading data...\n');
[~, ~, rawFull] = xlsread(INPUT_FILE);   % row 1 = headers, row 2+ = data

headers  = rawFull(1, :);
findCol  = @(name) find(strcmp(headers, name));

% Helper: extract a numeric column from DATA rows (rows 2:end of rawFull)
dataRows = rawFull(2:end, :);
col2num  = @(colIdx) cellfun( ...
    @(x) ternary(isnumeric(x) && isscalar(x), x, NaN), ...
    dataRows(:, colIdx));

% Group column
colGroup  = findCol('Group(HA=1 HTD=2 LA=3 LTD=4)');
if isempty(colGroup)
    error('Column "Group(HA=1 HTD=2 LA=3 LTD=4)" not found.');
end
groupVals = round(col2num(colGroup));

%% =========================================================
%  Color palette
% ---------------------------------------------------------
COLOR_HA  = [182/255, 208/255, 147/255];  % #b6d093  绿
COLOR_HTD = [245/255, 207/255, 168/255];  % #f5cfa8  橙
COLOR_LA  = [148/255, 206/255, 229/255];  % #94cee5  蓝
COLOR_LTD = [233/255, 162/255, 163/255];  % #e9a2a3  红

%% =========================================================
%  PART A — HA vs LA : Alpha_RawPower, Alpha_NormPower
%% =========================================================
% Note: RawPower is calculated by integrating the PSD over the band,
% so its physical unit is signal unit squared (e.g., \muV^2), not \muV^2/Hz.
metricsA     = {'Alpha_RawPower',         'Alpha_NormPower'};
titlesA      = {'Alpha Absolute Power',     'Relative Alpha Power'};
yUnitsA      = {'\muV^2',                 '%'};
groupCodesA  = [1, 3];
groupNamesA  = {'HA', 'LA'};
groupColorsA = [COLOR_HA; COLOR_LA];
figNamesA    = {'EEG_HA_LA_Alpha_AbsolutePower.png'};

fprintf('\n========== PART A: HA vs LA ==========\n');
[statsA, figsA] = runComparison( ...
    headers, dataRows, col2num, groupVals, ...
    metricsA, titlesA, yUnitsA, ...
    groupCodesA, groupNamesA, groupColorsA);

for k = 1:numel(figsA)
    outPath = fullfile(OUTPUT_DIR, figNamesA{k});
    exportgraphics(figsA(k), outPath, 'Resolution', 800);
    fprintf('Saved: %s\n', outPath);
end

%% =========================================================
%  PART B — HTD vs LTD : ThetaDelta_RawPower, ThetaDelta_NormPower
%% =========================================================
metricsB     = {'ThetaDelta_RawPower',        'ThetaDelta_NormPower'};
titlesB      = {'Theta + Delta Absolute Power', 'Relative Theta + Delta Power'};
yUnitsB      = {'\muV^2',                     '%'};
groupCodesB  = [2, 4];
groupNamesB  = {'HTD', 'LTD'};
groupColorsB = [COLOR_HTD; COLOR_LTD];
figNamesB    = {'EEG_HTD_LTD_ThetaDelta_AbsolutePower.png'};

fprintf('\n========== PART B: HTD vs LTD ==========\n');
[statsB, figsB] = runComparison( ...
    headers, dataRows, col2num, groupVals, ...
    metricsB, titlesB, yUnitsB, ...
    groupCodesB, groupNamesB, groupColorsB);

for k = 1:numel(figsB)
    outPath = fullfile(OUTPUT_DIR, figNamesB{k});
    exportgraphics(figsB(k), outPath, 'Resolution', 800);
    fprintf('Saved: %s\n', outPath);
end

%% =========================================================
%  Save statistics → 4 sheets in one Excel file
%% =========================================================
fprintf('\nWriting statistics...\n');
outputStats = fullfile(OUTPUT_DIR, OUTPUT_STATS);
writeStatsSheets(outputStats, statsA, metricsA, groupNamesA, 'HA_vs_LA');
writeStatsSheets(outputStats, statsB, metricsB, groupNamesB, 'HTD_vs_LTD');
fprintf('Stats saved: %s\n', outputStats);

fprintf('\nDone.\n');


%% =========================================================
%  FUNCTION  runComparison
%  -------------------------------------------------------
%  Inputs
%    headers      : 1×nCol cell array of column names
%    dataRows     : (nRow)×nCol cell array of data
%    col2num      : anonymous function col2num(colIdx) → numeric vector
%    groupValsAll : numeric vector of group codes for every data row
%    metricNames  : 1×nMetrics cell of column names to analyse
%    metricTitles : 1×nMetrics cell of descriptive names
%                   → ylabel  = metricTitle
%                   → title   = "metricTitle (unit)"
%    yUnits       : 1×nMetrics cell of unit strings, e.g. '\muV^2/Hz', '%'
%    groupCodes   : 1×2 numeric  [code_group1, code_group2]
%    groupNames   : 1×2 cell     {'Name1', 'Name2'}
%    groupColors  : 2×3 RGB matrix
%  Outputs
%    S    : stats struct
%    figs : figure handles for absolute-power metrics only
%% =========================================================
function [S, figs] = runComparison( ...
        headers, dataRows, col2num, groupValsAll, ...
        metricNames, metricTitles, yUnits, ...
        groupCodes, groupNames, groupColors)

    nMetrics = numel(metricNames);
    nGroups  = 2;

    % --- Locate metric columns (search in headers) ---
    metricCols = zeros(1, nMetrics);
    for mi = 1:nMetrics
        idx = find(strcmp(headers, metricNames{mi}));
        if isempty(idx)
            error('Column not found: "%s"', metricNames{mi});
        end
        metricCols(mi) = idx;
    end

    % --- Build numeric data matrix (all rows) ---
    nRows   = size(dataRows, 1);
    allData = zeros(nRows, nMetrics);
    for mi = 1:nMetrics
        allData(:, mi) = col2num(metricCols(mi));
    end

    % --- Filter: keep only the two target groups, drop NaN rows ---
    validMask = ismember(groupValsAll, groupCodes) & all(~isnan(allData), 2);
    gv = groupValsAll(validMask);
    ad = allData(validMask, :);

    fprintf('Valid records: %d\n', sum(validMask));
    groupIdx  = cell(1, nGroups);
    nPerGroup = zeros(1, nGroups);
    for g = 1:nGroups
        groupIdx{g}  = find(gv == groupCodes(g));
        nPerGroup(g) = numel(groupIdx{g});
        fprintf('  %s: n=%d\n', groupNames{g}, nPerGroup(g));
    end

    % --- Wilcoxon rank-sum tests ---
    MWp  = zeros(nMetrics, 1);
    MWz  = zeros(nMetrics, 1);
    Mdn  = zeros(nMetrics, nGroups);
    IQRv = zeros(nMetrics, nGroups);
    Mn   = zeros(nMetrics, nGroups);
    SDv  = zeros(nMetrics, nGroups);

    for mi = 1:nMetrics
        feat = ad(:, mi);
        for g = 1:nGroups
            gd         = feat(groupIdx{g});
            Mdn(mi,g)  = median(gd);
            IQRv(mi,g) = iqr(gd);
            Mn(mi,g)   = mean(gd);
            SDv(mi,g)  = std(gd);
        end
        [p_mw, ~, st] = ranksum(feat(groupIdx{1}), feat(groupIdx{2}));
        MWp(mi) = p_mw;
        MWz(mi) = ternary(isfield(st,'zval'), st.zval, NaN);
        fprintf('[%s]  z=%.3f, p=%.4f %s\n', ...
                metricNames{mi}, MWz(mi), p_mw, getSigStr(p_mw));
    end

    % --- Pack stats struct ---
    S.MWp       = MWp;
    S.MWz       = MWz;
    S.Median    = Mdn;
    S.IQR       = IQRv;
    S.Mean      = Mn;
    S.SD        = SDv;
    S.groupIdx  = groupIdx;
    S.nPerGroup = nPerGroup;
    S.allData   = ad;
    S.groupVals = gv;
    S.groupCodes= groupCodes;

    % --- Draw absolute-power figures only (skip relative/normalized power) ---
    plotMetricIdx = find(~contains(metricNames, '_NormPower'));
    figs = gobjects(1, numel(plotMetricIdx));

    for pi = 1:numel(plotMetricIdx)
        mi = plotMetricIdx(pi);
        feat = ad(:, mi);

        % ylabel = "Descriptive Name (unit)",  no title
        ylabelStr = sprintf('%s (%s)', metricTitles{mi}, yUnits{mi});

        fig = figure('Name', metricNames{mi}, 'Color', 'w', ...
                     'Position', [30 + (mi-1)*340, 50, 320, 420]);
        figs(pi) = fig;
        ax = axes('Parent', fig);
        hold(ax, 'on');

        for g = 1:nGroups
            gdata  = feat(groupIdx{g});
            q25    = prctile(gdata, 25);
            q50    = median(gdata);
            q75    = prctile(gdata, 75);
            iqrVal = q75 - q25;

            inFence  = gdata(gdata >= q25 - 1.5*iqrVal & gdata <= q75 + 1.5*iqrVal);
            wLow     = min(inFence);
            wHigh    = max(inFence);
            outliers = gdata(gdata < wLow | gdata > wHigh);

            bw = 0.35;
            c  = groupColors(g,:);

            % Box body
            fill(ax, [g-bw g+bw g+bw g-bw], [q25 q25 q75 q75], ...
                 c, 'FaceAlpha', 0.45, 'EdgeColor', c*0.75, 'LineWidth', 1.5);
            % Median line
            plot(ax, [g-bw g+bw], [q50 q50], '-', 'Color', c*0.65, 'LineWidth', 2.5);
            % Whisker lines
            plot(ax, [g g], [wLow  q25],  '-', 'Color', c*0.65, 'LineWidth', 1.2);
            plot(ax, [g g], [q75  wHigh], '-', 'Color', c*0.65, 'LineWidth', 1.2);
            % Whisker caps
            plot(ax, [g-bw/2 g+bw/2], [wLow  wLow],  '-', 'Color', c*0.65, 'LineWidth', 1.2);
            plot(ax, [g-bw/2 g+bw/2], [wHigh wHigh], '-', 'Color', c*0.65, 'LineWidth', 1.2);
            % Outliers
            if ~isempty(outliers)
                scatter(ax, repmat(g, numel(outliers), 1), outliers, 22, ...
                        c*0.75, 'o', 'LineWidth', 1.0, 'MarkerEdgeAlpha', 0.7);
            end
            % Jittered individual points
            jitter = (rand(numel(gdata), 1) - 0.5) * 0.28;
            scatter(ax, g + jitter, gdata, 20, c*0.8, 'filled', 'MarkerFaceAlpha', 0.55);
        end

        % Significance bracket
        yMax   = max(feat);
        yMin   = min(feat);
        yRange = yMax - yMin;
        yTop   = yMax + 0.08 * yRange;
        sigStr = getSigStr(MWp(mi));

        if ~isempty(sigStr)
            plot(ax, [1, 2], [yTop yTop], 'k-', 'LineWidth', 1.0);
            text(ax, 1.5, yTop + 0.03*yRange, sigStr, ...
                 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        end

        % Axis formatting
        nSig = double(~isempty(sigStr));
        yPad = 0.08*yRange + nSig * 0.16*yRange;
        ylim(ax, [yMin - 0.08*yRange, yMax + yPad]);
        xlim(ax, [0.5, 2.5]);
        set(ax, 'XTick', 1:2, 'XTickLabel', groupNames, ...
            'FontSize', 10, 'Box', 'off');   % Box off: removes all 4 borders
        ax.XAxis.Visible = 'on';
        ax.YAxis.Visible = 'on';
        % Keep only left and bottom borders
        set(ax, 'TickDir', 'out');
        xtickangle(ax, 0);
        ylabel(ax, ylabelStr, 'FontSize', 10, 'FontWeight', 'bold');
        % No title
        grid(ax, 'off');   % Remove horizontal and vertical grid lines
        hold(ax, 'off');
    end
end


%% =========================================================
%  FUNCTION  writeStatsSheets
%  Appends two sheets (<prefix>_Summary, <prefix>_RawData)
%  to the output Excel file.
%% =========================================================
function writeStatsSheets(outputStats, S, metricNames, groupNames, sheetPrefix)

    nMetrics = numel(metricNames);
    g1 = groupNames{1};
    g2 = groupNames{2};

    % ---- Summary sheet ----
    headerRow = {'Metric', ...
        [g1 '_Mean'], [g1 '_SD'], [g1 '_Median'], [g1 '_IQR'], ...
        [g2 '_Mean'], [g2 '_SD'], [g2 '_Median'], [g2 '_IQR'], ...
        'z', 'p', 'Significance'};

    dataRows = cell(nMetrics, numel(headerRow));
    for mi = 1:nMetrics
        dataRows(mi,:) = { ...
            metricNames{mi}, ...
            round(S.Mean(mi,1),4),   round(S.SD(mi,1),4), ...
            round(S.Median(mi,1),4), round(S.IQR(mi,1),4), ...
            round(S.Mean(mi,2),4),   round(S.SD(mi,2),4), ...
            round(S.Median(mi,2),4), round(S.IQR(mi,2),4), ...
            round(S.MWz(mi),4),      round(S.MWp(mi),6), ...
            getSigStr(S.MWp(mi)) };
    end
    xlswrite(outputStats, [headerRow; dataRows], [sheetPrefix '_Summary']);

    % ---- RawData sheet ----
    rawHeader = [{'Group'}, metricNames];
    nValid    = size(S.allData, 1);
    rawRows   = cell(nValid, numel(rawHeader));

    for i = 1:nValid
        lbl = '';
        for g = 1:2
            if S.groupVals(i) == S.groupCodes(g)
                lbl = groupNames{g};
                break;
            end
        end
        rawRows{i,1} = lbl;
        for mi = 1:nMetrics
            rawRows{i, mi+1} = S.allData(i, mi);
        end
    end
    xlswrite(outputStats, [rawHeader; rawRows], [sheetPrefix '_RawData']);
end


%% =========================================================
%  Local helper functions
%% =========================================================
function out = ternary(cond, valTrue, valFalse)
    if cond, out = valTrue; else, out = valFalse; end
end

function sig = getSigStr(p)
    if     p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,             sig = '';
    end
end
