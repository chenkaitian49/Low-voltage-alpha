%% EEG Visual Reading vs Algorithm Classification Consistency Analysis
% Comparing visual classification (Group) with algorithm classification (Mean_Freq_Hz)
% Group: HA=1, HTD=2, LA=3, LTD=4
%   Alpha rhythm: HA(1), LA(3); Non-Alpha rhythm: HTD(2), LTD(4)
% Mean_Freq_Hz: >=8 -> Alpha; <8 -> Non-Alpha; NaN -> excluded
%
% v4 modifications:
%   - Removed Section 9 (ROC curve + Youden Index + Bootstrap CI)
%   - Distribution figures: HA & HTD plotted in one panel; LA & LTD in another
%   - Group colours unchanged: #b6d093 HA, #f5cfa8 HTD, #94cee5 LA, #e9a2a3 LTD
%   - Section 7: detailed list of inconsistent patients (row, ID, group, freq, mismatch)
%   - Histogram bars: uniform FaceAlpha=0.75 with deep-colour edge borders

clear; clc; close all;

%% ========== Path Settings (modify as needed) ==========
input_file  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
output_dir  = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅\频率分类的一致性';
% =======================================================

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% 1. Load data
fprintf('Reading data: %s\n', input_file);
opts = detectImportOptions(input_file);
T = readtable(input_file, opts);

colNames = T.Properties.VariableNames;
groupCol = colNames{contains(colNames, 'Group', 'IgnoreCase', true)};
freqCol  = 'Mean_Freq_Hz';

Group        = T.(groupCol);
Mean_Freq_Hz = T.(freqCol);

fprintf('Total patients: %d\n', height(T));

%% 2. Exclude patients with missing Mean_Freq_Hz
validIdx = ~isnan(Mean_Freq_Hz);
fprintf('Excluded (missing Mean_Freq_Hz): %d\n', sum(~validIdx));
fprintf('Included in analysis: %d\n', sum(validIdx));

Group_v  = Group(validIdx);
FreqHz_v = Mean_Freq_Hz(validIdx);

%% 3. Binary classification labels
visual_alpha = (Group_v == 1) | (Group_v == 3);
algo_alpha   = FreqHz_v >= 8;

%% 4. Overall consistency rate
n_total     = length(Group_v);
n_agree     = sum(visual_alpha == algo_alpha);
overall_acc = n_agree / n_total * 100;

fprintf('\n========== Consistency Analysis Results ==========\n');
fprintf('Overall consistency: %.2f%% (%d / %d)\n', overall_acc, n_agree, n_total);

%% 5. Cohen's kappa
TP = sum( visual_alpha &  algo_alpha);
FN = sum( visual_alpha & ~algo_alpha);
FP = sum(~visual_alpha &  algo_alpha);
TN = sum(~visual_alpha & ~algo_alpha);
confMat = [TP, FN; FP, TN];

Po = (TP + TN) / n_total;
P_vis_alpha  = (TP + FN) / n_total;
P_vis_non    = (FP + TN) / n_total;
P_algo_alpha = (TP + FP) / n_total;
P_algo_non   = (FN + TN) / n_total;
Pe    = P_vis_alpha * P_algo_alpha + P_vis_non * P_algo_non;
kappa = (Po - Pe) / (1 - Pe);

fprintf('Confusion matrix (rows=visual, cols=algorithm, 1=Alpha, 0=NonAlpha):\n');
fprintf('              Algo:Alpha  Algo:NonAlpha\n');
fprintf('Visual:Alpha     %3d           %3d\n', TP, FN);
fprintf('Visual:NonAlpha  %3d           %3d\n', FP, TN);
fprintf('\nCohen''s kappa = %.4f\n', kappa);

if kappa < 0,        kappa_interp = 'Poor';
elseif kappa < 0.20, kappa_interp = 'Slight';
elseif kappa < 0.40, kappa_interp = 'Fair';
elseif kappa < 0.60, kappa_interp = 'Moderate';
elseif kappa < 0.80, kappa_interp = 'Substantial';
else,                kappa_interp = 'Almost Perfect';
end
fprintf('Kappa interpretation: %s\n', kappa_interp);

%% 6. Per-group statistics
groupLabels = {'HA', 'HTD', 'LA', 'LTD'};
groupIDs    = [1, 2, 3, 4];

fprintf('\n========== Per-group Consistency ==========\n');
fprintf('%-10s %6s %6s %8s %10s\n', 'Group','N','Agree','Disagree','Rate(%)');

grp_n        = zeros(4,1);
grp_agree    = zeros(4,1);
grp_disagree = zeros(4,1);
grp_acc      = zeros(4,1);

for i = 1:4
    idx   = Group_v == groupIDs(i);
    g_n   = sum(idx);
    g_agr = sum(visual_alpha(idx) == algo_alpha(idx));
    g_dis = g_n - g_agr;
    g_acc = g_agr / g_n * 100;
    grp_n(i) = g_n;  grp_agree(i) = g_agr;
    grp_disagree(i) = g_dis;  grp_acc(i) = g_acc;
    fprintf('%-10s %6d %6d %8d %9.1f%%\n', ...
        groupLabels{i}, g_n, g_agr, g_dis, g_acc);
end

%% =========================================================
%% 7. Inconsistent patients — detailed list
%% =========================================================
disagree_mask = (visual_alpha ~= algo_alpha);
T_valid       = T(validIdx, :);   % rows matching Group_v / FreqHz_v

% Patient name column (Name = English pinyin)
idCol = 'Name';
if ~any(strcmp(colNames, idCol))
    idCol = '';   % fall back to row number only if column absent
end

% Pre-compute mapping: valid-row index -> original row number in T
validOrigRows = find(validIdx);   % validOrigRows(r) = original row of r-th valid patient

dis_rows    = find(disagree_mask);
n_dis_total = length(dis_rows);

fprintf('\n========== Inconsistent Patients ==========\n');
if isempty(idCol)
    fprintf('%-6s %-8s %-10s %-12s %-22s\n', ...
        'Row','Group','Freq(Hz)','Visual','Algorithm');
else
    fprintf('%-6s %-20s %-8s %-10s %-12s %-22s\n', ...
        'Row','Name','Group','Freq(Hz)','Visual','Algorithm');
end

for k = 1:n_dis_total
    r       = dis_rows(k);
    grpID   = Group_v(r);
    freq    = FreqHz_v(r);
    origRow = validOrigRows(r);   % row number in original table T
    vis_str  = ternary_str(visual_alpha(r), 'Alpha',    'Non-Alpha');
    algo_str = ternary_str(algo_alpha(r),   'Alpha',    'Non-Alpha');

    if isempty(idCol)
        fprintf('%-6d %-8s %-10.3f %-12s %-22s\n', ...
            origRow, groupLabels{grpID}, freq, vis_str, algo_str);
    else
        pid = T_valid.(idCol)(r);
        if isnumeric(pid), pid_str = num2str(pid);
        else,              pid_str = char(string(pid)); end
        fprintf('%-6d %-20s %-8s %-10.3f %-12s %-22s\n', ...
            origRow, pid_str, groupLabels{grpID}, freq, vis_str, algo_str);
    end
end
fprintf('Total inconsistent: %d\n', n_dis_total);

%% =========================================================
%% 8. Dominant frequency distribution (histogram + KDE)
%%    consistent vs inconsistent patients, per visual group
%% =========================================================
% Palette: #b6d093(green) #f5cfa8(orange) #94cee5(blue) #e9a2a3(red)

palette_hex = {'#b6d093', '#f5cfa8', '#94cee5', '#e9a2a3'};
hex2rgb     = @(h) [hex2dec(h(2:3)), hex2dec(h(4:5)), hex2dec(h(6:7))] / 255;
palette     = cellfun(hex2rgb, palette_hex, 'UniformOutput', false);

% Frequency axis: integer steps
freqMin  = floor(min(FreqHz_v)) - 1;
freqMax  = ceil(max(FreqHz_v))  + 1;
nBins    = 20;
binEdges = linspace(freqMin, freqMax, nBins+1);
binWidth = binEdges(2) - binEdges(1);
xi       = linspace(freqMin, freqMax, 500);

figNames    = {'EEG_freq_density_HA_HTD.png', 'EEG_freq_density_LA_LTD.png'};
panelGroups = {[1, 2]; [3, 4]};   % Panel 1: HA+HTD;  Panel 2: LA+LTD

for p = 1:2
    fig = figure('Position', [100, 100, 1000, 500], 'Color', 'white');
    ax = axes(fig);
    hold on;

    legend_handles = [];
    legend_labels  = {};

    for ki = 1:numel(panelGroups{p})
        i        = panelGroups{p}(ki);   % group index (1-4) same as palette index
        clr      = palette{i};
        clrLight = clr * 0.5 + [1 1 1] * 0.5;

        idx_grp    = (Group_v == groupIDs(i));
        freq_grp   = FreqHz_v(idx_grp);
        agree_mask = (visual_alpha(idx_grp) == algo_alpha(idx_grp));

        freq_agree    = freq_grp( agree_mask);
        freq_disagree = freq_grp(~agree_mask);
        n_grp = length(freq_grp);

        clrEdge      = clr      * 0.55;   % dark edge for consistent bars
        clrEdgeLight = clrLight * 0.55;   % dark edge for inconsistent bars

        % ---- Histograms (normalised to probability density) ----
        if ~isempty(freq_disagree)
            counts_dis = histcounts(freq_disagree, binEdges) / (n_grp * binWidth);
            for b = 1:nBins
                if counts_dis(b) > 0
                    xb = [binEdges(b), binEdges(b+1), binEdges(b+1), binEdges(b)];
                    yb = [0, 0, counts_dis(b), counts_dis(b)];
                    patch(xb, yb, clrLight, ...
                        'FaceAlpha', 0.75, 'EdgeColor', clrEdgeLight, 'LineWidth', 1.2);
                end
            end
        end

        if ~isempty(freq_agree)
            counts_agr = histcounts(freq_agree, binEdges) / (n_grp * binWidth);
            for b = 1:nBins
                if counts_agr(b) > 0
                    xb = [binEdges(b), binEdges(b+1), binEdges(b+1), binEdges(b)];
                    yb = [0, 0, counts_agr(b), counts_agr(b)];
                    patch(xb, yb, clr, ...
                        'FaceAlpha', 0.75, 'EdgeColor', clrEdge, 'LineWidth', 1.2);
                end
            end
        end

        % ---- Overall KDE curve ----
        if n_grp >= 2
            [f_all, ~] = ksdensity(freq_grp, xi);
            plot(xi, f_all, '-', 'Color', clr * 0.65, 'LineWidth', 2.2);
        end

        % ---- Legend proxy patches ----
        h_agr = patch(nan, nan, clr,      'FaceAlpha', 0.75, 'EdgeColor', clrEdge,      'LineWidth', 1.2);
        h_dis = patch(nan, nan, clrLight, 'FaceAlpha', 0.75, 'EdgeColor', clrEdgeLight, 'LineWidth', 1.2);
        legend_handles = [legend_handles, h_agr, h_dis]; %#ok<AGROW>
        legend_labels  = [legend_labels, ...
            {sprintf('%s Consistent (n=%d)',   groupLabels{i}, sum( agree_mask))}, ...
            {sprintf('%s Inconsistent (n=%d)', groupLabels{i}, sum(~agree_mask))}]; %#ok<AGROW>
    end

    % ---- 8 Hz boundary line ----
    xline(8, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.3);

    % ---- Axes formatting ----
    set(ax, 'Color', [0.95 0.96 0.98], 'FontSize', 10, 'Box', 'off', ...
        'GridColor', [1 1 1], 'GridAlpha', 1, 'XGrid', 'on', 'YGrid', 'on', ...
        'XTick', freqMin:1:freqMax, ...
        'TickDir', 'out');

    xlim([freqMin, freqMax]);
    xlabel('Dominant Frequency (Hz)', 'FontSize', 10);
    ylabel('Probability Density', 'FontSize', 10);

    % Title: consistency rates only
    gi1 = panelGroups{p}(1);
    gi2 = panelGroups{p}(2);
    title(sprintf('%s: %.1f%%  |  %s: %.1f%%', ...
        groupLabels{gi1}, grp_acc(gi1), ...
        groupLabels{gi2}, grp_acc(gi2)), ...
        'FontSize', 11, 'FontWeight', 'bold');

    % ---- Fix bottom margin ----
    drawnow;
    orig_pos = ax.Position;
    bottom_shift = 0.06;
    ax.Position = [orig_pos(1), orig_pos(2) + bottom_shift, ...
                   orig_pos(3), orig_pos(4) - bottom_shift];

    % ---- Legend inside axes (northeast) ----
    lgd = legend(legend_handles, legend_labels, ...
        'Location', 'northeast', 'FontSize', 9, 'Box', 'off');

    % ---- Align legend right edge to axes right ----
    drawnow;
    ax_pos = ax.Position;
    ti     = ax.TightInset;
    plot_right = ax_pos(1) + ax_pos(3) + ti(3);
    lgd.Units  = 'normalized';
    lgd.Position(1) = plot_right - lgd.Position(3) - 0.002;

    hold off;

    % ---- Save individual figure ----
    fig_path = fullfile(output_dir, figNames{p});
    exportgraphics(fig, fig_path, 'Resolution', 300);
    fprintf('\nFigure saved: %s\n', fig_path);
end

%% 8. Save results to text file
result_path = fullfile(output_dir, 'consistency_results.txt');
fid = fopen(result_path, 'w', 'n', 'UTF-8');
fprintf(fid, '========== EEG Alpha Rhythm Classification Consistency ==========\r\n\r\n');
fprintf(fid, 'Analysis date: %s\r\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Input file: %s\r\n\r\n', input_file);
fprintf(fid, 'Total patients: %d\r\n', height(T));
fprintf(fid, 'Excluded (missing Mean_Freq_Hz): %d\r\n', sum(~validIdx));
fprintf(fid, 'Included in analysis: %d\r\n\r\n', n_total);

fprintf(fid, 'Classification rules:\r\n');
fprintf(fid, '  Visual: HA(1)/LA(3) = Alpha; HTD(2)/LTD(4) = Non-Alpha\r\n');
fprintf(fid, '  Algorithm: Mean_Freq_Hz >= 8 = Alpha; < 8 = Non-Alpha\r\n\r\n');

fprintf(fid, 'Confusion matrix (rows=visual, cols=algorithm):\r\n');
fprintf(fid, '                  Algo:Alpha  Algo:NonAlpha\r\n');
fprintf(fid, 'Visual:Alpha          %3d           %3d\r\n', TP, FN);
fprintf(fid, 'Visual:NonAlpha       %3d           %3d\r\n\r\n', FP, TN);

fprintf(fid, 'Overall consistency: %.2f%% (%d / %d)\r\n', overall_acc, n_agree, n_total);
fprintf(fid, 'Cohen''s kappa = %.4f  [%s]\r\n\r\n', kappa, kappa_interp);

fprintf(fid, 'Per-group details:\r\n');
fprintf(fid, '%-12s %5s %8s %10s %12s\r\n', 'Group','N','Agree','Disagree','Rate(%)');
for i = 1:4
    fprintf(fid, '%-12s %5d %8d %10d %11.1f%%\r\n', ...
        groupLabels{i}, grp_n(i), grp_agree(i), grp_disagree(i), grp_acc(i));
end

fprintf(fid, '\r\n========== Inconsistent Patients ==========\r\n');
if isempty(idCol)
    fprintf(fid, '%-6s %-8s %-10s %-12s %-22s\r\n', ...
        'Row','Group','Freq(Hz)','Visual','Algorithm');
else
    fprintf(fid, '%-6s %-20s %-8s %-10s %-12s %-22s\r\n', ...
        'Row','Name','Group','Freq(Hz)','Visual','Algorithm');
end
for k = 1:n_dis_total
    r       = dis_rows(k);
    grpID   = Group_v(r);
    freq    = FreqHz_v(r);
    origRow = validOrigRows(r);
    vis_str  = ternary_str(visual_alpha(r), 'Alpha',    'Non-Alpha');
    algo_str = ternary_str(algo_alpha(r),   'Alpha',    'Non-Alpha');
    if isempty(idCol)
        fprintf(fid, '%-6d %-8s %-10.3f %-12s %-22s\r\n', ...
            origRow, groupLabels{grpID}, freq, vis_str, algo_str);
    else
        pid = T_valid.(idCol)(r);
        if isnumeric(pid), pid_str = num2str(pid);
        else,              pid_str = char(string(pid)); end
        fprintf(fid, '%-6d %-20s %-8s %-10.3f %-12s %-22s\r\n', ...
            origRow, pid_str, groupLabels{grpID}, freq, vis_str, algo_str);
    end
end
fprintf(fid, 'Total inconsistent: %d\r\n', n_dis_total);

fclose(fid);
fprintf('Results saved: %s\n', result_path);

fprintf('\nAnalysis complete.\n');

%% =========================================================
%% Helper function
%% =========================================================
function s = ternary_str(cond, a, b)
% Return string a if cond is true, else b
    if cond, s = a; else, s = b; end
end