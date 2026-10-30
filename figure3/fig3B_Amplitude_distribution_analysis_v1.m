%% EEG Broadband Amplitude Distribution Analysis
% Visualising the distribution of Broadband_Amplitude_Mean across visual groups
% Group: HA=1, HTD=2, LA=3, LTD=4
%
% Modifications from Freq_consistency_analysis_v4.m:
%   - Target variable changed from Mean_Freq_Hz to Broadband_Amplitude_Mean
%   - Removed all consistency / Cohen's-kappa analysis (no visual reference)
%   - Removed 8 Hz boundary line
%   - Each panel shows the full distribution per group (no consistent/inconsistent split)
%   - Panel 1: HA & LA  (Alpha groups);  Panel 2: HTD & LTD  (Non-Alpha groups)

clear; clc; close all;

%% ========== Path Settings (modify as needed) ==========
input_file  = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
output_dir  = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅\波幅分类的一致性';
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
ampCol   = 'Broadband_Amplitude_Mean';

Group = T.(groupCol);
Amp   = T.(ampCol);

fprintf('Total patients: %d\n', height(T));

%% 2. Exclude patients with missing Broadband_Amplitude_Mean
validIdx = ~isnan(Amp);
fprintf('Excluded (missing Broadband_Amplitude_Mean): %d\n', sum(~validIdx));
fprintf('Included in analysis: %d\n', sum(validIdx));

Group_v = Group(validIdx);
Amp_v   = Amp(validIdx);

%% 3. Per-group basic statistics
groupLabels = {'HA', 'HTD', 'LA', 'LTD'};
groupIDs    = [1, 2, 3, 4];

fprintf('\n========== Per-group Broadband Amplitude Statistics ==========\n');
fprintf('%-10s %6s %10s %10s %10s %10s\n', ...
    'Group', 'N', 'Mean', 'Median', 'SD', 'IQR');

grp_n      = zeros(4,1);
grp_mean   = zeros(4,1);
grp_median = zeros(4,1);
grp_sd     = zeros(4,1);
grp_iqr    = zeros(4,1);

for i = 1:4
    idx          = (Group_v == groupIDs(i));
    amp_grp      = Amp_v(idx);
    grp_n(i)     = sum(idx);
    grp_mean(i)  = mean(amp_grp);
    grp_median(i)= median(amp_grp);
    grp_sd(i)    = std(amp_grp);
    grp_iqr(i)   = iqr(amp_grp);
    fprintf('%-10s %6d %10.3f %10.3f %10.3f %10.3f\n', ...
        groupLabels{i}, grp_n(i), grp_mean(i), grp_median(i), grp_sd(i), grp_iqr(i));
end

%% 4. Distribution figures (histogram + KDE), two panels
% Palette: #b6d093(green/HA) #f5cfa8(orange/HTD) #94cee5(blue/LA) #e9a2a3(red/LTD)

palette_hex = {'#b6d093', '#f5cfa8', '#94cee5', '#e9a2a3'};
hex2rgb     = @(h) [hex2dec(h(2:3)), hex2dec(h(4:5)), hex2dec(h(6:7))] / 255;
palette     = cellfun(hex2rgb, palette_hex, 'UniformOutput', false);

ampMin   = floor(min(Amp_v)) - 1;
ampMax   = ceil(max(Amp_v))  + 1;
nBins    = 20;
binEdges = linspace(ampMin, ampMax, nBins+1);
binWidth = binEdges(2) - binEdges(1);
xi       = linspace(ampMin, ampMax, 500);

figNames    = {'EEG_amplitude_density_HA_LA.png', 'EEG_amplitude_density_HTD_LTD.png'};
panelGroups = {[1, 3]; [2, 4]};   % Panel 1: HA+LA (Alpha);  Panel 2: HTD+LTD (Non-Alpha)

for p = 1:2
    fig = figure('Position', [100, 100, 1000, 500], 'Color', 'white');
    ax  = axes(fig);
    hold on;

    legend_handles = [];
    legend_labels  = {};

    for ki = 1:numel(panelGroups{p})
        i    = panelGroups{p}(ki);
        clr  = palette{i};

        idx_grp = (Group_v == groupIDs(i));
        amp_grp = Amp_v(idx_grp);
        n_grp   = length(amp_grp);

        clrEdge = clr * 0.55;

        % ---- Histogram (normalised to probability density) ----
        counts = histcounts(amp_grp, binEdges) / (n_grp * binWidth);
        for b = 1:nBins
            if counts(b) > 0
                xb = [binEdges(b), binEdges(b+1), binEdges(b+1), binEdges(b)];
                yb = [0, 0, counts(b), counts(b)];
                patch(xb, yb, clr, ...
                    'FaceAlpha', 0.75, 'EdgeColor', clrEdge, 'LineWidth', 1.2);
            end
        end

        % ---- KDE curve ----
        if n_grp >= 2
            [f_all, ~] = ksdensity(amp_grp, xi);
            plot(xi, f_all, '-', 'Color', clr * 0.65, 'LineWidth', 2.2);
        end

        % ---- Legend proxy patch ----
        h = patch(nan, nan, clr, 'FaceAlpha', 0.75, 'EdgeColor', clrEdge, 'LineWidth', 1.2);
        legend_handles = [legend_handles, h]; %#ok<AGROW>
        legend_labels  = [legend_labels, ...
            {sprintf('%s  (n=%d, mean=%.2f)', groupLabels{i}, n_grp, grp_mean(i))}]; %#ok<AGROW>
    end

    % ---- Axes formatting ----
    set(ax, 'Color', [0.95 0.96 0.98], 'FontSize', 10, 'Box', 'off', ...
        'GridColor', [1 1 1], 'GridAlpha', 1, 'XGrid', 'on', 'YGrid', 'on', ...
        'TickDir', 'out');

    xlim([ampMin, ampMax]);
    xlabel('Broadband envelope amplitude (μV)', 'FontSize', 10);
    ylabel('Probability Density',       'FontSize', 10);

    gi1 = panelGroups{p}(1);
    gi2 = panelGroups{p}(2);
    title(sprintf('%s  |  %s', groupLabels{gi1}, groupLabels{gi2}), ...
        'FontSize', 11, 'FontWeight', 'bold');

    % ---- Fix bottom margin ----
    drawnow;
    orig_pos = ax.Position;
    bottom_shift = 0.06;
    ax.Position = [orig_pos(1), orig_pos(2) + bottom_shift, ...
                   orig_pos(3), orig_pos(4) - bottom_shift];

    % ---- Legend ----
    lgd = legend(legend_handles, legend_labels, ...
        'Location', 'northeast', 'FontSize', 9, 'Box', 'off');

    drawnow;
    ax_pos     = ax.Position;
    ti         = ax.TightInset;
    plot_right = ax_pos(1) + ax_pos(3) + ti(3);
    lgd.Units  = 'normalized';
    lgd.Position(1) = plot_right - lgd.Position(3) - 0.002;

    hold off;

    fig_path = fullfile(output_dir, figNames{p});
    exportgraphics(fig, fig_path, 'Resolution', 300);
    fprintf('\nFigure saved: %s\n', fig_path);
end

%% 5. Save statistics to text file
result_path = fullfile(output_dir, 'amplitude_distribution_results.txt');
fid = fopen(result_path, 'w', 'n', 'UTF-8');
fprintf(fid, '========== EEG Broadband Amplitude Distribution ==========\r\n\r\n');
fprintf(fid, 'Analysis date: %s\r\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Input file: %s\r\n\r\n', input_file);
fprintf(fid, 'Total patients: %d\r\n', height(T));
fprintf(fid, 'Excluded (missing Broadband_Amplitude_Mean): %d\r\n', sum(~validIdx));
fprintf(fid, 'Included in analysis: %d\r\n\r\n', sum(validIdx));

fprintf(fid, '%-10s %6s %10s %10s %10s %10s\r\n', ...
    'Group', 'N', 'Mean', 'Median', 'SD', 'IQR');
for i = 1:4
    fprintf(fid, '%-10s %6d %10.3f %10.3f %10.3f %10.3f\r\n', ...
        groupLabels{i}, grp_n(i), grp_mean(i), grp_median(i), grp_sd(i), grp_iqr(i));
end

fclose(fid);
fprintf('Results saved: %s\n', result_path);

fprintf('\nAnalysis complete.\n');