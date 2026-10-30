%% FOOOF_TwoPatients_Plot.m
% 功能：指定2个患者的 .set 文件，固定 min_peak_height=0.2
% 在同一坐标轴绘制两组原始功率谱 + 拟合曲线，并标注检测到的峰值
% 患者1：HA Group（高波幅 alpha）
% 患者2：LTD Group（低波幅 ThetaDelta）

clear; clc; close all;

%% ==================== 设置 ====================

% --- 患者1：HA Group ---
file_path_1 = 'F:\EEG\LowVoltageResearch\PreprocessedData\highAlpha_ICA\SunEnQiang_pre_rejected_processed.set';
group_label_1 = 'HA Group';

% --- 患者2：LTD Group ---
file_path_2 = 'F:\EEG\LowVoltageResearch\PreprocessedData\lowThetaDelta_ICA\GaoXingLong_pre_rejected_processed.set';
group_label_2 = 'LTD Group';

% 输出文件夹
output_folder = 'F:\EEG\LowVoltageResearch\OutputResult\V2\频率+波幅示意图\频率示意图';

% 要分析的通道列表（多个通道取平均）
target_channels = {'P3', 'P4', 'O1', 'O2', 'T5', 'T6'};

% 固定 min_peak_height = 0.2
min_peak_height = 0.2;

% ==================== 色系定义 ====================
% Original Spectrum  → 患者1 绿色 #b6d093 / 患者2 粉红色 #e9a2a3
% Aperiodic Fit      → 两组统一  浅黄 #f5cfa8
% Full Model Fit     → 两组统一  浅蓝 #94cee5
% Selected Peak      → 两组统一  浅紫 #c8b4d8
% Unselected Peak    → 两组统一  浅灰 #c8c8c8

COLOR_APERIODIC_SHARED = [245, 207, 168] / 255;   % #f5cfa8  浅黄（两组共用）
COLOR_FIT_SHARED       = [148, 206, 229] / 255;   % #94cee5  浅蓝（两组共用）
COLOR_SELECTED_SHARED  = [200, 180, 216] / 255;   % #c8b4d8  浅紫（两组共用）
COLOR_UNSEL_SHARED     = [200, 200, 200] / 255;   % #c8c8c8  浅灰（两组共用）

% ---- 患者1（HA Group） ----
COLOR_RAW_1       = [182, 208, 147] / 255;   % #b6d093  绿色
COLOR_APERIODIC_1 = COLOR_APERIODIC_SHARED;
COLOR_FIT_1       = COLOR_FIT_SHARED;
COLOR_SELECTED_1  = COLOR_SELECTED_SHARED;
COLOR_UNSEL_1     = COLOR_UNSEL_SHARED;

% ---- 患者2（LTD Group） ----
COLOR_RAW_2       = [233, 162, 163] / 255;   % #e9a2a3  粉红色
COLOR_APERIODIC_2 = COLOR_APERIODIC_SHARED;
COLOR_FIT_2       = COLOR_FIT_SHARED;
COLOR_SELECTED_2  = COLOR_SELECTED_SHARED;
COLOR_UNSEL_2     = COLOR_UNSEL_SHARED;

% 添加工具箱
addpath('F:\matlab_downloadtoolbox\EEGToolbox\eeglab2023.0');
eeglab nogui;

%% ==================== 确保输出文件夹存在 ====================
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% ==================== 打包两组参数，循环处理 ====================
file_paths    = {file_path_1,    file_path_2};
group_labels  = {group_label_1,  group_label_2};

COLOR_RAW      = {COLOR_RAW_1,      COLOR_RAW_2};
COLOR_APERIODIC= {COLOR_APERIODIC_1, COLOR_APERIODIC_2};
COLOR_FIT      = {COLOR_FIT_1,      COLOR_FIT_2};
COLOR_SELECTED = {COLOR_SELECTED_1,  COLOR_SELECTED_2};
COLOR_UNSEL    = {COLOR_UNSEL_1,     COLOR_UNSEL_2};

% 存储每组结果（用于绘图）
results_all     = cell(1, 2);
filenames_all   = cell(1, 2);
valid_group     = false(1, 2);

for g = 1:2
    filepath = file_paths{g};
    [input_folder, filename_base, ~] = fileparts(filepath);
    datafile = [filename_base, '.set'];

    fprintf('处理 [%s]: %s\n', group_labels{g}, datafile);
    fprintf('%s\n', repmat('-', 1, 50));

    try
        %% 加载数据
        EEG = pop_loadset('filename', datafile, 'filepath', input_folder);
        fs  = EEG.srate;

        if ndims(EEG.data) == 3
            eeg_data = reshape(EEG.data, size(EEG.data,1), []);
        else
            eeg_data = EEG.data;
        end
        eeg_data = double(eeg_data);

        fprintf('  数据: %d 通道, %.1f 秒\n', size(eeg_data,1), size(eeg_data,2)/fs);

        %% 查找目标通道
        chan_indices    = [];
        chan_names_found = {};

        if ~isempty(EEG.chanlocs)
            all_labels = {EEG.chanlocs.labels};
            for i = 1:length(target_channels)
                idx = find(strcmpi(all_labels, target_channels{i}));
                if ~isempty(idx)
                    chan_indices(end+1)   = idx(1);
                    chan_names_found{end+1} = target_channels{i};
                end
            end
        end

        if isempty(chan_indices)
            chan_indices = 1;
            if ~isempty(EEG.chanlocs)
                chan_names_found = {EEG.chanlocs(1).labels};
            else
                chan_names_found = {'Ch1'};
            end
            fprintf('  警告: 未找到指定通道，使用 %s\n', chan_names_found{1});
        end

        fprintf('  使用通道: %s\n', strjoin(chan_names_found, ', '));

        %% 计算多通道平均功率谱
        nfft    = 2^nextpow2(fs * 2);
        psd_all = [];

        for i = 1:length(chan_indices)
            [psd_temp, freqs] = pwelch(eeg_data(chan_indices(i),:)', ...
                hanning(nfft), nfft/2, nfft, fs);
            psd_all(:, i) = psd_temp;
        end

        psd      = mean(psd_all, 2);
        freq_idx = freqs >= 1 & freqs <= 40;
        freqs_plot = freqs(freq_idx);
        psd_plot   = psd(freq_idx);

        %% FOOOF 拟合
        settings = struct();
        settings.freq_range        = [1, 40];
        settings.peak_width_limits = [0.5, 12];
        settings.max_n_peaks       = 10;
        settings.min_peak_height   = min_peak_height;
        settings.aperiodic_mode    = 'fixed';

        fooof_result = fooof_fit_simple(freqs_plot, psd_plot, settings);

        n_peaks = size(fooof_result.peak_params, 1);
        fprintf('  检测到峰值数: %d\n', n_peaks);
        if n_peaks > 0
            for pk = 1:n_peaks
                fprintf('    峰值 %d: CF=%.2f Hz, Power=%.3f, BW=%.3f\n', ...
                    pk, fooof_result.peak_params(pk,1), ...
                    fooof_result.peak_params(pk,2), ...
                    fooof_result.peak_params(pk,3));
            end
        end

        %% 确定"选中峰"：1-13Hz 中中心频率最大的峰
        selected_pk = NaN;
        if n_peaks > 0
            cfs      = fooof_result.peak_params(:, 1);
            in_range = cfs >= 1 & cfs <= 13;
            if any(in_range)
                cfs_in_range = cfs(in_range);
                [~, max_rel] = max(cfs_in_range);
                cfs_all_idx  = find(in_range);
                selected_pk  = cfs_all_idx(max_rel);
            end
        end

        results_all{g}   = fooof_result;
        results_all{g}.selected_pk = selected_pk;
        filenames_all{g} = filename_base;
        valid_group(g)   = true;

    catch ME
        fprintf('  错误: %s\n', ME.message);
        fprintf('  位置: %s (行 %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

%% ==================== 绘图（同一坐标轴） ====================
fig = figure('Position', [100, 100, 900, 520], 'Visible', 'off', 'Color', 'white');
hold on;

% 用于收集各类句柄（循环内填充，循环后组装图例）
h_raw_all  = cell(1, 2);
h_aper_all = cell(1, 2);
h_fit_all  = cell(1, 2);

for g = 1:2
    if ~valid_group(g)
        continue;
    end

    fooof_result = results_all{g};
    selected_pk  = fooof_result.selected_pk;
    n_peaks      = size(fooof_result.peak_params, 1);

    x           = fooof_result.freqs;
    y_raw       = fooof_result.power_spectrum;
    y_aperiodic = fooof_result.aperiodic_fit;
    y_fit       = fooof_result.fooofed_spectrum;

    % --- 原始功率谱 ---
    h_raw  = plot(x, y_raw, '-', ...
        'Color', COLOR_RAW{g}, ...
        'LineWidth', 1.5);

    % --- 非周期拟合 ---
    h_aper = plot(x, y_aperiodic, '--', ...
        'Color', COLOR_APERIODIC{g}, ...
        'LineWidth', 1.6);

    % --- 完整拟合（最突出） ---
    h_fit  = plot(x, y_fit, '-', ...
        'Color', COLOR_FIT{g}, ...
        'LineWidth', 2.4);

    % 收集各类句柄，图例在循环结束后统一组装
    h_raw_all{g}  = h_raw;
    h_aper_all{g} = h_aper;
    h_fit_all{g}  = h_fit;

    % --- 标记 1 Hz 处 offset 交点（黄色虚线左端点） ---
    offset_val   = fooof_result.aperiodic_params(1);   % = aperiodic_fit 在1Hz处的值
    COLOR_OFFSET = COLOR_APERIODIC{g};                 % 与 Aperiodic fit 同色
    h_offset_dot = plot(1, offset_val, 'o', ...
        'MarkerSize',      9, ...
        'MarkerFaceColor', COLOR_OFFSET, ...
        'MarkerEdgeColor', COLOR_OFFSET * 0.75, ...
        'LineWidth', 1.2);

    % 全局只保存第一次出现的句柄用于图例
    if ~exist('h_offset_legend', 'var')
        h_offset_legend = h_offset_dot;
    end

    %% 标注峰值
    h_sel_g   = [];
    h_unsel_g = [];

    if n_peaks > 0
        y_range = range(y_fit);

        for pk = 1:n_peaks
            cf = fooof_result.peak_params(pk, 1);
            [~, cf_idx] = min(abs(x - cf));
            y_on_fit    = y_fit(cf_idx);

            is_selected = (pk == selected_pk);

            if is_selected
                mk_color   = COLOR_SELECTED{g};
                mk_edge    = mk_color * 0.75;
                mk_size    = 9;
                txt_color  = mk_color * 0.85;
                txt_weight = 'bold';
            else
                mk_color   = COLOR_UNSEL{g};
                mk_edge    = mk_color * 0.7;
                mk_size    = 8;
                txt_color  = [0.5, 0.5, 0.5];
                txt_weight = 'normal';
            end

            h_dot = plot(cf, y_on_fit, 'o', ...
                'MarkerSize',      mk_size, ...
                'MarkerFaceColor', mk_color, ...
                'MarkerEdgeColor', mk_edge, ...
                'LineWidth', 1.2);

            % 记录图例句柄（每组只记录第一次出现）
            if is_selected && isempty(h_sel_g)
                h_sel_g = h_dot;
            elseif ~is_selected && isempty(h_unsel_g)
                h_unsel_g = h_dot;
            end

            % 文字标注：统一在峰值上方，偏移量小；
            % 患者2额外向右微移 0.4 Hz，避免与患者1标注文字重叠
            y_offset  = y_range * 0.02;
            x_offset  = 0.4 * (g - 1);   % g=1: 0, g=2: +0.4 Hz
            text(cf + x_offset, y_on_fit + y_offset, sprintf('%.1f Hz', cf), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'bottom', ...
                'FontSize',   8.5, ...
                'FontWeight', txt_weight, ...
                'Color',      txt_color);
        end

        % Peak 类型：全局首次出现时保存句柄
        if ~isempty(h_unsel_g) && ~exist('h_unsel_legend', 'var')
            h_unsel_legend = h_unsel_g;
        end
        if ~isempty(h_sel_g) && ~exist('h_sel_legend', 'var')
            h_sel_legend = h_sel_g;
        end
    end
end

%% ==================== 坐标轴美化 ====================
xlabel('Frequency (Hz)', 'FontSize', 12, 'Color', [0.2 0.2 0.2]);
ylabel('log_{10}(\muV^{2}/Hz)', 'FontSize', 12, 'Color', [0.2 0.2 0.2]);

xlim([1, 40]);
xticks([1, 5, 10, 15, 20, 25, 30, 35, 40]);   % 显式标出 1 Hz

ax = gca;
ax.FontSize  = 11;
ax.Color     = 'white';
ax.LineWidth = 0.8;
ax.Box       = 'off';
ax.XAxis.Color = [0.3 0.3 0.3];
ax.YAxis.Color = [0.3 0.3 0.3];
ax.TickDir    = 'out';
ax.TickLength = [0.012, 0.025];
grid off;

%% ==================== 图例组装（顺序：两条Original → Aperiodic → Full Fit → Unselected Peak → Selected Peak） ====================
legend_handles = [];
legend_labels  = {};

% 两条 Original Spectrum（绿、粉，放在一起）
for g = 1:2
    if valid_group(g)
        legend_handles(end+1) = h_raw_all{g};
        legend_labels{end+1}  = sprintf('Original Spectrum (%s)', group_labels{g});
    end
end

% Aperiodic Fit（两组同色，只显示一条）
for g = 1:2
    if valid_group(g)
        legend_handles(end+1) = h_aper_all{g};
        legend_labels{end+1}  = 'Aperiodic Fit';
        break;   % 只加一次
    end
end

% Full Model Fit（两组同色，只显示一条）
for g = 1:2
    if valid_group(g)
        legend_handles(end+1) = h_fit_all{g};
        legend_labels{end+1}  = 'Full Model Fit';
        break;
    end
end

% Unselected Peak（两组同色，只显示一条）
if exist('h_unsel_legend', 'var')
    legend_handles(end+1) = h_unsel_legend;
    legend_labels{end+1}  = 'Unselected Peak';
end

% Selected Peak（两组同色，只显示一条）
if exist('h_sel_legend', 'var')
    legend_handles(end+1) = h_sel_legend;
    legend_labels{end+1}  = 'Selected Peak';
end

% Offset at 1 Hz
if exist('h_offset_legend', 'var')
    legend_handles(end+1) = h_offset_legend;
    legend_labels{end+1}  = 'Aperiodic Offset';
end

if ~isempty(legend_handles)
    lg = legend(legend_handles, legend_labels, ...
        'Location', 'northeast', ...
        'FontSize', 9, ...
        'Box', 'off');
    lg.ItemTokenSize = [20, 10];
end

%% ==================== 保存图片 ====================
output_filename = 'TwoGroups_HA_vs_LTD_FOOOF.png';
output_fullpath = fullfile(output_folder, output_filename);
exportgraphics(fig, output_fullpath, 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);

fprintf('\n图片已保存: %s\n', output_fullpath);
fprintf('========================================\n');
fprintf('处理完成！\n');
fprintf('========================================\n');


%% ==================== 简化的 FOOOF 拟合函数 ====================
function results = fooof_fit_simple(freqs, power_spectrum, settings)

    freqs          = freqs(:);
    power_spectrum = power_spectrum(:);

    if any(power_spectrum <= 0)
        min_pos = min(power_spectrum(power_spectrum > 0));
        power_spectrum(power_spectrum <= 0) = min_pos / 10;
    end

    log_freqs  = log10(freqs);
    log_powers = log10(power_spectrum);

    valid = isfinite(log_freqs) & isfinite(log_powers);
    p = polyfit(log_freqs(valid), log_powers(valid), 1);
    aperiodic_params = [p(2), -p(1)];

    for iter = 1:5
        aperiodic_fit = aperiodic_params(1) - aperiodic_params(2) * log_freqs;
        flat_spectrum = log_powers - aperiodic_fit;

        [peak_params, peak_fit] = fit_peaks_simple(freqs, flat_spectrum, settings);

        spectrum_no_peaks = log_powers - peak_fit;
        valid = isfinite(log_freqs) & isfinite(spectrum_no_peaks);
        if sum(valid) < 3, break; end

        p = polyfit(log_freqs(valid), spectrum_no_peaks(valid), 1);
        aperiodic_params_new = [p(2), -p(1)];

        if max(abs(aperiodic_params_new - aperiodic_params)) < 0.001
            break;
        end
        aperiodic_params = aperiodic_params_new;
    end

    aperiodic_fit = aperiodic_params(1) - aperiodic_params(2) * log_freqs;
    flat_spectrum = log_powers - aperiodic_fit;
    [peak_params, peak_fit] = fit_peaks_simple(freqs, flat_spectrum, settings);

    full_fit  = aperiodic_fit + peak_fit;
    ss_res    = sum((log_powers - full_fit).^2);
    ss_tot    = sum((log_powers - mean(log_powers)).^2);
    r_squared = 1 - ss_res / ss_tot;

    results.aperiodic_params = aperiodic_params;
    results.peak_params      = peak_params;
    results.r_squared        = r_squared;
    results.aperiodic_fit    = aperiodic_fit;
    results.peak_fit         = peak_fit;
    results.fooofed_spectrum = full_fit;
    results.power_spectrum   = log_powers;
    results.freqs            = freqs;
end


function [peak_params, peak_fit] = fit_peaks_simple(freqs, flat_spectrum, settings)

    peak_params   = [];
    peak_fit      = zeros(size(flat_spectrum));
    flat_spectrum = flat_spectrum(:);
    freqs         = freqs(:);

    threshold = settings.min_peak_height;
    residual  = flat_spectrum;
    n_peaks   = 0;

    while n_peaks < settings.max_n_peaks
        [max_val, max_idx] = max(residual);

        if max_val < threshold, break; end

        center_freq = freqs(max_idx);
        half_max    = max_val / 2;

        left_idx = max_idx;
        while left_idx > 1 && residual(left_idx) > half_max
            left_idx = left_idx - 1;
        end

        right_idx = max_idx;
        while right_idx < length(residual) && residual(right_idx) > half_max
            right_idx = right_idx + 1;
        end

        fwhm      = freqs(right_idx) - freqs(left_idx);
        bandwidth = max(fwhm / 2.355, 0.5);

        if bandwidth < settings.peak_width_limits(1) || bandwidth > settings.peak_width_limits(2)
            residual(max_idx) = 0;
            continue;
        end

        gaussian = max_val * exp(-((freqs - center_freq).^2) / (2 * bandwidth^2));

        n_peaks = n_peaks + 1;
        peak_params(n_peaks, :) = [center_freq, max_val, bandwidth];

        peak_fit = peak_fit + gaussian;
        residual = residual - gaussian;
    end

    if ~isempty(peak_params)
        [~, sort_idx] = sort(peak_params(:, 1));
        peak_params = peak_params(sort_idx, :);
    end
end