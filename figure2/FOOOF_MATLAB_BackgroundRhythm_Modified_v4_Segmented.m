%% FOOOF_MATLAB_BackgroundRhythm_Modified_v4_Segmented.m
% 功能：纯 MATLAB 实现 FOOOF (specparam) 算法 - 分段版本
% 
% v4 修改：
% - 将数据按1分钟分段，每段单独进行FOOOF分析
% - 最后不足1分钟的数据：若<30秒则合并到前一段，否则单独分析
% - Excel输出精简：移除 has_alpha/theta/delta 列
%
% 功能说明：
% 1. 按频段检测峰值（使用阈值0.15）
%    - Delta: 1.5-4 Hz (排除1Hz，避免DC干扰)
%    - Theta: 4-8 Hz
%    - Alpha: 8-13 Hz
%    - Beta: 13-30 Hz
% 2. 背景节律判定规则：
%    - 在Delta/Theta/Alpha三个频段中检测是否存在超过阈值的峰值
%    - 如果多个频段都有峰值，选择频率最快的作为背景节律
%    - 优先级：Alpha > Theta > Delta
%    - 可能没有背景节律（None）
% 3. Beta单独检测和报告
% 4. 每个被试的每个时间段单独分析并输出
% 5. 同时输出全脑结果和后头部电极结果

clear; 
clc;
close all;

%% ==================== 路径设置 ====================
% ★★★ 请修改为你的数据路径 ★★★
filepath = 'F:\EEG\LowVoltageResearch\PreprocessedData\TestingData\RunICA';
outpath = 'F:\EEG\LowVoltageResearch\OutputResult\外部验证数据\FOOOF_Segmented_T20';

%% ==================== 分段设置 ====================
segment_duration_sec = 60;  % 每段时长（秒）
min_last_segment_sec = 30;  % 最后一段的最小时长，小于此值则合并到前一段

%% ==================== 通道选择设置 ====================
% 后头部电极列表
posterior_channels = {'P3', 'P4', 'O1', 'O2', 'T5', 'T6', 'Pz' ,'PO3', 'PO4', 'PO7', 'PO8', 'Oz'};

%% ==================== FOOOF 参数设置 ====================
fooof_settings = struct();
fooof_settings.freq_range = [1, 40];         % 分析频率范围 (Hz)
fooof_settings.peak_width_limits = [0.5, 12]; % 峰值带宽限制 (Hz)
fooof_settings.max_n_peaks = 10;             % 最大峰值数量
fooof_settings.min_peak_height = 0.2;       % 最小峰值高度阈值（对数功率单位）  修改此处
fooof_settings.aperiodic_mode = 'fixed';     % 'fixed' 或 'knee'

%% ==================== 频段定义 ====================
band_settings = struct();
band_settings.delta_range = [1.5, 4];    % Delta 范围 (排除1Hz)
band_settings.theta_range = [4, 8];      % Theta 范围
band_settings.alpha_range = [8, 13];     % Alpha 范围
band_settings.beta_range = [13, 30];     % Beta 范围
band_settings.peak_threshold = 0.15;     % 峰值检测阈值

%% ==================== 图片保存设置 ====================
save_figures = true;           % 是否保存每个被试的图片
max_segments_per_figure = 6;   % 每张图最多显示的分段数（超过则分多张图）

%% ==================== 创建输出文件夹 ====================
if ~exist(outpath, 'dir')
    mkdir(outpath);
end

% 创建图片子文件夹
fig_folder = fullfile(outpath, 'Figures');
if save_figures && ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

%% ==================== 添加工具箱路径 ====================
% ★★★ 请修改为你的EEGLAB路径 ★★★
addpath('F:\matlab_downloadtoolbox\EEGToolbox\eeglab2023.0');
eeglab nogui;

%% ==================== 检查数据文件 ====================
fprintf('\n========== FOOOF 背景节律分析 (v4 - 分段版本) ==========\n');
fprintf('数据路径: %s\n', filepath);
fprintf('分段设置: 每段 %d 秒，最后一段最小 %d 秒\n', segment_duration_sec, min_last_segment_sec);
fprintf('后头部电极: %s\n', strjoin(posterior_channels, ', '));

if ~exist(filepath, 'dir')
    error('错误: 文件夹不存在 - %s', filepath);
end

set_files = dir(fullfile(filepath, '*.set'));

if isempty(set_files)
    error('错误: 没有找到 .set 文件');
end

num_subjects = length(set_files);
fprintf('找到 %d 个被试\n\n', num_subjects);

%% ==================== 初始化结果存储（动态数组）====================
% 使用cell数组，因为每个被试可能有不同数量的分段
results_global_all = {};
results_posterior_all = {};

%% ==================== 主循环：处理每个被试 ====================
for f = 1:num_subjects
    datafile = set_files(f).name;
    [~, subject_id, ~] = fileparts(datafile);
    
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    fprintf('[%d/%d] 被试: %s\n', f, num_subjects, subject_id);
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    
    %% 加载数据
    try
        EEG = pop_loadset('filename', datafile, 'filepath', filepath);
    catch ME
        fprintf('  ✗ 加载失败: %s\n\n', ME.message);
        % 添加错误记录
        error_result = create_error_result(subject_id, 1, 'Error');
        results_global_all = [results_global_all; error_result];
        results_posterior_all = [results_posterior_all; error_result];
        continue;
    end
    
    fs = EEG.srate;
    fprintf('  采样率: %d Hz\n', fs);
    
    %% 处理数据维度
    if ndims(EEG.data) == 3
        eeg_data = reshape(EEG.data, size(EEG.data,1), []);
        fprintf('  数据类型: Epoched → 拼接后 %.1f 秒\n', size(eeg_data,2)/fs);
    else
        eeg_data = EEG.data;
        fprintf('  数据类型: Continuous, %.1f 秒\n', size(eeg_data,2)/fs);
    end
    eeg_data = double(eeg_data);
    
    total_samples = size(eeg_data, 2);
    total_duration_sec = total_samples / fs;
    
    %% 计算分段
    segment_samples = segment_duration_sec * fs;
    min_last_samples = min_last_segment_sec * fs;
    
    % 计算分段边界
    segment_boundaries = compute_segment_boundaries(total_samples, segment_samples, min_last_samples);
    num_segments = size(segment_boundaries, 1);
    
    fprintf('  总时长: %.1f 秒 → 分为 %d 段\n', total_duration_sec, num_segments);
    for seg = 1:num_segments
        seg_start_sec = (segment_boundaries(seg, 1) - 1) / fs;
        seg_end_sec = segment_boundaries(seg, 2) / fs;
        fprintf('    段 %d: %.1f - %.1f 秒 (%.1f 秒)\n', seg, seg_start_sec, seg_end_sec, seg_end_sec - seg_start_sec);
    end
    
    %% 获取通道标签
    if ~isempty(EEG.chanlocs)
        all_labels = {EEG.chanlocs.labels};
    else
        all_labels = arrayfun(@(x) sprintf('Ch%d', x), 1:size(eeg_data,1), 'UniformOutput', false);
    end
    
    %% 找到后头部电极的索引
    posterior_idx = [];
    posterior_found = {};
    for ch_name = posterior_channels
        ch_idx = find(strcmpi(all_labels, ch_name{1}));
        if ~isempty(ch_idx)
            posterior_idx = [posterior_idx, ch_idx(1)];
            posterior_found = [posterior_found, ch_name{1}];
        end
    end
    
    %% 对每个分段进行分析
    segment_results_global = {};  % 存储该被试所有分段的全脑结果
    segment_results_posterior = {};  % 存储该被试所有分段的后头部结果
    segment_fooof_global = {};  % 存储FOOOF详细结果用于画图
    segment_fooof_posterior = {};
    
    for seg = 1:num_segments
        seg_start = segment_boundaries(seg, 1);
        seg_end = segment_boundaries(seg, 2);
        seg_data = eeg_data(:, seg_start:seg_end);
        
        seg_start_sec = (seg_start - 1) / fs;
        seg_end_sec = seg_end / fs;
        seg_label = sprintf('%d-%ds', round(seg_start_sec), round(seg_end_sec));
        
        fprintf('\n  ▶ 分析段 %d/%d (%s)\n', seg, num_segments, seg_label);
        
        %% 全脑分析
        fprintf('    【全脑】(%d 通道)\n', size(seg_data, 1));
        
        [result_global, fooof_detail_global] = analyze_segment_with_details(seg_data, 1:size(seg_data,1), all_labels, ...
            fs, fooof_settings, band_settings, subject_id, seg, seg_label);
        
        results_global_all = [results_global_all; result_global];
        segment_results_global{seg} = result_global;
        segment_fooof_global{seg} = fooof_detail_global;
        print_segment_results(result_global, '全脑');
        
        %% 后头部电极分析
        if isempty(posterior_idx)
            fprintf('    【后头部】✗ 未找到指定电极\n');
            result_posterior = create_error_result(subject_id, seg, 'NoChannels');
            result_posterior.segment_label = seg_label;
            fooof_detail_posterior = [];
        else
            fprintf('    【后头部】(%d 通道: %s)\n', length(posterior_idx), strjoin(posterior_found, ', '));
            
            [result_posterior, fooof_detail_posterior] = analyze_segment_with_details(seg_data, posterior_idx, all_labels, ...
                fs, fooof_settings, band_settings, subject_id, seg, seg_label);
            
            print_segment_results(result_posterior, '后头部');
        end
        
        results_posterior_all = [results_posterior_all; result_posterior];
        segment_results_posterior{seg} = result_posterior;
        segment_fooof_posterior{seg} = fooof_detail_posterior;
    end
    
    %% 为该被试生成图片
    if save_figures && num_segments > 0
        % 全脑图片
        plot_subject_segments(subject_id, segment_results_global, segment_fooof_global, ...
            'Global', fig_folder, max_segments_per_figure, band_settings);
        
        % 后头部图片（如果有有效数据）
        if ~isempty(posterior_idx)
            plot_subject_segments(subject_id, segment_results_posterior, segment_fooof_posterior, ...
                'Posterior', fig_folder, max_segments_per_figure, band_settings);
        end
        
        fprintf('  ✓ 图片已保存到: %s\n', fig_folder);
    end
    end
    
    fprintf('\n');
% end

%% ==================== 保存结果到 Excel ====================
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('保存结果...\n');

excel_file = fullfile(outpath, 'FOOOF_Results_Segmented.xlsx');

% 保存全脑结果
T_global = create_results_table_simplified(results_global_all);
writetable(T_global, excel_file, 'Sheet', 'Global_AllChannels');
fprintf('✓ 已保存全脑结果到 Sheet: Global_AllChannels\n');

% 保存后头部结果
T_posterior = create_results_table_simplified(results_posterior_all);
writetable(T_posterior, excel_file, 'Sheet', 'Posterior');
fprintf('✓ 已保存后头部结果到 Sheet: Posterior\n');

fprintf('✓ Excel文件: %s\n', excel_file);

%% ==================== 汇总统计 ====================
fprintf('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('【汇总统计】(共 %d 名被试，%d 个分段)\n', num_subjects, length(results_global_all));
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

% 全脑统计
print_summary_stats(results_global_all, '全脑');

% 后头部统计
print_summary_stats(results_posterior_all, '后头部');

fprintf('\n========== 分析完成！ ==========\n');

%% ==================== 辅助函数 ====================

function boundaries = compute_segment_boundaries(total_samples, segment_samples, min_last_samples)
    % 计算分段边界
    % 输入:
    %   total_samples - 总样本数
    %   segment_samples - 每段样本数
    %   min_last_samples - 最后一段最小样本数
    % 输出:
    %   boundaries - N x 2 矩阵，每行是 [起始样本, 结束样本]
    
    if total_samples <= segment_samples
        % 数据不足一段，作为整体处理
        boundaries = [1, total_samples];
        return;
    end
    
    num_full_segments = floor(total_samples / segment_samples);
    remaining_samples = mod(total_samples, segment_samples);
    
    boundaries = [];
    
    for i = 1:num_full_segments
        start_sample = (i - 1) * segment_samples + 1;
        end_sample = i * segment_samples;
        boundaries = [boundaries; start_sample, end_sample];
    end
    
    % 处理剩余样本
    if remaining_samples > 0
        if remaining_samples >= min_last_samples
            % 剩余部分足够长，单独成段
            start_sample = num_full_segments * segment_samples + 1;
            end_sample = total_samples;
            boundaries = [boundaries; start_sample, end_sample];
        else
            % 剩余部分太短，合并到最后一段
            if ~isempty(boundaries)
                boundaries(end, 2) = total_samples;
            else
                boundaries = [1, total_samples];
            end
        end
    end
end


function result = create_error_result(subject_id, segment_num, error_type)
    % 创建错误结果结构体
    
    result = struct();
    result.subject_id = subject_id;
    result.segment_num = segment_num;
    result.segment_label = '';
    result.background_rhythm = error_type;
    result.background_peak_freq = NaN;
    result.background_peak_power = NaN;
    result.delta_peak_freq = NaN;
    result.delta_peak_power = NaN;
    result.theta_peak_freq = NaN;
    result.theta_peak_power = NaN;
    result.alpha_peak_freq = NaN;
    result.alpha_peak_power = NaN;
    result.beta_peak_freq = NaN;
    result.beta_peak_power = NaN;
    result.aperiodic_exp = NaN;
    result.aperiodic_offset = NaN;
    result.r_squared = NaN;
    result.all_peaks_str = 'None';
end


function [result, fooof_detail] = analyze_segment_with_details(seg_data, channel_idx, all_labels, fs, fooof_settings, band_settings, subject_id, segment_num, segment_label)
    % 对单个数据段进行FOOOF分析，同时返回FOOOF详细结果用于画图
    %
    % 输出:
    %   result - 简化的结果结构体
    %   fooof_detail - 包含freqs, psd, fooof_results的详细结果
    
    result = struct();
    result.subject_id = subject_id;
    result.segment_num = segment_num;
    result.segment_label = segment_label;
    
    fooof_detail = struct();
    
    num_channels = length(channel_idx);
    
    %% 计算每个通道的功率谱密度，然后取平均
    nfft = 2^nextpow2(fs * 2);
    
    % 计算每个通道的PSD
    psd_all = [];
    for i = 1:num_channels
        [psd_temp, freqs] = pwelch(seg_data(channel_idx(i),:)', hanning(nfft), nfft/2, nfft, fs);
        psd_all(:, i) = psd_temp;
    end
    
    % 对所有通道的PSD取平均
    psd_avg = mean(psd_all, 2);
    
    %% 限制频率范围
    freq_idx = freqs >= fooof_settings.freq_range(1) & freqs <= fooof_settings.freq_range(2);
    freqs = freqs(freq_idx);
    psd_avg = psd_avg(freq_idx);
    
    % 保存用于画图
    fooof_detail.freqs = freqs;
    fooof_detail.psd = psd_avg;
    
    %% 对平均后的PSD做一次FOOOF
    fooof_results = fooof_fit_simple(freqs, psd_avg, fooof_settings);
    fooof_detail.fooof_results = fooof_results;
    
    if isempty(fooof_results)
        result.background_rhythm = 'Error';
        result.background_peak_freq = NaN;
        result.background_peak_power = NaN;
        result.delta_peak_freq = NaN;
        result.delta_peak_power = NaN;
        result.theta_peak_freq = NaN;
        result.theta_peak_power = NaN;
        result.alpha_peak_freq = NaN;
        result.alpha_peak_power = NaN;
        result.beta_peak_freq = NaN;
        result.beta_peak_power = NaN;
        result.aperiodic_exp = NaN;
        result.aperiodic_offset = NaN;
        result.r_squared = NaN;
        result.all_peaks_str = 'None';
        return;
    end
    
    %% 记录FOOOF基本参数
    result.aperiodic_offset = fooof_results.aperiodic_params(1);
    result.aperiodic_exp = fooof_results.aperiodic_params(2);
    result.r_squared = fooof_results.r_squared;
    
    %% 按频段检测峰值
    peak_params = fooof_results.peak_params;
    
    [has_delta, result.delta_peak_freq, result.delta_peak_power] = ...
        detect_band_peak(peak_params, band_settings.delta_range, band_settings.peak_threshold);
    
    [has_theta, result.theta_peak_freq, result.theta_peak_power] = ...
        detect_band_peak(peak_params, band_settings.theta_range, band_settings.peak_threshold);
    
    [has_alpha, result.alpha_peak_freq, result.alpha_peak_power] = ...
        detect_band_peak(peak_params, band_settings.alpha_range, band_settings.peak_threshold);
    
    [has_beta, result.beta_peak_freq, result.beta_peak_power] = ...
        detect_band_peak(peak_params, band_settings.beta_range, band_settings.peak_threshold);
    
    %% 确定最终的背景节律（优先选择最快的）
    result.background_rhythm = determine_background_rhythm(has_alpha, has_theta, has_delta);
    
    % 获取背景节律的峰值参数
    switch result.background_rhythm
        case 'Alpha'
            result.background_peak_freq = result.alpha_peak_freq;
            result.background_peak_power = result.alpha_peak_power;
        case 'Theta'
            result.background_peak_freq = result.theta_peak_freq;
            result.background_peak_power = result.theta_peak_power;
        case 'Delta'
            result.background_peak_freq = result.delta_peak_freq;
            result.background_peak_power = result.delta_peak_power;
        otherwise
            result.background_peak_freq = NaN;
            result.background_peak_power = NaN;
    end
    
    %% 记录所有峰值（字符串格式）
    if ~isempty(peak_params)
        peak_freqs = peak_params(:, 1)';
        result.all_peaks_str = sprintf('%.1f ', peak_freqs);
    else
        result.all_peaks_str = 'None';
    end
end


function plot_subject_segments(subject_id, segment_results, segment_fooof, region_name, fig_folder, max_per_fig, band_settings)
    % 为单个被试绘制所有分段的FOOOF结果
    %
    % 输入:
    %   subject_id - 被试ID
    %   segment_results - cell数组，每个元素是一个分段的结果
    %   segment_fooof - cell数组，每个元素是FOOOF详细结果
    %   region_name - 'Global' 或 'Posterior'
    %   fig_folder - 图片保存文件夹
    %   max_per_fig - 每张图最多显示的分段数
    %   band_settings - 频段设置
    
    num_segments = length(segment_results);
    
    if num_segments == 0
        return;
    end
    
    % 计算需要多少张图
    num_figures = ceil(num_segments / max_per_fig);
    
    % 定义频段颜色
    band_colors = struct();
    band_colors.delta = [0.6, 0.2, 0.8];  % 紫色
    band_colors.theta = [0.2, 0.6, 0.2];  % 绿色
    band_colors.alpha = [0.8, 0.2, 0.2];  % 红色
    band_colors.beta = [0.2, 0.4, 0.8];   % 蓝色
    
    for fig_num = 1:num_figures
        % 计算当前图包含的分段
        start_seg = (fig_num - 1) * max_per_fig + 1;
        end_seg = min(fig_num * max_per_fig, num_segments);
        segs_in_fig = start_seg:end_seg;
        n_segs = length(segs_in_fig);
        
        % 创建图形
        fig = figure('Position', [50, 50, 1600, 250 * n_segs], 'Visible', 'off');
        
        for i = 1:n_segs
            seg_idx = segs_in_fig(i);
            result = segment_results{seg_idx};
            fooof_detail = segment_fooof{seg_idx};
            
            if isempty(fooof_detail) || isempty(fooof_detail.fooof_results)
                continue;
            end
            
            fooof_results = fooof_detail.fooof_results;
            freqs = fooof_detail.freqs;
            
            %% 左图：原始功率谱 + 拟合
            subplot(n_segs, 2, (i-1)*2 + 1);
            
            plot(freqs, fooof_results.power_spectrum, 'b-', 'LineWidth', 1.5); hold on;
            plot(freqs, fooof_results.aperiodic_fit, 'r--', 'LineWidth', 1.2);
            plot(freqs, fooof_results.fooofed_spectrum, 'Color', [0, 0.7, 0], 'LineWidth', 1.5);
            
            % 添加频段背景色
            ylims = ylim;
            fill([band_settings.delta_range(1), band_settings.delta_range(2), band_settings.delta_range(2), band_settings.delta_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.delta, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([band_settings.theta_range(1), band_settings.theta_range(2), band_settings.theta_range(2), band_settings.theta_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.theta, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([band_settings.alpha_range(1), band_settings.alpha_range(2), band_settings.alpha_range(2), band_settings.alpha_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.alpha, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([band_settings.beta_range(1), band_settings.beta_range(2), band_settings.beta_range(2), band_settings.beta_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.beta, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            
            xlabel('Frequency (Hz)');
            ylabel('Log Power');
            title(sprintf('Segment %d (%s) - Power Spectrum', seg_idx, result.segment_label), 'Interpreter', 'none');
            legend('Original', 'Aperiodic', 'Full Fit', 'Location', 'northeast');
            xlim([1 40]);
            grid on;
            
            %% 右图：扁平化功率谱 + 峰值检测
            subplot(n_segs, 2, (i-1)*2 + 2);
            
            flat_spectrum = fooof_results.power_spectrum - fooof_results.aperiodic_fit;
            plot(freqs, flat_spectrum, 'k-', 'LineWidth', 1.2); hold on;
            plot(freqs, fooof_results.peak_fit, 'm-', 'LineWidth', 1.5);
            
            % 添加频段背景色
            ylims = ylim;
            fill([band_settings.delta_range(1), band_settings.delta_range(2), band_settings.delta_range(2), band_settings.delta_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.delta, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([band_settings.theta_range(1), band_settings.theta_range(2), band_settings.theta_range(2), band_settings.theta_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.theta, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([band_settings.alpha_range(1), band_settings.alpha_range(2), band_settings.alpha_range(2), band_settings.alpha_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.alpha, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            fill([band_settings.beta_range(1), band_settings.beta_range(2), band_settings.beta_range(2), band_settings.beta_range(1)], ...
                 [ylims(1), ylims(1), ylims(2), ylims(2)], band_colors.beta, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            
            % 标记峰值位置
            if ~isempty(fooof_results.peak_params)
                for pk = 1:size(fooof_results.peak_params, 1)
                    cf = fooof_results.peak_params(pk, 1);
                    pw = fooof_results.peak_params(pk, 2);
                    
                    % 根据频率确定颜色
                    if cf >= band_settings.delta_range(1) && cf <= band_settings.delta_range(2)
                        pk_color = band_colors.delta;
                        band_name = 'δ';
                    elseif cf >= band_settings.theta_range(1) && cf <= band_settings.theta_range(2)
                        pk_color = band_colors.theta;
                        band_name = 'θ';
                    elseif cf >= band_settings.alpha_range(1) && cf <= band_settings.alpha_range(2)
                        pk_color = band_colors.alpha;
                        band_name = 'α';
                    elseif cf >= band_settings.beta_range(1) && cf <= band_settings.beta_range(2)
                        pk_color = band_colors.beta;
                        band_name = 'β';
                    else
                        pk_color = [0.5, 0.5, 0.5];
                        band_name = '';
                    end
                    
                    plot(cf, pw, 'o', 'MarkerSize', 10, 'MarkerFaceColor', pk_color, 'MarkerEdgeColor', 'k');
                    text(cf, pw + 0.08, sprintf('%s%.1fHz', band_name, cf), 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
                end
            end
            
            % 绘制阈值线
            yline(band_settings.peak_threshold, 'r--', sprintf('Threshold=%.2f', band_settings.peak_threshold), 'LineWidth', 1);
            
            xlabel('Frequency (Hz)');
            ylabel('Residual Power');
            
            % 构建标题信息
            title_str = sprintf('BG: %s', result.background_rhythm);
            if ~isnan(result.background_peak_freq)
                title_str = sprintf('%s (%.1fHz)', title_str, result.background_peak_freq);
            end
            title_str = sprintf('%s | R²=%.3f | exp=%.2f', title_str, result.r_squared, result.aperiodic_exp);
            title(title_str);
            
            xlim([1 40]);
            grid on;
        end
        
        % 添加总标题
        if num_figures > 1
            sgtitle(sprintf('%s - %s Analysis (Page %d/%d)', subject_id, region_name, fig_num, num_figures), ...
                'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
        else
            sgtitle(sprintf('%s - %s Analysis (%d segments)', subject_id, region_name, num_segments), ...
                'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        % 保存图片
        if num_figures > 1
            output_filename = sprintf('%s_%s_Segments_Page%d.png', subject_id, region_name, fig_num);
        else
            output_filename = sprintf('%s_%s_Segments.png', subject_id, region_name);
        end
        output_fullpath = fullfile(fig_folder, output_filename);
        saveas(fig, output_fullpath);
        close(fig);
    end
end


function print_segment_results(result, region_name)
    % 打印单个分段的结果
    
    fprintf('      Aperiodic: exp=%.3f, offset=%.3f, R²=%.3f\n', ...
        result.aperiodic_exp, result.aperiodic_offset, result.r_squared);
    
    % 显示各频段峰值
    fprintf('      ');
    if ~isnan(result.delta_peak_freq)
        fprintf('Delta:%.1fHz ', result.delta_peak_freq);
    end
    if ~isnan(result.theta_peak_freq)
        fprintf('Theta:%.1fHz ', result.theta_peak_freq);
    end
    if ~isnan(result.alpha_peak_freq)
        fprintf('Alpha:%.1fHz ', result.alpha_peak_freq);
    end
    if ~isnan(result.beta_peak_freq)
        fprintf('Beta:%.1fHz ', result.beta_peak_freq);
    end
    fprintf('\n');
    
    % 背景节律
    fprintf('      ★ 背景节律: %s', result.background_rhythm);
    if ~strcmp(result.background_rhythm, 'None') && ~strcmp(result.background_rhythm, 'Error') && ~strcmp(result.background_rhythm, 'NoChannels')
        fprintf(' (%.2f Hz)', result.background_peak_freq);
    end
    fprintf('\n');
end


function print_summary_stats(results_all, region_name)
    % 打印汇总统计
    
    fprintf('\n【%s统计】\n', region_name);
    
    num_segments = length(results_all);
    
    % 提取背景节律
    bg_rhythms = cellfun(@(x) x.background_rhythm, results_all, 'UniformOutput', false);
    
    % 背景节律分布
    n_alpha = sum(strcmp(bg_rhythms, 'Alpha'));
    n_theta = sum(strcmp(bg_rhythms, 'Theta'));
    n_delta = sum(strcmp(bg_rhythms, 'Delta'));
    n_none = sum(strcmp(bg_rhythms, 'None'));
    
    fprintf('  背景节律分布 (共 %d 个分段):\n', num_segments);
    fprintf('    Alpha: %d (%.1f%%)\n', n_alpha, n_alpha/num_segments*100);
    fprintf('    Theta: %d (%.1f%%)\n', n_theta, n_theta/num_segments*100);
    fprintf('    Delta: %d (%.1f%%)\n', n_delta, n_delta/num_segments*100);
    fprintf('    None:  %d (%.1f%%)\n', n_none, n_none/num_segments*100);
    
    % Beta检测率
    n_beta = sum(cellfun(@(x) ~isnan(x.beta_peak_freq), results_all));
    fprintf('  Beta检测率: %d/%d (%.1f%%)\n', n_beta, num_segments, n_beta/num_segments*100);
end


function T = create_results_table_simplified(results_all)
    % 创建精简的结果表格（移除has_alpha等布尔列）
    
    num_rows = length(results_all);
    
    % 预分配
    subject_id = cell(num_rows, 1);
    segment_num = zeros(num_rows, 1);
    segment_label = cell(num_rows, 1);
    background_rhythm = cell(num_rows, 1);
    background_peak_freq = nan(num_rows, 1);
    background_peak_power = nan(num_rows, 1);
    delta_peak_freq = nan(num_rows, 1);
    delta_peak_power = nan(num_rows, 1);
    theta_peak_freq = nan(num_rows, 1);
    theta_peak_power = nan(num_rows, 1);
    alpha_peak_freq = nan(num_rows, 1);
    alpha_peak_power = nan(num_rows, 1);
    beta_peak_freq = nan(num_rows, 1);
    beta_peak_power = nan(num_rows, 1);
    aperiodic_exp = nan(num_rows, 1);
    aperiodic_offset = nan(num_rows, 1);
    r_squared = nan(num_rows, 1);
    all_peaks_str = cell(num_rows, 1);
    
    for i = 1:num_rows
        r = results_all{i};
        subject_id{i} = r.subject_id;
        segment_num(i) = r.segment_num;
        segment_label{i} = r.segment_label;
        background_rhythm{i} = r.background_rhythm;
        background_peak_freq(i) = r.background_peak_freq;
        background_peak_power(i) = r.background_peak_power;
        delta_peak_freq(i) = r.delta_peak_freq;
        delta_peak_power(i) = r.delta_peak_power;
        theta_peak_freq(i) = r.theta_peak_freq;
        theta_peak_power(i) = r.theta_peak_power;
        alpha_peak_freq(i) = r.alpha_peak_freq;
        alpha_peak_power(i) = r.alpha_peak_power;
        beta_peak_freq(i) = r.beta_peak_freq;
        beta_peak_power(i) = r.beta_peak_power;
        aperiodic_exp(i) = r.aperiodic_exp;
        aperiodic_offset(i) = r.aperiodic_offset;
        r_squared(i) = r.r_squared;
        all_peaks_str{i} = r.all_peaks_str;
    end
    
    T = table(subject_id, segment_num, segment_label, ...
              background_rhythm, background_peak_freq, background_peak_power, ...
              delta_peak_freq, delta_peak_power, ...
              theta_peak_freq, theta_peak_power, ...
              alpha_peak_freq, alpha_peak_power, ...
              beta_peak_freq, beta_peak_power, ...
              aperiodic_exp, aperiodic_offset, r_squared, ...
              all_peaks_str, ...
              'VariableNames', {...
                  'Subject', 'Segment_Num', 'Segment_Time', ...
                  'Background_Rhythm', 'BG_Peak_Freq_Hz', 'BG_Peak_Power', ...
                  'Delta_Peak_Freq_Hz', 'Delta_Peak_Power', ...
                  'Theta_Peak_Freq_Hz', 'Theta_Peak_Power', ...
                  'Alpha_Peak_Freq_Hz', 'Alpha_Peak_Power', ...
                  'Beta_Peak_Freq_Hz', 'Beta_Peak_Power', ...
                  'Aperiodic_Exponent', 'Aperiodic_Offset', 'R_Squared', ...
                  'All_Peaks_Hz'});
end


function [has_peak, peak_freq, peak_power] = detect_band_peak(peak_params, freq_range, threshold)
    % 在指定频段内检测是否存在超过阈值的峰值
    
    has_peak = false;
    peak_freq = NaN;
    peak_power = NaN;
    
    if isempty(peak_params)
        return;
    end
    
    % 筛选在指定频段内的峰值
    in_range = peak_params(:,1) >= freq_range(1) & peak_params(:,1) <= freq_range(2);
    range_peaks = peak_params(in_range, :);
    
    if isempty(range_peaks)
        return;
    end
    
    % 找功率最高的峰值
    [max_power, max_idx] = max(range_peaks(:, 2));
    
    % 检查是否超过阈值
    if max_power >= threshold
        has_peak = true;
        peak_freq = range_peaks(max_idx, 1);
        peak_power = max_power;
    end
end


function rhythm = determine_background_rhythm(has_alpha, has_theta, has_delta)
    % 根据检测结果确定背景节律类型
    % 规则：如果多个频段都有峰值，选择频率最快的
    % 优先级：Alpha > Theta > Delta
    
    if has_alpha
        rhythm = 'Alpha';
    elseif has_theta
        rhythm = 'Theta';
    elseif has_delta
        rhythm = 'Delta';
    else
        rhythm = 'None';
    end
end


%% ==================== FOOOF 核心函数 ====================

function results = fooof_fit_simple(freqs, power_spectrum, settings)
    % FOOOF 拟合主函数
    
    results = [];
    
    freqs = freqs(:);
    power_spectrum = power_spectrum(:);
    
    % 处理非正值
    if any(power_spectrum <= 0)
        power_spectrum(power_spectrum <= 0) = min(power_spectrum(power_spectrum > 0)) / 10;
    end
    
    % 对数变换
    log_freqs = log10(freqs);
    log_powers = log10(power_spectrum);
    
    % 拟合非周期性成分
    valid = isfinite(log_freqs) & isfinite(log_powers);
    p = polyfit(log_freqs(valid), log_powers(valid), 1);
    aperiodic_params = [p(2), -p(1)];  % [offset, exponent]
    
    % 迭代拟合（5次）
    for iter = 1:5
        aperiodic_fit = aperiodic_params(1) - aperiodic_params(2) * log_freqs;
        flat_spectrum = log_powers - aperiodic_fit;
        
        [peak_params, peak_fit] = fit_peaks_simple(freqs, flat_spectrum, settings);
        
        spectrum_no_peaks = log_powers - peak_fit;
        valid = isfinite(log_freqs) & isfinite(spectrum_no_peaks);
        if sum(valid) < 3
            break;
        end
        p = polyfit(log_freqs(valid), spectrum_no_peaks(valid), 1);
        aperiodic_params_new = [p(2), -p(1)];
        
        if max(abs(aperiodic_params_new - aperiodic_params)) < 0.001
            break;
        end
        aperiodic_params = aperiodic_params_new;
    end
    
    % 最终拟合
    aperiodic_fit = aperiodic_params(1) - aperiodic_params(2) * log_freqs;
    flat_spectrum = log_powers - aperiodic_fit;
    [peak_params, peak_fit] = fit_peaks_simple(freqs, flat_spectrum, settings);
    
    full_fit = aperiodic_fit + peak_fit;
    
    % 计算 R²
    ss_res = sum((log_powers - full_fit).^2);
    ss_tot = sum((log_powers - mean(log_powers)).^2);
    r_squared = 1 - ss_res / ss_tot;
    
    results.aperiodic_params = aperiodic_params;
    results.peak_params = peak_params;
    results.r_squared = r_squared;
    results.aperiodic_fit = aperiodic_fit;
    results.peak_fit = peak_fit;
    results.fooofed_spectrum = full_fit;
    results.power_spectrum = log_powers;
    results.freqs = freqs;
end


function [peak_params, peak_fit] = fit_peaks_simple(freqs, flat_spectrum, settings)
    % 检测并拟合峰值
    
    peak_params = [];
    peak_fit = zeros(size(flat_spectrum));
    
    flat_spectrum = flat_spectrum(:);
    freqs = freqs(:);
    
    threshold = settings.min_peak_height;
    residual = flat_spectrum;
    n_peaks = 0;
    
    while n_peaks < settings.max_n_peaks
        [max_val, max_idx] = max(residual);
        
        if max_val < threshold
            break;
        end
        
        center_freq = freqs(max_idx);
        
        half_max = max_val / 2;
        
        left_idx = max_idx;
        while left_idx > 1 && residual(left_idx) > half_max
            left_idx = left_idx - 1;
        end
        
        right_idx = max_idx;
        while right_idx < length(residual) && residual(right_idx) > half_max
            right_idx = right_idx + 1;
        end
        
        fwhm = freqs(right_idx) - freqs(left_idx);
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
