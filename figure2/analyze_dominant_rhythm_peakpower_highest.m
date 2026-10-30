%% 分析EEG背景频率主导节律 (修改版)
% 该脚本读取Excel文件中的Posterior表，判断每个患者总体的主导频率
% 方法：在1-13Hz范围内（Alpha、Theta、Delta三个频段），选择Peak_Power最大的频段作为主导
% 只考虑Alpha、Theta、Delta三个频段，不考虑Beta
% 修改内容：
%   1. 增加R_Squared > 0.9的筛选条件
%   2. 增加Aperiodic_Exponent和Aperiodic_Offset的平均值输出
%   3. Mean_Freq_Hz计算方式：每个segment取Delta/Theta/Alpha中Peak_Power最大者
%      对应的频率作为该segment的背景频率，再对所有segment求均值
%      （若最大功率为NaN或对应频率为NaN，则排除该segment）
%   4. 主导节律判断（Dominant_Rhythm_Freq）改为直接根据Mean_Freq_Hz判断
%      （1-4Hz为Delta，4-8Hz为Theta，8-13Hz为Alpha）
%   5. 移除Background_Rhythm筛选条件，只保留R_Squared > 0.9
% 输出结果到新的Excel文件

clear; clc;

%% 设置参数
inputFile = 'lowThetaDelta_FOOOF_Results_Segmented.xlsx';  % 输入文件
sheetName = 'Posterior';                                % 分析的表
outputFile = 'lowThetaDelta_DominantRhythm_HighestPeakPower_Results.xlsx';            % 输出文件

% 频段定义 (Hz) - 只考虑Delta, Theta, Alpha
freqBands.Delta = [0.5, 4];
freqBands.Theta = [4, 8];
freqBands.Alpha = [8, 13];

% R_Squared阈值
R_Squared_Threshold = 0.9;

%% 读取数据
fprintf('正在读取数据: %s\n', inputFile);
data = readtable(inputFile, 'Sheet', sheetName);

% 获取所有唯一的患者
subjects = unique(data.Subject);
numSubjects = length(subjects);
fprintf('共发现 %d 位患者\n', numSubjects);
fprintf('只考虑Alpha、Theta、Delta三个频段\n');
fprintf('Mean_Freq_Hz计算方式：每个segment取Peak_Power最大频段对应的频率，再求均值\n');
fprintf('主导节律判断方式：根据Mean_Freq_Hz直接判断（1-4Hz:Delta, 4-8Hz:Theta, 8-13Hz:Alpha）\n');
fprintf('R_Squared阈值: > %.2f\n\n', R_Squared_Threshold);

%% 初始化结果表格
results = table();
results.Subject = subjects;
results.Total_Segments = zeros(numSubjects, 1);
results.Valid_Segments = zeros(numSubjects, 1);
results.Excluded_By_R_Squared = zeros(numSubjects, 1);  % 新增：因R_Squared被排除的segment数

% 频率统计
results.Mean_Freq_Hz = zeros(numSubjects, 1);
results.Std_Freq_Hz = zeros(numSubjects, 1);

% Power统计
results.Alpha_Mean_Power = zeros(numSubjects, 1);
results.Theta_Mean_Power = zeros(numSubjects, 1);
results.Delta_Mean_Power = zeros(numSubjects, 1);
results.Beta_Mean_Power  = zeros(numSubjects, 1);  % 新增：Beta频段平均Peak_Power（仅统计，不参与主导判断）
results.Power_Ratio_Alpha = zeros(numSubjects, 1);
results.Power_Ratio_Theta = zeros(numSubjects, 1);
results.Power_Ratio_Delta = zeros(numSubjects, 1);

% Aperiodic参数统计（新增）
results.Mean_Aperiodic_Exponent = zeros(numSubjects, 1);
results.Std_Aperiodic_Exponent = zeros(numSubjects, 1);
results.Mean_Aperiodic_Offset = zeros(numSubjects, 1);
results.Std_Aperiodic_Offset = zeros(numSubjects, 1);

% 主导节律判断
results.Dominant_Rhythm_Freq = cell(numSubjects, 1);    % 基于Mean_Freq_Hz范围判断
results.Dominant_Rhythm_Power = cell(numSubjects, 1);   % 保留列名兼容

%% 分析每个患者
fprintf('%-45s | %-6s | %-8s | %-10s | %-10s | %-8s | %-8s\n', ...
    'Subject', 'Valid', 'MeanHz', 'Freq判断', 'Power判断', 'Exponent', 'Offset');
fprintf('%s\n', repmat('-', 1, 110));

for i = 1:numSubjects
    subjectName = subjects{i};
    
    % 提取该患者的所有数据
    subjectData = data(strcmp(data.Subject, subjectName), :);
    totalSegments = height(subjectData);
    
    % 筛选条件：只保留 R_Squared > 0.9（移除 Background_Rhythm 筛选）
    rSquared = subjectData.R_Squared;
    rSquaredMask = rSquared > R_Squared_Threshold;
    excludedByRSquared = sum(~rSquaredMask);

    validData = subjectData(rSquaredMask, :);
    validSegments = height(validData);

    %% 新逻辑：计算 Mean_Freq_Hz
    % 每个segment取Delta/Theta/Alpha中Peak_Power最大者对应的频率，再求均值
    segFreqList = [];  % 收集每个 segment 的代表频率

    if validSegments > 0
        deltaPower = validData.Delta_Peak_Power;
        thetaPower = validData.Theta_Peak_Power;
        alphaPower = validData.Alpha_Peak_Power;
        betaPower  = validData.Beta_Peak_Power;

        deltaFreq  = validData.Delta_Peak_Freq_Hz;
        thetaFreq  = validData.Theta_Peak_Freq_Hz;
        alphaFreq  = validData.Alpha_Peak_Freq_Hz;

        powerMat = [deltaPower, thetaPower, alphaPower];  % N×3，列顺序：Delta/Theta/Alpha
        freqMat  = [deltaFreq,  thetaFreq,  alphaFreq ];  % N×3，与 powerMat 列对应

        for s = 1:validSegments
            rowPower = powerMat(s, :);

            % 三个功率全为 NaN → 跳过
            if all(isnan(rowPower))
                continue
            end

            % 找最大功率对应的频段索引（忽略 NaN）
            [maxPow, maxIdx] = max(rowPower, [], 'omitnan');

            % 最大功率 ≤ 0 或对应频率为 NaN → 跳过
            if maxPow <= 0 || isnan(freqMat(s, maxIdx))
                continue
            end

            segFreqList(end+1) = freqMat(s, maxIdx);  %#ok<AGROW>
        end

        if ~isempty(segFreqList)
            meanFreq = mean(segFreqList);
            stdFreq  = std(segFreqList);
        else
            meanFreq = NaN;
            stdFreq  = NaN;
        end

        %% Peak Power统计（用于功率比例输出）
        alphaPower(isnan(alphaPower)) = 0;
        thetaPower(isnan(thetaPower)) = 0;
        deltaPower(isnan(deltaPower)) = 0;
        betaPower(isnan(betaPower))   = 0;

        alphaMeanPower = mean(alphaPower);
        thetaMeanPower = mean(thetaPower);
        deltaMeanPower = mean(deltaPower);
        betaMeanPower  = mean(betaPower);

        % 计算功率占比
        totalPower = sum(alphaPower) + sum(thetaPower) + sum(deltaPower);
        if totalPower > 0
            alphaRatio = sum(alphaPower) / totalPower * 100;
            thetaRatio = sum(thetaPower) / totalPower * 100;
            deltaRatio = sum(deltaPower) / totalPower * 100;
        else
            alphaRatio = 0; thetaRatio = 0; deltaRatio = 0;
        end

        %% 主导节律判断：直接根据 Mean_Freq_Hz 判断
        if ~isnan(meanFreq)
            if meanFreq >= 1 && meanFreq < 4
                dominantFreq  = 'Delta';
            elseif meanFreq >= 4 && meanFreq < 8
                dominantFreq  = 'Theta';
            elseif meanFreq >= 8 && meanFreq <= 13
                dominantFreq  = 'Alpha';
            else
                dominantFreq  = 'Undetermined';  % 超出 1-13Hz 范围
            end
            dominantPower = dominantFreq;  % 保持列名兼容
        else
            dominantFreq  = 'Undetermined';
            dominantPower = 'Undetermined';
        end

    else
        alphaMeanPower = 0; thetaMeanPower = 0; deltaMeanPower = 0; betaMeanPower = 0;
        alphaRatio = 0; thetaRatio = 0; deltaRatio = 0;
        dominantFreq  = 'Undetermined';
        dominantPower = 'Undetermined';
        meanFreq = NaN;
        stdFreq  = NaN;
    end
    
    %% Aperiodic参数统计（新增，只考虑有效segment）
    if validSegments > 0
        aperiodicExp = validData.Aperiodic_Exponent;
        aperiodicOff = validData.Aperiodic_Offset;
        
        % 排除NaN值进行统计
        aperiodicExp = aperiodicExp(~isnan(aperiodicExp));
        aperiodicOff = aperiodicOff(~isnan(aperiodicOff));
        
        if ~isempty(aperiodicExp)
            meanAperiodicExp = mean(aperiodicExp);
            stdAperiodicExp = std(aperiodicExp);
        else
            meanAperiodicExp = NaN;
            stdAperiodicExp = NaN;
        end
        
        if ~isempty(aperiodicOff)
            meanAperiodicOff = mean(aperiodicOff);
            stdAperiodicOff = std(aperiodicOff);
        else
            meanAperiodicOff = NaN;
            stdAperiodicOff = NaN;
        end
    else
        meanAperiodicExp = NaN;
        stdAperiodicExp = NaN;
        meanAperiodicOff = NaN;
        stdAperiodicOff = NaN;
    end
    
    %% 保存结果
    results.Total_Segments(i) = totalSegments;
    results.Valid_Segments(i) = validSegments;
    results.Excluded_By_R_Squared(i) = excludedByRSquared;
    results.Mean_Freq_Hz(i) = round(meanFreq, 2);
    results.Std_Freq_Hz(i) = round(stdFreq, 2);
    results.Alpha_Mean_Power(i) = round(alphaMeanPower, 4);
    results.Theta_Mean_Power(i) = round(thetaMeanPower, 4);
    results.Delta_Mean_Power(i) = round(deltaMeanPower, 4);
    results.Beta_Mean_Power(i)  = round(betaMeanPower,  4);
    results.Power_Ratio_Alpha(i) = round(alphaRatio, 1);
    results.Power_Ratio_Theta(i) = round(thetaRatio, 1);
    results.Power_Ratio_Delta(i) = round(deltaRatio, 1);
    results.Mean_Aperiodic_Exponent(i) = round(meanAperiodicExp, 4);
    results.Std_Aperiodic_Exponent(i) = round(stdAperiodicExp, 4);
    results.Mean_Aperiodic_Offset(i) = round(meanAperiodicOff, 4);
    results.Std_Aperiodic_Offset(i) = round(stdAperiodicOff, 4);
    results.Dominant_Rhythm_Freq{i} = dominantFreq;
    results.Dominant_Rhythm_Power{i} = dominantPower;
    
    % 打印进度
    fprintf('%-45s | %-6d | %-8.2f | %-10s | %-10s | %-8.4f | %-8.4f\n', ...
        subjectName, validSegments, meanFreq, dominantFreq, dominantPower, ...
        meanAperiodicExp, meanAperiodicOff);
end

%% 输出汇总统计
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('汇总统计\n');
fprintf('%s\n', repmat('=', 1, 70));

fprintf('\nR_Squared筛选信息:\n');
fprintf('  阈值: > %.2f\n', R_Squared_Threshold);
totalExcluded = sum(results.Excluded_By_R_Squared);
totalSegs = sum(results.Total_Segments);
fprintf('  因R_Squared被排除的segment总数: %d / %d (%.1f%%)\n', ...
    totalExcluded, totalSegs, totalExcluded/totalSegs*100);

fprintf('\n基于Mean_Freq_Hz的主导节律分布（1-4Hz:Delta, 4-8Hz:Theta, 8-13Hz:Alpha）:\n');
alphaFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Alpha'));
thetaFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Theta'));
deltaFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Delta'));
undeterFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Undetermined'));
fprintf('  Alpha: %d (%.1f%%)\n', alphaFreq, alphaFreq/numSubjects*100);
fprintf('  Theta: %d (%.1f%%)\n', thetaFreq, thetaFreq/numSubjects*100);
fprintf('  Delta: %d (%.1f%%)\n', deltaFreq, deltaFreq/numSubjects*100);
fprintf('  Undetermined: %d (%.1f%%)\n', undeterFreq, undeterFreq/numSubjects*100);

fprintf('\nAperiodic参数整体统计（所有患者平均值的均值）:\n');
fprintf('  Mean Aperiodic Exponent: %.4f ± %.4f\n', ...
    mean(results.Mean_Aperiodic_Exponent, 'omitnan'), ...
    std(results.Mean_Aperiodic_Exponent, 'omitnan'));
fprintf('  Mean Aperiodic Offset: %.4f ± %.4f\n', ...
    mean(results.Mean_Aperiodic_Offset, 'omitnan'), ...
    std(results.Mean_Aperiodic_Offset, 'omitnan'));

%% 保存结果到Excel
fprintf('\n正在保存结果到: %s\n', outputFile);
writetable(results, outputFile, 'Sheet', 'Results');

% 创建频段定义说明表
bandInfo = table();
bandInfo.Band = {'Delta'; 'Theta'; 'Alpha'};
bandInfo.Freq_Range_Hz = {'1-4'; '4-8'; '8-13'};
writetable(bandInfo, outputFile, 'Sheet', 'Frequency_Bands');

% 创建分析参数说明表
analysisParams = table();
analysisParams.Parameter = {'R_Squared_Threshold'; 'Included_Rhythms'; 'NaN_Handling_Power'; 'Mean_Freq_Hz_Method'; 'Dominant_Rhythm_Method'};
analysisParams.Value = {sprintf('> %.2f', R_Squared_Threshold); 'Alpha, Theta, Delta'; 'Segments with all-NaN power or NaN freq excluded'; 'Per-segment max Peak_Power band freq, then mean across segments'; '1-4Hz:Delta, 4-8Hz:Theta, 8-13Hz:Alpha'};
analysisParams.Description = {'只纳入R_Squared大于此阈值的数据段'; '只考虑这三个频段的Peak_Power和Peak_Freq'; '功率全为NaN或对应频率为NaN的segment被排除'; '每个segment取Delta/Theta/Alpha中Peak_Power最大者对应的频率，再对所有segment求均值'; '直接根据Mean_Freq_Hz的范围判断主导节律'};
writetable(analysisParams, outputFile, 'Sheet', 'Analysis_Parameters');

fprintf('分析完成！\n');
