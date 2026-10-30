%% 分析EEG背景频率主导节律 (修改版)
% 该脚本读取Excel文件中的Posterior表，判断每个患者总体的主导频率
% 方法：基于平均频率判断所属频段
% 只考虑Alpha、Theta、Delta三个频段，不考虑Beta
% 修改内容：
%   1. 增加R_Squared > 0.9的筛选条件
%   2. 增加Aperiodic_Exponent和Aperiodic_Offset的平均值输出
% 输出结果到新的Excel文件

clear; clc;

%% 设置参数
inputFile = 'lowThetaDelta_FOOOF_Results_Segmented.xlsx';  % 输入文件
sheetName = 'Posterior';                                % 分析的表
outputFile = 'lowThetaDelta_DominantRhythm_Fast_Results_v3.xlsx';            % 输出文件

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
results.Power_Ratio_Alpha = zeros(numSubjects, 1);
results.Power_Ratio_Theta = zeros(numSubjects, 1);
results.Power_Ratio_Delta = zeros(numSubjects, 1);

% Aperiodic参数统计（新增）
results.Mean_Aperiodic_Exponent = zeros(numSubjects, 1);
results.Std_Aperiodic_Exponent = zeros(numSubjects, 1);
results.Mean_Aperiodic_Offset = zeros(numSubjects, 1);
results.Std_Aperiodic_Offset = zeros(numSubjects, 1);

% 主导节律判断
results.Dominant_Rhythm_Freq = cell(numSubjects, 1);    % 基于平均频率
results.Dominant_Rhythm_Power = cell(numSubjects, 1);   % 基于平均功率

%% 分析每个患者
fprintf('%-45s | %-6s | %-8s | %-10s | %-10s | %-8s | %-8s\n', ...
    'Subject', 'Valid', 'MeanHz', 'Freq判断', 'Power判断', 'Exponent', 'Offset');
fprintf('%s\n', repmat('-', 1, 110));

for i = 1:numSubjects
    subjectName = subjects{i};
    
    % 提取该患者的所有数据
    subjectData = data(strcmp(data.Subject, subjectName), :);
    totalSegments = height(subjectData);
    
    % 获取背景节律列、频率列和R_Squared列
    bgRhythm = subjectData.Background_Rhythm;
    bgFreq = subjectData.BG_Peak_Freq_Hz;
    rSquared = subjectData.R_Squared;
    
    % 筛选条件1：R_Squared > 0.9
    rSquaredMask = rSquared > R_Squared_Threshold;
    excludedByRSquared = sum(~rSquaredMask);
    
    % 筛选条件2：Background_Rhythm为Alpha/Theta/Delta
    if iscell(bgRhythm)
        rhythmMask = strcmpi(bgRhythm, 'Alpha') | strcmpi(bgRhythm, 'Theta') | strcmpi(bgRhythm, 'Delta');
    else
        rhythmMask = ismember(bgRhythm, {'Alpha', 'Theta', 'Delta'});
    end
    
    % 综合筛选：同时满足两个条件
    validMask = rSquaredMask & rhythmMask;
    validData = subjectData(validMask, :);
    validSegments = height(validData);
    
    %% 频率统计（只考虑有效segment）
    validFreq = validData.BG_Peak_Freq_Hz;
    validFreq = validFreq(~isnan(validFreq));
    if ~isempty(validFreq)
        meanFreq = mean(validFreq);
        stdFreq = std(validFreq);
    else
        meanFreq = NaN;
        stdFreq = NaN;
    end
    
    %% 基于平均频率判断主导节律
    if ~isnan(meanFreq)
        dominantFreq = classifyFreqBand(meanFreq, freqBands);
    else
        dominantFreq = 'Undetermined';
    end
    
    %% Peak Power统计（只考虑有效segment，NaN视为0）
    if validSegments > 0
        alphaPower = validData.Alpha_Peak_Power;
        thetaPower = validData.Theta_Peak_Power;
        deltaPower = validData.Delta_Peak_Power;
        
        % NaN视为0
        alphaPower(isnan(alphaPower)) = 0;
        thetaPower(isnan(thetaPower)) = 0;
        deltaPower(isnan(deltaPower)) = 0;
        
        alphaMeanPower = mean(alphaPower);
        thetaMeanPower = mean(thetaPower);
        deltaMeanPower = mean(deltaPower);
        
        % 计算功率占比
        totalPower = sum(alphaPower) + sum(thetaPower) + sum(deltaPower);
        if totalPower > 0
            alphaRatio = sum(alphaPower) / totalPower * 100;
            thetaRatio = sum(thetaPower) / totalPower * 100;
            deltaRatio = sum(deltaPower) / totalPower * 100;
        else
            alphaRatio = 0; thetaRatio = 0; deltaRatio = 0;
        end
        
        % 基于平均功率判断主导频段
        rhythmNames = {'Alpha', 'Theta', 'Delta'};
        powerValues = [alphaMeanPower, thetaMeanPower, deltaMeanPower];
        [~, maxPowerIdx] = max(powerValues);
        if max(powerValues) > 0
            dominantPower = rhythmNames{maxPowerIdx};
        else
            dominantPower = 'Undetermined';
        end
    else
        alphaMeanPower = 0; thetaMeanPower = 0; deltaMeanPower = 0;
        alphaRatio = 0; thetaRatio = 0; deltaRatio = 0;
        dominantPower = 'Undetermined';
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

fprintf('\n基于平均频率的主导节律分布:\n');
alphaFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Alpha'));
thetaFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Theta'));
deltaFreq = sum(strcmpi(results.Dominant_Rhythm_Freq, 'Delta'));
fprintf('  Alpha: %d (%.1f%%)\n', alphaFreq, alphaFreq/numSubjects*100);
fprintf('  Theta: %d (%.1f%%)\n', thetaFreq, thetaFreq/numSubjects*100);
fprintf('  Delta: %d (%.1f%%)\n', deltaFreq, deltaFreq/numSubjects*100);

fprintf('\n基于平均功率的主导节律分布:\n');
alphaPwr = sum(strcmpi(results.Dominant_Rhythm_Power, 'Alpha'));
thetaPwr = sum(strcmpi(results.Dominant_Rhythm_Power, 'Theta'));
deltaPwr = sum(strcmpi(results.Dominant_Rhythm_Power, 'Delta'));
fprintf('  Alpha: %d (%.1f%%)\n', alphaPwr, alphaPwr/numSubjects*100);
fprintf('  Theta: %d (%.1f%%)\n', thetaPwr, thetaPwr/numSubjects*100);
fprintf('  Delta: %d (%.1f%%)\n', deltaPwr, deltaPwr/numSubjects*100);

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
bandInfo.Freq_Range_Hz = {'0.5-4'; '4-8'; '8-13'};
writetable(bandInfo, outputFile, 'Sheet', 'Frequency_Bands');

% 创建分析参数说明表（新增）
analysisParams = table();
analysisParams.Parameter = {'R_Squared_Threshold'; 'Included_Rhythms'; 'NaN_Handling_Power'};
analysisParams.Value = {sprintf('> %.2f', R_Squared_Threshold); 'Alpha, Theta, Delta'; 'Treated as 0'};
analysisParams.Description = {'只纳入R_Squared大于此阈值的数据段'; '只考虑这三个频段的背景节律'; '功率值为NaN时视为0处理'};
writetable(analysisParams, outputFile, 'Sheet', 'Analysis_Parameters');

fprintf('分析完成！\n');

%% 辅助函数：根据频率分类频段
function band = classifyFreqBand(freq, freqBands)
    if freq >= freqBands.Delta(1) && freq < freqBands.Theta(1)
        band = 'Delta';
    elseif freq >= freqBands.Theta(1) && freq < freqBands.Alpha(1)
        band = 'Theta';
    elseif freq >= freqBands.Alpha(1) && freq < 13
        band = 'Alpha';
    else
        band = 'OutOfRange';
    end
end
