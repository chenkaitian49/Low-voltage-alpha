%% 分析EEG背景频率主导节律 (v3 修改版)
% 该脚本读取Excel文件中的Posterior表，判断每个患者总体的主导频率
%
% 【主导节律判断逻辑】
%   1. 对R_Squared > 0.9的有效segment，分别统计Alpha/Theta/Delta/Beta
%      各频段在所有有效segment中出现的比例
%   2. 直接取Alpha/Theta/Delta中占比最高的频段作为主导节律
%   3. 若三者均为0（无任何有效节律），主导节律标记为 'Undetermined'
%
% 【平均频率计算逻辑（新）】
%   仅对含有主导节律的segment计算峰频率均值（而非全部有效segment）
%
% 【Beta节律（新）】
%   对Beta频段（13-30 Hz）同样做比例>=50%的判断，
%   若满足则额外输出 Beta_Mean_Power（基于含Beta的segment），
%   Beta结果独立于主导节律判断（不参与Alpha>Theta>Delta优先级竞争）
%
% 其他保留逻辑：
%   - R_Squared > 0.9 筛选
%   - Aperiodic_Exponent / Aperiodic_Offset 统计
%   - Alpha/Theta/Delta Peak Power 统计（基于各自含有该节律的segment）
%   - 输出结果到新的Excel文件

clear; clc;

%% ========== 参数设置 ==========
inputFile  = 'Test_FOOOF_Results_Segmented.xlsx';
sheetName  = 'Posterior';
outputFile = 'Test_DominantRhythm_PeakPower_Results_v3.xlsx';

% 频段定义 (Hz)
freqBands.Delta = [0.5,  4];
freqBands.Theta = [4,    8];
freqBands.Alpha = [8,   13];
freqBands.Beta  = [13,  30];

% R_Squared 质量阈值
R_Squared_Threshold = 0.9;

% Beta节律判断阈值（Beta_Peak_Power列，与主导节律判断无关）
PROPORTION_THRESHOLD = 0.50;   % 50%（仅用于Beta达标判断）

%% ========== 读取数据 ==========
fprintf('正在读取数据: %s (Sheet: %s)\n', inputFile, sheetName);
data = readtable(inputFile, 'Sheet', sheetName);

subjects    = unique(data.Subject);
numSubjects = length(subjects);
fprintf('共发现 %d 位患者\n', numSubjects);
fprintf('主导节律判断: 取Alpha/Theta/Delta中占比最高者\n');
fprintf('Beta达标阈值: >= %.0f%%\n', PROPORTION_THRESHOLD * 100);
fprintf('R_Squared 阈值: > %.2f\n\n', R_Squared_Threshold);

%% ========== 初始化结果表格 ==========
results = table();
results.Subject = subjects;

% Segment 计数
results.Total_Segments          = zeros(numSubjects, 1);
results.Valid_Segments          = zeros(numSubjects, 1);   % 通过R²筛选的segment数
results.Excluded_By_R_Squared   = zeros(numSubjects, 1);

% 各频段在有效segment中的出现比例（%）
results.Proportion_Alpha_pct    = zeros(numSubjects, 1);
results.Proportion_Theta_pct    = zeros(numSubjects, 1);
results.Proportion_Delta_pct    = zeros(numSubjects, 1);
results.Proportion_Beta_pct     = zeros(numSubjects, 1);

% 主导节律（基于比例+优先级）
results.Dominant_Rhythm         = cell(numSubjects, 1);

% 主导节律对应segment数
results.Dominant_Rhythm_Seg_N   = zeros(numSubjects, 1);

% 峰频率（仅基于含主导节律的segment）
results.Mean_Freq_Hz            = nan(numSubjects, 1);
results.Std_Freq_Hz             = nan(numSubjects, 1);

% Peak Power（各频段，基于含该节律的segment，NaN视为0）
results.Alpha_Mean_Power        = zeros(numSubjects, 1);
results.Theta_Mean_Power        = zeros(numSubjects, 1);
results.Delta_Mean_Power        = zeros(numSubjects, 1);
results.Beta_Mean_Power         = nan(numSubjects, 1);    % 仅比例>=50%时有值

% Beta 是否达标标记
results.Beta_Above_Threshold    = false(numSubjects, 1);

% Aperiodic 参数（基于所有通过R²筛选的有效segment）
results.Mean_Aperiodic_Exponent = nan(numSubjects, 1);
results.Std_Aperiodic_Exponent  = nan(numSubjects, 1);
results.Mean_Aperiodic_Offset   = nan(numSubjects, 1);
results.Std_Aperiodic_Offset    = nan(numSubjects, 1);

%% ========== 逐患者分析 ==========
fprintf('%-45s | %-6s | %-6s | %-6s | %-6s | %-6s | %-12s | %-8s | %-8s\n', ...
    'Subject', 'Valid', 'α%', 'θ%', 'δ%', 'β%', 'Dominant', 'MeanHz', 'BetaPwr');
fprintf('%s\n', repmat('-', 1, 125));

for i = 1:numSubjects
    subjectName = subjects{i};
    subjectData = data(strcmp(data.Subject, subjectName), :);
    totalSegments = height(subjectData);

    % ---- 步骤1：R_Squared 筛选 ----
    rSquared       = subjectData.R_Squared;
    rSquaredMask   = rSquared > R_Squared_Threshold;
    excludedByRSq  = sum(~rSquaredMask);
    validData      = subjectData(rSquaredMask, :);
    validSegments  = height(validData);

    % ---- 步骤2：统计各频段在有效segment中的出现比例 ----
    %
    % Alpha/Theta/Delta：依据 Background_Rhythm 列的标签判断
    %   （每个segment只有一个主标签，三者互斥）
    %
    % Beta：完全独立判断，依据 Beta_Peak_Power 列是否有实际峰值
    %   （非NaN且>0即视为该segment含有Beta节律）
    %   这样即使该segment的Background_Rhythm被标为Alpha，
    %   只要检测到Beta峰，就计入Beta的出现次数
    if validSegments > 0
        bgRhythm = validData.Background_Rhythm;
        if ~iscell(bgRhythm)
            bgRhythm = cellstr(bgRhythm);
        end

        % Alpha/Theta/Delta 基于 Background_Rhythm 标签
        nAlpha = sum(strcmpi(bgRhythm, 'Alpha'));
        nTheta = sum(strcmpi(bgRhythm, 'Theta'));
        nDelta = sum(strcmpi(bgRhythm, 'Delta'));

        propAlpha = nAlpha / validSegments;
        propTheta = nTheta / validSegments;
        propDelta = nDelta / validSegments;

        % Beta 独立判断：Beta_Peak_Power 有实际值（非NaN且>0）的segment数
        betaPowerCol = validData.Beta_Peak_Power;
        betaHasPeak  = ~isnan(betaPowerCol) & betaPowerCol > 0;
        nBeta        = sum(betaHasPeak);
        propBeta     = nBeta / validSegments;
    else
        nAlpha = 0; nTheta = 0; nDelta = 0; nBeta = 0;
        propAlpha = 0; propTheta = 0; propDelta = 0; propBeta = 0;
        betaHasPeak = false(0,1);
    end

    % ---- 步骤3：确定主导节律 ----
    % 直接取Alpha/Theta/Delta中占比最高者作为主导节律
    % 若三者均为0（无任何有效节律），保持 'Undetermined'
    dominantRhythm = 'Undetermined';
    if validSegments > 0
        maxProp = max([propAlpha, propTheta, propDelta]);
        if maxProp > 0
            [~, maxIdx] = max([propAlpha, propTheta, propDelta]);
            rhythmNames = {'Alpha', 'Theta', 'Delta'};
            dominantRhythm = rhythmNames{maxIdx};
        end
        % maxProp == 0：三者全为0，保持 'Undetermined'
    end
    isFallback = false;  % 保留变量以兼容控制台打印逻辑（不再使用）

    % ---- 步骤4：仅对含主导节律的segment计算峰频率 ----
    if ~strcmp(dominantRhythm, 'Undetermined') && validSegments > 0
        dominantMask   = strcmpi(bgRhythm, dominantRhythm);
        dominantData   = validData(dominantMask, :);
        dominantSegN   = height(dominantData);

        domFreq = dominantData.BG_Peak_Freq_Hz;
        domFreq = domFreq(~isnan(domFreq));
        if ~isempty(domFreq)
            meanFreq = mean(domFreq);
            stdFreq  = std(domFreq);
        else
            meanFreq = NaN; stdFreq = NaN;
        end
    else
        dominantSegN = 0;
        meanFreq = NaN; stdFreq = NaN;
    end

    % ---- 步骤5：Peak Power 统计 ----
    % Alpha/Theta/Delta：基于含该节律的segment，NaN视为0
    % Beta：仅当比例>=50%时计算，否则输出NaN
    if validSegments > 0
        % Alpha power（含Alpha的segment）
        alphaMask = strcmpi(bgRhythm, 'Alpha');
        if any(alphaMask)
            ap = validData.Alpha_Peak_Power(alphaMask);
            ap(isnan(ap)) = 0;
            alphaMeanPower = mean(ap);
        else
            alphaMeanPower = 0;
        end

        % Theta power（含Theta的segment）
        thetaMask = strcmpi(bgRhythm, 'Theta');
        if any(thetaMask)
            tp = validData.Theta_Peak_Power(thetaMask);
            tp(isnan(tp)) = 0;
            thetaMeanPower = mean(tp);
        else
            thetaMeanPower = 0;
        end

        % Delta power（含Delta的segment）
        deltaMask = strcmpi(bgRhythm, 'Delta');
        if any(deltaMask)
            dp = validData.Delta_Peak_Power(deltaMask);
            dp(isnan(dp)) = 0;
            deltaMeanPower = mean(dp);
        else
            deltaMeanPower = 0;
        end

        % Beta power（独立判断，仅当Beta出现比例>=50%时计算）
        % 使用 betaHasPeak 掩码（基于Beta_Peak_Power列），与主导节律无关
        betaAboveThreshold = propBeta >= PROPORTION_THRESHOLD;
        if betaAboveThreshold
            bp = validData.Beta_Peak_Power(betaHasPeak);
            bp(isnan(bp)) = 0;
            betaMeanPower = mean(bp);
        else
            betaMeanPower = NaN;
        end
    else
        alphaMeanPower = 0; thetaMeanPower = 0; deltaMeanPower = 0;
        betaMeanPower = NaN; betaAboveThreshold = false;
    end

    % ---- 步骤6：Aperiodic 参数（所有通过R²筛选的有效segment）----
    if validSegments > 0
        expVals = validData.Aperiodic_Exponent;
        offVals = validData.Aperiodic_Offset;
        expVals = expVals(~isnan(expVals));
        offVals = offVals(~isnan(offVals));

        if ~isempty(expVals)
            meanAperExp = mean(expVals); stdAperExp = std(expVals);
        else
            meanAperExp = NaN; stdAperExp = NaN;
        end
        if ~isempty(offVals)
            meanAperOff = mean(offVals); stdAperOff = std(offVals);
        else
            meanAperOff = NaN; stdAperOff = NaN;
        end
    else
        meanAperExp = NaN; stdAperExp = NaN;
        meanAperOff = NaN; stdAperOff = NaN;
    end

    % ---- 保存结果 ----
    results.Total_Segments(i)          = totalSegments;
    results.Valid_Segments(i)          = validSegments;
    results.Excluded_By_R_Squared(i)   = excludedByRSq;
    results.Proportion_Alpha_pct(i)    = round(propAlpha * 100, 1);
    results.Proportion_Theta_pct(i)    = round(propTheta * 100, 1);
    results.Proportion_Delta_pct(i)    = round(propDelta * 100, 1);
    results.Proportion_Beta_pct(i)     = round(propBeta  * 100, 1);
    results.Dominant_Rhythm{i}         = dominantRhythm;
    results.Dominant_Rhythm_Seg_N(i)   = dominantSegN;
    results.Mean_Freq_Hz(i)            = round(meanFreq, 2);
    results.Std_Freq_Hz(i)             = round(stdFreq,  2);
    results.Alpha_Mean_Power(i)        = round(alphaMeanPower, 4);
    results.Theta_Mean_Power(i)        = round(thetaMeanPower, 4);
    results.Delta_Mean_Power(i)        = round(deltaMeanPower, 4);
    results.Beta_Mean_Power(i)         = round(betaMeanPower,  4);
    results.Beta_Above_Threshold(i)    = betaAboveThreshold;
    results.Mean_Aperiodic_Exponent(i) = round(meanAperExp, 4);
    results.Std_Aperiodic_Exponent(i)  = round(stdAperExp,  4);
    results.Mean_Aperiodic_Offset(i)   = round(meanAperOff, 4);
    results.Std_Aperiodic_Offset(i)    = round(stdAperOff,  4);

    % 控制台打印
    betaPwrStr = 'N/A';
    if betaAboveThreshold && ~isnan(betaMeanPower)
        betaPwrStr = sprintf('%.4f', betaMeanPower);
    end
    % 控制台显示时，fallback情况加括号提示，Excel中不含此标记
    displayRhythm = dominantRhythm;
    if isFallback
        displayRhythm = [dominantRhythm, '(fb)'];
    end
    fprintf('%-45s | %-6d | %-6.1f | %-6.1f | %-6.1f | %-6.1f | %-14s | %-8.2f | %s\n', ...
        subjectName, validSegments, ...
        propAlpha*100, propTheta*100, propDelta*100, propBeta*100, ...
        displayRhythm, meanFreq, betaPwrStr);
end

%% ========== 汇总统计 ==========
fprintf('\n%s\n', repmat('=', 1, 75));
fprintf('汇总统计\n');
fprintf('%s\n', repmat('=', 1, 75));

totalSegs    = sum(results.Total_Segments);
totalExclude = sum(results.Excluded_By_R_Squared);
fprintf('\nR_Squared 筛选:\n');
fprintf('  阈值: > %.2f\n', R_Squared_Threshold);
fprintf('  排除 segment: %d / %d (%.1f%%)\n', ...
    totalExclude, totalSegs, totalExclude/totalSegs*100);

fprintf('\n主导节律分布（取Alpha/Theta/Delta中占比最高者）:\n');
fprintf('  （三者均为0时标记为Undetermined）\n');
rhythmList = {'Alpha','Theta','Delta','Undetermined'};
for r = 1:length(rhythmList)
    n = sum(strcmp(results.Dominant_Rhythm, rhythmList{r}));
    fprintf('  %-15s: %d (%.1f%%)\n', rhythmList{r}, n, n/numSubjects*100);
end

fprintf('\nBeta节律（比例>=%.0f%%）:\n', PROPORTION_THRESHOLD*100);
nBetaAbove = sum(results.Beta_Above_Threshold);
fprintf('  Beta达标患者数: %d / %d (%.1f%%)\n', ...
    nBetaAbove, numSubjects, nBetaAbove/numSubjects*100);

fprintf('\nAperiodic 参数整体统计（所有患者均值的均值±SD）:\n');
fprintf('  Exponent: %.4f ± %.4f\n', ...
    mean(results.Mean_Aperiodic_Exponent, 'omitnan'), ...
    std(results.Mean_Aperiodic_Exponent, 'omitnan'));
fprintf('  Offset:   %.4f ± %.4f\n', ...
    mean(results.Mean_Aperiodic_Offset, 'omitnan'), ...
    std(results.Mean_Aperiodic_Offset, 'omitnan'));

%% ========== 保存结果到 Excel ==========
fprintf('\n正在保存结果到: %s\n', outputFile);
writetable(results, outputFile, 'Sheet', 'Results');

% 频段定义说明
bandInfo = table( ...
    {'Delta';'Theta';'Alpha';'Beta'}, ...
    {'0.5-4';'4-8';'8-13';'13-30'}, ...
    'VariableNames', {'Band','Freq_Range_Hz'});
writetable(bandInfo, outputFile, 'Sheet', 'Frequency_Bands');

% 分析参数说明
analysisParams = table( ...
    {'R_Squared_Threshold'; 'Dominant_Logic'; ...
     'Freq_Calc_Scope';    'Power_Calc_Scope';      'NaN_Handling_Power'; 'Beta_Logic'}, ...
    {sprintf('> %.2f', R_Squared_Threshold); ...
     '取Alpha/Theta/Delta中占比最高者'; ...
     '仅含主导节律的segment'; ...
     '含该频段节律的segment（NaN→0）'; ...
     '功率NaN视为0'; ...
     '独立判断，不参与优先级竞争'}, ...
    {'只纳入R²大于阈值的segment'; ...
     '无阈值限制，直接选最高占比频段；三者均0则标记Undetermined'; ...
     '峰频率均值只对含主导节律的segment计算'; ...
     'Power均值只对含该频段节律的segment计算'; ...
     '功率缺失视为该段无峰'; ...
     sprintf('Beta比例>=%d%%时独立输出Beta_Mean_Power', round(PROPORTION_THRESHOLD*100))}, ...
    'VariableNames', {'Parameter','Value','Description'});
writetable(analysisParams, outputFile, 'Sheet', 'Analysis_Parameters');

fprintf('分析完成！\n');