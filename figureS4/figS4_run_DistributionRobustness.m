%% run_DistributionRobustness.m  (script version)
% Distribution shape comparison between two candidate amplitude features
% for DoC EEG quantification:
%   FullBand_RawPower   vs   Broadband_Amplitude_Mean
%
% For each feature, on the raw distribution:
%   - Skewness                (Fisher; 2000-bootstrap 95% CI)
%   - Dynamic range           max/min ratio, p95/p5 ratio
%   - SD/MAD ratio            std / (1.4826*MAD)  (1 for Gaussian; >1 = heavy tails/outliers)
%
% Just edit the CONFIG block and press Run.
% Requires: Statistics and Machine Learning Toolbox.

clear; clc; close all;

%% ============================ CONFIG ===================================
inputFile = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';   % <-- EDIT
outputDir = 'F:\EEG\LowVoltageResearch\OutputResult\V3\频率波幅方法对比\分别稳健性';    % <-- EDIT
N_BOOT    = 2000;
rng(42);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

featCols = {'FullBand_RawPower','Broadband_Amplitude_Mean'};
featLab  = {'Broadband absolute power(μV²)','Broadband envelope amplitude(μV)'};

% Plot colors: green (#b6d093) / blue (#94cee5)
PALETTE = [182 208 147; 148 206 229]/255;

%% ============================ LOAD =====================================
try
    opts = detectImportOptions(inputFile, 'VariableNamingRule','preserve');
    T = readtable(inputFile, opts);
catch
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    try
        opts = detectImportOptions(inputFile);
        T = readtable(inputFile, opts);
    catch
        T = readtable(inputFile);
    end
end
vn = T.Properties.VariableNames;
featCols = resolveCols(vn, featCols);

%% ============================ LOG ======================================
logFile = fullfile(outputDir,'Distribution_Log.txt');
fid = fopen(logFile,'w');
logb = @(varargin) logBoth(fid, varargin{:});
logb('Distribution analysis\n');
logb('Input : %s\n', inputFile);
logb('N rows: %d   bootstraps: %d\n\n', height(T), N_BOOT);

%% ===================== PER-FEATURE STATISTICS =========================
Res = struct([]);
DATA = cell(numel(featCols),1);   % store cleaned vectors for figures

for i = 1:numel(featCols)
    x  = double(T.(featCols{i}));
    x  = x(~isnan(x));
    DATA{i} = x;
    n  = numel(x);

    % --- distribution shape
    sk = skewness(x, 1);                       % 1 = no bias correction (Fisher)

    % bootstrap CI for skewness
    sb = nan(N_BOOT,1);
    for b = 1:N_BOOT
        idx   = randi(n, n, 1);
        sb(b) = skewness(x(idx), 1);
    end
    skCI = prctile(sb, [2.5 97.5]);

    % --- dynamic range
    mn = min(x); mx = max(x);
    p5  = prctile(x, 5);  p95 = prctile(x, 95);
    if mn > 0, ratio_maxmin = mx/mn; else, ratio_maxmin = NaN; end
    if p5 > 0, ratio_p95p5 = p95/p5; else, ratio_p95p5 = NaN; end

    % --- SD/MAD ratio
    sdv  = std(x);
    madv = 1.4826 * mad(x, 1);                  % 1 = median absolute deviation, scaled
    sd_mad = sdv / madv;

    % --- store
    Res(end+1).feature = featLab{i};            %#ok<SAGROW>
    Res(end).n            = n;
    Res(end).mean         = mean(x);
    Res(end).median       = median(x);
    Res(end).SD           = sdv;
    Res(end).MAD_scaled   = madv;
    Res(end).skewness     = sk;
    Res(end).skew_CI_lo   = skCI(1);
    Res(end).skew_CI_hi   = skCI(2);
    Res(end).min          = mn;
    Res(end).max          = mx;
    Res(end).p5           = p5;
    Res(end).p95          = p95;
    Res(end).ratio_maxmin = ratio_maxmin;
    Res(end).ratio_p95p5  = ratio_p95p5;
    Res(end).SD_MAD       = sd_mad;

    logb('======================================================\n');
    logb('%s  (n=%d)\n', featLab{i}, n);
    logb('======================================================\n');
    logb('  mean=%.4g  median=%.4g  SD=%.4g  MAD(scaled)=%.4g\n', ...
        mean(x), median(x), sdv, madv);
    logb('  Skewness   = %.3f  [95%% CI %.3f, %.3f]\n', sk, skCI(1), skCI(2));
    logb('  Dynamic range:\n');
    logb('    max/min  = %.3g\n', ratio_maxmin);
    logb('    p95/p5   = %.3g\n', ratio_p95p5);
    logb('  SD/MAD ratio = %.3f   (1.0 = Gaussian; higher = heavier tails)\n\n', sd_mad);
end

%% ============================ FIGURE: distributions ===================
% Create and save one independent figure for each feature
for i = 1:numel(featCols)
    x = DATA{i};
    c = PALETTE(i,:);

    % 每个指标单独建立一个图窗
    f1 = figure('Color','w','Position',[100 100 500 380]);
    ax = axes(f1);
    hold(ax,'on');

    histogram(ax, x, ...
        'FaceColor',c, ...
        'EdgeColor',c*0.55, ...
        'LineWidth',1.6, ...
        'FaceAlpha',0.85);

    xline(ax, mean(x), '-', ...
        'Color',[.2 .2 .2], ...
        'LineWidth',1.2, ...
        'Label','mean');

    xlabel(ax, featLab{i},'FontSize',9);
    ylabel(ax,'Number of patients','FontSize',9);
    box(ax,'off');
    hold(ax,'off');

    % 将指标名称转换为合法文件名
    fileTag = regexprep(featLab{i}, '[^A-Za-z0-9_-]', '_');
    baseName = fullfile(outputDir, ...
        ['Fig_Distribution_' fileTag]);

    saveas(f1, [baseName '.fig']);
    print(f1, [baseName '.png'], '-dpng', '-r300');

    % 保存后关闭当前图窗
    close(f1);
end
%% ============================ EXPORT XLSX =============================
xlsxFile = fullfile(outputDir,'Distribution_Results.xlsx');
S1 = struct2table(Res);
writetable(S1, xlsxFile, 'Sheet','Distribution');

logb('\nAll outputs written to: %s\n', outputDir);
fclose(fid);
fprintf('Done. Results in %s\n', outputDir);


%% ======================== LOCAL FUNCTIONS =============================
function logBoth(fid, varargin)
fprintf(varargin{:});
if fid > 0, fprintf(fid, varargin{:}); end
end

function cols = resolveCols(vn, wanted)
cols = cell(size(wanted));
for i = 1:numel(wanted)
    key = regexprep(wanted{i},'[^A-Za-z0-9]','');
    hit = find(strcmpi(regexprep(vn,'[^A-Za-z0-9]',''), key),1);
    if isempty(hit)
        tok = key(1:min(8,numel(key)));
        hit = find(contains(lower(regexprep(vn,'[^A-Za-z0-9]','')),lower(tok)),1);
    end
    if isempty(hit)
        error('Column not found for "%s". Available: %s', wanted{i}, strjoin(vn,', '));
    end
    cols{i} = vn{hit};
end
end
