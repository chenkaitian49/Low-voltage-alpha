%% draw_venn_pie.m
% Venn + embedded pie chart for good vs bad prognosis
% F = Frequency (>=8Hz), A = Amplitude (>median), S = Stability (>median)
% Usage: set INPUT_FILE and OUTPUT_DIR, then run.

clear; clc; close all;

%% ========== User config ==========
INPUT_FILE = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
OUTPUT_DIR = 'F:\EEG\LowVoltageResearch\OutputResult\V2\预后韦恩图\调色版本';
%% ==================================

OUTPUT_DIR = strtrim(OUTPUT_DIR);
if OUTPUT_DIR(end)=='\' || OUTPUT_DIR(end)=='/'
    OUTPUT_DIR = OUTPUT_DIR(1:end-1);
end
if ~exist(OUTPUT_DIR,'dir')
    [s,m] = mkdir(OUTPUT_DIR);
    if ~s, error('Cannot create output dir: %s\n%s',OUTPUT_DIR,m); end
end

%% ---------- 1. Load data ----------
if ~exist(INPUT_FILE,'file'), error('File not found: %s',INPUT_FILE); end
fprintf('Reading: %s\n',INPUT_FILE);
[~,~,raw] = xlsread(INPUT_FILE);
headers  = raw(1,:);
findCol  = @(n) find(strcmp(headers,n));

colDiag     = findCol('post Diagnosis ( VS=1 MCS-=2 MCS+=3 EMCS=4)');
colFreq     = findCol('Mean_Freq_Hz');
colAmp      = findCol('Broadband_Amplitude_Mean');
colAlphaSta = findCol('Alpha_Stability');
colDTSta    = findCol('DeltaTheta_Stability');

for ci = 1:5
    cols = {colDiag,colFreq,colAmp,colAlphaSta,colDTSta};
    nms  = {'post Diagnosis','Mean_Freq_Hz','Broadband_Amplitude_Mean','Alpha_Stability','DeltaTheta_Stability'};
    if isempty(cols{ci}), error('Column not found: %s',nms{ci}); end
end

dataRaw   = raw(2:end,:);
safeNum   = @(c) cellfun(@(x) ternary(isnumeric(x)&&isscalar(x),x,NaN), dataRaw(:,c));
diagVals  = safeNum(colDiag);
freqHz    = safeNum(colFreq);
amplitude = safeNum(colAmp);
alphaStab = safeNum(colAlphaSta);
dtStab    = safeNum(colDTSta);
fprintf('Loaded %d rows.\n',size(dataRaw,1));

%% ---------- 2. Stability & filter ----------
nAll = length(freqHz);
stability = NaN(nAll,1);
for i = 1:nAll
    if ~isnan(freqHz(i)) && freqHz(i)>=8
        stability(i) = alphaStab(i);
    else
        stability(i) = dtStab(i);
    end
end
valid     = ~isnan(diagVals)&~isnan(freqHz)&~isnan(amplitude)&~isnan(stability);
diagVals  = diagVals(valid);  freqHz    = freqHz(valid);
amplitude = amplitude(valid); stability = stability(valid);
fprintf('Valid samples: %d\n',sum(valid));

%% ---------- 3. Groups ----------
idxGood = diagVals==3 | diagVals==4;
idxBad  = diagVals==1 | diagVals==2;
fprintf('Good (MCS+/EMCS): %d\nBad (UWS/MCS-): %d\n',sum(idxGood),sum(idxBad));

%% ---------- 3b. Global cutoffs (all valid samples, good+bad combined) ----------
% 与 analyze_8subgroup3.m 保持一致：
%   Freq    → 固定 8 Hz（神经生理学界值）
%   Amp     → 全样本中位数
%   Stability → 全样本中位数
globalMedAmp  = median(amplitude, 'omitnan');
globalMedStab = median(stability,  'omitnan');
fprintf('\n--- Global cutoffs (all valid samples, N=%d) ---\n', sum(valid));
fprintf('  Freq cutoff      : 8 Hz (fixed)\n');
fprintf('  Amplitude median : %.4f\n', globalMedAmp);
fprintf('  Stability median : %.4f\n', globalMedStab);

%% ---------- 4. Stats ----------
fprintf('\n--- Good ---\n');
[vG,piG] = calcGroupStats(freqHz(idxGood), amplitude(idxGood), stability(idxGood), ...
                          globalMedAmp, globalMedStab);
fprintf('\n--- Bad ---\n');
[vB,piB] = calcGroupStats(freqHz(idxBad),  amplitude(idxBad),  stability(idxBad), ...
                          globalMedAmp, globalMedStab);

%% ---------- 5. Venn raster geometry ----------
R      = 0.36;   % circle radius (slightly smaller to leave room for labels & pie)
spread = 0.21;   % center offset

cenF = spread*[cosd(90),  sind(90) ];   % top
cenA = spread*[cosd(210), sind(210)];   % bottom-left
cenS = spread*[cosd(330), sind(330)];   % bottom-right

res = 1000;
lim = 1.0;
[gx,gy] = meshgrid(linspace(-lim,lim,res), linspace(-lim,lim,res));
gx = gx(:); gy = gy(:);

inF = (gx-cenF(1)).^2+(gy-cenF(2)).^2 <= R^2;
inA = (gx-cenA(1)).^2+(gy-cenA(2)).^2 <= R^2;
inS = (gx-cenS(1)).^2+(gy-cenS(2)).^2 <= R^2;

rMask.F    =  inF&~inA&~inS;
rMask.A    = ~inF& inA&~inS;
rMask.S    = ~inF&~inA& inS;
rMask.FA   =  inF& inA&~inS;
rMask.FS   =  inF&~inA& inS;
rMask.AS   = ~inF& inA& inS;
rMask.FAS  =  inF& inA& inS;
rMask.none = ~inF&~inA&~inS;

regionNames = {'none','F','A','S','FA','FS','AS','FAS'};
% maxDark: upper bound for the darkest shade (0=white, 1=pure black).
% Set to 0.80 so the deepest region renders as 80% gray, not pure black.
maxDark  = 0.80;
pct2gray = @(p) [1 1 1]*(1 - p/100 * maxDark);

%% ---------- 6. Draw combined Venn+Pie figures ----------
groupLabels  = {'Prognosis (MCS+/EMCS)', 'Prognosis (UWS/MCS-)'};
groupN       = {sum(idxGood), sum(idxBad)};
vennData     = {vG, vB};
piData       = {piG, piB};
bgColorsRGB  = {[233 162 163]/255, [148 206 229]/255};

% Pie color palettes — harmonize with background
% Good (warm reddish bg): warm browns/oranges like reference figure
goodPieC = {[0.94 0.94 0.92];   % FN=0 light cream
            [0.80 0.62 0.45];   % FN=1 tan
            [0.65 0.35 0.18];   % FN=2 brown-orange
            [0.38 0.18 0.08]};  % FN=3 dark brown
% Bad  (cool bluish bg): cool blue-greys
badPieC  = {[0.92 0.94 0.96];   % FN=0 near white
            [0.60 0.72 0.84];   % FN=1 light blue
            [0.28 0.46 0.65];   % FN=2 mid blue
            [0.12 0.22 0.40]};  % FN=3 dark navy
pieColorSets = {goodPieC, badPieC};

vennFnames = {'venn_good.png', 'venn_bad.png'};

% Circle outline points
tOut  = linspace(0,2*pi,500);

% F/A/S label positions: just outside each circle along its radial direction
dimOff    = R + 0.1;   % radial offset from origin
dimAngles = [90, 210, 330];
dimLabels = {'Freq', 'Amp', 'Stab'};

% Pie inset parameters (in axis units)
pieCenter = [0.58, 0.72];    % upper-right corner
pieR      = 0.18;            % pie radius
fnLabOff  = pieR + 0.05;     % label radial offset from pie center

% ---- Global max percentage for normalization (across BOTH groups) ----
% All 14 region values (7 regions × 2 groups) share one normalization denominator
vennRegions7 = {'F','A','S','FA','FS','AS','FAS'};
allPctsGlobal = [];
for g = 1:2
    vd = vennData{g};
    for ri = 1:7
        allPctsGlobal(end+1) = vd.(vennRegions7{ri}); %#ok<AGROW>
    end
end
globalMaxPct = max(allPctsGlobal);
fprintf('Global max pct (both groups combined) = %.2f%%\n', globalMaxPct);

for g = 1:2
    vd  = vennData{g};
    pd  = piData{g};
    bgC = bgColorsRGB{g};
    cols= pieColorSets{g};
    maxPct = globalMaxPct;   % global normalization denominator (both groups combined)

    fig = figure('Color',bgC,'Position',[100 100 560 600]); % taller for bottom margin

    % Axes with bottom margin: leave ~8% at bottom
    ax = axes('Position',[0 0.08 1 0.92]);
    hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
    xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);

    % --- Solid background ---
    patch(ax,[-lim lim lim -lim],[-lim -lim lim lim],bgC,'EdgeColor','none');

    % --- Build RGB image for Venn (normalized per group) ---
    % pct2gray maps 0%->white(1,1,1), 100%->black(0,0,0)
    % After normalization: region/maxPct * 100 so max region = 100% = black
    imgR = ones(res,res); imgG = ones(res,res); imgB = ones(res,res);
    for ri = 1:length(regionNames)
        rn  = regionNames{ri};
        m   = rMask.(rn);
        pct = vd.(rn);
        % Normalize: scale so that maxPct maps to 100% (darkest)
        if strcmp(rn, 'none')
            normPct = 0;   % background outside circles stays white
        else
            normPct = pct / maxPct * 100;
        end
        c   = pct2gray(normPct);
        mImg = reshape(m,res,res);
        imgR(mImg) = c(1); imgG(mImg) = c(2); imgB(mImg) = c(3);
    end
    % background pixels get panel bg color
    nImg = reshape(rMask.none,res,res);
    imgR(nImg)=bgC(1); imgG(nImg)=bgC(2); imgB(nImg)=bgC(3);

    rgbImg = cat(3,imgR,imgG,imgB);
    image(ax, linspace(-lim,lim,res), linspace(-lim,lim,res), rgbImg);
    axis(ax,'xy');

    % --- Circle outlines ---
    cens = {cenF,cenA,cenS};
    for i = 1:3
        plot(ax, R*cos(tOut)+cens{i}(1), R*sin(tOut)+cens{i}(2), ...
            'Color',[0.15 0.15 0.15],'LineWidth',1.5);
    end

    % --- F / A / S labels: outside each circle at its corner ---
    for i = 1:3
        lx = dimOff*cosd(dimAngles(i));
        ly = dimOff*sind(dimAngles(i));
        text(ax,lx,ly,dimLabels{i}, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',15,'FontWeight','bold','Color',[0.08 0.08 0.08]);
    end

    % --- Group label in upper-left (inside panel, over background) ---
    text(ax,-0.90,0.92,groupLabels{g}, ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',14,'FontWeight','bold','Color',[0.12 0.12 0.12]);

    % --- Embedded pie chart in lower-right ---
    total    = sum(pd);
    pcts     = pd/total*100;
    startAng = pi/2;   % start at top, clockwise

    for i = 1:4
        if pd(i)==0, continue; end
        sweep   = pcts(i)/100*2*pi;
        endAng  = startAng - sweep;
        theta   = linspace(startAng,endAng,200);
        xPoly   = [pieCenter(1), pieR*cos(theta)+pieCenter(1), pieCenter(1)];
        yPoly   = [pieCenter(2), pieR*sin(theta)+pieCenter(2), pieCenter(2)];
        patch(ax,xPoly,yPoly,cols{i},'EdgeColor','white','LineWidth',0.8);

        % FN label outside slice
        midAng = (startAng+endAng)/2;
        lx = pieCenter(1) + fnLabOff*cos(midAng);
        ly = pieCenter(2) + fnLabOff*sin(midAng);
        text(ax,lx,ly,sprintf('FN=%d',i-1), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',7.5,'Color',[0.10 0.10 0.10]);

        startAng = endAng;
    end

    % Thin border around pie
    plot(ax, pieR*cos(tOut)+pieCenter(1), pieR*sin(tOut)+pieCenter(2), ...
        'Color',[0.3 0.3 0.3],'LineWidth',0.8);

    % --- Save ---
    outPath = fullfile(OUTPUT_DIR,vennFnames{g});
    exportgraphics(fig,outPath,'Resolution',300);
    fprintf('Saved: %s\n',outPath);
    close(fig);
end

%% ---------- 7. Grayscale legend (single, shared by both groups) ----------
% X-axis max = globalMaxPct, so the darkest shade always corresponds to
% the actual proportion of the single densest block across both groups.
nTicks   = 6;
tickVals = linspace(0, globalMaxPct, nTicks);

fig = figure('Color','white','Position',[300 300 1120 160]);
ax  = axes('Position',[0.04 0.45 0.92 0.15]);
gradImg = repmat(linspace(1, 1-maxDark, 512), 20, 1);  % cap matches maxDark
imagesc(ax, linspace(0,globalMaxPct,512), 1, gradImg);
colormap(ax, gray(512));
ax.XLim = [0 globalMaxPct];
ax.XTick = tickVals;
ax.XTickLabel = arrayfun(@(x) sprintf('%.1f%%',x), tickVals, 'UniformOutput', false);
ax.XAxisLocation='bottom'; ax.YTick=[]; ax.FontSize=10; ax.Box='on';
xlabel(ax, sprintf('Proportion of patients in group (%%)', globalMaxPct), 'FontSize',11);
outPath = fullfile(OUTPUT_DIR, 'legend_grayscale.png');
exportgraphics(fig, outPath, 'Resolution',300);
fprintf('Saved: %s\n', outPath);
close(fig);

fprintf('\nDone! 3 figures saved to: %s\n',OUTPUT_DIR);

%% ---------- 8. Export statistics to Excel ----------
xlsPath = fullfile(OUTPUT_DIR, 'venn_statistics.xlsx');

groupLabelsXLS = {'Good Prognosis (MCS+/EMCS)', 'Bad Prognosis (UWS/MCS-)'};
groupNums      = {sum(idxGood), sum(idxBad)};
vennDataAll    = {vG, vB};
piDataAll      = {piG, piB};
regionOrder    = {'F','A','S','FA','FS','AS','FAS','none'};

% Build data table as cell array
allRows = {'Group', 'n', 'Region', 'Count', 'Percentage (%)'};

% 在表头后记录全局分界值
allRows(end+1,:) = {'[Global cutoffs]', '', 'Freq >= 8 Hz (fixed)',      '', ''};
allRows(end+1,:) = {'[Global cutoffs]', '', 'Amplitude median',  round(globalMedAmp,4),  ''};
allRows(end+1,:) = {'[Global cutoffs]', '', 'Stability median',  round(globalMedStab,4), ''};
allRows(end+1,:) = {'','','','',''};  % blank separator

for g = 1:2
    gLabel = groupLabelsXLS{g};
    gN     = groupNums{g};
    vd     = vennDataAll{g};
    pd     = piDataAll{g};

    % Venn regions
    for ri = 1:length(regionOrder)
        rn  = regionOrder{ri};
        pct = vd.(rn);
        cnt = round(pct/100 * gN);
        allRows(end+1,:) = {gLabel, gN, ['Venn_' rn], cnt, round(pct,2)};
    end

    % FN distribution
    total = sum(pd);
    for i = 1:4
        allRows(end+1,:) = {gLabel, gN, sprintf('FN=%d',i-1), pd(i), round(pd(i)/total*100,2)};
    end

    % Medians
    if g == 1
        freq_g = freqHz(idxGood); amp_g = amplitude(idxGood); stab_g = stability(idxGood);
    else
        freq_g = freqHz(idxBad);  amp_g = amplitude(idxBad);  stab_g = stability(idxBad);
    end
    allRows(end+1,:) = {gLabel, gN, 'Median_Amplitude(global)',  round(globalMedAmp,  4), ''};
    allRows(end+1,:) = {gLabel, gN, 'Median_Stability(global)',  round(globalMedStab, 4), ''};
    allRows(end+1,:) = {gLabel, gN, 'N_HighFreq(>=8Hz)', sum(freq_g>=8), round(sum(freq_g>=8)/gN*100,2)};
    allRows(end+1,:) = {'','','','',''};  % blank separator row
end

xlswrite(xlsPath, allRows, 'Statistics');
fprintf('Statistics saved: %s\n', xlsPath);

%% ========== Local functions ==========

function result = ternary(cond,a,b)
    if cond; result=a; else; result=b; end
end

function [vennPct,pieCnt] = calcGroupStats(freq,amp,stab,medAmp,medStab)
    % medAmp / medStab 由外部传入（全样本中位数），不在组内重新计算
    n       = length(freq);
    fprintf('  Amplitude cutoff=%.4f  Stability cutoff=%.4f  (global medians)\n',medAmp,medStab);
    HF = freq>=8;
    HA = amp >medAmp;
    HS = stab>medStab;
    vennPct.F    = sum( HF&~HA&~HS)/n*100;
    vennPct.A    = sum(~HF& HA&~HS)/n*100;
    vennPct.S    = sum(~HF&~HA& HS)/n*100;
    vennPct.FA   = sum( HF& HA&~HS)/n*100;
    vennPct.FS   = sum( HF&~HA& HS)/n*100;
    vennPct.AS   = sum(~HF& HA& HS)/n*100;
    vennPct.FAS  = sum( HF& HA& HS)/n*100;
    vennPct.none = sum(~HF&~HA&~HS)/n*100;
    FN = double(HF)+double(HA)+double(HS);
    pieCnt = [sum(FN==0),sum(FN==1),sum(FN==2),sum(FN==3)];
end