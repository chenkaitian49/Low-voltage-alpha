%% run_AgreementAnalysis.m  (script version)
% Individual-level agreement between algorithmic and visual EEG
% classification (amplitude & frequency) for DoC patients.
%
% Reference standard = visual grading:
%   Amplitude : high = HA(1)+HTD(2), low = LA(3)+LTD(4)
%   Frequency : high = HA(1)+LA(3),  low = HTD(2)+LTD(4)
%
% For each algorithm: ROC + Youden cut-point -> dichotomize; then Cohen's
% kappa (95% CI by 2000 bootstraps), Po, sensitivity, specificity, PABAK,
% and 2x2 confusion matrix against visual judgement at the Youden cut-point.
%
% Just edit the CONFIG block and press Run.
% Requires: Statistics and Machine Learning Toolbox (perfcurve).

clear; clc; close all;

%% ============================ CONFIG ===================================
inputFile = 'F:\EEG\LowVoltageResearch\OutputResult\V3\频率波幅方法对比\4种分类方法总结果.xlsx';   % <-- EDIT
outputDir = 'F:\EEG\LowVoltageResearch\OutputResult\V3\频率波幅方法对比\结果';                          % <-- EDIT
N_BOOT    = 2000;
rng(42);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

ampCols  = {'Broadband_Amplitude_Mean','FullBand_RawPower','Mean_Aperiodic_Offset'};
freqCols = {'Maximun_Mean_Freq_Hz','Highest_Mean_Freq_Hz','Fastest_Mean_Freq_Hz'};
ampLab   = {'Broadband Envelope Amplitude','Broadband Absolute Power','Aperiodic Offset'};
freqLab  = {'Maximum Freq','Highest-power Freq','Fastest Freq'};

% Plot colors (shared by both figures): blue / green / purple
c_clin  = [148 206 229]/255;   % #94cee5
c_eeg   = [182 208 147]/255;   % #b6d093
c_joint = [200 180 216]/255;   % #c8b4d8
PALETTE = [c_clin; c_eeg; c_joint];

% Confusion-matrix colors mapped to the reference figure, but lighter
% for paper-friendly display.
c_tp = [205 229 241]/255;   % lighter #94cee5
c_fn = [244 210 211]/255;   % lighter #e9a2a3
c_fp = [250 229 205]/255;   % lighter #f5cfa8
c_tn = [220 233 207]/255;   % lighter #b6d093
CM_COLORS = [c_tp; c_fn; c_fp; c_tn];

%% ============================ LOAD =====================================
% Version-robust read (older MATLAB lacks 'VariableNamingRule').
try
    % Newer MATLAB (R2019b+): keep original-ish names, suppress warnings
    opts = detectImportOptions(inputFile, 'VariableNamingRule','preserve');
    T = readtable(inputFile, opts);
catch
    % Older MATLAB: fall back. Suppress the auto-rename warning.
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    try
        opts = detectImportOptions(inputFile);
        T = readtable(inputFile, opts);
    catch
        T = readtable(inputFile);   % last-resort default reader
    end
end
vn = T.Properties.VariableNames;

% Resolve the Group column (first var whose name contains 'HA' and '1',
% robust to header sanitization). Falls back to first numeric column whose
% values are all in {1,2,3,4}.
gIdx = find(contains(vn,'HA') & contains(vn,'1'),1);
if isempty(gIdx)
    for j = 1:width(T)
        col = T{:,j};
        if isnumeric(col)
            u = unique(col(~isnan(col)));
            if ~isempty(u) && all(ismember(u,[1 2 3 4]))
                gIdx = j; break;
            end
        end
    end
end
if isempty(gIdx); gIdx = 1; end
grp = double(T{:,gIdx});

% Resolve feature columns (fuzzy-match against possibly-sanitized headers)
ampCols  = resolveCols(vn, ampCols);
freqCols = resolveCols(vn, freqCols);

% Visual binary labels
ampHigh  = double(ismember(grp,[1 2]));
freqHigh = double(ismember(grp,[1 3]));

%% ============================ LOG ======================================
logFile = fullfile(outputDir,'Agreement_Log.txt');
fid = fopen(logFile,'w');
logb = @(varargin) logBoth(fid, varargin{:});
logb('Individual-level agreement analysis\n');
logb('Input : %s\n', inputFile);
logb('N rows: %d   bootstraps: %d\n\n', height(T), N_BOOT);

%% ===================== ALGORITHM vs VISUAL =============================
domains = {'AMPLITUDE', ampCols, ampLab, ampHigh; ...
           'FREQUENCY', freqCols, freqLab, freqHigh};
Res = struct([]);       % flat results for export
ConfMat = struct([]);   % 2x2 confusion matrices for each method

for dd = 1:size(domains,1)
    name = domains{dd,1}; cols = domains{dd,2};
    labs = domains{dd,3}; yvis = domains{dd,4};
    logb('======================================================\n');
    logb('%s  -- algorithm (dichotomized) vs visual\n', name);
    logb('======================================================\n');

    for i = 1:numel(cols)
        x = double(T.(cols{i}));
        ok = ~isnan(x);
        s = x(ok); y = yvis(ok);

        % orient so higher score -> high class
        [~,~,~,auc] = perfcurve(y, s, 1);
        flip = auc < 0.5;
        if flip, s = -s; end

        % Youden cut-point
        [fpr,tpr,thr] = perfcurve(y, s, 1);
        J = tpr - fpr; [~,jx] = max(J);
        cut = thr(jx);
        pred = double(s >= cut);
        cutReport = cut; if flip, cutReport = -cut; end

        % kappa & friends
        [k,po,sens,spec] = kappaStats(y,pred);
        pabak = 2*po - 1;

        % 2x2 confusion matrix against visual judgement at this Youden cut-point.
        % Rows = visual judgement, columns = quantitative judgement:
        %   [Visual high & Quant high, Visual high & Quant low;
        %    Visual low  & Quant high, Visual low  & Quant low ]
        TP = sum(y==1 & pred==1);
        FN = sum(y==1 & pred==0);
        FP = sum(y==0 & pred==1);
        TN = sum(y==0 & pred==0);
        C2 = [TP FN; FP TN];

        % bootstrap CI
        ks = nan(N_BOOT,1);
        for b = 1:N_BOOT
            idx = randi(numel(y),numel(y),1);
            yb = y(idx); sb = s(idx);
            if numel(unique(yb))<2, continue; end
            [f2,t2,th2] = perfcurve(yb,sb,1);
            [~,j2] = max(t2-f2);
            pb = double(sb >= th2(j2));
            if numel(unique(pb))<2, ks(b)=0; continue; end
            ks(b) = kappaStats(yb,pb);
        end
        ci = prctile(ks(~isnan(ks)),[2.5 97.5]);
        interp = kappaLabel(k);

        Res(end+1).domain = name;            %#ok<SAGROW>
        Res(end).method = labs{i};
        Res(end).cutoff = cutReport;
        Res(end).kappa  = k;
        Res(end).ci_lo  = ci(1);
        Res(end).ci_hi  = ci(2);
        Res(end).po     = po*100;
        Res(end).sens   = sens*100;
        Res(end).spec   = spec*100;
        Res(end).pabak  = pabak;
        Res(end).interp = interp;
        Res(end).n      = numel(y);
        Res(end).TP     = TP;
        Res(end).FN     = FN;
        Res(end).FP     = FP;
        Res(end).TN     = TN;

        ConfMat(end+1).domain = name;          %#ok<SAGROW>
        ConfMat(end).method   = labs{i};
        ConfMat(end).cutoff   = cutReport;
        ConfMat(end).matrix   = C2;
        ConfMat(end).n        = numel(y);
        ConfMat(end).kappa    = k;
        ConfMat(end).po       = po*100;
        ConfMat(end).sens     = sens*100;
        ConfMat(end).spec     = spec*100;
        ConfMat(end).youden   = sens + spec - 1;

        logb('\n  %s (n=%d)  cutoff=%.3f\n', labs{i}, numel(y), cutReport);
        logb('    Cohen kappa = %.3f [95%% CI %.3f-%.3f]  (%s)\n', k, ci(1), ci(2), interp);
        logb('    Po = %.1f%%   PABAK = %.3f\n', po*100, pabak);
        logb('    sensitivity = %.1f%%   specificity = %.1f%%\n', sens*100, spec*100);
        logb('    confusion matrix [visual rows: High Low; quantitative cols: High Low] = [%d %d; %d %d]\n', ...
            C2(1,1), C2(1,2), C2(2,1), C2(2,2));
    end
    logb('\n');
end

%% ============================ FIGURE 1: ROC ============================
f1 = figure('Color','w','Position',[100 100 1000 430]);
colors = PALETTE;
panels = {freqCols, freqLab, freqHigh, 'Frequency'; ...
          ampCols, ampLab, ampHigh, 'Amplitude'};
for pn = 1:2
    cols = panels{pn,1}; labs = panels{pn,2}; yvis = panels{pn,3};
    ax = subplot(1,2,pn); hold(ax,'on');
    hh = gobjects(numel(cols),1); leg = cell(numel(cols),1);
    for i = 1:numel(cols)
        x = double(T.(cols{i})); ok = ~isnan(x);
        s = x(ok); y = yvis(ok);
        [~,~,~,auc] = perfcurve(y,s,1);
        if auc<0.5, s=-s; [~,~,~,auc]=perfcurve(y,s,1); end
        [fpr,tpr,~] = perfcurve(y,s,1);
        hh(i) = plot(ax,fpr,tpr,'Color',colors(i,:),'LineWidth',2.2);
        [~,jx] = max(tpr-fpr);
        plot(ax,fpr(jx),tpr(jx),'o','Color',colors(i,:), ...
            'MarkerFaceColor',colors(i,:),'MarkerSize',7,'MarkerEdgeColor',colors(i,:)*0.6);
        st = lookupRes(Res, panels{pn,4}, labs{i});
        % kappa label above the Youden point, same color as the curve
        text(ax, fpr(jx), tpr(jx)+0.035, sprintf('\\kappa=%.2f', st.kappa), ...
            'Color',colors(i,:)*0.75,'FontSize',8,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','bottom');
        leg{i} = sprintf('%s (AUC=%.2f)', labs{i}, auc);
    end
    plot(ax,[0 1],[0 1],'--','Color',[.6 .6 .6],'LineWidth',0.8);
    xlim(ax,[0 1]); ylim(ax,[0 1]); axis(ax,'square'); box(ax,'off');
    xlabel(ax,'1 - Specificity'); ylabel(ax,'Sensitivity');
    title(ax,panels{pn,4},'FontWeight','bold');
    legend(ax,hh,leg,'Location','southeast','Box','off','FontSize',8);
    hold(ax,'off');
end
saveas(f1, fullfile(outputDir,'Fig1_ROC.fig'));
print(f1, fullfile(outputDir,'Fig1_ROC.png'),'-dpng','-r300');

%% ====================== FIGURE 2: CONFUSION MATRICES =====================
% Six 2x2 confusion matrices at each method's own Youden cut-point.
% Rows = visual judgement; columns = quantitative judgement.
f2 = figure('Color','w','Position',[80 60 980 560]);
plotOrder = {'FREQUENCY','FREQUENCY','FREQUENCY','AMPLITUDE','AMPLITUDE','AMPLITUDE'};
plotMethods = {freqLab{2}, freqLab{3}, freqLab{1}, ...
               ampLab{3},  ampLab{2},  ampLab{1}};
plotTitles = {'Highest Power', 'Fastest Freq', 'Maximum Freq', ...
              'Aperiodic Offset', 'Absolute Power', 'Broadband Envelope Amplitude'};

tl = tiledlayout(f2,2,3,'TileSpacing','none','Padding','none');
for pp = 1:6
    ax = nexttile(tl,pp);
    st = lookupConfMat(ConfMat, plotOrder{pp}, plotMethods{pp});
    plotConfusionCard(ax, st.matrix, CM_COLORS, plotTitles{pp});
end
saveas(f2, fullfile(outputDir,'Fig2_ConfusionMatrices_Youden.fig'));
print(f2, fullfile(outputDir,'Fig2_ConfusionMatrices_Youden.png'),'-dpng','-r300');

%% ============================ EXPORT XLSX ==============================
xlsxFile = fullfile(outputDir,'Agreement_Results.xlsx');

S1 = struct2table(Res);
S1.Properties.VariableNames = {'Domain','Method','Cutoff','Kappa','CI_low', ...
    'CI_high','Po_pct','Sens_pct','Spec_pct','PABAK','Interpretation','N', ...
    'TP_VisualHigh_QuantHigh','FN_VisualHigh_QuantLow', ...
    'FP_VisualLow_QuantHigh','TN_VisualLow_QuantLow'};
writetable(S1, xlsxFile, 'Sheet','Algorithm_vs_Visual');

% Export the six 2x2 matrices in a long, easy-to-copy table.
Crows = {};
for q = 1:numel(ConfMat)
    C = ConfMat(q).matrix;
    Crows(end+1,:) = {ConfMat(q).domain, ConfMat(q).method, ConfMat(q).cutoff, ...
        'Visual high', 'Quant high', C(1,1)}; %#ok<SAGROW>
    Crows(end+1,:) = {ConfMat(q).domain, ConfMat(q).method, ConfMat(q).cutoff, ...
        'Visual high', 'Quant low',  C(1,2)}; %#ok<SAGROW>
    Crows(end+1,:) = {ConfMat(q).domain, ConfMat(q).method, ConfMat(q).cutoff, ...
        'Visual low',  'Quant high', C(2,1)}; %#ok<SAGROW>
    Crows(end+1,:) = {ConfMat(q).domain, ConfMat(q).method, ConfMat(q).cutoff, ...
        'Visual low',  'Quant low',  C(2,2)}; %#ok<SAGROW>
end
S2 = cell2table(Crows, 'VariableNames', ...
    {'Domain','Method','Cutoff','Visual_judgement','Quantitative_judgement','Count'});
writetable(S2, xlsxFile, 'Sheet','Confusion_2x2_Youden');

logb('\nAll outputs written to: %s\n', outputDir);
fclose(fid);
fprintf('Done. Results in %s\n', outputDir);


%% ======================== LOCAL FUNCTIONS =============================
function st = lookupRes(Res, domainTitle, method)
% domainTitle like '(A) Amplitude' -> match 'AMPLITUDE'
if contains(upper(domainTitle),'AMPL'), dom='AMPLITUDE'; else, dom='FREQUENCY'; end
st = struct('kappa',NaN,'ci_lo',NaN,'ci_hi',NaN,'sens',NaN,'spec',NaN,'po',NaN);
for q = 1:numel(Res)
    if strcmp(Res(q).domain,dom) && strcmp(Res(q).method,method)
        st.kappa=Res(q).kappa; st.ci_lo=Res(q).ci_lo; st.ci_hi=Res(q).ci_hi;
        st.sens=Res(q).sens; st.spec=Res(q).spec; st.po=Res(q).po; return;
    end
end
end

function st = lookupConfMat(ConfMat, domain, method)
st = struct('matrix',nan(2,2),'cutoff',NaN,'kappa',NaN,'po',NaN,'n',NaN, ...
    'sens',NaN,'spec',NaN,'youden',NaN);
for q = 1:numel(ConfMat)
    if strcmp(ConfMat(q).domain,domain) && strcmp(ConfMat(q).method,method)
        st.matrix = ConfMat(q).matrix;
        st.cutoff = ConfMat(q).cutoff;
        st.kappa  = ConfMat(q).kappa;
        st.po     = ConfMat(q).po;
        st.n      = ConfMat(q).n;
        st.sens   = ConfMat(q).sens;
        st.spec   = ConfMat(q).spec;
        st.youden = ConfMat(q).youden;
        return;
    end
end
end

function plotConfusionCard(ax, C, cmColors, ttl)
% Draw a denser publication-style confusion-matrix card.
% C = [TP FN; FP TN], rows visual High/Low, columns quantitative High/Low.
axes(ax); %#ok<LAXES>
cla(ax); hold(ax,'on');
axis(ax,[0 1 0 1]);
axis(ax,'off');

N = sum(C(:));
pct = 100 * C / N;

% Title.
text(ax,0.5,0.955,ttl,'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',10,'Color',[0.1 0.1 0.1]);

% Actual number of valid subjects for this method.  The confusion matrix
% contains one entry per subject, so summing all four cells gives the exact
% sample size used for this card.  Place it between the title and the matrix
% to keep the annotation centered and visually separate from cell values.
text(ax,0.5,0.865,sprintf('n=%d',N), ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',8.5,'FontWeight','normal','Color',[0.30 0.30 0.30]);

% Top header labels: move closer to the boxes.
text(ax,0.315,0.745,'Computed high','HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontSize',7.8,'Color',[0.25 0.25 0.25]);
text(ax,0.685,0.745,'Computed low','HorizontalAlignment','center', ...
    'VerticalAlignment','middle','FontSize',7.8,'Color',[0.25 0.25 0.25]);

% Left labels: vertical arrangement and closer to the matrix.
text(ax,0.060,0.555,'Visual high','Rotation',90, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',7.8,'Color',[0.25 0.25 0.25]);
text(ax,0.060,0.275,'Visual low','Rotation',90, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',7.8,'Color',[0.25 0.25 0.25]);

% Main 2x2 rounded boxes: occupy more of each tile.
boxW = 0.345; boxH = 0.255;
x1 = 0.11; x2 = 0.475;
yTop = 0.41; yBot = 0.125;
vals = [C(1,1), C(1,2), C(2,1), C(2,2)];
cellPct = [pct(1,1), pct(1,2), pct(2,1), pct(2,2)];
countTextColors = [0.12 0.39 0.71; 0.77 0.18 0.18; 0.70 0.40 0.05; 0.22 0.55 0.20];
positions = [x1 yTop boxW boxH; x2 yTop boxW boxH; x1 yBot boxW boxH; x2 yBot boxW boxH];
for k = 1:4
    pos = positions(k,:);
    rectangle(ax,'Position',pos,'Curvature',0.10, ...
        'FaceColor',cmColors(k,:),'EdgeColor','none');
    cx = pos(1) + pos(3)/2;
    cy = pos(2) + pos(4)/2;
    text(ax,cx,cy+0.030,sprintf('%d',vals(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',18.5,'FontWeight','bold','Color',countTextColors(k,:));
    text(ax,cx,pos(2)+0.036,sprintf('%.1f%%',cellPct(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',8,'Color',countTextColors(k,:));
end
hold(ax,'off');
end

function logBoth(fid, varargin)
fprintf(varargin{:});
if fid > 0
    fprintf(fid, varargin{:});
end
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

function [k,po,sens,spec] = kappaStats(y,pred)
TP = sum(pred==1 & y==1); TN = sum(pred==0 & y==0);
FP = sum(pred==1 & y==0); FN = sum(pred==0 & y==1);
n = numel(y);
po = (TP+TN)/n;
pe = ((TP+FP)*(TP+FN) + (FN+TN)*(FP+TN)) / n^2;
k  = (po - pe) / (1 - pe);
sens = TP/(TP+FN);
spec = TN/(TN+FP);
end

function lab = kappaLabel(k)
if k < 0,        lab = 'poor';
elseif k < 0.20, lab = 'slight';
elseif k < 0.40, lab = 'fair';
elseif k < 0.60, lab = 'moderate';
elseif k < 0.80, lab = 'substantial';
else,            lab = 'almost perfect';
end
end

function icc = icc_A1(Y)
% ICC(A,1): two-way random, absolute agreement, single measurement
[n,k] = size(Y);
gm  = mean(Y(:));
MSr = k * sum((mean(Y,2)-gm).^2) / (n-1);
MSc = n * sum((mean(Y,1)-gm).^2) / (k-1);
SSt = sum((Y(:)-gm).^2);
SSe = SSt - (n-1)*MSr - (k-1)*MSc;
MSe = SSe / ((n-1)*(k-1));
icc = (MSr - MSe) / (MSr + (k-1)*MSe + k*(MSc-MSe)/n);
end
