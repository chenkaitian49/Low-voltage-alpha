%% =========================================================
%  run_ModelEvaluation_AllInOne.m
%
%  两模型框架（M_eeg / M_joint）完整评估，一个脚本搞定：
%    ★ bootstrap 只跑一次 ★ —— 同一组重采样里同时算出
%       AUC 的 CI、NRI/IDI 的 CI、校准 optimism 校正，
%    然后分别绘制【三张独立的图】（各自单独窗口、单独存盘，不拼大图）：
%       图1 区分度   discrimination_2models.png/.pdf
%       图2 临床效用 dca_nri_idi_2models.png/.pdf
%       图3 校准     calibration_2models.png/.pdf
%    并各自写出统计日志 txt。
%    另外导出 EEG 模型 OR/95%CI/Wald P（每1 SD及原始单位）csv + txt。
%
%  队列：全程只用剔除基线 MCS+/EMCS（pre Diagnosis>=3）后的患者（87 例）。
%
%  依赖：threeModelCommon.m（须在同目录或 path 中）
%  ★ 仅需修改「用户配置区」★
%% =========================================================

clc; clear; close all;
C = threeModelCommon();

%% ── 用户配置区 ───────────────────────────────────────────
cfg.INPUT_FILE = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
OUTPUT_DIR     = 'F:\EEG\LowVoltageResearch\OutputResult\V3\3模型对比\2model_result';

cfg.COL_POST_DIAG = 'post Diagnosis ( VS=1 MCS-=2 MCS+=3 EMCS=4)';
cfg.COL_PRE_DIAG  = 'pre Diagnosis ( VS=1 MCS-=2 MCS+=3 EMCS=4)';
cfg.COL_AGE       = 'Age';
cfg.COL_ETIOLOGY  = 'Etiology (trauma=1 stroke=2 ischemia=3)';
cfg.COL_DURATION  = 'Duration';
cfg.COL_PRE_CRS   = 'pre CRS-R';
cfg.COL_VAR1      = 'Mean_Freq_Hz';
cfg.COL_VAR2      = 'Broadband_Amplitude_Mean';
cfg.COL_ALPHA_STAB= 'Alpha_Stability';
cfg.COL_DELTA_STAB= 'DeltaTheta_Stability';
cfg.COL_NAME_EN   = 'Name';

cfg.GRP_THR  = 3;        % post Diagnosis >= 3 → 阳性(MCS+/EMCS)
PRE_EXCL_THR = 3;        % pre Diagnosis >= 3 被剔除
N_BOOT       = 2000;     % ★ 唯一的 bootstrap 次数（三图共用）
N_BINS       = 5;        % 校准/HL 分箱数
NRI_THR      = [0.2 0.4 0.6];   % 分类 NRI 风险阈值

% ── 统一视觉风格 ─────────────────────────────────────────
c_eeg   = [148 206 229]/255;   % #94cee5 蓝
c_joint = [182 208 147]/255;   % #b6d093 绿
mcol = {c_eeg, c_joint};
mnames = {'M_{eeg}','M_{joint}'};
mlog   = {'M_eeg','M_joint'};
%% ─────────────────────────────────────────────────────────

if ~exist(OUTPUT_DIR,'dir'); mkdir(OUTPUT_DIR); end

%% =========================================================
%  1. 数据 + 主拟合（只做一次）
%% =========================================================
D = C.loadData(cfg.INPUT_FILE, cfg);
idx = ~(D.pre >= PRE_EXCL_THR);
y   = D.grp(idx);
N   = sum(idx);
fprintf('队列（剔除基线 MCS+/EMCS）: n=%d (阳性=%d, 阴性=%d)\n', N, sum(y==1), sum(y==0));

[~, Xe, Xj, dinfo] = C.buildDesign(D, idx);
Xsets = {Xe, Xj};
mfeat = {dinfo.eegNames, dinfo.jointNames};

fits = cell(1,2); P = cell(1,2);
for m = 1:2
    fits{m} = C.fitLogit(Xsets{m}, y);
    P{m}    = fits{m}.prob;
end
pE = P{1}; pJ = P{2};

%% =========================================================
%  1.1 新增：导出 EEG 模型 OR / 95%CI / Wald P（量纲核查）
%  用于核对 Supplementary Table 4；不改变模型、预测值或后续绘图。
%  threeModelCommon.fitLogit 对 EEG 变量先做 z-score 标准化。
%  统一使用正态近似 Wald 95%CI：beta +/- norminv(0.975)*SE，再取 exp。
%  这与 fitglm 的 Wald P 使用同一个标准正态近似，避免 coefCI 的 t 临界值造成不一致。
%  另外把同一模型的结果还原为“每原始单位”，用于查旧表的量纲。
%  输出文件：eeg_model_OR_95CI.csv / eeg_model_OR_95CI.txt
%% =========================================================
mdlEEG  = fits{1}.mdl;                         % M_eeg: 三项 EEG 特征同时入模
coefEEG = mdlEEG.Coefficients;
alphaEEG = 0.05;                               % 双侧95%置信区间
zCritEEG = norminv(1 - alphaEEG/2);             % ≈1.959964，正态分布临界值
nEEG    = numel(dinfo.eegNames);
rowsEEG = (2 : nEEG+1)';                       % 第1行为截距，之后对应三项 EEG 特征

assert(size(Xe,2) == nEEG, 'EEG 特征列数与名称不一致，请检查 buildDesign。');
assert(numel(fits{1}.sg) == nEEG, '标准化 SD 数量与 EEG 特征不一致。');

metricEEG = string(dinfo.eegNames(:));
sdEEG     = fits{1}.sg(:);                     % 训练/开发建模队列使用的标准差
betaEEG   = coefEEG.Estimate(rowsEEG);
seEEG     = coefEEG.SE(rowsEEG);
pEEG      = coefEEG.pValue(rowsEEG);           % fitglm 输出的 Wald P

% 正态近似 Wald 95%CI：与双侧正态 Wald P 使用一致的临界值
ciBetaEEG = [betaEEG - zCritEEG .* seEEG, ...
             betaEEG + zCritEEG .* seEEG];

% 以 z=beta/SE 重新计算双侧 P，仅核对 fitglm 的 P 是否采用相同口径
waldZEEG   = betaEEG ./ seEEG;
pEEGcheck  = 2 .* normcdf(-abs(waldZEEG));
if any(abs(pEEGcheck - pEEG) > 1e-6)
    warning('fitglm 给出的 Wald P 与标准正态近似不一致，请检查模型与CI计算方法。');
end

% 每增加 1 SD 的效应 —— 若表注写“per one-standard-deviation increase”，用这些数值
orSD   = exp(betaEEG);
ciLoSD = exp(ciBetaEEG(:,1));
ciHiSD = exp(ciBetaEEG(:,2));

% 同一模型换算为每增加 1 个原始单位 —— 仅用于核对旧版 OR 量纲
orUnit   = exp(betaEEG ./ sdEEG);
ciLoUnit = exp(ciBetaEEG(:,1) ./ sdEEG);
ciHiUnit = exp(ciBetaEEG(:,2) ./ sdEEG);

OR_table = table(metricEEG, sdEEG, betaEEG, seEEG, orSD, ciLoSD, ciHiSD, ...
    pEEG, orUnit, ciLoUnit, ciHiUnit, ...
    'VariableNames', {'EEG_Metric','SD_in_Modelling_Cohort', ...
    'Beta_per1SD','SE_per1SD','OR_per1SD','CI95_L_per1SD','CI95_U_per1SD', ...
    'Wald_P','OR_perRawUnit','CI95_L_perRawUnit','CI95_U_perRawUnit'});

orCsvPath = fullfile(OUTPUT_DIR, 'eeg_model_OR_95CI.csv');
orTxtPath = fullfile(OUTPUT_DIR, 'eeg_model_OR_95CI.txt');
writetable(OR_table, orCsvPath);

fidOR = fopen(orTxtPath, 'w');
if fidOR < 0
    error('无法创建 OR 输出日志: %s', orTxtPath);
end
fprintf(fidOR, 'M_eeg: multivariable logistic regression (n=%d)\n', N);
fprintf(fidOR, 'Outcome: follow-up MCS+/EMCS versus UWS/MCS-\n');
fprintf(fidOR, 'All three EEG features are in the same model.\n');
fprintf(fidOR, '95%% CI method: two-sided normal-approximation Wald, z critical = %.8f.\n', zCritEEG);
fprintf(fidOR, 'Per 1 SD: OR = exp(beta_z), CI = exp(beta_z +/- zcrit*SE_z).\n');
fprintf(fidOR, 'Per raw unit: OR = exp(beta_z/SD), CI = exp((beta_z +/- zcrit*SE_z)/SD).\n');
fprintf(fidOR, 'Wald P: fitglm coefficient P, cross-checked against two-sided normal z=beta_z/SE_z.\n\n');
for k = 1:nEEG
    fprintf(fidOR, '%s\n', char(metricEEG(k)));
    fprintf(fidOR, '  SD (modelling cohort) = %.6f\n', sdEEG(k));
    fprintf(fidOR, '  beta_z = %.6f; SE_z = %.6f; Wald z = %.6f\n', ...
        betaEEG(k), seEEG(k), waldZEEG(k));
    fprintf(fidOR, '  Wald P (fitglm) = %.6f; P (normal check) = %.6f\n', ...
        pEEG(k), pEEGcheck(k));
    fprintf(fidOR, '  OR per 1 SD = %.6f (%.6f-%.6f)\n', ...
        orSD(k), ciLoSD(k), ciHiSD(k));
    fprintf(fidOR, '  OR per raw unit = %.6f (%.6f-%.6f)\n\n', ...
        orUnit(k), ciLoUnit(k), ciHiUnit(k));
end
fclose(fidOR);

fprintf('\n===== M_eeg: OR and 95%%CI, per 1 SD and per raw unit =====\n');
disp(OR_table);
fprintf('OR 明细已保存：\n  %s\n  %s\n', orCsvPath, orTxtPath);

% 主拟合点估计
auc0 = zeros(1,2);
for m = 1:2
    [~,~,~,auc0(m)] = C.computeROC(P{m}, y);
end
R_joint = C.nri_idi(pE, pJ, y, NRI_THR);
slope_app = zeros(1,2); int_app = zeros(1,2);
for m = 1:2
    [int_app(m), slope_app(m)] = calSlopeIntercept(P{m}, y);
end
HL = zeros(1,2); HLdof = zeros(1,2); HLp = zeros(1,2);
for m = 1:2
    [HL(m),HLdof(m),HLp(m)] = hosmerLemeshow(P{m}, y, N_BINS);
end
brier = cellfun(@(p) mean((p-y).^2), P);

% DeLong（点估计，无需 bootstrap）
[pJE,~,dJE] = C.delongTest(pJ, pE, y);

%% =========================================================
%  2. 唯一的 bootstrap 循环（三类指标一起算）
%% =========================================================
fprintf('运行 bootstrap（N=%d，三图共用）...\n', N_BOOT);
rowsAll = find(idx); Nn = numel(rowsAll);

% 预分配
auc_b   = nan(N_BOOT, 2);     % AUC：对样本内分数重采样（标准做法）
nri_b   = nan(N_BOOT, 1);     % joint vs eeg
idi_b   = nan(N_BOOT, 1);
opt_s   = nan(N_BOOT, 2);     % 校准 slope optimism
opt_i   = nan(N_BOOT, 2);     % 校准 intercept optimism
% 校准曲线 CI 用的网格
xgrid   = linspace(0,1,50);
calYB   = nan(N_BOOT, 50, 2); % 每模型一层

rng(42);
for b = 1:N_BOOT
    sel = randsample(Nn, Nn, true);
    rb  = rowsAll(sel);
    yb  = D.grp(rb);
    if numel(unique(yb)) < 2; continue; end

    % --- (a) AUC CI：对主拟合分数按同一 sel 重采样（不重拟合，标准）---
    for m = 1:2
        sb = P{m}(sel);
        [~,~,~,auc_b(b,m)] = C.computeROC(sb, yb);
    end

    % --- 以下需要重拟合：NRI/IDI、optimism、校准曲线 ---
    [~, Xeb, Xjb] = C.buildDesign(D, rb);
    Xb = {Xeb, Xjb};
    % 列数需与主设计一致（某 etiology 类别缺失会少列）→ 跳过该次
    okwidth = size(Xeb,2)==size(Xe,2) && size(Xjb,2)==size(Xj,2);
    if ~okwidth; continue; end

    fb = cell(1,2); pbb = cell(1,2); pbo = cell(1,2);
    good = true;
    for m = 1:2
        try
            fb{m}  = C.fitLogit(Xb{m}, yb);
            pbb{m} = fb{m}.prob;                                  % boot 预测 boot
            Xz_o   = (Xsets{m} - fb{m}.mu) ./ fb{m}.sg;
            pbo{m} = predict(fb{m}.mdl, Xz_o);                    % boot 预测 orig
        catch
            good = false; break;
        end
    end
    if ~good; continue; end

    % --- (b) NRI / IDI（boot 重拟合，参照 EEG）---
    Rj = C.nri_idi(pbb{1}, pbb{2}, yb, NRI_THR);
    nri_b(b) = Rj.NRI_cont;
    idi_b(b) = Rj.IDI;

    % --- (c) 校准 optimism（slope & intercept）---
    for m = 1:2
        [ib_,sb_] = calSlopeIntercept(pbb{m}, yb);   % apparent on boot
        [io_,so_] = calSlopeIntercept(pbo{m}, y);    % boot-model on orig
        opt_s(b,m) = sb_ - so_;
        opt_i(b,m) = ib_ - io_;
    end

    % --- (d) 校准曲线 CI（对主拟合分数重采样后核回归）---
    for m = 1:2
        calYB(b,:,m) = kreg(P{m}(sel), yb, xgrid);
    end
end

% 汇总 CI
ci_auc = nan(2,2);
for m=1:2; v=auc_b(~isnan(auc_b(:,m)),m); ci_auc(m,:)=quantile(v,[.025 .975]); end
v=nri_b(~isnan(nri_b)); ci_nri=quantile(v,[.025 .975]);
v=idi_b(~isnan(idi_b)); ci_idi=quantile(v,[.025 .975]);
slope_cor = slope_app - mean(opt_s,1,'omitnan');
int_cor   = int_app   - mean(opt_i,1,'omitnan');

% 校准曲线带
calMean = nan(2,50); calLo = nan(2,50); calHi = nan(2,50);
for m=1:2
    layer = calYB(:,:,m);
    calMean(m,:) = mean(layer,1,'omitnan');
    calLo(m,:)   = quantile(layer,0.025,1);
    calHi(m,:)   = quantile(layer,0.975,1);
end

fprintf('bootstrap 完成。\n\n');

%% =========================================================
%  2.5 两模型诊断分类指标（Youden 最优阈值）
%       灵敏度 / 特异度 / 准确率 / Youden / PPV / NPV ...
%% =========================================================
for m = 1:2
    d0.name = mlog{m};
    DX(m) = diagMetrics(P{m}, y, C, d0);   % 计算并填充
    clear d0;
    fprintf('[%s] Sen=%.3f  Spe=%.3f  Acc=%.3f  Youden=%.3f  (cut=%.4f)\n', ...
        mlog{m}, DX(m).Sen, DX(m).Spe, DX(m).Acc, DX(m).Youden, DX(m).cutoff);
end
fprintf('\n');

% 写出独立汇总表（txt + csv）
writeDiagMetrics(OUTPUT_DIR, N, y, DX);

%% =========================================================
%  3. 写统计日志（三份）
%% =========================================================
writeDiscLog(OUTPUT_DIR, PRE_EXCL_THR, N, y, mlog, mfeat, fits, ...
    auc0, ci_auc, dJE, pJE, DX);
writeDCALog(OUTPUT_DIR, scen(), N, y, P, R_joint, ...
    ci_nri, ci_idi, NRI_THR);
writeCalLog(OUTPUT_DIR, N, mlog, slope_app, int_app, slope_cor, int_cor, ...
    HL, HLdof, HLp, brier);

%% =========================================================
%  4. 图1：区分度（ROC + AUC 森林）—— 独立窗口
%% =========================================================
fig1 = figure('Name','Discrimination','Position',[120 380 470 430],'Color','w');
% ROC
ax = axes(fig1); hold(ax,'on');
plot(ax,[0 1],[0 1],'--','Color',[.6 .6 .6],'LineWidth',1.1);
hh = gobjects(1,2);
for m=1:2
    [fpr,tpr] = C.computeROC(P{m}, y);
    hh(m)=plot(ax,fpr,tpr,'-','Color',mcol{m}*0.85,'LineWidth',(m==2)*2.4+(m<2)*2.0);
end
xlim(ax,[0 1]);ylim(ax,[0 1]);axis(ax,'square');
xlabel(ax,'1 - Specificity','FontSize',10);ylabel(ax,'Sensitivity','FontSize',10);
rocLegendNames = arrayfun(@(m) sprintf('%s (AUC = %.3f, Accuracy = %.3f)', ...
    mnames{m}, auc0(m), DX(m).Acc), ...
    1:2, 'UniformOutput', false);
legend(ax,hh,rocLegendNames,'Location','southeast','FontSize',8.5,'Box','off');
set(ax,'TickDir','out','Box','off','FontSize',9);
exportgraphics(fig1,fullfile(OUTPUT_DIR,'discrimination_2models.png'),'Resolution',300);
exportgraphics(fig1,fullfile(OUTPUT_DIR,'discrimination_2models.pdf'),'ContentType','vector');
fprintf('图1 区分度已保存\n');

%% =========================================================
%  5. 图2：临床效用（DCA + NRI/IDI 森林）—— 独立窗口
%% =========================================================
pt=linspace(0.01,0.80,100); prev=mean(y);
nb=cell(1,2); for m=1:2; nb{m}=netBenefit(P{m},y,pt); end
nb_all=prev-(1-prev).*(pt./(1-pt)); nb_none=zeros(size(pt));

fig2 = figure('Name','ClinicalUtility','Position',[150 120 500 440],'Color','w');
% DCA
ax = axes(fig2); hold(ax,'on');
plot(ax,pt,nb_none,'-','Color',[.55 .55 .55],'LineWidth',1.1);
plot(ax,pt,nb_all,'--','Color',[.35 .35 .35],'LineWidth',1.1);
for m=1:2; plot(ax,pt,nb{m},'-','Color',mcol{m}*0.85,'LineWidth',(m==2)*2.4+(m<2)*2.0); end
ymax=max([nb{:}])*1.15;
xlim(ax,[0 0.8]);ylim(ax,[max(-0.05,-prev*0.5),max(0.05,ymax)]);
xlabel(ax,'Threshold probability','FontSize',10);ylabel(ax,'Net benefit','FontSize',10);
legend(ax,{'None','All',mnames{:}},'Location','northeast','FontSize',8.5,'Box','off');
set(ax,'TickDir','out','Box','off','FontSize',9);
exportgraphics(fig2,fullfile(OUTPUT_DIR,'dca_nri_idi_2models.png'),'Resolution',300);
exportgraphics(fig2,fullfile(OUTPUT_DIR,'dca_nri_idi_2models.pdf'),'ContentType','vector');
fprintf('图2 临床效用已保存\n');

%% =========================================================
%  6. 图3：校准（曲线 + slope 对比）—— 独立窗口
%% =========================================================
fig3 = figure('Name','Calibration','Position',[180 -100 490 460],'Color','w');
% 校准曲线
ax = axes(fig3); hold(ax,'on');
plot(ax,[0 1],[0 1],'--','Color',[.55 .55 .55],'LineWidth',1.1);
hh = gobjects(1,2);
for m=1:2
    xs=xgrid; ys=calMean(m,:); yl=calLo(m,:); yu=calHi(m,:);
    good=~isnan(ys);
    fill(ax,[xs(good) fliplr(xs(good))],[yl(good) fliplr(yu(good))],mcol{m}, ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    hh(m)=plot(ax,xs(good),min(max(ys(good),0),1),'-','Color',mcol{m}*0.85, ...
        'LineWidth',(m==2)*2.4+(m<2)*2.0);
    [bx,by]=binCalib(P{m},y,N_BINS);
    scatter(ax,bx,by,40,mcol{m},'filled','MarkerEdgeColor','w','LineWidth',0.6,'HandleVisibility','off');
end
% rug（joint）
yr=-0.04; pr=P{2}(:)';
for ii=1:numel(pr); line(ax,[pr(ii) pr(ii)],[yr-0.012 yr+0.012],'Color',c_joint*0.7,'LineWidth',0.8,'HandleVisibility','off'); end
xlim(ax,[0 1]);ylim(ax,[-0.06 1.0]);axis(ax,'square');
xlabel(ax,'Predicted probability','FontSize',10);ylabel(ax,'Observed proportion','FontSize',10);
legend(ax,hh,mnames,'Location','northwest','FontSize',8.5,'Box','off');
set(ax,'TickDir','out','Box','off','FontSize',9);
exportgraphics(fig3,fullfile(OUTPUT_DIR,'calibration_2models.png'),'Resolution',300);
exportgraphics(fig3,fullfile(OUTPUT_DIR,'calibration_2models.pdf'),'ContentType','vector');
fprintf('图3 校准已保存\n');

fprintf('\n全部完成。三张图 + 三份统计日志已输出至:\n  %s\n', OUTPUT_DIR);

%% =========================================================
%  局部函数
%% =========================================================
function s = scen(); s = 'Excl. baseline MCS+/EMCS'; end

function s = pstr(p)
    if p < 0.001; s='<0.001'; else; s=sprintf('%.3f',p); end
end

% ── 诊断分类指标（基于 Youden 最优阈值）─────────────────────
function d = diagMetrics(p, y, C, d)
    p = p(:); y = y(:);
    [fpr,tpr,thr,~] = C.computeROC(p, y);
    yj = tpr - fpr; [~,iy] = max(yj);     % Youden 最大处
    cut = thr(iy);
    pred = double(p >= cut);
    TP=sum(pred==1 & y==1); TN=sum(pred==0 & y==0);
    FP=sum(pred==1 & y==0); FN=sum(pred==0 & y==1);
    N = numel(y);
    d.cutoff = cut;
    d.TP=TP; d.TN=TN; d.FP=FP; d.FN=FN;
    d.Sen    = TP / max(TP+FN,1);          % 灵敏度
    d.Spe    = TN / max(TN+FP,1);          % 特异度
    d.Acc    = (TP+TN) / N;                % 准确率
    d.Youden = d.Sen + d.Spe - 1;          % Youden 指数
    d.PPV    = TP / max(TP+FP,1);
    d.NPV    = TN / max(TN+FN,1);
    d.PLR    = d.Sen / max(1-d.Spe,eps);   % 阳性似然比
    d.NLR    = (1-d.Sen) / max(d.Spe,eps); % 阴性似然比
    d.F1     = 2*TP / max(2*TP+FP+FN,1);
end

function writeDiagMetrics(OUT, N, y, DX)
    % --- txt 汇总 ---
    fid=fopen(fullfile(OUT,'diagnostic_metrics_summary.txt'),'w');
    L=@(varargin) fprintf2(fid,varargin{:});
    sep=repmat('=',1,64);
    L('%s\n  两模型诊断分类指标（Youden 最优阈值）\n%s\n',sep,sep);
    L('生成时间: %s\n',datestr(now));
    L('样本: n=%d (阳性=%d, 阴性=%d)\n\n',N,sum(y==1),sum(y==0));
    L('%-12s %8s %8s %8s %8s %8s %8s %8s\n', ...
        'Model','Cutoff','Sen','Spe','Acc','Youden','PPV','NPV');
    L('%s\n',repmat('-',1,64));
    for m=1:numel(DX)
        d=DX(m);
        L('%-12s %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            d.name, d.cutoff, d.Sen, d.Spe, d.Acc, d.Youden, d.PPV, d.NPV);
    end
    L('\n混淆矩阵与附加指标:\n%s\n',repmat('-',1,64));
    for m=1:numel(DX)
        d=DX(m);
        L('[%s]\n',d.name);
        L('  TP=%d  TN=%d  FP=%d  FN=%d\n',d.TP,d.TN,d.FP,d.FN);
        L('  PLR=%.3f  NLR=%.3f  F1=%.3f\n\n',d.PLR,d.NLR,d.F1);
    end
    fclose(fid);

    % --- csv（便于直接放入论文/Excel）---
    fid=fopen(fullfile(OUT,'diagnostic_metrics.csv'),'w');
    fprintf(fid,'Model,Cutoff,Sensitivity,Specificity,Accuracy,Youden,PPV,NPV,PLR,NLR,F1,TP,TN,FP,FN\n');
    for m=1:numel(DX)
        d=DX(m);
        fprintf(fid,'%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d\n', ...
            d.name,d.cutoff,d.Sen,d.Spe,d.Acc,d.Youden,d.PPV,d.NPV,d.PLR,d.NLR,d.F1, ...
            d.TP,d.TN,d.FP,d.FN);
    end
    fclose(fid);
    fprintf('诊断指标汇总已保存 (diagnostic_metrics_summary.txt / .csv)\n');
end

function forestRow(ax, val, ci, col, yy)
    plot(ax,[ci(1) ci(2)],[yy yy],'-','Color',col*0.75,'LineWidth',2.4);
    plot(ax,[ci(1) ci(1)],[yy-0.1 yy+0.1],'-','Color',col*0.75,'LineWidth',2.4);
    plot(ax,[ci(2) ci(2)],[yy-0.1 yy+0.1],'-','Color',col*0.75,'LineWidth',2.4);
    scatter(ax,val,yy,95,col,'filled','MarkerEdgeColor','w','LineWidth',1.0);
    text(ax,ci(2)+0.012,yy,sprintf('%.2f [%.2f, %.2f]',val,ci(1),ci(2)), ...
        'FontSize',8.5,'VerticalAlignment','middle');
end

function nb = netBenefit(p, y, pt)
    N=numel(y); nb=zeros(size(pt));
    for i=1:numel(pt)
        thr=pt(i); pred=p>=thr;
        TP=sum(pred&y==1); FP=sum(pred&y==0);
        nb(i)=TP/N-(FP/N)*(thr/(1-thr));
    end
end

function [a,b] = calSlopeIntercept(p, y)
    eps=1e-6; pp=min(max(p,eps),1-eps); lp=log(pp./(1-pp));
    mdl=fitglm(lp,y,'Distribution','binomial','Link','logit');
    a=mdl.Coefficients.Estimate(1); b=mdl.Coefficients.Estimate(2);
end

function [HL,dof,pval] = hosmerLemeshow(p, y, g)
    [ps,ord]=sort(p); ys=y(ord);
    edges=quantile(ps,linspace(0,1,g+1)); edges(1)=edges(1)-1e-9; edges(end)=edges(end)+1e-9;
    HL=0; used=0;
    for i=1:g
        mask=ps>edges(i)&ps<=edges(i+1); n=sum(mask);
        if n==0; continue; end
        obs=sum(ys(mask)); exp=sum(ps(mask));
        if exp<=0||exp>=n; continue; end
        HL=HL+(obs-exp)^2/(exp*(1-exp/n)); used=used+1;
    end
    dof=max(used-2,1); pval=1-chi2cdf(HL,dof);
end

function [bx,by] = binCalib(p, y, g)
    [ps,ord]=sort(p); ys=y(ord);
    edges=quantile(ps,linspace(0,1,g+1)); edges(1)=edges(1)-1e-9; edges(end)=edges(end)+1e-9;
    bx=nan(g,1); by=nan(g,1);
    for i=1:g
        mask=ps>edges(i)&ps<=edges(i+1);
        if sum(mask)==0; continue; end
        bx(i)=mean(ps(mask)); by(i)=mean(ys(mask));
    end
    keep=~isnan(bx); bx=bx(keep); by=by(keep);
end

function yhat = kreg(x, y, xq)
    x=x(:); y=y(:); xq=xq(:)';
    h=0.9*min(std(x),iqr(x)/1.34)*numel(x)^(-1/5);
    if h<=0||isnan(h); h=0.1; end
    yhat=nan(1,numel(xq));
    for k=1:numel(xq)
        w=exp(-0.5*((x-xq(k))/h).^2); sw=sum(w);
        if sw>0; yhat(k)=sum(w.*y)/sw; end
    end
end

function fprintf2(fid, varargin)
    fprintf(varargin{:});
    if ~isempty(fid)&&fid>0; fprintf(fid,varargin{:}); end
end

function writeDiscLog(OUT, PRE_EXCL_THR, N, y, mlog, mfeat, fits, ...
        auc0, ci_auc, dJE, pJE, DX)
    fid=fopen(fullfile(OUT,'discrimination_summary.txt'),'w');
    L=@(varargin) fprintf2(fid,varargin{:});
    sep=repmat('=',1,60);
    L('%s\n  区分度分析（ROC / AUC / DeLong）\n%s\n',sep,sep);
    L('生成时间: %s\n',datestr(now));
    L('队列: 剔除基线 MCS+/EMCS（pre Diagnosis>=%d）\n',PRE_EXCL_THR);
    L('样本: n=%d (阳性=%d, 阴性=%d)\n\n',N,sum(y==1),sum(y==0));
    L('%s\n  Logistic 系数（特征已标准化）\n%s\n',repmat('-',1,60),repmat('-',1,60));
    for m=1:2
        co=fits{m}.mdl.Coefficients;
        L('[%s]\n',mlog{m}); L('  (Intercept) = %+.4f\n',co.Estimate(1));
        for k=1:numel(mfeat{m})
            L('  %-14s = %+.4f  (p=%.4f)\n',mfeat{m}{k},co.Estimate(k+1),co.pValue(k+1));
        end
        L('\n');
    end
    L('%s\n  AUC + Youden 最优点分类指标\n%s\n',repmat('-',1,60),repmat('-',1,60));
    for m=1:2
        d=DX(m);
        L('[%s]  AUC=%.3f [%.3f, %.3f]\n',mlog{m},auc0(m),ci_auc(m,1),ci_auc(m,2));
        L('  最优阈值=%.4f  Sen=%.3f  Spe=%.3f  Youden=%.3f\n',d.cutoff,d.Sen,d.Spe,d.Youden);
        L('  TP=%d TN=%d FP=%d FN=%d  Acc=%.3f  PPV=%.3f  NPV=%.3f\n\n', ...
            d.TP,d.TN,d.FP,d.FN,d.Acc,d.PPV,d.NPV);
    end
    L('%s\n  DeLong 配对检验\n%s\n',repmat('-',1,60),repmat('-',1,60));
    L('  M_joint vs M_eeg  : dAUC=%+.3f, p=%.4f\n',dJE,pJE);
    fclose(fid);
end

function writeDCALog(OUT, scen_label, N, y, P, R_joint, ci_nri, ci_idi, NRI_THR)
    fid=fopen(fullfile(OUT,'dca_nri_idi_summary.txt'),'w');
    L=@(varargin) fprintf2(fid,varargin{:});
    sep=repmat('=',1,60); prev=mean(y);
    L('%s\n  临床效用分析（DCA / NRI / IDI）\n%s\n',sep,sep);
    L('生成时间: %s\n',datestr(now));
    L('队列: %s | n=%d (阳性=%d, 阴性=%d)\n',scen_label,N,sum(y==1),sum(y==0));
    L('事件率 = %.3f\n\n',prev);
    L('%s\n  决策曲线：代表性阈值净获益\n%s\n',repmat('-',1,60),repmat('-',1,60));
    L('  %-6s %-10s %-10s %-10s\n','pt','M_eeg','M_joint','TreatAll');
    for q=[0.10 0.20 0.30 0.40 0.50]
        nbe=nbAt(P{1},y,q); nbj=nbAt(P{2},y,q);
        L('  %-6.2f %-10.4f %-10.4f %-10.4f\n',q,nbe,nbj,prev-(1-prev)*(q/(1-q)));
    end
    L('\n%s\n  增量价值（参照 M_eeg）\n%s\n',repmat('-',1,60),repmat('-',1,60));
    L('[M_joint vs M_eeg]\n');
    L('  连续 NRI = %+.3f [%.3f, %.3f]  (事件组=%+.3f, 非事件组=%+.3f)\n', ...
        R_joint.NRI_cont,ci_nri(1),ci_nri(2),R_joint.NRI_event,R_joint.NRI_nonevent);
    L('  IDI      = %+.3f [%.3f, %.3f]\n',R_joint.IDI,ci_idi(1),ci_idi(2));
    L('  分类 NRI@[%.2f %.2f %.2f] = %+.3f\n\n',NRI_THR,R_joint.NRI_cat);
    fclose(fid);
end

function nb = nbAt(p, y, thr)
    N=numel(y); pred=p>=thr;
    TP=sum(pred&y==1); FP=sum(pred&y==0);
    nb=TP/N-(FP/N)*(thr/(1-thr));
end

function writeCalLog(OUT, N, mlog, slope_app, int_app, slope_cor, int_cor, HL, HLdof, HLp, brier)
    fid=fopen(fullfile(OUT,'calibration_summary.txt'),'w');
    L=@(varargin) fprintf2(fid,varargin{:});
    sep=repmat('=',1,60);
    L('%s\n  校准分析\n%s\n',sep,sep);
    L('生成时间: %s\n  n=%d\n\n',datestr(now),N);
    for m=1:2
        L('[%s]\n',mlog{m});
        L('  Apparent       : slope=%.3f, intercept=%+.3f\n',slope_app(m),int_app(m));
        L('  Optimism-corr. : slope=%.3f, intercept=%+.3f\n',slope_cor(m),int_cor(m));
        L('  Hosmer-Lemeshow: chi2=%.2f (df=%d), p=%.3f\n',HL(m),HLdof(m),HLp(m));
        L('  Brier score    : %.3f\n\n',brier(m));
    end
    fclose(fid);
end
