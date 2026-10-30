%% =========================================================
%  export_JointModel.m
%
%  从开发队列（87 例，剔除基线 MCS+/EMCS 后）拟合并【冻结】联合模型
%  M_joint，导出做外部验证所需的全部要素：
%     - logistic 回归系数（在“标准化特征”空间）
%     - 标准化参数 mu / sg（★ 外部数据必须沿用开发集的 mu/sg）
%     - 固定的特征列清单 jointNames（★ 含具体哪几个 Etiology 哑变量及顺序）
%     - 开发集固定 cutoff 值（★ 外部验证计算 sensitivity/specificity/accuracy/Youden 时必须沿用）
%     - 关键配置（结局阈值、基线剔除阈值、频率门控等）
%
%  产出：
%     joint_model_frozen.mat      —— 供 apply_JointModel_External.m 加载
%     joint_model_coefficients.csv/.txt —— 人类可读系数表（放论文补充材料）
%     joint_model_cutoff.csv      —— 开发集 cutoff 及对应指标
%
%  依赖：threeModelCommon.m（同目录或在 path 中）
%
%  ★★ 只需修改「用户配置区」，与主评估脚本保持完全一致 ★★
%% =========================================================

clc; clear; close all;
C = threeModelCommon();

%% ── 用户配置区（须与 run_ModelEvaluation_AllInOne2.m 完全一致）──
cfg.INPUT_FILE = 'E:\WorkPlace\低波幅快波机器学习excel-最高占比\合并结果_v4_filled.xlsx';
OUTPUT_DIR     = 'F:\EEG\LowVoltageResearch\OutputResult\V3\3模型对比\模型';

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
PRE_EXCL_THR = 3;        % pre Diagnosis >= 3 被剔除（与主分析一致）
%% ─────────────────────────────────────────────────────────

if ~exist(OUTPUT_DIR,'dir'); mkdir(OUTPUT_DIR); end

%% 1. 载入数据 + 同样的基线剔除
D   = C.loadData(cfg.INPUT_FILE, cfg);
idx = ~(D.pre >= PRE_EXCL_THR);
y   = D.grp(idx);
N   = sum(idx);
fprintf('开发队列（剔除基线 MCS+/EMCS）: n=%d (阳性=%d, 阴性=%d)\n', ...
        N, sum(y==1), sum(y==0));

%% 2. 构造 M_joint 设计矩阵并拟合（与主脚本同一路径）
[~, ~, Xj, dinfo] = C.buildDesign(D, idx);
fitJ = C.fitLogit(Xj, y);          % 内部已做列标准化（mu/sg 存于 fitJ）

% 主拟合 apparent AUC（自检用）
[~,~,~,aucJ] = C.computeROC(fitJ.prob, y);
fprintf('M_joint apparent AUC = %.3f\n', aucJ);

% 开发集固定 cutoff：按 Youden index 最大原则确定
% ★ 外部验证时必须直接使用该 cutoff，不能在测试集上重新寻优。
pDev = fitJ.prob(:);
yDev = (y(:) == 1);

assert(any(yDev==1) && any(yDev==0), ...
    '开发集必须同时包含阳性和阴性，才能计算 sensitivity/specificity/cutoff。');

thrList = unique(pDev(~isnan(pDev)));
thrList = sort(thrList(:), 'descend');
assert(~isempty(thrList), '预测概率为空或全为 NaN，无法计算 cutoff。');

sensList   = nan(numel(thrList),1);
specList   = nan(numel(thrList),1);
accList    = nan(numel(thrList),1);
youdenList = nan(numel(thrList),1);

for t = 1:numel(thrList)
    yhat = pDev >= thrList(t);

    TP = sum( yhat &  yDev);
    TN = sum(~yhat & ~yDev);
    FP = sum( yhat & ~yDev);
    FN = sum(~yhat &  yDev);

    sensList(t)   = TP / (TP + FN);
    specList(t)   = TN / (TN + FP);
    accList(t)    = (TP + TN) / numel(yDev);
    youdenList(t) = sensList(t) + specList(t) - 1;
end

% 若多个 cutoff 的 Youden 相同：先选 accuracy 最高者；若仍并列，因 thrList 为降序，取更高 cutoff。
maxYouden = max(youdenList);
candIdx = find(abs(youdenList - maxYouden) <= 1e-12);
[~, bestRelIdx] = max(accList(candIdx));
bestIdx = candIdx(bestRelIdx);

cutoffJ = thrList(bestIdx);
sensJ   = sensList(bestIdx);
specJ   = specList(bestIdx);
accJ    = accList(bestIdx);
youdenJ = youdenList(bestIdx);

fprintf(['M_joint 开发集固定 cutoff = %.6f ', ...
         '(Sensitivity=%.3f, Specificity=%.3f, Accuracy=%.3f, Youden=%.3f)\n'], ...
        cutoffJ, sensJ, specJ, accJ, youdenJ);

%% 3. 打包冻结模型
model = struct();
model.type          = 'M_joint';
model.created       = datestr(now);
model.N_dev         = N;
model.npos_dev      = sum(y==1);
model.nneg_dev      = sum(y==0);

% --- 核心：系数 + 标准化参数 + 固定列清单 ---
co = fitJ.mdl.Coefficients;
model.intercept     = co.Estimate(1);               % 标准化空间截距
model.beta          = co.Estimate(2:end);           % 标准化空间斜率（按列顺序）
model.featNames     = dinfo.jointNames(:);          % ★ 固定列名与顺序
model.mu            = fitJ.mu(:);                    % ★ 各列训练均值
model.sg            = fitJ.sg(:);                    % ★ 各列训练标准差
model.mdl           = fitJ.mdl;                      % 完整 GLM 对象（备查）

% --- 复现特征所需的配置（★ 外部必须沿用）---
model.cfg           = cfg;
model.PRE_EXCL_THR  = PRE_EXCL_THR;
model.GRP_THR       = cfg.GRP_THR;
model.freqGate      = 8;     % Dominant Stability 频率门控：freq>=8 取 Alpha 否则 DeltaTheta
model.etioRef       = 1;     % Etiology 参照组 = trauma
model.apparentAUC   = aucJ;

% --- 开发集固定 cutoff 及其开发集表观指标（★ 外部验证必须沿用 cutoff）---
model.cutoff        = cutoffJ;       % 推荐外部脚本直接使用 model.cutoff
model.cutoff_dev    = cutoffJ;       % 明确标记：来自开发集
model.cutoffMethod  = 'Max Youden index on development cohort';
model.cutoffTieRule = 'If multiple thresholds have the same Youden, choose the highest accuracy; if still tied, choose the highest cutoff.';
model.sens_dev_at_cutoff   = sensJ;
model.spec_dev_at_cutoff   = specJ;
model.acc_dev_at_cutoff    = accJ;
model.youden_dev_at_cutoff = youdenJ;
model.cutoffTable_dev = table(thrList, sensList, specList, accList, youdenList, ...
    'VariableNames', {'Cutoff','Sensitivity','Specificity','Accuracy','Youden'});

% 一致性自检：系数个数应等于特征列数
assert(numel(model.beta)==numel(model.featNames), ...
    '系数个数与特征列名个数不一致，导出中止。');

save(fullfile(OUTPUT_DIR,'joint_model_frozen.mat'), 'model');
fprintf('已保存冻结模型: joint_model_frozen.mat\n');

%% 4. 人类可读系数表（论文补充材料）
fid = fopen(fullfile(OUTPUT_DIR,'joint_model_coefficients.csv'),'w');
fprintf(fid,'Term,Beta_std,Mu_train,Sg_train,pValue\n');
fprintf(fid,'(Intercept),%.6f,,,%.4g\n', model.intercept, co.pValue(1));
for k = 1:numel(model.featNames)
    fprintf(fid,'%s,%.6f,%.6f,%.6f,%.4g\n', ...
        model.featNames{k}, model.beta(k), model.mu(k), model.sg(k), co.pValue(k+1));
end
fclose(fid);

fid = fopen(fullfile(OUTPUT_DIR,'joint_model_coefficients.txt'),'w');
fprintf(fid,'冻结 M_joint 模型（标准化特征空间）\n');
fprintf(fid,'生成时间: %s\n', model.created);
fprintf(fid,'开发队列 n=%d (阳性=%d, 阴性=%d), apparent AUC=%.3f\n', ...
    N, sum(y==1), sum(y==0), aucJ);
fprintf(fid,'开发集固定 cutoff = %.6f（Max Youden index）\n', model.cutoff);
fprintf(fid,'cutoff 对应开发集指标: Sensitivity=%.3f, Specificity=%.3f, Accuracy=%.3f, Youden=%.3f\n\n', ...
    model.sens_dev_at_cutoff, model.spec_dev_at_cutoff, ...
    model.acc_dev_at_cutoff, model.youden_dev_at_cutoff);
fprintf(fid,'%-16s %12s %12s %12s\n','Term','Beta(std)','Mu','Sg');
fprintf(fid,'%s\n',repmat('-',1,54));
fprintf(fid,'%-16s %12.4f\n','(Intercept)', model.intercept);
for k = 1:numel(model.featNames)
    fprintf(fid,'%-16s %12.4f %12.4f %12.4f\n', ...
        model.featNames{k}, model.beta(k), model.mu(k), model.sg(k));
end
fprintf(fid,'\n说明：预测时先用 (x-Mu)/Sg 标准化，再 logit = Intercept + sum(Beta.*x_std)。\n');
fprintf(fid,'Etiology 哑变量参照组 = trauma；外部缺某病因类别时该列补 0。\n');
fprintf(fid,'外部验证计算分类指标时使用 yhat = prob >= cutoff，不能用测试集重新计算 cutoff。\n');
fclose(fid);

%% 5. 开发集 cutoff 表（方便论文补充材料或外部脚本核对）
fid = fopen(fullfile(OUTPUT_DIR,'joint_model_cutoff.csv'),'w');
fprintf(fid,'Cutoff,Method,Sensitivity_dev,Specificity_dev,Accuracy_dev,Youden_dev,AUC_dev\n');
fprintf(fid,'%.6f,"%s",%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
    model.cutoff, model.cutoffMethod, ...
    model.sens_dev_at_cutoff, model.spec_dev_at_cutoff, ...
    model.acc_dev_at_cutoff, model.youden_dev_at_cutoff, model.apparentAUC);
fclose(fid);

fprintf('已保存系数表: joint_model_coefficients.csv / .txt\n');
fprintf('已保存开发集 cutoff: joint_model_cutoff.csv\n');
fprintf('\n导出完成。下一步在外部数据上运行 apply_JointModel_External.m\n');
