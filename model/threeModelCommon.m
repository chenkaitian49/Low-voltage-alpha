function C = threeModelCommon()
%% =========================================================
%  threeModelCommon.m
%
%  三模型框架（M_clin / M_eeg / M_joint）的共享工具函数集合。
%  通过 C = threeModelCommon(); 取得一个结构体，里面每个字段是一个
%  函数句柄，供两个主脚本（敏感性分析 forest、DCA+NRI/IDI）调用，
%  确保两份图共用同一套数据加载、结局定义、模型拟合与统计方法。
%
%  设计要点（与 NAPI 方法学决策一致）：
%   - 三个模型在同一队列、同一结局上拟合 → DeLong 配对检验、NRI/IDI 成立。
%   - M_eeg 沿用原 analyze_JointROC 的三变量（Dominant Freq / Broadband Amp
%     / Dominant Stability），Dominant Stability 按 Mean_Freq>=8 自动取
%     Alpha_Stability 否则 DeltaTheta_Stability。
%   - M_clin 使用 Age, Etiology(类别), Duration, pre CRS-R, pre Diagnosis。
%     Etiology 为名义类别（trauma/stroke/ischemia），做哑变量编码，
%     不作为连续 1/2/3 进入模型。
%   - M_joint = M_eeg 特征 + M_clin 特征 全部并入同一 logistic 回归。
%   - 结局 grp：post Diagnosis >= 3 → 1 (MCS+/EMCS)，否则 0 (UWS/MCS-)。
%% =========================================================

C.loadData      = @loadData;
C.buildDesign   = @buildDesign;
C.fitLogit      = @fitLogit;
C.computeROC    = @computeROC;
C.aucCI_boot    = @aucCI_boot;
C.delongTest    = @delongTest;
C.nri_idi       = @nri_idi;
C.zcolsTrain    = @zcolsTrain;
C.zcolsApply    = @zcolsApply;

end

%% ---------------------------------------------------------
function D = loadData(INPUT_FILE, cfg)
%LOADDATA  读取 Excel、计算 Dominant Stability、构建结局与各模型原始特征。
%  返回结构体 D，字段含全部有效样本（EEG+临床均非缺失）的行。
%  不在此处做基线剔除——剔除在主脚本按场景控制。

    if ~isfile(INPUT_FILE)
        error('找不到输入文件: %s', INPUT_FILE);
    end
    T = readTablePreserve(INPUT_FILE);

    getcol = @(nm) getColumn(T, nm);

    post_diag = double(getcol(cfg.COL_POST_DIAG));
    pre_diag  = double(getcol(cfg.COL_PRE_DIAG));
    age       = double(getcol(cfg.COL_AGE));
    etiology  = double(getcol(cfg.COL_ETIOLOGY));
    duration  = double(getcol(cfg.COL_DURATION));
    pre_crs   = double(getcol(cfg.COL_PRE_CRS));
    freq      = double(getcol(cfg.COL_VAR1));
    amp       = double(getcol(cfg.COL_VAR2));
    alpha_s   = double(getcol(cfg.COL_ALPHA_STAB));
    delta_s   = double(getcol(cfg.COL_DELTA_STAB));

    try
        name_en = getcol(cfg.COL_NAME_EN);
        if ~iscell(name_en); name_en = cellstr(string(name_en)); end
    catch
        name_en = repmat({'N/A'}, height(T), 1);
    end

    % Dominant Stability：Mean_Freq>=8 取 Alpha，否则取 DeltaTheta
    dom_stab = zeros(size(freq));
    is_alpha = freq >= 8;
    dom_stab( is_alpha) = alpha_s( is_alpha);
    dom_stab(~is_alpha) = delta_s(~is_alpha);

    % 有效样本：EEG 三特征 + 临床五变量 + 结局 均非缺失
    valid = ~isnan(freq) & ~isnan(amp) & ~isnan(dom_stab) ...
          & ~isnan(age) & ~isnan(etiology) & ~isnan(duration) ...
          & ~isnan(pre_crs) & ~isnan(pre_diag) & ~isnan(post_diag);

    D.name     = name_en(valid);
    D.post     = post_diag(valid);
    D.pre      = pre_diag(valid);
    D.age      = age(valid);
    D.etiology = etiology(valid);
    D.duration = duration(valid);
    D.preCRS   = pre_crs(valid);
    D.freq     = freq(valid);
    D.amp      = amp(valid);
    D.domstab  = dom_stab(valid);
    D.isAlpha  = is_alpha(valid);

    % 结局：post Diagnosis >= 3 → 1 (MCS+/EMCS)
    D.grp = double(D.post >= cfg.GRP_THR);
    D.N   = numel(D.grp);
end

%% ---------------------------------------------------------
function T = readTablePreserve(INPUT_FILE)
%READTABLEPRESERVE  读取 Excel 并尽力保留原始列标题（跨 MATLAB 版本兼容）。
%   不同版本设置原始列名的方式不同，依次尝试，全失败再退回默认读取。
    T = [];
    % 方式1：readtable(...,'VariableNamingRule','preserve')  (R2019b+)
    if isempty(T)
        try
            T = readtable(INPUT_FILE, 'VariableNamingRule', 'preserve');
        catch
        end
    end
    % 方式2：detectImportOptions + VariableNamingRule
    if isempty(T)
        try
            opts = detectImportOptions(INPUT_FILE);
            opts.VariableNamingRule = 'preserve';
            T = readtable(INPUT_FILE, opts);
        catch
        end
    end
    % 方式3：旧参数 'PreserveVariableNames'  (R2019a)
    if isempty(T)
        try
            T = readtable(INPUT_FILE, 'PreserveVariableNames', true);
        catch
        end
    end
    % 方式4：默认（列名会被改成合法标识符；getColumn 有兜底匹配）
    if isempty(T)
        T = readtable(INPUT_FILE);
    end
end

%% ---------------------------------------------------------
function v = getColumn(T, nm)
%GETCOLUMN  按列名取列；若原始名未保留（被 MATLAB 改成合法标识符），
%   则用 matlab.lang.makeValidName 把目标名做同样清洗后再匹配。
    vars = T.Properties.VariableNames;
    % 1) 直接匹配原始名
    if ismember(nm, vars)
        v = T.(nm); return;
    end
    % 2) 兜底：把目标名与表中各列名都做"合法标识符"清洗后匹配
    try
        nm_clean = matlab.lang.makeValidName(nm);
    catch
        nm_clean = regexprep(nm, '[^A-Za-z0-9]', '_');
    end
    if ismember(nm_clean, vars)
        v = T.(nm_clean); return;
    end
    % 3) 再兜底：清洗表中列名后比对（应对 makeValidName 细节差异）
    vars_clean = cell(size(vars));
    for i = 1:numel(vars)
        try
            vars_clean{i} = matlab.lang.makeValidName(vars{i});
        catch
            vars_clean{i} = regexprep(vars{i}, '[^A-Za-z0-9]', '_');
        end
    end
    hit = find(strcmp(nm_clean, vars_clean), 1);
    if ~isempty(hit)
        v = T.(vars{hit}); return;
    end
    % 4) 仍未找到 → 报错并列出可用列名，便于排查
    error(['未找到列 "%s"（清洗后 "%s"）。\n表中可用列名为:\n  %s\n' ...
           '请核对 Excel 表头与配置区列名是否一致。'], ...
           nm, nm_clean, strjoin(vars, '\n  '));
end

%% ---------------------------------------------------------
function [Xc, Xe, Xj, info] = buildDesign(D, idx)
%BUILDDESIGN  为给定行索引 idx 构造三模型的设计矩阵（未含截距，fitglm 自动加）。
%   Xc = M_clin 特征：Age, Duration, preCRS, preDiag(序数), Etiology哑变量(2列)
%   Xe = M_eeg  特征：Freq, Amp, DomStability
%   Xj = M_joint 特征：[Xc, Xe]
%
%   Etiology 哑变量以 trauma(=1) 为参照组，生成 stroke / ischemia 两个指示列。
%   注意：哑变量是否生成某一列取决于该子集是否出现该类别。info 记录列名。

    age   = D.age(idx);
    dur   = D.duration(idx);
    crs   = D.preCRS(idx);
    pred  = D.pre(idx);          % 初始诊断作为序数预测子
    etio  = D.etiology(idx);

    freq  = D.freq(idx);
    amp   = D.amp(idx);
    stab  = D.domstab(idx);

    % --- Etiology 哑变量（参照 = trauma=1）---
    etio_levels = [2 3];             % stroke, ischemia
    etio_names  = {'Etio_stroke','Etio_ischemia'};
    keep = false(1, numel(etio_levels));
    Edum = [];
    for k = 1:numel(etio_levels)
        col = double(etio == etio_levels(k));
        if any(col) && ~all(col)     % 该类别在子集中存在且非全部，才纳入
            Edum = [Edum, col]; %#ok<AGROW>
            keep(k) = true;
        end
    end
    used_etio_names = etio_names(keep);

    Xc = [age, dur, crs, pred, Edum];
    Xe = [freq, amp, stab];
    Xj = [Xc, Xe];

    info.clinNames  = [{'Age','Duration','preCRS','preDiag'}, used_etio_names];
    info.eegNames   = {'DomFreq','BroadAmp','DomStability'};
    info.jointNames = [info.clinNames, info.eegNames];
end

%% ---------------------------------------------------------
function [mu, sg] = zcolsTrain(X)
%ZCOLSTRAIN  返回各列均值/标准差（标准差为0时置1，避免除零）。
    mu = mean(X, 1);
    sg = std(X, 0, 1);
    sg(sg == 0) = 1;
end

function Xz = zcolsApply(X, mu, sg)
%ZCOLSAPPLY  用给定 mu/sg 标准化（连续列居中缩放；哑变量0/1列会被一并缩放，
%   对 logistic 拟合无碍，仅改变系数尺度，不改变预测概率与 ROC）。
    Xz = (X - mu) ./ sg;
end

%% ---------------------------------------------------------
function out = fitLogit(X, y)
%FITLOGIT  拟合 logistic 回归，返回模型与样本内预测概率。
%   对特征做列标准化以改善小样本数值稳定性（不影响 AUC / 概率）。
    [mu, sg] = zcolsTrain(X);
    Xz = zcolsApply(X, mu, sg);
    mdl = fitglm(Xz, y, 'Distribution','binomial','Link','logit');
    out.mdl  = mdl;
    out.mu   = mu;
    out.sg   = sg;
    out.prob = predict(mdl, Xz);
end

%% ---------------------------------------------------------
function [fpr, tpr, thr, auc] = computeROC(scores, labels)
%COMPUTEROC  ROC（封装 perfcurve，含备用实现），正类 = 1。
    try
        [fpr, tpr, thr, auc] = perfcurve(labels, scores, 1);
    catch
        thr = [max(scores)+1; sort(unique(scores),'descend')];
        tpr = zeros(numel(thr),1); fpr = zeros(numel(thr),1);
        pos = labels==1; neg = labels==0;
        for i = 1:numel(thr)
            pred = scores >= thr(i);
            tpr(i) = sum(pred & pos) / max(sum(pos),1);
            fpr(i) = sum(pred & neg) / max(sum(neg),1);
        end
        auc = abs(trapz(fpr, tpr));
        if auc < 0.5; auc = 1-auc; tpr = 1-tpr; fpr = 1-fpr; end
    end
end

%% ---------------------------------------------------------
function [auc, ci, auc_boot] = aucCI_boot(scores, labels, nboot, seed)
%AUCCI_BOOT  对给定一组（样本内）预测分数做 bootstrap，返回 AUC 点估计与
%   分位数 95% CI。重抽样在 (scores,labels) 配对上进行（评估端 bootstrap）。
    if nargin < 4; seed = 42; end
    [~,~,~,auc] = computeROC(scores, labels);
    N = numel(labels);
    auc_boot = NaN(nboot,1);
    rng(seed);
    for b = 1:nboot
        ib = randsample(N, N, true);
        yb = labels(ib);
        if numel(unique(yb)) < 2; continue; end
        [~,~,~,ab] = computeROC(scores(ib), yb);
        auc_boot(b) = ab;
    end
    auc_boot = auc_boot(~isnan(auc_boot));
    ci = quantile(auc_boot, [0.025 0.975]);
end

%% ---------------------------------------------------------
function [p, z, dAUC] = delongTest(s1, s2, labels)
%DELONGTEST  DeLong 配对检验比较两条相关 ROC 的 AUC（同一组受试者）。
%   s1, s2 = 两个模型对同一批受试者的预测分数；labels = 0/1 结局。
%   返回双侧 p 值、z 统计量、AUC 差 (AUC1 - AUC2)。
    pos = labels == 1; neg = labels == 0;
    X = s1(pos); Y = s1(neg);   % model 1
    U = s2(pos); V = s2(neg);   % model 2
    m = numel(X); n = numel(Y);

    [a1, v10_1, v01_1] = structComponents(X, Y);
    [a2, v10_2, v01_2] = structComponents(U, V);

    % 协方差结构（DeLong 1988）
    S10 = cov([v10_1, v10_2]);     % m×2 → 2×2
    S01 = cov([v01_1, v01_2]);     % n×2 → 2×2
    S   = S10 / m + S01 / n;

    L = [1; -1];
    var_d = L' * S * L;
    dAUC = a1 - a2;
    if var_d <= 0
        z = 0; p = 1; return;
    end
    z = dAUC / sqrt(var_d);
    p = 2 * (1 - normcdf(abs(z)));
end

function [auc, v10, v01] = structComponents(Xp, Yn)
%STRUCTCOMPONENTS  DeLong 结构分量（placement values / midrank 处理 ties）。
    m = numel(Xp); n = numel(Yn);
    v10 = zeros(m,1); v01 = zeros(n,1);
    for i = 1:m
        v10(i) = mean(psi(Xp(i), Yn));
    end
    for j = 1:n
        v01(j) = mean(psi(Xp, Yn(j)));
    end
    auc = mean(v10);
end

function val = psi(a, b)
%PSI  比较核：a>b →1, a==b →0.5, a<b →0（支持向量广播）。
    val = double(a > b) + 0.5 * double(a == b);
end

%% ---------------------------------------------------------
function R = nri_idi(p_old, p_new, labels, thresholds)
%NRI_IDI  计算分类改善指标，比较旧模型(p_old)与新模型(p_new)。
%   - 连续 NRI（无阈值）：事件组中 p 上升比例 - 下降比例，加上
%     非事件组中 p 下降比例 - 上升比例。
%   - 分类 NRI（给定阈值 thresholds，可为标量或向量）：按风险分层的净重分类。
%   - IDI：事件组平均(p_new-p_old) - 非事件组平均(p_new-p_old)。
%   返回结构体 R，含点估计与（对 IDI/连续NRI 的）解析或自助说明。
    if nargin < 4 || isempty(thresholds); thresholds = []; end

    ev  = labels == 1;
    nev = labels == 0;
    d   = p_new - p_old;

    % --- 连续 NRI ---
    up_ev   = mean(d(ev)  > 0);
    down_ev = mean(d(ev)  < 0);
    up_nev  = mean(d(nev) > 0);
    down_nev= mean(d(nev) < 0);
    R.NRI_event    = up_ev - down_ev;
    R.NRI_nonevent = down_nev - up_nev;
    R.NRI_cont     = R.NRI_event + R.NRI_nonevent;

    % --- IDI ---
    R.IS_new = mean(p_new(ev));   R.IS_old = mean(p_old(ev));
    R.IP_new = mean(p_new(nev));  R.IP_old = mean(p_old(nev));
    R.IDI    = (R.IS_new - R.IS_old) - (R.IP_new - R.IP_old);

    % --- 分类 NRI（可选）---
    R.thresholds = thresholds;
    if ~isempty(thresholds)
        cat_old = sum(p_old(:) >= thresholds(:)', 2);  % 风险类别序号
        cat_new = sum(p_new(:) >= thresholds(:)', 2);
        moved_up   = cat_new > cat_old;
        moved_down = cat_new < cat_old;
        R.NRI_cat_event    = mean(moved_up(ev))   - mean(moved_down(ev));
        R.NRI_cat_nonevent = mean(moved_down(nev))- mean(moved_up(nev));
        R.NRI_cat          = R.NRI_cat_event + R.NRI_cat_nonevent;
    end
end
