% Program for CPM of sbMCN on EF outcomes
% Adds: (1) LOOCV in training (plots + CSV) (2) CSV saves for train/test predictions
% Robust to NaNs and rank deficiency; includes heatmaps and mask saving.
% Yuyao Zhao, Updated: Aug 20, 2025

clear; clc;

%% -------- Paths --------
trainPath    = '.../3_sbMCN/sbMCN_CT_train/Outputs_Y2to10/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y2to10';  % training imaging
testPath     = '.../3_sbMCN/sbMCN_CT_test/Outputs_Y2to10/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y2to10';   % test imaging
BehavPath    = '.../0_data_behavior';                                      % behavior data
outputFolder = '.../3_sbMCN/sbMCN_CT_CPM/Outputs_Y2to10/reg_gender_ROIave_no_MVM_WB_AAL_Y2to10_CPM_092425';
figureFolder = '.../3_sbMCN/sbMCN_CT_CPM/Figures_Y2to10/reg_gender_ROIave_no_MVM_WB_AAL_Y2to10_CPM_092425';

%% -------- Parameters --------
thresholds = [0.05, 0.01];  % p-value thresholds for edge selection
all_cols   = 1:78;

% Ensure output directories exist
if ~exist(fullfile(outputFolder,'sig_masks'),'dir'),    mkdir(fullfile(outputFolder,'sig_masks'));    end
if ~exist(fullfile(outputFolder,'nonsig_masks'),'dir'), mkdir(fullfile(outputFolder,'nonsig_masks')); end
if ~exist(fullfile(outputFolder,'stats'),'dir'),        mkdir(fullfile(outputFolder,'stats'));        end
if ~exist(fullfile(outputFolder,'predictions'),'dir'),  mkdir(fullfile(outputFolder,'predictions'));  end
if ~exist(figureFolder,'dir'),                          mkdir(figureFolder);                          end

%% -------- Network define (Yeo-7) --------
networks = struct( ...
    'Name',  {'Vis','Som','Lim','Pos','Neg'}, ...
    'Cols',  {39:52,[1:2,17:18,20,53:54,65:70],[5:6,21:22,27:28,71:72,75:78],[7:14,19,29:30,33:34,55:60],[3:4,15:16,23:26,31:32,35:36,61:64,73:74]} ...
    );

nNet      = numel(networks);
net_names = {networks(1:nNet).Name};

%% -------- Load imaging data --------
trainData = load(fullfile(trainPath,'workspace_z_sorted_no_MVM.mat'), 'MCM_subjects_norm','AveragedResults');
testData  = load(fullfile(testPath, 'workspace_z_sorted_no_MVM.mat'), 'MCM_subjects_norm','AveragedResults');
MCM_train = trainData.MCM_subjects_norm;
MCM_test  = testData.MCM_subjects_norm;
AveragedResults_tr = trainData.AveragedResults;
AveragedResults_te = testData.AveragedResults;

%% -------- Load behaviors/covariates --------
T_tr = readtable(fullfile(BehavPath,'behav_10oldest_train_merged_scaled_no_MVM.csv'));
T_te = readtable(fullfile(BehavPath,'behav_10oldest_test_merged_scaled_no_MVM.csv'));
behavNames = {'nt_fic_as_scaled','nt_dccs_as_scaled','nt_ls_as_scaled','nt_ps_as_scaled','ssp_length_scaled','composite_score'};

mapTr = containers.Map(T_tr.SubjectID, 1:height(T_tr));
mapTe = containers.Map(T_te.SubjectID, 1:height(T_te));

N_tr = numel(AveragedResults_tr);
N_te = numel(AveragedResults_te);

Ytr_all = nan(N_tr, numel(behavNames));
Yte_all = nan(N_te, numel(behavNames));
Sex_tr  = nan(N_tr,1);
Year_tr = nan(N_tr,1);
Sex_te  = nan(N_te,1);
Year_te = nan(N_te,1);
SID_tr  = strings(N_tr,1);
SID_te  = strings(N_te,1);

for s = 1:N_tr
    sid = AveragedResults_tr{s}{1};
    SID_tr(s) = string(sid);
    if isKey(mapTr, sid)
        Ytr_all(s,:) = T_tr{mapTr(sid), behavNames};
        sval = T_tr{mapTr(sid), "Sex"};
        if iscell(sval), sval = sval{1}; end
        if isstring(sval), sval = char(sval); end
        Sex_tr(s)  = double(strcmpi(sval,'M')); % 1=M else 0
        Year_tr(s) = double(T_tr{mapTr(sid), "Year"});
    end
end

for s = 1:N_te
    sid = AveragedResults_te{s}{1};
    SID_te(s) = string(sid);
    if isKey(mapTe, sid)
        Yte_all(s,:) = T_te{mapTe(sid), behavNames};
        sval = T_te{mapTe(sid), "Sex"};
        if iscell(sval), sval = sval{1}; end
        if isstring(sval), sval = char(sval); end
        Sex_te(s)  = double(strcmpi(sval,'M'));
        Year_te(s) = double(T_te{mapTe(sid), "Year"});
    end
end

%% -------- Precompute edge indexing & subject×edge matrices --------
nC = numel(all_cols);
idx_mask = triu(true(nC),1);
[ii, jj]  = find(idx_mask);
nEdges    = numel(ii);

Xtr = nan(N_tr, nEdges);
for i = 1:N_tr
    M = MCM_train{i}(all_cols,all_cols);
    Xtr(i,:) = M(idx_mask);
end
Xte = nan(N_te, nEdges);
for i = 1:N_te
    M = MCM_test{i}(all_cols,all_cols);
    Xte(i,:) = M(idx_mask);
end


Cov_tr = [Sex_tr Year_tr];
Cov_te = [Sex_te Year_te];

%% -------- Main loop over behaviors --------
for b = 1:numel(behavNames)
    Ytr = Ytr_all(:,b);
    Yte = Yte_all(:,b);
    fprintf('\n=== CPM for behavior: %s ===\n', behavNames{b});


    R_train   = nan(numel(thresholds),1);
    P_train   = nan(numel(thresholds),1);
    R_train_pos= nan(numel(thresholds),1);
    P_train_pos= nan(numel(thresholds),1);
    R_train_neg= nan(numel(thresholds),1);
    P_train_neg= nan(numel(thresholds),1);
    R_test    = nan(numel(thresholds),1);
    P_test    = nan(numel(thresholds),1);
    R_test_pos= nan(numel(thresholds),1);
    P_test_pos= nan(numel(thresholds),1);
    R_test_neg= nan(numel(thresholds),1);
    P_test_neg= nan(numel(thresholds),1);

    % Containers for LOOCV predictions
    LOOCV_pred_comb = nan(N_tr, numel(thresholds));
    LOOCV_pred_pos  = nan(N_tr, numel(thresholds));
    LOOCV_pred_neg  = nan(N_tr, numel(thresholds));
    n_pos_mask  = nan(N_tr, numel(thresholds));
    n_neg_mask  = nan(N_tr, numel(thresholds));


    for tIdx = 1:numel(thresholds)
        thr = thresholds(tIdx);
        figFolderB = fullfile(figureFolder, behavNames{b}, sprintf('thr%.3f',thr));
        if ~exist(figFolderB, 'dir'), mkdir(figFolderB); end

        % --- counters to track per-edge selection across all LOOCV folds
        pos_count = zeros(nEdges,1);
        neg_count = zeros(nEdges,1);

        % -----------------------
        % (A) LOOCV in training
        % -----------------------
        for k = 1:N_tr
            train_idx = true(N_tr,1); train_idx(k) = false;    % N-1 train
            test_idx_k = ~train_idx;                            % left-out

            % Partial corr edge selection on N-1
            [r_vec_k, p_vec_k] = partialcorr(Xtr(train_idx,:), Ytr(train_idx), Cov_tr(train_idx,:), ...
                'Rows','complete','Type','Spearman');

            pos_mask_k = r_vec_k > 0 & p_vec_k < thr;
            neg_mask_k = r_vec_k < 0 & p_vec_k < thr;



            % --- accumulate how often each edge is selected (sign-specific)
            pos_count = pos_count + pos_mask_k(:);
            neg_count = neg_count + neg_mask_k(:);

            if ~any(pos_mask_k | neg_mask_k)
                LOOCV_pred_comb(k,tIdx) = NaN;
                LOOCV_pred_pos(k,tIdx)  = NaN;
                LOOCV_pred_neg(k,tIdx)  = NaN;
                n_pos_mask(k,tIdx) = NaN;
                n_neg_mask(k,tIdx) = NaN;
                continue;
            end

            % Build features for N-1 and left-out
            sumpos_tr_k = sum(Xtr(train_idx, pos_mask_k), 2);
            sumneg_tr_k = sum(Xtr(train_idx, neg_mask_k),  2);
            sumpos_ko   = sum(Xtr(test_idx_k,  pos_mask_k), 2);
            sumneg_ko   = sum(Xtr(test_idx_k,  neg_mask_k), 2);

            X_comb = [sumpos_tr_k, sumneg_tr_k, Cov_tr(train_idx,:), ones(sum(train_idx),1)];
            X_pos  = [sumpos_tr_k,              Cov_tr(train_idx,:), ones(sum(train_idx),1)];
            X_neg  = [             sumneg_tr_k, Cov_tr(train_idx,:), ones(sum(train_idx),1)];
            y_tr_k = Ytr(train_idx);

            % Drop rows with any NaN
            valid_rows_comb = all(~isnan(X_comb),2) & ~isnan(y_tr_k);
            valid_rows_pos  = all(~isnan(X_pos ),2) & ~isnan(y_tr_k);
            valid_rows_neg  = all(~isnan(X_neg ),2) & ~isnan(y_tr_k);

            if sum(valid_rows_comb) < 5
                LOOCV_pred_comb(k,tIdx) = NaN;
                LOOCV_pred_pos(k,tIdx)  = NaN;
                LOOCV_pred_neg(k,tIdx)  = NaN;
                n_pos_mask(k,tIdx) = NaN;
                n_neg_mask(k,tIdx) = NaN;
                continue;
            end

            Xc = X_comb(valid_rows_comb,:); yc = y_tr_k(valid_rows_comb);
            Xp = X_pos(valid_rows_pos,:);   yp = y_tr_k(valid_rows_pos);
            Xn = X_neg(valid_rows_neg,:);   yn = y_tr_k(valid_rows_neg);

            % Fit with backslash (robust to rank deficiency)
            b_vec_k = Xc \ yc;
            b_pos_k = Xp \ yp;
            b_neg_k = Xn \ yn;

            cov_ko = Cov_tr(test_idx_k,:);
            if any(isnan([sumpos_ko, sumneg_ko, cov_ko]))
                LOOCV_pred_comb(k,tIdx) = NaN;
                LOOCV_pred_pos(k,tIdx)  = NaN;
                LOOCV_pred_neg(k,tIdx)  = NaN;
                n_pos_mask(k,tIdx) = NaN;
                n_neg_mask(k,tIdx) = NaN;
            else
                LOOCV_pred_comb(k,tIdx) = b_vec_k(1)*sumpos_ko + b_vec_k(2)*sumneg_ko + b_vec_k(3)*cov_ko(1) + b_vec_k(4)*cov_ko(2) + b_vec_k(5);
                LOOCV_pred_pos(k,tIdx)  = b_pos_k(1)*sumpos_ko + b_pos_k(2)*cov_ko(1) + b_pos_k(3)*cov_ko(2) + b_pos_k(4);
                LOOCV_pred_neg(k,tIdx)  = b_neg_k(1)*sumneg_ko + b_neg_k(2)*cov_ko(1) + b_neg_k(3)*cov_ko(2) + b_neg_k(4);
                n_pos_mask(k,tIdx) = sum(pos_mask_k);
                n_neg_mask(k,tIdx) = sum(neg_mask_k);
            end
        end

        % Training LOOCV performance (combined)
        valid_loocv = ~isnan(LOOCV_pred_comb(:,tIdx)) & ~isnan(Ytr);
        [R_tr, p_tr] = corr(LOOCV_pred_comb(valid_loocv,tIdx), Ytr(valid_loocv), 'Type','Spearman','Rows','complete');
        [R_tr_pos, p_tr_pos] = corr(LOOCV_pred_pos(valid_loocv,tIdx), Ytr(valid_loocv), 'Type','Spearman','Rows','complete');
        [R_tr_neg, p_tr_neg] = corr(LOOCV_pred_neg(valid_loocv,tIdx), Ytr(valid_loocv), 'Type','Spearman','Rows','complete');

        R_train(tIdx)   = R_tr;
        P_train(tIdx)   = p_tr;
        R_train_pos(tIdx)   = R_tr_pos;
        P_train_pos(tIdx)   = p_tr_pos;
        R_train_neg(tIdx)   = R_tr_neg;
        P_train_neg(tIdx)   = p_tr_neg;

        % Plot LOOCV observed vs predicted (combined)
        scatter_with_fit(LOOCV_pred_comb(valid_loocv,tIdx), Ytr(valid_loocv), ...
            sprintf('%s LOOCV (combined) thr=%.3f', behavNames{b}, thr), ...
            fullfile(figFolderB, sprintf('%s_thr%.3f_train_LOOCV_scatter.png', behavNames{b}, thr)), 'k');

        % Save LOOCV CSV for R
        T_loocv = table( ...
            SID_tr, repmat("Train-LOOCV",N_tr,1), Ytr, ...
            LOOCV_pred_comb(:,tIdx), LOOCV_pred_pos(:,tIdx), LOOCV_pred_neg(:,tIdx), ...
            n_pos_mask(:,tIdx), n_neg_mask(:,tIdx), Sex_tr, Year_tr, ...
            'VariableNames', {'SubjectID','Set','Observed','Pred_Combined','Pred_Pos','Pred_Neg','Pos_mask_size', 'Neg_mask_size', 'Sex','Year'});
        writetable(T_loocv, fullfile(outputFolder,'predictions', ...
            sprintf('train_LOOCV_%s_thr%.3f.csv', behavNames{b}, thr)));

        % -----------------------
        % (B) Train fit → TEST
        % -----------------------
        % --- intersection masks = edges selected in ALL LOOCV folds
        pos_mask_inter = (pos_count > 0.8*N_tr);
        neg_mask_inter = (neg_count > 0.8*N_tr);

        % 78x78 symmetric masks
        pos_mask_mat_inter = false(nC,nC); neg_mask_mat_inter = false(nC,nC);
        pos_mask_mat_inter(idx_mask) = pos_mask_inter; pos_mask_mat_inter = pos_mask_mat_inter | pos_mask_mat_inter';
        neg_mask_mat_inter(idx_mask) = neg_mask_inter; neg_mask_mat_inter = neg_mask_mat_inter | neg_mask_mat_inter';


        % --- full-train masks = dge-behavior partial corr on full training
        [r_vec_full, p_vec_full] = partialcorr(Xtr, Ytr, Cov_tr, 'Rows','complete','Type','Spearman');
        %pos_mask = r_vec > 0 & p_vec < thr;
        %neg_mask = r_vec < 0 & p_vec < thr;

        % 78x78 symmetric masks
        %pos_mask_mat = false(nC, nC); neg_mask_mat = false(nC, nC);
        %pos_mask_mat(idx_mask) = pos_mask; pos_mask_mat = pos_mask_mat | pos_mask_mat';
        %neg_mask_mat(idx_mask) = neg_mask; neg_mask_mat = neg_mask_mat | neg_mask_mat';

        

        % Sum features using the intersection masks
        sumpos_tr = sum(Xtr(:, pos_mask_inter), 2);
        sumneg_tr = sum(Xtr(:, neg_mask_inter), 2);
        sumpos_te = sum(Xte(:, pos_mask_inter), 2);
        sumneg_te = sum(Xte(:, neg_mask_inter), 2);

        % Define valid rows for full-train fit (no NaNs anywhere)
        valid_tr_full = all(~isnan([sumpos_tr, sumneg_tr, Cov_tr, ones(N_tr,1)]),2) & ~isnan(Ytr);

        X_comb_full = [sumpos_tr(valid_tr_full), sumneg_tr(valid_tr_full), Cov_tr(valid_tr_full,:), ones(sum(valid_tr_full),1)];
        X_pos_full  = [sumpos_tr(valid_tr_full),                              Cov_tr(valid_tr_full,:), ones(sum(valid_tr_full),1)];
        X_neg_full  = [                          sumneg_tr(valid_tr_full),   Cov_tr(valid_tr_full,:), ones(sum(valid_tr_full),1)];
        y_full      = Ytr(valid_tr_full);

        % Fit with backslash
        if ~isempty(X_comb_full)
            b_vec = X_comb_full \ y_full;
        else
            b_vec = [NaN;NaN;NaN;NaN;NaN];
        end
        if ~isempty(X_pos_full)
            b_pos = X_pos_full  \ y_full;
        else
            b_pos = [NaN;NaN;NaN;NaN];
        end
        if ~isempty(X_neg_full)
            b_neg = X_neg_full  \ y_full;
        else
            b_neg = [NaN;NaN;NaN;NaN];
        end

        % Predict TRAIN (info) and TEST
        Ypred_tr = b_vec(1)*sumpos_tr + b_vec(2)*sumneg_tr + b_vec(3)*Cov_tr(:,1) + b_vec(4)*Cov_tr(:,2) + b_vec(5);
        Ypred_te = b_vec(1)*sumpos_te + b_vec(2)*sumneg_te + b_vec(3)*Cov_te(:,1) + b_vec(4)*Cov_te(:,2) + b_vec(5);
        Ypred_pos= b_pos(1)*sumpos_te + b_pos(2)*Cov_te(:,1) + b_pos(3)*Cov_te(:,2) + b_pos(4);
        Ypred_neg= b_neg(1)*sumneg_te + b_neg(2)*Cov_te(:,1) + b_neg(3)*Cov_te(:,2) + b_neg(4);

        valid_te = ~isnan(Ypred_te) & ~isnan(Yte);

        [R_pos, P_pos] = corr(Ypred_pos(valid_te), Yte(valid_te), 'rows','complete','Type','Spearman');
        [R_neg, P_neg] = corr(Ypred_neg(valid_te), Yte(valid_te), 'rows','complete','Type','Spearman');
        [R_te,  p_te ] = corr(Ypred_te(valid_te),  Yte(valid_te),  'rows','complete','Type','Spearman');


        R_test(tIdx)      = R_te;     P_test(tIdx)    = p_te;
        R_test_pos(tIdx)  = R_pos;    P_test_pos(tIdx)= P_pos;
        R_test_neg(tIdx)  = R_neg;    P_test_neg(tIdx)= P_neg;

        fprintf('thr=%.3f| train edges: pos %d-%d, neg %d-%d| test edges: pos %d, neg %d| Train LOOCV Rpos=%.3f (p=%.3g), Rneg=%.3f (p=%.3g), R=%.3f (p=%.3g) | Test Rpos=%.3f (p=%.3g), Rneg=%.3f (p=%.3g), R=%.3f (p=%.3g)\n', ...
            thr, min(n_pos_mask(:,tIdx)), max(n_pos_mask(:,tIdx)), min(n_neg_mask(:,tIdx)), max(n_neg_mask(:,tIdx)), sum(pos_mask_inter), sum(neg_mask_inter), R_train_pos(tIdx), P_train_pos(tIdx), R_train_neg(tIdx), P_train_neg(tIdx), R_train(tIdx), P_train(tIdx), R_test_pos(tIdx), P_test_pos(tIdx), R_test_neg(tIdx), P_test_neg(tIdx), R_test(tIdx), P_test(tIdx));

        % Plot TEST observed vs predicted (combined)
        scatter_with_fit(Ypred_te(valid_te), Yte(valid_te), ...
            sprintf('%s Test (combined) thr=%.3f', behavNames{b}, thr), ...
            fullfile(figFolderB, sprintf('%s_thr%.3f_test_scatter.png', behavNames{b}, thr)), 'k');

        % Save per-threshold stats
        stat_savedir = fullfile(outputFolder,'stats',sprintf('%s_thr%.3f.mat', behavNames{b}, thr));
        save(stat_savedir, 'R_train','P_train','R_train_pos','P_train_pos','R_train_neg','P_train_neg',...
            'R_test','P_test', 'R_test_pos','P_test_pos','R_test_neg','P_test_neg');

        % Save TEST CSV for R
        T_test = table( ...
            SID_te, repmat("Test",N_te,1), Yte, ...
            Ypred_te, Ypred_pos, Ypred_neg, ...
            Sex_te, Year_te, ...
            'VariableNames', {'SubjectID','Set','Observed','Pred_Combined','Pred_Pos','Pred_Neg','Sex','Year'});
        writetable(T_test, fullfile(outputFolder,'predictions', ...
            sprintf('test_FULL_%s_thr%.3f.csv', behavNames{b}, thr)));

        % -------- NEW SAVES: r_vec/p_vec of selected edges, and sum-feature correlations --------
        % r/p for selected positive/negative edges (from full-train r_vec_full/p_vec_full)
        sel_pos_idx = find(pos_mask_inter);
        sel_neg_idx = find(neg_mask_inter);

        % Edge (i,j) pairs for saving
        edge_i = ii(:); edge_j = jj(:);

        % Tables for r/p of selected edges
        T_rpos = table(edge_i(sel_pos_idx), edge_j(sel_pos_idx), ...
            r_vec_full(sel_pos_idx), p_vec_full(sel_pos_idx), ...
            'VariableNames', {'i','j','r','p'});
        T_rneg = table(edge_i(sel_neg_idx), edge_j(sel_neg_idx), ...
            r_vec_full(sel_neg_idx), p_vec_full(sel_neg_idx), ...
            'VariableNames', {'i','j','r','p'});

        % Write CSVs
        writetable(T_rpos, fullfile(outputFolder,'stats', ...
            sprintf('rvec_SELECTED_POS_INTERSECT_%s_thr%.3f.csv', behavNames{b}, thr)));
        writetable(T_rneg, fullfile(outputFolder,'stats', ...
            sprintf('rvec_SELECTED_NEG_INTERSECT_%s_thr%.3f.csv', behavNames{b}, thr)));

        % Correlations of network strengths (sum of features) with behavior
        [Rt_pos_tr, Pt_pos_tr] = corr(sumpos_tr, Ytr, 'Type','Spearman','Rows','pairwise');
        [Rt_neg_tr, Pt_neg_tr] = corr(sumneg_tr, Ytr, 'Type','Spearman','Rows','pairwise');
        [Rt_pos_te, Pt_pos_te] = corr(sumpos_te, Yte, 'Type','Spearman','Rows','pairwise');
        [Rt_neg_te, Pt_neg_te] = corr(sumneg_te, Yte, 'Type','Spearman','Rows','pairwise');

        T_sumcorr = table( ...
            repmat(behavNames(b),4,1), repmat(thr,4,1), ...
            {'Train','Train','Test','Test'}', ...
            {'Pos','Neg','Pos','Neg'}', ...
            [Rt_pos_tr; Rt_neg_tr; Rt_pos_te; Rt_neg_te], ...
            [Pt_pos_tr; Pt_neg_tr; Pt_pos_te; Pt_neg_te], ...
            'VariableNames', {'Behavior','Thr','Set','Sign','SpearmanR','p'});
        writetable(T_sumcorr, fullfile(outputFolder,'stats', ...
            sprintf('sumStrength_CORRs_INTERSECT_%s_thr%.3f.csv', behavNames{b}, thr)));

        % -------- Significance-dependent figures & mask saving --------
        if P_test(tIdx)<0.05 || P_test_pos(tIdx)<0.05 || P_test_neg(tIdx)<0.05
            % Network edge counts
            for iNet = 1:nNet
                idx_i = networks(iNet).Cols;
                is_involved = ismember(ii, idx_i) | ismember(jj, idx_i);
                total_pos = sum(pos_mask_inter & is_involved);
                total_neg = sum(neg_mask_inter & is_involved);
                fprintf('  %-12s: Pos=%3d, Neg=%3d\n', networks(iNet).Name, total_pos, total_neg);
            end

            % Pos-only scatter (TEST)
            scatter_with_fit(Ypred_pos(valid_te), Yte(valid_te), ...
                sprintf('%s Test (positive mask) thr=%.3f', behavNames{b}, thr), ...
                fullfile(figFolderB, sprintf('%s_thr%.3f_pos_scatter.png', behavNames{b}, thr)), 'r');

            % Neg-only scatter (TEST)
            scatter_with_fit(Ypred_neg(valid_te), Yte(valid_te), ...
                sprintf('%s Test (negative mask) thr=%.3f', behavNames{b}, thr), ...
                fullfile(figFolderB, sprintf('%s_thr%.3f_neg_scatter.png', behavNames{b}, thr)), 'b');

            % Heatmaps
            total_edges  = zeros(nNet, nNet);
            pos_edges_mat= zeros(nNet, nNet);
            neg_edges_mat= zeros(nNet, nNet);
            for iNet = 1:nNet
                idx_i = networks(iNet).Cols;
                for jNet = 1:nNet
                    idx_j = networks(jNet).Cols;
                    is_pair = (ismember(ii, idx_i) & ismember(jj, idx_j)) | ...
                        (ismember(ii, idx_j) & ismember(jj, idx_i));
                    total_edges(iNet,jNet)   = sum(is_pair);
                    pos_edges_mat(iNet,jNet) = sum(pos_mask_inter & is_pair);
                    neg_edges_mat(iNet,jNet) = sum(neg_mask_inter & is_pair);
                end
            end

            prop_signed = (pos_edges_mat - neg_edges_mat) ./ total_edges;
            mask_upper  = triu(true(size(prop_signed)), 1);
            prop_signed(total_edges==0) = NaN;
            prop_plot = prop_signed * 100; prop_plot(mask_upper) = NaN;

            figure('Color','w'); imagesc(prop_plot, 'AlphaData', ~isnan(prop_plot));
            colormap(redblue(45)); cb=colorbar;
            xticks(1:nNet); yticks(1:nNet); xticklabels(net_names); yticklabels(net_names);
            set(gca,'XTickLabelRotation',45); xlabel('Network'); ylabel('Network');
            title('Selected Edges% (High - Low)');
            set(cb,'TickLabels', strcat(string(get(cb,'Ticks')),'%')); set(gca,'Color','w','TickLength',[0 0]); box off;
            saveas(gcf, fullfile(figFolderB, sprintf('%s_thr%.3f_net_heatmap.png', behavNames{b}, thr))); close(gcf);

            pos_signed = pos_edges_mat ./ total_edges; pos_signed(total_edges==0)=NaN;
            pos_plot = pos_signed * 100; pos_plot(mask_upper) = NaN;
            figure('Color','w'); imagesc(pos_plot, 'AlphaData', ~isnan(pos_plot));
            colormap(redblue(45)); cb=colorbar;
            xticks(1:nNet); yticks(1:nNet); xticklabels(net_names); yticklabels(net_names);
            set(gca,'XTickLabelRotation',45); xlabel('Network'); ylabel('Network');
            title('Selected Edges% (High)');
            set(cb,'TickLabels', strcat(string(get(cb,'Ticks')),'%')); set(gca,'Color','w','TickLength',[0 0]); box off;
            saveas(gcf, fullfile(figFolderB, sprintf('%s_thr%.3f_net_pos_heatmap.png', behavNames{b}, thr))); close(gcf);

            neg_signed = neg_edges_mat ./ total_edges; neg_signed(total_edges==0)=NaN;
            neg_plot = neg_signed * 100; neg_plot(mask_upper) = NaN;
            figure('Color','w'); imagesc(neg_plot, 'AlphaData', ~isnan(neg_plot));
            colormap(redblue(45)); cb=colorbar;
            xticks(1:nNet); yticks(1:nNet); xticklabels(net_names); yticklabels(net_names);
            set(gca,'XTickLabelRotation',45); xlabel('Network'); ylabel('Network');
            title('Selected Edges% (Low)');
            set(cb,'TickLabels', strcat(string(get(cb,'Ticks')),'%')); set(gca,'Color','w','TickLength',[0 0]); box off;
            saveas(gcf, fullfile(figFolderB, sprintf('%s_thr%.3f_net_neg_heatmap.png', behavNames{b}, thr))); close(gcf);

            % Save masks (significant)
            save(fullfile(outputFolder,'sig_masks',sprintf('%s_thr%.3f_pos_mask_INTERSECT.mat', behavNames{b}, thr)), ...
                'pos_mask_mat_inter','pos_mask_inter','idx_mask','ii','jj');
            save(fullfile(outputFolder,'sig_masks',sprintf('%s_thr%.3f_neg_mask_INTERSECT.mat', behavNames{b}, thr)), ...
                'neg_mask_mat_inter','neg_mask_inter','idx_mask','ii','jj');
        else
            % Save masks (non-significant)
            save(fullfile(outputFolder,'nonsig_masks',sprintf('%s_thr%.3f_pos_mask_INTERSECT.mat', behavNames{b}, thr)), ...
                'pos_mask_mat_inter','pos_mask_inter','idx_mask','ii','jj');
            save(fullfile(outputFolder,'nonsig_masks',sprintf('%s_thr%.3f_neg_mask_INTERSECT.mat', behavNames{b}, thr)), ...
                'neg_mask_mat_inter','neg_mask_inter','idx_mask','ii','jj');
        end
        % --- NEW: summary counts/means for manuscript text
        mean_r_pos = mean(r_vec_full(sel_pos_idx), 'omitnan');
        mean_r_neg = mean(r_vec_full(sel_neg_idx), 'omitnan');
        T_summary = table( ...
            behavNames(b), thr, ...
            sum(pos_mask_inter), sum(neg_mask_inter), ...
            mean_r_pos, mean_r_neg, ...
            (sum(pos_mask_inter)+sum(neg_mask_inter)) / nEdges * 100, ...
            'VariableNames', {'Behavior','Thr','N_Pos','N_Neg','Mean_r_Pos','Mean_r_Neg','PctOfAllEdges'});
        writetable(T_summary, fullfile(outputFolder,'stats', ...
            sprintf('edgeSUMMARY_INTERSECT_%s_thr%.3f.csv', behavNames{b}, thr)));

        % -------- NEW: Compare network-strength vs strongest single edge (Williams/Steiger test) --------
        % We'll do this on TRAIN (same sample) because the two correlations must share the same Y.
        % Use rank-transforms so Pearson == Spearman.

        % Helper for safe rank transform with NaNs
        rankz = @(v) tiedrank(v);  % returns doubles; NaNs propagate

        % Common pieces
        n_obs = N_tr;
        edge_i = ii(:); edge_j = jj(:);

        % Container for a small results table to append
        rows = {};

        % ---------- POSITIVE: network (sumpos_tr) vs strongest + edge ----------
        if any(pos_mask_inter)
            % Strongest positive edge by full-train r
            [~, pos_max_idx_all] = max(r_vec_full);              % strongest r>0 overall
            % If you prefer "strongest among selected edges", use:
            % [~, kpos] = max(r_vec_full(sel_pos_idx)); pos_max_idx_all = sel_pos_idx(kpos);

            x_net  = sumpos_tr;               % network strength (pos)
            x_edge = Xtr(:, pos_max_idx_all); % single strongest + edge
            y      = Ytr;

            % Valid rows (no NaNs)
            V = ~isnan(x_net) & ~isnan(x_edge) & ~isnan(y);
            xn = rankz(x_net(V)); xe = rankz(x_edge(V)); yr = rankz(y(V));
            nv = sum(V);

            if nv >= 6
                % Pearson on ranks (equivalent to Spearman)
                r12 = corr(xn, yr, 'Type','Pearson');
                r13 = corr(xe, yr, 'Type','Pearson');
                r23 = corr(xn, xe, 'Type','Pearson');

                [t_w, p_w] = williams_test_dep_corr(r12, r13, r23, nv);

                % Also save plain Spearman rs (without rank Pearson, just for reporting)
                r_net_spear  = corr(x_net(V), y(V), 'Type','Spearman');
                r_edge_spear = corr(x_edge(V), y(V), 'Type','Spearman');

                rows(end+1,:) = {'Pos', behavNames{b}, thr, nv, ...
                    r_net_spear, r_edge_spear, r23, ...
                    t_w, p_w, ...
                    edge_i(pos_max_idx_all), edge_j(pos_max_idx_all)};
            end
        end

        % ---------- NEGATIVE: network (sumneg_tr) vs strongest – edge ----------
        if any(neg_mask_inter)
            % Most negative edge by full-train r
            [~, neg_min_idx_all] = min(r_vec_full);              % most negative overall
            % If you prefer "strongest among selected edges", use:
            % [~, kneg] = min(r_vec_full(sel_neg_idx)); neg_min_idx_all = sel_neg_idx(kneg);

            x_net  = sumneg_tr;               % network strength (neg)
            x_edge = Xtr(:, neg_min_idx_all); % single strongest – edge
            y      = Ytr;

            V = ~isnan(x_net) & ~isnan(x_edge) & ~isnan(y);
            xn = rankz(x_net(V)); xe = rankz(x_edge(V)); yr = rankz(y(V));
            nv = sum(V);

            if nv >= 6
                r12 = corr(xn, yr, 'Type','Pearson');
                r13 = corr(xe, yr, 'Type','Pearson');
                r23 = corr(xn, xe, 'Type','Pearson');

                [t_w, p_w] = williams_test_dep_corr(r12, r13, r23, nv);

                r_net_spear  = corr(x_net(V), y(V), 'Type','Spearman');
                r_edge_spear = corr(x_edge(V), y(V), 'Type','Spearman');

                rows(end+1,:) = {'Neg', behavNames{b}, thr, nv, ...
                    r_net_spear, r_edge_spear, r23, ...
                    t_w, p_w, ...
                    edge_i(neg_min_idx_all), edge_j(neg_min_idx_all)};
            end
        end

        % Write CSV (append or create)
        if ~isempty(rows)
            T_steiger = cell2table(rows, 'VariableNames', ...
                {'Sign','Behavior','Thr','N', ...
                'r_net_spearman','r_edge_spearman','r_net_edge_pearson_on_ranks', ...
                't_williams','p_williams','edge_i','edge_j'});

            outcsv = fullfile(outputFolder,'stats', ...
                sprintf('steiger_COMPARE_net_vs_edge_INTERSECT_%s_thr%.3f.csv', behavNames{b}, thr));

            if exist(outcsv,'file')
                writetable(T_steiger, outcsv, 'WriteMode','append');
            else
                writetable(T_steiger, outcsv);
            end
        end

    end % thresholds
end % behaviors

disp('Done with LOOCV, test predictions, plots, masks, and CSV exports.');

%% ================= Local helper functions =================

function scatter_with_fit(x, y, ttl, savepath, markerEdge)
% Scatter plot with LS line; Spearman R/p in title; saved to file.
f = figure('Color','w');
scatter(x, y, 50, 'o', 'MarkerEdgeColor', markerEdge, 'LineWidth', 1.2);
hold on; lsline;
xlabel('Predicted'); ylabel('Observed');
[Rtmp, Ptmp] = corr(x(:), y(:), 'Type','Spearman','Rows','complete');
title(sprintf('%s (Spearman R=%.3f, p=%.3g)', ttl, Rtmp, Ptmp));
set(gca,'Box','off');
outdir = fileparts(savepath);
if ~isempty(outdir) && ~exist(outdir, 'dir'); mkdir(outdir); end
saveas(f, savepath); close(f);
end


function [tval, pval] = williams_test_dep_corr(r12, r13, r23, n)
% Williams/Steiger test for comparing two dependent correlations
% that share one variable (Steiger, 1980; also known as Williams' t).
% r12 = corr(X1,Y), r13 = corr(X2,Y), r23 = corr(X1,X2), n = sample size.
% Returns t (df = n-3) and two-sided p-value.
%
% Notes:
% - We recommend computing r's as Pearson on rank-transformed data when
%   your primary analysis uses Spearman (as done in the caller).

    if n < 6
        tval = NaN; pval = NaN; return;
    end

    % Guard numerical issues
    r12 = max(min(r12, 0.999999), -0.999999);
    r13 = max(min(r13, 0.999999), -0.999999);
    r23 = max(min(r23, 0.999999), -0.999999);

    % Williams test statistic (Steiger, 1980; dependent/overlapping case)
    % Reference implementation follows common formulas used in cocor.
    % Compute the denominator pieces
    num = (r12 - r13) * sqrt((n - 3) * (1 + r23));
    den_inner = 2 * (1 - r23) * (1 + r23 - r12^2 - r13^2 - 2*r12*r13*r23);
    den = sqrt(den_inner);

    % Fallback if den is tiny (avoid division by ~0)
    if den <= 0
        tval = NaN; pval = NaN; return;
    end

    tval = num / den;
    df = n - 3;

    % Two-sided p from t approx (for moderate/large n this matches z closely)
    pval = 2 * (1 - tcdf(abs(tval), df));
end
