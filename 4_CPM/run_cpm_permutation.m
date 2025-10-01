% run_cpm_permutation.m
% Permutation tests (default 1000) for CPM:
%  - LOOCV (train): shuffle Y_tr, re-run LOOCV within each iteration, record Spearman r's
%  - TEST (optional): shuffle Y_tr, rebuild masks on train, fit model, predict test (Y_te is NOT shuffled)
%
% per-fold partialcorr feature selection,
% intersection masks (>= prop_intersect of folds), linear regression with covariates.
%
% OUTPUTS (under outputFolder/stats):
%   perm_LOOCV_<Behavior>_thr0.xxx.csv    (iteration-wise r_pos/r_neg/r_comb + mask sizes)
%   perm_TEST_<Behavior>_thr0.xxx.csv     (if DO_TEST_PERM = true)
%   perm_SUMMARY_<Behavior>_thr0.xxx.csv  (true r's, permutation p-values)
%
% 2025-09-25

clear; clc;

%% ===================== CONFIG =====================
% --- Match your latest CPM run ---
trainPath    = '.../3_sbMCN/sbMCN_CT_train/Outputs_Y0to2/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y0to2';
testPath     = '.../3_sbMCN/sbMCN_CT_test/Outputs_Y0to2/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y0to2';
BehavPath    = '.../0_data_behavior';
outputFolder = '.../3_sbMCN/sbMCN_CT_CPM/Outputs_Y0to2/reg_gender_ROIave_no_MVM_WB_AAL_Y0to2_CPM_092425';

if ~exist(fullfile(outputFolder,'stats'),'dir'), mkdir(fullfile(outputFolder,'stats')); end

% Behaviors & thresholds to analyze
behavNames  = {'composite_score'};
% options: 'nt_fic_as_scaled','nt_dccs_as_scaled','nt_ls_as_scaled','nt_ps_as_scaled','ssp_length_scaled',
thresholds  = 0.05;

% Node list used in your run
all_cols    = 1:78;

% Permutation settings
N_PERM             = 1000;
DO_TEST_PERM       = true;   % also compute test-set permutation distribution
PROP_INTERSECT     = 0.80;   % intersection threshold (>= 80% folds)
MIN_TRAIN_ROWS_FIT = 5;      % guard for regression stability
RNG_SEED           = 12345;  % reproducibility

%% ===================== LOAD DATA =====================
rng(RNG_SEED);

trainData = load(fullfile(trainPath,'workspace_z_sorted_no_MVM.mat'), 'MCM_subjects_norm','AveragedResults');
testData  = load(fullfile(testPath, 'workspace_z_sorted_no_MVM.mat'), 'MCM_subjects_norm','AveragedResults');

MCM_train = trainData.MCM_subjects_norm;
MCM_test  = testData.MCM_subjects_norm;
AveragedResults_tr = trainData.AveragedResults;
AveragedResults_te = testData.AveragedResults;

T_tr = readtable(fullfile(BehavPath,'behav_10oldest_train_merged_scaled_no_MVM.csv'));
T_te = readtable(fullfile(BehavPath,'behav_10oldest_test_merged_scaled_no_MVM.csv'));

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
        sval = T_tr{mapTr(sid), "Sex"}; if iscell(sval), sval = sval{1}; end; if isstring(sval), sval = char(sval); end
        Sex_tr(s)  = double(strcmpi(sval,'M'));
        Year_tr(s) = double(T_tr{mapTr(sid), "Year"});
    end
end
for s = 1:N_te
    sid = AveragedResults_te{s}{1};
    SID_te(s) = string(sid);
    if isKey(mapTe, sid)
        Yte_all(s,:) = T_te{mapTe(sid), behavNames};
        sval = T_te{mapTe(sid), "Sex"}; if iscell(sval), sval = sval{1}; end; if isstring(sval), sval = char(sval); end
        Sex_te(s)  = double(strcmpi(sval,'M'));
        Year_te(s) = double(T_te{mapTe(sid), "Year"});
    end
end

% Edge indexing & matrices
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

%% ===================== HELPERS =====================
spearman = @(a,b) corr(a,b,'Type','Spearman','Rows','complete');

function [loocv_pos, loocv_neg, loocv_comb, pos_count, neg_count] = do_loocv_pred(X, y, Cov, thr, MIN_ROWS)
    % Returns LOOCV predictions (pos/neg/combined) and per-edge selection counts
    [N, nEdges] = size(X);
    loocv_pos  = nan(N,1);
    loocv_neg  = nan(N,1);
    loocv_comb = nan(N,1);
    pos_count  = zeros(nEdges,1);
    neg_count  = zeros(nEdges,1);

    for k = 1:N
        train_idx = true(N,1); train_idx(k) = false;
        test_idx  = ~train_idx;

        [r_vec_k, p_vec_k] = partialcorr(X(train_idx,:), y(train_idx), Cov(train_idx,:), 'Rows','complete','Type','Spearman');
        pos_mask_k = (r_vec_k > 0) & (p_vec_k < thr);
        neg_mask_k = (r_vec_k < 0) & (p_vec_k < thr);

        pos_count = pos_count + pos_mask_k(:);
        neg_count = neg_count + neg_mask_k(:);

        if ~any(pos_mask_k | neg_mask_k)
            continue;
        end

        sumpos_tr_k = sum(X(train_idx, pos_mask_k), 2);
        sumneg_tr_k = sum(X(train_idx, neg_mask_k), 2);
        sumpos_ko   = sum(X(test_idx,  pos_mask_k), 2);
        sumneg_ko   = sum(X(test_idx,  neg_mask_k), 2);

        Xc = [sumpos_tr_k, sumneg_tr_k, Cov(train_idx,:), ones(sum(train_idx),1)];
        Xp = [sumpos_tr_k,              Cov(train_idx,:), ones(sum(train_idx),1)];
        Xn = [             sumneg_tr_k, Cov(train_idx,:), ones(sum(train_idx),1)];
        yk = y(train_idx);

        valid_c = all(~isnan(Xc),2) & ~isnan(yk);
        valid_p = all(~isnan(Xp),2) & ~isnan(yk);
        valid_n = all(~isnan(Xn),2) & ~isnan(yk);

        if sum(valid_c) >= MIN_ROWS
            bc = Xc(valid_c,:) \ yk(valid_c);
            if ~any(isnan([sumpos_ko, sumneg_ko, Cov(test_idx,:)]))
                loocv_comb(k) = bc(1)*sumpos_ko + bc(2)*sumneg_ko + bc(3)*Cov(test_idx,1) + bc(4)*Cov(test_idx,2) + bc(5);
            end
        end
        if sum(valid_p) >= MIN_ROWS
            bp = Xp(valid_p,:) \ yk(valid_p);
            if ~any(isnan([sumpos_ko, Cov(test_idx,:)]))
                loocv_pos(k) = bp(1)*sumpos_ko + bp(2)*Cov(test_idx,1) + bp(3)*Cov(test_idx,2) + bp(4);
            end
        end
        if sum(valid_n) >= MIN_ROWS
            bn = Xn(valid_n,:) \ yk(valid_n);
            if ~any(isnan([sumneg_ko, Cov(test_idx,:)]))
                loocv_neg(k) = bn(1)*sumneg_ko + bn(2)*Cov(test_idx,1) + bn(3)*Cov(test_idx,2) + bn(4);
            end
        end
    end
end

function [Ypred_pos, Ypred_neg, Ypred_comb, pos_mask_inter, neg_mask_inter] = fit_train_predict_test(Xtr, ytr, Cov_tr, Xte, Cov_te, thr, prop_intersect, MIN_ROWS)
    % Build intersection masks from LOOCV (under current labels), then fit on full train and predict test
    [~, ~, ~, pos_count, neg_count] = do_loocv_pred(Xtr, ytr, Cov_tr, thr, MIN_ROWS);
    N = size(Xtr,1);
    pos_mask_inter = (pos_count >= ceil(prop_intersect*N));
    neg_mask_inter = (neg_count >= ceil(prop_intersect*N));

    sumpos_tr = sum(Xtr(:, pos_mask_inter), 2);
    sumneg_tr = sum(Xtr(:, neg_mask_inter), 2);
    sumpos_te = sum(Xte(:, pos_mask_inter), 2);
    sumneg_te = sum(Xte(:, neg_mask_inter), 2);

    valid_tr_full = all(~isnan([sumpos_tr, sumneg_tr, Cov_tr, ones(N,1)]),2) & ~isnan(ytr);

    Xc = [sumpos_tr(valid_tr_full), sumneg_tr(valid_tr_full), Cov_tr(valid_tr_full,:), ones(sum(valid_tr_full),1)];
    Xp = [sumpos_tr(valid_tr_full),                              Cov_tr(valid_tr_full,:), ones(sum(valid_tr_full),1)];
    Xn = [                          sumneg_tr(valid_tr_full),   Cov_tr(valid_tr_full,:), ones(sum(valid_tr_full),1)];
    y  = ytr(valid_tr_full);

    if ~isempty(Xc), bc = Xc \ y; else, bc = nan(5,1); end
    if ~isempty(Xp), bp = Xp \ y; else, bp = nan(4,1); end
    if ~isempty(Xn), bn = Xn \ y; else, bn = nan(4,1); end

    Ypred_comb = bc(1)*sumpos_te + bc(2)*sumneg_te + bc(3)*Cov_te(:,1) + bc(4)*Cov_te(:,2) + bc(5);
    Ypred_pos  = bp(1)*sumpos_te + bp(2)*Cov_te(:,1) + bp(3)*Cov_te(:,2) + bp(4);
    Ypred_neg  = bn(1)*sumneg_te + bn(2)*Cov_te(:,1) + bn(3)*Cov_te(:,2) + bn(4);
end

%% ===================== PERMUTATION LOOPS =====================
for b = 1:numel(behavNames)
    Ytr = Ytr_all(:,b);
    Yte = Yte_all(:,b);
    beh = behavNames{b};

    for tIdx = 1:numel(thresholds)
        thr = thresholds(tIdx);
        fprintf('\n=== Permutation for %s, thr=%.3f ===\n', beh, thr);

        % ---------- TRUE statistics ----------
        % LOOCV (true)
        [loocv_pos_true, loocv_neg_true, loocv_comb_true, pos_count_true, neg_count_true] = ...
            do_loocv_pred(Xtr, Ytr, Cov_tr, thr, MIN_TRAIN_ROWS_FIT);

        Vtr_pos  = ~isnan(loocv_pos_true)  & ~isnan(Ytr);
        Vtr_neg  = ~isnan(loocv_neg_true)  & ~isnan(Ytr);
        Vtr_comb = ~isnan(loocv_comb_true) & ~isnan(Ytr);

        r_tr_pos_true  = spearman(loocv_pos_true(Vtr_pos),   Ytr(Vtr_pos));
        r_tr_neg_true  = spearman(loocv_neg_true(Vtr_neg),   Ytr(Vtr_neg));
        r_tr_comb_true = spearman(loocv_comb_true(Vtr_comb), Ytr(Vtr_comb));

        % TEST (true; using intersection under true labels)
        [Ypred_pos_true, Ypred_neg_true, Ypred_comb_true, ~, ~] = ...
            fit_train_predict_test(Xtr, Ytr, Cov_tr, Xte, Cov_te, thr, PROP_INTERSECT, MIN_TRAIN_ROWS_FIT);

        Vte_pos  = ~isnan(Ypred_pos_true)  & ~isnan(Yte);
        Vte_neg  = ~isnan(Ypred_neg_true)  & ~isnan(Yte);
        Vte_comb = ~isnan(Ypred_comb_true) & ~isnan(Yte);

        r_te_pos_true  = spearman(Ypred_pos_true(Vte_pos),   Yte(Vte_pos));
        r_te_neg_true  = spearman(Ypred_neg_true(Vte_neg),   Yte(Vte_neg));
        r_te_comb_true = spearman(Ypred_comb_true(Vte_comb), Yte(Vte_comb));

        % ---------- PERM distributions ----------
        permLOO = nan(N_PERM, 6);  % [r_pos, r_neg, r_comb, npos_mask, nneg_mask, iterNvalid]
        permTES = nan(N_PERM, 3);  % [r_pos, r_neg, r_comb]

        % store true on row 1 if you want (optional)
        % (we'll keep distributions purely from permutations)
        for it = 1:N_PERM
            if mod(it,100)==0, fprintf('  perm %d/%d\n', it, N_PERM); end

            % Shuffle TRAIN labels
            yperm = Ytr(randperm(N_tr));

            % LOOCV under permuted labels
            [loocv_pos, loocv_neg, loocv_comb, pos_count, neg_count] = ...
                do_loocv_pred(Xtr, yperm, Cov_tr, thr, MIN_TRAIN_ROWS_FIT);

            Vp_pos  = ~isnan(loocv_pos)  & ~isnan(yperm);
            Vp_neg  = ~isnan(loocv_neg)  & ~isnan(yperm);
            Vp_comb = ~isnan(loocv_comb) & ~isnan(yperm);

            rp_pos  = spearman(loocv_pos(Vp_pos),   yperm(Vp_pos));
            rp_neg  = spearman(loocv_neg(Vp_neg),   yperm(Vp_neg));
            rp_comb = spearman(loocv_comb(Vp_comb), yperm(Vp_comb));

            permLOO(it,1) = rp_pos;
            permLOO(it,2) = rp_neg;
            permLOO(it,3) = rp_comb;
            permLOO(it,4) = sum(pos_count >= ceil(PROP_INTERSECT*N_tr)); % size of INTERSECTION pos
            permLOO(it,5) = sum(neg_count >= ceil(PROP_INTERSECT*N_tr)); % size of INTERSECTION neg
            permLOO(it,6) = sum(Vp_comb); % N valid rows used for combined

            if DO_TEST_PERM
                % TEST under permuted train labels:
                [Ypred_pos_p, Ypred_neg_p, Ypred_comb_p, ~, ~] = ...
                    fit_train_predict_test(Xtr, yperm, Cov_tr, Xte, Cov_te, thr, PROP_INTERSECT, MIN_TRAIN_ROWS_FIT);

                Vtp_pos  = ~isnan(Ypred_pos_p)  & ~isnan(Yte);
                Vtp_neg  = ~isnan(Ypred_neg_p)  & ~isnan(Yte);
                Vtp_comb = ~isnan(Ypred_comb_p) & ~isnan(Yte);

                permTES(it,1) = spearman(Ypred_pos_p(Vtp_pos),   Yte(Vtp_pos));
                permTES(it,2) = spearman(Ypred_neg_p(Vtp_neg),   Yte(Vtp_neg));
                permTES(it,3) = spearman(Ypred_comb_p(Vtp_comb), Yte(Vtp_comb));
            end
        end

        % ---------- p-values (right-tailed, as in your exemplar) ----------
        % LOOCV
        p_tr_pos  = mean(permLOO(:,1) >= r_tr_pos_true,  'omitnan');
        p_tr_neg  = mean(permLOO(:,2) >= r_tr_neg_true,  'omitnan');
        p_tr_comb = mean(permLOO(:,3) >= r_tr_comb_true, 'omitnan');

        % TEST (if enabled)
        if DO_TEST_PERM
            p_te_pos  = mean(permTES(:,1) >= r_te_pos_true,  'omitnan');
            p_te_neg  = mean(permTES(:,2) >= r_te_neg_true,  'omitnan');
            p_te_comb = mean(permTES(:,3) >= r_te_comb_true, 'omitnan');
        else
            [p_te_pos, p_te_neg, p_te_comb] = deal(NaN);
        end

        % ---------- Save distributions ----------
        T_permLOO = table((1:N_PERM)', permLOO(:,1), permLOO(:,2), permLOO(:,3), permLOO(:,4), permLOO(:,5), permLOO(:,6), ...
            'VariableNames', {'Iter','r_pos','r_neg','r_comb','N_pos_inter','N_neg_inter','N_valid_comb'});
        writetable(T_permLOO, fullfile(outputFolder,'stats', ...
            sprintf('perm_LOOCV_%s_thr%.3f.csv', beh, thr)));

        if DO_TEST_PERM
            T_permTES = table((1:N_PERM)', permTES(:,1), permTES(:,2), permTES(:,3), ...
                'VariableNames', {'Iter','r_pos','r_neg','r_comb'});
            writetable(T_permTES, fullfile(outputFolder,'stats', ...
                sprintf('perm_TEST_%s_thr%.3f.csv', beh, thr)));
        end

        % ---------- Save summary ----------
        T_sum = table( ...
            string(beh), thr, ...
            r_tr_pos_true,  p_tr_pos, ...
            r_tr_neg_true,  p_tr_neg, ...
            r_tr_comb_true, p_tr_comb, ...
            r_te_pos_true,  p_te_pos, ...
            r_te_neg_true,  p_te_neg, ...
            r_te_comb_true, p_te_comb, ...
            'VariableNames', {'Behavior','Thr', ...
                              'Train_r_Pos','Train_p_Pos_perm', ...
                              'Train_r_Neg','Train_p_Neg_perm', ...
                              'Train_r_Comb','Train_p_Comb_perm', ...
                              'Test_r_Pos','Test_p_Pos_perm', ...
                              'Test_r_Neg','Test_p_Neg_perm', ...
                              'Test_r_Comb','Test_p_Comb_perm'});

        writetable(T_sum, fullfile(outputFolder,'stats', ...
            sprintf('perm_SUMMARY_%s_thr%.3f.csv', beh, thr)));

        fprintf(' Done: %s thr=%.3f | LOOCV p (pos/neg/comb)= %.4f/%.4f/%.4f | TEST p= %.4f/%.4f/%.4f\n', ...
            beh, thr, p_tr_pos, p_tr_neg, p_tr_comb, p_te_pos, p_te_neg, p_te_comb);
    end
end

disp('Permutation testing complete.');
