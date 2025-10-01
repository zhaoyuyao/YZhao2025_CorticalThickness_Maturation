% run_williams_from_outputs.m
% Read CPM output CSVs and compute Williams/Steiger tests comparing
% correlations of Observed with Pred_Pos, Pred_Neg, and Pred_Combined
% for both TRAIN (LOOCV) and TEST sets.
%
% Yuyao, 2025-09-25

clear; clc;

%% --------- CONFIG (edit paths as needed) ----------
outputFolder = '.../3_sbMCN/sbMCN_CT_CPM/Outputs_Y2to10/reg_gender_ROIave_no_MVM_WB_AAL_Y2to10_CPM_092425';
predDir = fullfile(outputFolder, 'predictions');
statsDir = fullfile(outputFolder, 'stats');
if ~exist(statsDir,'dir'), mkdir(statsDir); end

% Regex to parse behavior and thr from filenames we already write
reTrain = '^train_LOOCV_(?<Behav>.+)_thr(?<Thr>\d+\.\d{3})\.csv$';
reTest  = '^test_FULL_(?<Behav>.+)_thr(?<Thr>\d+\.\d{3})\.csv$';

%% --------- LIST FILES ----------
trainFiles = dir(fullfile(predDir, 'train_LOOCV_*_thr*.csv'));
testFiles  = dir(fullfile(predDir, 'test_FULL_*_thr*.csv'));

if isempty(trainFiles) && isempty(testFiles)
    error('No prediction CSVs found in %s. Run the CPM script first.', predDir);
end

%% --------- UTILITIES ----------
rankz = @(v) tiedrank(v);
safe_corr_pearson = @(a,b) corr(a,b,'Type','Pearson','Rows','pairwise');
safe_corr_spear   = @(a,b) corr(a,b,'Type','Spearman','Rows','pairwise');
do_steiger = @(y,ypA,ypB) steiger_pair(y,ypA,ypB,rankz,safe_corr_pearson);

%% --------- PROCESS TRAIN FILES ----------
for f = 1:numel(trainFiles)
    fn = trainFiles(f).name;
    tok = regexp(fn, reTrain, 'names');
    if isempty(tok), continue; end
    behav = string(tok.Behav);
    thr   = str2double(tok.Thr);

    T = readtable(fullfile(predDir, fn));
    if ~all(ismember({'Observed','Pred_Pos','Pred_Neg','Pred_Combined'}, T.Properties.VariableNames))
        warning('Missing required columns in %s; skipping.', fn);
        continue;
    end

    Y   = T.Observed;
    Ppos= T.Pred_Pos;
    Pneg= T.Pred_Neg;
    Pcmb= T.Pred_Combined;

    % Spearman performances (on all pairwise-available rows)
    [Rs_pos, Ps_pos]   = safe_corr_spear(Ppos, Y);
    [Rs_neg, Ps_neg]   = safe_corr_spear(Pneg, Y);
    [Rs_comb, Ps_comb] = safe_corr_spear(Pcmb, Y);

    % Pairwise Williams/Steiger on rank-Pearson
    [N_pc, rA_pc, rB_pc, rAB_pc, t_pc, p_pc] = do_steiger(Y, Ppos, Pcmb);
    [N_nc, rA_nc, rB_nc, rAB_nc, t_nc, p_nc] = do_steiger(Y, Pneg, Pcmb);
    [N_pn, rA_pn, rB_pn, rAB_pn, t_pn, p_pn] = do_steiger(Y, Ppos, Pneg);

    % Save TRAIN Steiger
    rows = {
        'pos_vs_comb', N_pc, rA_pc, rB_pc, rAB_pc, t_pc, p_pc
        'neg_vs_comb', N_nc, rA_nc, rB_nc, rAB_nc, t_nc, p_nc
        'pos_vs_neg',  N_pn, rA_pn, rB_pn, rAB_pn, t_pn, p_pn
    };
    T_pairs_tr = cell2table(rows, 'VariableNames', ...
        {'Comparison','N','r_y_modelA_rankPearson','r_y_modelB_rankPearson','r_modelA_modelB_rankPearson','t_williams','p_williams'});

    T_pairs_tr.Behavior = repmat(behav, height(T_pairs_tr),1);
    T_pairs_tr.Thr      = repmat(thr,   height(T_pairs_tr),1);
    T_pairs_tr.Set      = repmat("Train", height(T_pairs_tr),1);
    T_pairs_tr = movevars(T_pairs_tr, {'Behavior','Thr','Set'}, 'Before',1);

    out_paircsv_tr = fullfile(statsDir, sprintf('modelCOMPARE_steiger_TRAIN_%s_thr%.3f.csv', behav, thr));
    writetable(T_pairs_tr, out_paircsv_tr);

    % Save unified performance summary (TRAIN-only row here)
    T_perf_tr = table(behav, thr, Rs_pos, Ps_pos, Rs_neg, Ps_neg, Rs_comb, Ps_comb, ...
        'VariableNames', {'Behavior','Thr','Train_Rs_Pos','Train_p_Pos','Train_Rs_Neg','Train_p_Neg','Train_Rs_Comb','Train_p_Comb'});
    out_perf = fullfile(statsDir, sprintf('modelCOMPARE_perf_%s_thr%.3f.csv', behav, thr));
    if exist(out_perf,'file')
        % append columns if exist? We'll merge on TEST pass; (re)write train-only
        writetable(T_perf_tr, out_perf);  % write fresh; TEST pass will merge
    else
        writetable(T_perf_tr, out_perf);
    end
end

%% --------- PROCESS TEST FILES ----------
for f = 1:numel(testFiles)
    fn = testFiles(f).name;
    tok = regexp(fn, reTest, 'names');
    if isempty(tok), continue; end
    behav = string(tok.Behav);
    thr   = str2double(tok.Thr);

    T = readtable(fullfile(predDir, fn));
    if ~all(ismember({'Observed','Pred_Pos','Pred_Neg','Pred_Combined'}, T.Properties.VariableNames))
        warning('Missing required columns in %s; skipping.', fn);
        continue;
    end

    Y   = T.Observed;
    Ppos= T.Pred_Pos;
    Pneg= T.Pred_Neg;
    Pcmb= T.Pred_Combined;

    % Spearman performances (TEST)
    [Rs_pos, Ps_pos]   = safe_corr_spear(Ppos, Y);
    [Rs_neg, Ps_neg]   = safe_corr_spear(Pneg, Y);
    [Rs_comb, Ps_comb] = safe_corr_spear(Pcmb, Y);

    % Pairwise Williams/Steiger on rank-Pearson
    [N_pc, rA_pc, rB_pc, rAB_pc, t_pc, p_pc] = do_steiger(Y, Ppos, Pcmb);
    [N_nc, rA_nc, rB_nc, rAB_nc, t_nc, p_nc] = do_steiger(Y, Pneg, Pcmb);
    [N_pn, rA_pn, rB_pn, rAB_pn, t_pn, p_pn] = do_steiger(Y, Ppos, Pneg);

    % Save TEST Steiger
    rows = {
        'pos_vs_comb', N_pc, rA_pc, rB_pc, rAB_pc, t_pc, p_pc
        'neg_vs_comb', N_nc, rA_nc, rB_nc, rAB_nc, t_nc, p_nc
        'pos_vs_neg',  N_pn, rA_pn, rB_pn, rAB_pn, t_pn, p_pn
    };
    T_pairs_te = cell2table(rows, 'VariableNames', ...
        {'Comparison','N','r_y_modelA_rankPearson','r_y_modelB_rankPearson','r_modelA_modelB_rankPearson','t_williams','p_williams'});

    T_pairs_te.Behavior = repmat(behav, height(T_pairs_te),1);
    T_pairs_te.Thr      = repmat(thr,   height(T_pairs_te),1);
    T_pairs_te.Set      = repmat("Test", height(T_pairs_te),1);
    T_pairs_te = movevars(T_pairs_te, {'Behavior','Thr','Set'}, 'Before',1);

    out_paircsv_te = fullfile(statsDir, sprintf('modelCOMPARE_steiger_TEST_%s_thr%.3f.csv', behav, thr));
    writetable(T_pairs_te, out_paircsv_te);

    % Merge TEST perf with TRAIN perf if exists
    out_perf = fullfile(statsDir, sprintf('modelCOMPARE_perf_%s_thr%.3f.csv', behav, thr));
    T_perf_test = table(behav, thr, Rs_pos, Ps_pos, Rs_neg, Ps_neg, Rs_comb, Ps_comb, ...
        'VariableNames', {'Behavior','Thr','Test_Rs_Pos','Test_p_Pos','Test_Rs_Neg','Test_p_Neg','Test_Rs_Comb','Test_p_Comb'});

    if exist(out_perf,'file')
        T_train_only = readtable(out_perf);
        % Defensive merge (Behavior, Thr)
        T_merged = outerjoin(T_train_only, T_perf_test, 'Keys', {'Behavior','Thr'}, 'MergeKeys', true);
        writetable(T_merged, out_perf);
    else
        % If train not processed yet, just write test; train pass will merge later if rerun
        writetable(T_perf_test, out_perf);
    end
end

disp('Done: Williams/Steiger comparisons and performance summaries written to stats/.');

%% ====== Local helpers ======
function [N, rYA, rYB, rAB, tval, pval] = steiger_pair(Y, A, B, rankz, safe_corr_pearson)
    V = ~isnan(Y) & ~isnan(A) & ~isnan(B);
    N = sum(V);
    if N < 6
        rYA = NaN; rYB = NaN; rAB = NaN; tval = NaN; pval = NaN; return;
    end
    Yr = rankz(Y(V));
    Ar = rankz(A(V));
    Br = rankz(B(V));

    rYA = safe_corr_pearson(Yr, Ar);
    rYB = safe_corr_pearson(Yr, Br);
    rAB = safe_corr_pearson(Ar, Br);

    [tval, pval] = williams_test_dep_corr(rYA, rYB, rAB, N);
end
