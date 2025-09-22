%% load processed perm results
BootPath = '.../Outputs/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/GlobalNodalMetrics_Yeo';
path = '.../Outputs/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/ttest_Yeo';

load(fullfile(path,'permutation_ttest_results_allMeasures.mat'), ...
    'real_ge','pval_ge','real_deg','pval_deg','real_seg','pval_seg',...
    'timepointPairs');
load(fullfile(BootPath,'metrics_by_network.mat'), ...
        'obs_ge','obs_deg','obs_seg', ...
        'ge_lobes','deg_lobes','seg_lobes');
% === SETTINGS ===
networks = struct( ...
    'Name',  {'Vis','Som','Lim','Pos','Neg','WB'});
years          = [0,1,2,4,6,8,10]; 

output_csv = fullfile(pwd, 'SCN_metrics_long_globZ2_train.csv');

measures = { ...
    struct('name','GE',       'data',ge_lobes,   'pvals',pval_ge), ...
    struct('name','Degree',   'data',deg_lobes,  'pvals',pval_deg), ...
    struct('name','Seg','data',seg_lobes,'pvals',pval_seg)
};

all_rows = {};

% Loop over each measure
for mIdx = 1:numel(measures)
    mName  = measures{mIdx}.name;
    mData  = measures{mIdx}.data;   % dims: years × networks × subjects
    mPvals = measures{mIdx}.pvals;  % dims: pairs × networks
    
    if isempty(mData)
        continue;
    end

    % Correct p-values (FDR) per measure
    valid_pvals = mPvals(~isnan(mPvals));
    if ~isempty(valid_pvals)
        [~, ~, ~, adj_valid_p] = fdr_bh(valid_pvals, 0.05, 'pdep', 'no');
        adj_p = NaN(size(mPvals));
        adj_p(~isnan(mPvals)) = adj_valid_p;
    else
        adj_p = mPvals;
    end
    
    for netIdx = 1:numel(networks)
        netName = networks(netIdx).Name;

        for yIdx = 1:numel(years)
            vals = squeeze(mData(yIdx, netIdx, :));
            mean_val = mean(vals, 'omitnan');
            std_val  = std(vals, 'omitnan');
            
            % p-value from previous year (if exists)
            if yIdx < numel(years)
                % pair index is same as yIdx for consecutive pairs
                pval_here = adj_p(yIdx, netIdx);
            else
                pval_here = NaN;
            end
            
            % Store row
            all_rows(end+1,:) = { ...
                mName, ...         % Measure
                years(yIdx), ...   % Year
                netName, ...       % Network
                mean_val, ...      % Mean
                std_val, ...       % STD
                pval_here ...      % p_value
            };
        end
    end
end

% Convert to table
T = cell2table(all_rows, ...
    'VariableNames', {'Measure','Year','Network','Mean','STD','p_value'});

% Write CSV
writetable(T, output_csv);
fprintf('Exported SCN metrics to: %s\n', output_csv);
