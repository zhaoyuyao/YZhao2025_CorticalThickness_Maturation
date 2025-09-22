%% Plot SCN by year (LH first, then RH; networks within each hemisphere)
folderPath   = '.../3_WB_matprep_regoverallave_globZ2_no_MVM_AAL';
outputFolder = '.../Figures/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/SCN_by_year_new';
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end

fileNames = {'T1_WB_reg_gender_WBCT.csv','T2_WB_reg_gender_WBCT.csv','T3_WB_reg_gender_WBCT.csv', ...
             'T4_WB_reg_gender_WBCT.csv','T5_WB_reg_gender_WBCT.csv','T6_WB_reg_gender_WBCT.csv','T7_WB_reg_gender_WBCT.csv'};

% --- 0) Order: all Left (odd) by networks, then all Right (even) by networks
networks = struct( ...
    'Name',  {'Vis','Som','Lim','Pos','Neg'}, ...
    'Cols',  {39:52,[1:2,17:18,20,53:54,65:70],[5:6,21:22,27:28,71:72,75:78],[7:14,19,29:30,33:34,55:60],[3:4,15:16,23:26,31:32,35:36,61:64,73:74]} ...
);
reorder  = [];
for n = 1:numel(networks)
    cols_n  = networks(n).Cols(:)'; 
    reorder  = [reorder,  sort(cols_n)]; %#ok<AGROW>
end
reorder_idx = unique(reorder, 'stable');

% --- 1) Global color limits (using reordered matrices)
globalMin = Inf; globalMax = -Inf;
for i = 1:numel(fileNames)
    dataTable = readtable(fullfile(folderPath, fileNames{i}));
    res_mat   = table2array(dataTable(:,2:end));
    valid_idx = reorder_idx(reorder_idx >= 1 & reorder_idx <= size(res_mat,1));
    res_ord   = res_mat(valid_idx, valid_idx);
    globalMin = min(globalMin, min(res_ord(:)));
    globalMax = max(globalMax, max(res_ord(:)));
end

% --- 2) Plot (square + transparent figure/axes as you had)
for i = 1:numel(fileNames)
    dataTable = readtable(fullfile(folderPath, fileNames{i}));
    res_mat   = table2array(dataTable(:,2:end));
    valid_idx = reorder_idx(reorder_idx >= 1 & reorder_idx <= size(res_mat,1));
    res_ord   = res_mat(valid_idx, valid_idx);

    fig = figure('Units','inches','Position',[1 1 6 6], ...
                 'Color','none','InvertHardcopy','off');
    ax = axes(fig,'Position',[0.1 0.1 0.8 0.8]);
    imagesc(ax, res_ord);
    colormap(ax, jet);
    colorbar(ax);
    caxis(ax, [globalMin, globalMax]);
    axis(ax, 'square');
    set(ax,'Color','none');

    outName = [fileNames{i}(1:2), '_SCN.png'];
    exportgraphics(fig, fullfile(outputFolder, outName), ...
                   'BackgroundColor','none','ContentType','image');
    close(fig);
end
