
% Program for plot group averaged sbMCN matrix (after sbMCN_EF_EBDS_average_Y0to2/Y2to10.m)
%% Yuyao Zhao, Aug 19, 2025



clear;clc;
outputFolder = '.../Outputs_Y0to10/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y0to10';
figureFolder = '.../Figures_Y0to10/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y0to10';

load(fullfile(outputFolder,'workspace_z_sorted_no_MVM.mat'),'MCM_av_normalized');


%% Figure of group netowrk
networks = struct( ...
    'Name',  {'Vis','Som','Lim','Pos','Neg'}, ...
    'Cols',  {39:52,[1:2,17:18,20,53:54,65:70],[5:6,21:22,27:28,71:72,75:78],[7:14,19,29:30,33:34,55:60],[3:4,15:16,23:26,31:32,35:36,61:64,73:74]} ...
);
left_order  = []; right_order = [];
for n = 1:numel(networks)
    cols_n  = networks(n).Cols(:)'; 
    left_order  = [left_order,  sort(cols_n(mod(cols_n,2)==1))]; %#ok<AGROW>
    right_order = [right_order, sort(cols_n(mod(cols_n,2)==0))]; %#ok<AGROW>
end
reorder_idx = unique([left_order, right_order], 'stable');

% Reorder both rows and columns
MCM_av_normalized_plot=MCM_av_normalized(reorder_idx, reorder_idx);

fig = figure('Units','inches','Position',[1 1 6 6], ...
                 'Color','none','InvertHardcopy','off');
ax = axes(fig,'Position',[0.1 0.1 0.8 0.8]);
imagesc(ax, MCM_av_normalized_plot);
colormap(ax, jet);
colorbar(ax);
caxis(ax, [0, 1]);
axis(ax, 'square');
set(ax,'Color','none');

outName = 'group_averaged_sbMCN_z_sort_no_MVM.png';
exportgraphics(fig, fullfile(figureFolder, outName), ...
                   'BackgroundColor','none','ContentType','image');
close(fig);