
% Program for plot heatmaps and connectograms for CPM selected edges
%% Yuyao Zhao, Aug 20, 2025

clear; clc;

%% Paths
mask_dir = '.../3_sbMCN/sbMCN_CT_CPM/Outputs_Y0to2/reg_gender_ROIave_no_MVM_WB_AAL_Y0to2_CPM_092425/sig_masks'; % Path to saved mask .mat file
pos_mask_file = fullfile(mask_dir,'composite_score_thr0.050_pos_mask_INTERSECT.mat');
neg_mask_file = fullfile(mask_dir,'composite_score_thr0.050_neg_mask_INTERSECT.mat');
figFolder = '.../3_sbMCN/sbMCN_CT_CPM/R_fig';
behav_name = 'composite_score';
thr = 0.05;
all_cols = 1:78;

% Load mask
load(pos_mask_file);
load(neg_mask_file);





%% Parameters
baseNetworks = struct( ...
    'Name',  {'Vis','Som','Lim','Pos','Neg'}, ...
    'Cols',  {39:52,[1:2,17:18,20,53:54,65:70],[5:6,21:22,27:28,71:72,75:78],[7:14,19,29:30,33:34,55:60],[3:4,15:16,23:26,31:32,35:36,61:64,73:74]} ...
);


% Network colors
net_colors = [
    30 144 255;   % Vis = dodgerblue
     0 139   0;   % Som = chartreuse4
   255 127  80;   % Lim = coral
   205  38  38;   % Pos = firebrick3
    24 116 205;   % Neg = dodgerblue3
] / 255;





% Split every network into Left (odd roi index) and Right (even roi index)
nBase = numel(baseNetworks);
networksLR = struct('Name',{},'Cols',{});
for k = 1:nBase
    cols = baseNetworks(k).Cols(:)';
    Lcols = cols(mod(cols,2)==1);
    Rcols = cols(mod(cols,2)==0);
    networksLR(end+1).Name = [baseNetworks(k).Name '-L']; networksLR(end).Cols = Lcols; %#ok<SAGROW>
    networksLR(end+1).Name = [baseNetworks(k).Name '-R']; networksLR(end).Cols = Rcols; %#ok<SAGROW>
end

% Custom order: first all L, then all R
orderLR = [1:2:numel(networksLR), 2:2:numel(networksLR)];
networksLR = networksLR(orderLR);
grpNames   = {networksLR.Name};
nGrp       = numel(networksLR);




%% heatmaps (collapsed to 5 networks, no L/R split)

% Use the baseNetworks directly (Vis, Som, Lim, Pos, Neg)
grpNames5 = {baseNetworks.Name};
nGrp5     = numel(baseNetworks);

% Build the upper-tri index across all 78 ROIs once
nC = numel(all_cols);
idx_mask = triu(true(nC),1);
[ii, jj] = find(idx_mask);

% Count totals / pos / neg edges per network pair (5x5)
total_edges5   = zeros(nGrp5, nGrp5);
pos_edges_mat5 = zeros(nGrp5, nGrp5);
neg_edges_mat5 = zeros(nGrp5, nGrp5);

for i = 1:nGrp5
    idx_i = baseNetworks(i).Cols;
    for j = 1:nGrp5
        idx_j = baseNetworks(j).Cols;

        % Any edge whose (ii,jj) is in (i-group, j-group) or (j-group, i-group)
        is_pair = (ismember(ii, idx_i) & ismember(jj, idx_j)) | ...
                  (ismember(ii, idx_j) & ismember(jj, idx_i));

        total_edges5(i,j)   = sum(is_pair);
        pos_edges_mat5(i,j) = sum(pos_mask_inter & is_pair);
        neg_edges_mat5(i,j) = sum(neg_mask_inter & is_pair);
    end
end

% ----- Signed proportion heatmap (High-Low)
prop_signed5 = (pos_edges_mat5 - neg_edges_mat5) ./ total_edges5;
prop_signed5(total_edges5==0) = NaN;

% Plot lower triangle only
mask_upper5 = triu(true(size(prop_signed5)), 1);
prop_plot5 = prop_signed5 * 100; 
prop_plot5(mask_upper5) = NaN;

fig1 = figure;
imagesc(prop_plot5, 'AlphaData', ~isnan(prop_plot5));
axis square;             % force axes to equal size (square cells)
axis tight;              % trim whitespace around data

colormap(redblue(45)); colorbar_handle = colorbar;
xticks(1:nGrp5); yticks(1:nGrp5);
xticklabels(grpNames5); yticklabels(grpNames5);
xlabel('Networks'); ylabel('Networks');
title('Selected Edges% (High EF - Low EF)'); set(gca,'XTickLabelRotation',45);
caxis_vals = caxis;
set(colorbar_handle, 'Ticks', caxis_vals, ...
                     'TickLabels', {num2str(caxis_vals(1),'%.0f'), num2str(caxis_vals(2),'%.0f')});

set(gca,'Color','w','TickLength',[0 0]); box off;
set(gca,'FontSize',14)

saveas(fig1, fullfile(figFolder, sprintf('%s_thr%.3f_net5_heatmap.png', behav_name, thr)));

% ----- Pos heatmap (High)
pos_signed5 = pos_edges_mat5 ./ total_edges5; 
pos_signed5(total_edges5==0) = NaN;
pos_plot5 = pos_signed5 * 100; 
pos_plot5(triu(true(size(pos_plot5)),1)) = NaN;

fig2 = figure;
imagesc(pos_plot5, 'AlphaData', ~isnan(pos_plot5));
axis square;             % force axes to equal size (square cells)
axis tight;              % trim whitespace around data
colormap(slanCM('Reds')); colorbar_handle = colorbar;
xticks(1:nGrp5); yticks(1:nGrp5);
xticklabels(grpNames5); yticklabels(grpNames5);
xlabel('Networks'); ylabel('Networks');
title('Selected Edges% (High EF)'); set(gca,'XTickLabelRotation',45);
caxis_vals = caxis;
set(colorbar_handle, 'Ticks', caxis_vals, ...
                     'TickLabels', {num2str(caxis_vals(1),'%.0f'), num2str(caxis_vals(2),'%.0f')});

set(gca,'Color','w','TickLength',[0 0]); box off;
set(gca,'FontSize',20)

saveas(fig2, fullfile(figFolder, sprintf('%s_thr%.3f_net5_pos_heatmap.png', behav_name, thr)));

% ----- Neg heatmap (Low)
neg_signed5 = neg_edges_mat5 ./ total_edges5; 
neg_signed5(total_edges5==0) = NaN;
neg_plot5 = neg_signed5 * 100; 
neg_plot5(triu(true(size(neg_plot5)),1)) = NaN;

fig3 = figure;
imagesc(neg_plot5, 'AlphaData', ~isnan(neg_plot5));
axis square;             % force axes to equal size (square cells)
axis tight;              % trim whitespace around data

colormap(slanCM('Blues')); colorbar_handle = colorbar;
xticks(1:nGrp5); yticks(1:nGrp5);
xticklabels(grpNames5); yticklabels(grpNames5);
xlabel('Networks'); ylabel('Networks');
title('Selected Edges% (Low EF)'); set(gca,'XTickLabelRotation',45);
caxis_vals = caxis;
set(colorbar_handle, 'Ticks', caxis_vals, ...
                     'TickLabels', {num2str(caxis_vals(1),'%.0f'), num2str(caxis_vals(2),'%.0f')});

set(gca,'Color','w','TickLength',[0 0]); box off;
set(gca,'FontSize',20)

saveas(fig3, fullfile(figFolder, sprintf('%s_thr%.3f_net5_neg_heatmap.png', behav_name, thr)));



%% ==================== 4) CONNECTOGRAMS (pos / neg) ====================
% Build adjacency from masks
N  = numel(all_cols);
A_pos = zeros(N,N); A_neg = zeros(N,N);
A_pos(sub2ind([N N], ii(pos_mask_inter), jj(pos_mask_inter))) = 1;  A_pos = A_pos + A_pos.';
A_neg(sub2ind([N N], ii(neg_mask_inter), jj(neg_mask_inter))) = 1;  A_neg = A_neg + A_neg.';

% ---- Order: all Left groups (Vis-L..Neg-L), then all Right groups (Vis-R..Neg-R)
order = []; 
groupIdx = zeros(N,1);         % 1..nGrp

L_groups = 1:nBase;                    % [Vis-L, Som-L, Lim-L, Pos-L, Neg-L]
R_groups = (nBase+1):nGrp;             % [Vis-R, Som-R, Lim-R, Pos-R, Neg-R]

for k = L_groups                      % left hemisphere first
    c = networksLR(k).Cols(:)'; 
    order = [order c]; 
    groupIdx(c) = k;
end
for k = R_groups                      % then right hemisphere
    c = networksLR(k).Cols(:)'; 
    order = [order c]; 
    groupIdx(c) = k;
end

A_pos = A_pos(order, order);
A_neg = A_neg(order, order);
groupIdx = groupIdx(order);

% color index by base network (same colors for L/R of the same network)
nodeBaseIdx = [1:nBase, 1:nBase];      % map group -> base network
nodeBaseIdx = nodeBaseIdx(groupIdx);




% Helper: plot one connectogram
plot_connectogram = @(A, edge_color, title_str, save_name) ...
    local_plot_connectogram_LR(A, groupIdx, nodeBaseIdx, net_colors, {networksLR.Name}, ...
                               title_str, edge_color, fullfile(figFolder, save_name), nBase);

plot_connectogram(A_pos, [0.80 0.15 0.15], '', ...
                  sprintf('%s_thr%.3f_connectogram_pos_LR.png', behav_name, thr));
plot_connectogram(A_neg, [0.20 0.25 0.80], '', ...
                  sprintf('%s_thr%.3f_connectogram_neg_LR.png', behav_name, thr));


%% ==================== Binary mask heatmaps (78x78, white=selected) ====================
% Rebuild A_pos/A_neg if not present
if ~exist('A_pos','var') || ~exist('A_neg','var')
    N = numel(all_cols);
    idx_mask = triu(true(N),1);
    [ii, jj] = find(idx_mask);
    A_pos = zeros(N,N); A_neg = zeros(N,N);
    A_pos(sub2ind([N N], ii(pos_mask), jj(pos_mask))) = 1;  A_pos = A_pos + A_pos.';
    A_neg(sub2ind([N N], ii(neg_mask), jj(neg_mask))) = 1;  A_neg = A_neg + A_neg.';
end
% zero diagonal for clarity
A_pos(1:end+1:end) = 0;
A_neg(1:end+1:end) = 0;

% ---- helper function
function save_bw_mask(M, outpath)
    f = figure('Color','w','Visible','off','MenuBar','none','ToolBar','none');
    ax = axes('Parent',f,'Position',[0 0 1 1]);
    imagesc(M,'Parent',ax);
    axis(ax,'square'); axis(ax,'off'); axis(ax,'tight');
    colormap(ax, gray(2)); caxis(ax,[0 1]);
    % export the axes only (no toolbar, no margins)
    exportgraphics(ax, outpath, 'Resolution', 300);
    close(f);
end

% ---- save positive & negative
save_bw_mask(A_pos, fullfile(figFolder, sprintf('%s_thr%.3f_pos_mask_bw_78x78.png', behav_name, thr)));
save_bw_mask(A_neg, fullfile(figFolder, sprintf('%s_thr%.3f_neg_mask_bw_78x78.png', behav_name, thr)));

%% ==================== LOCAL FUNCTION ====================
function local_plot_connectogram_LR(A, groupIdx, nodeBaseIdx, net_colors, group_names, ttl, edge_rgb, savepath, nBase)
    % groupIdx: 1..(2*nBase), L groups are odd indices, R groups are even indices
    % nodeBaseIdx: 1..nBase (colors), same for L/R

    G = graph(A);
    n = size(A,1);



% ----- ANGLES (mirror; UNIFORM nodes with hemisphere separation)
n  = size(A,1);
nL = n/2;

% gaps between hemispheres (radians). tweak ~0.04–0.12
hemi_gap_top    = 0.08;   % gap around the top (between Vis-L and Vis-R)
hemi_gap_bottom = -0.08;   % gap around the bottom (between Neg-L and Neg-R)

% Left half spans from top+gap/2 to bottom-gap/2 (top -> bottom)
startL =  pi/2 + hemi_gap_top/2;
endL   = 3*pi/2 - hemi_gap_bottom/2;

% Right half spans from top-gap/2 to bottom+gap/2 (top -> bottom)
startR =  pi/2 - hemi_gap_top/2;
endR   = -pi/2 + hemi_gap_bottom/2;

thetaL = linspace(startL, endL, nL+1);   thetaL(end) = [];
thetaR = linspace(startR, endR, n-nL+1); thetaR(end) = [];

theta  = [thetaL, thetaR];
X = cos(theta); 
Y = sin(theta);


    % ----- figure/axes
    f = figure('Color','w'); 
    ax = axes(f); 
    hold(ax,'on');

    % ----- edges only
    p = plot(G, 'XData', X, 'YData', Y, 'NodeLabel', [], ...
                'LineWidth', 0.8, 'EdgeColor', edge_rgb, 'Parent', ax);
    p.EdgeAlpha = 0.6; 
    p.MarkerSize = 4; 
    p.NodeColor = 'none';

    % ----- nodes with explicit RGB
    scatter(ax, X, Y, 28, net_colors(nodeBaseIdx,:), 'filled', 'MarkerEdgeColor','none');


% ----- outer ring (one arc per L/R group), length = node span
r_in  = 1.10; 
r_out = 1.18;
nGrp  = numel(group_names);

for g = 1:nGrp
    idx = find(groupIdx==g);
    if isempty(idx), continue; end

    th_seg = theta(idx);                 % angles for nodes in this group (monotonic)
    % optional small inset so the patch doesn't sit on top of the dots:
    if numel(th_seg) > 1
        pad = 0.25 * min(abs(diff(th_seg)));   % a quarter of the smallest node spacing
    else
        % single-node group; use a fixed tiny arc
        pad = 0.01;
    end
    % Start and end angles for the arc, staying inside the gap
    th_start = th_seg(1)   + sign(th_seg(end)-th_seg(1)) * pad;
    th_end   = th_seg(end) - sign(th_seg(end)-th_seg(1)) * pad;

    th_dense = linspace(th_start, th_end, 120);
    x1 = r_in  * cos(th_dense);  y1 = r_in  * sin(th_dense);
    x2 = r_out * cos(fliplr(th_dense)); 
    y2 = r_out * sin(fliplr(th_dense));

    base = g; if g > nBase, base = g - nBase; end
    patch([x1 x2],[y1 y2], net_colors(base,:), ...
          'EdgeColor','none', 'FaceAlpha',0.9, 'Parent', ax);

    % label at the middle of the node span
    th_mid = 0.5*(th_start + th_end);
    text(1.26*cos(th_mid), 1.26*sin(th_mid), group_names{g}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Rotation', rad2deg(th_mid)-90, 'FontSize',20, 'FontWeight','bold', ...
        'Color', net_colors(base,:), 'Parent', ax);
end


    % ----- clean axes
    axis(ax,'equal'); 
    axis(ax,'off'); 
    set(ax, 'Position', [0.05 0.05 0.90 0.90]); % leave 5% margins

    
if ~isempty(ttl)
    title(ax, ttl, 'FontWeight','bold');
end

    % no legend
    legend(ax,'off');

    % ----- save
    if ~exist(fileparts(savepath),'dir'); mkdir(fileparts(savepath)); end
    exportgraphics(ax, savepath, 'Resolution', 300);
    close(f);
end

