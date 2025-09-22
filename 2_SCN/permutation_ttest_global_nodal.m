%% permutation_ttest_glaobal_nodal.m


%
% Author: Yuyao Zhao (Mar 2025)
% Adjust as needed (file names, variable names, path to BCT, etc.)
%% ==================== 1) CONFIGURATIONS =====================
clear; clc;
addpath(genpath('/Applications/MATLAB_R2025a.app/toolbox/BCT'));  % Adjust BCT path
folderPath = '.../Outputs/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/Permutation_1000';
BootPath = '.../Outputs/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/GlobalNodalMetrics_Yeo';
OutputPath = '.../Outputs/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/ttest_Yeo';
figureFolder = '.../Figures/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/ttest_Yeo';
if ~exist(figureFolder, 'dir'), mkdir(figureFolder); end
if ~exist(OutputPath, 'dir'), mkdir(OutputPath); end
% a) Define your networks

% --- AAL_Yeo ------
networks = struct( ...
    'Name',  {'Vis','Som','Lim','Pos','Neg','WB'}, ...
    'Cols',  {39:52,[1:2,17:18,20,53:54,65:70],[5:6,21:22,27:28,71:72,75:78],[7:14,19,29:30,33:34,55:60],[3:4,15:16,23:26,31:32,35:36,61:64,73:74],1:78}, ...
    'Color', {[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.5 0.5 0.5]} ...
);

% --- Destrieux_Yeo ------
% networks = struct( ...
%     'Name',  {'Vis','Som','Lim','Pos','Neg','WB'}, ...
%     'Cols',  {[2  11  19  20  21  22  42  44  51  57  58  59  61  65  76  85  93  94  95  96 116 118 125 131 132 133 135 139], ...
%     [3   4   8  28  29  33  34  36  41  45  74  77  78  82 102 103 107 108 110 115 119 148], ...
%     [1  23  24  31  32  35  37  43  50  63  75  97  98 105 106 109 111 117 124 137], ...
%     [7  12  15  17  18  26  27  39  40  46  47  48  49  52  53  56  60  62  64  66  67 68  69  81  86  89  91  92 100 101 113 114 120 121 122 123 126 127 130 134 136 138 140 141 142 143], ...
%     [5   6   9  10  13  14  16  25  30  38  54  55  70  71  72  73  79  80  83  84  87 88  90  99 104 112 128 129 144 145 146 147],...
%     1:148}, ...
%     'Color', {[0.3010 0.7450 0.9330],[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840],[0 0.4470 0.7410],[0.5 0.5 0.5]} ...
% );

% b) Timepoints & pairs
years          = [0,1,2,4,6,8,10]; 
timepointPairs = [0 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7];
nPairs         = size(timepointPairs,1);
numNets        = numel(networks);


%% =========== 2) LOAD PERMTEST DATA ========================
%
%   We do so by: 
%       a) For each timepoint pair (T1 vs T2):
%          - Real difference = measure(T2) - measure(T1)
%          - Null distribution from permutations
%          - p-value = fraction(|nullDiff| >= |realDiff|)
%

% a) load the raw adjacency permutations (permtest_age) 
%    => #Rows x #Permutations cell array. 
%    Typically, row1,2 => T1 vs T2, row3,4 => T2 vs T3, etc.
load(fullfile(folderPath,'permtest_age.mat'));  % variable permtest_age

% Load bootstraped data
if exist(fullfile(BootPath,'metrics_by_network.mat'),'file')
    load(fullfile(BootPath,'metrics_by_network.mat'), ...
        'obs_ge','obs_mod','obs_deg','obs_betw','obs_seg');
else
    warning('No nodal_metrics_by_network.mat found in Path.');
    obs_ge = []; obs_deg = []; obs_betw = []; obs_mod = [];obs_seg = [];
end

%% ============= 3) COMPUTE PERMUTATION-BASED PVALUES ===========

% We’ll store the “real” differences and p-values in large arrays:
real_ge  = zeros(nPairs,numNets); 
pval_ge  = zeros(nPairs,numNets);
real_mod   = zeros(nPairs,numNets); 
pval_mod  = zeros(nPairs,numNets); 

real_deg = zeros(nPairs,numNets); 
pval_deg = zeros(nPairs,numNets);
real_betw= zeros(nPairs,numNets);
pval_betw= zeros(nPairs,numNets);
real_seg   = zeros(nPairs,numNets); 
pval_seg  = zeros(nPairs,numNets);

% -- MAIN LOOP --
for pp = 1:nPairs-1
    t1 = timepointPairs(pp,1);
    t2 = timepointPairs(pp,2);
    
    % row indices in permtest_age
    rowA = 2*(pp-1) + 1;
    rowB = rowA + 1;
    
    permsetA = permtest_age(rowA, :);   % 1 x #Perm
    permsetB = permtest_age(rowB, :);   % 1 x #Perm
    numPerm  = length(permsetA);



    % +++++++ Build the null distribution +++++++ 
    geDist   = zeros(numPerm,numNets);
    modDist  = zeros(numPerm,numNets);
    degDist  = zeros(numPerm,numNets);
    betwDist = zeros(numPerm,numNets);
    segDist  = zeros(numPerm,numNets);

    for iPerm = 1:numPerm
        dataA = permsetA{iPerm};
        dataB = permsetB{iPerm};
        
        M_A = computeAllMeasures(dataA, networks);
        M_B = computeAllMeasures(dataB, networks);

        geDist(iPerm,:)   = M_B.ge    - M_A.ge;
        modDist(iPerm,:)  = M_B.mod  - M_A.mod;
        degDist(iPerm,:)  = M_B.deg   - M_A.deg;
        betwDist(iPerm,:) = M_B.betw  - M_A.betw;
        segDist(iPerm,:)  = M_B.seg  - M_A.seg;
    end

    % +++++++ One-sided p-value +++++++ 
    for iNet = 1:numNets
            % real difference = B - A
        real_ge(pp,iNet)   = obs_ge(pp+1,iNet)    - obs_ge(pp,iNet);
        real_mod(pp,iNet) = obs_mod(pp+1,iNet) - obs_mod(pp,iNet);
        real_deg(pp,iNet)  = obs_deg(pp+1,iNet)    - obs_deg(pp,iNet);
        real_betw(pp,iNet) = obs_betw(pp+1,iNet)    - obs_betw(pp,iNet);
        real_seg(pp,iNet) = obs_seg(pp+1,iNet) - obs_seg(pp,iNet);
        
        
        realVal = real_ge(pp,iNet);
        nullVal = geDist(:,iNet);
        pval_ge(pp,iNet) = mean(abs(nullVal) >= abs(realVal));

        realVal = real_mod(pp,iNet);
        nullVal = modDist(:,iNet);
        pval_mod(pp,iNet) = mean(abs(nullVal) >= abs(realVal));

        realVal = real_deg(pp,iNet);
        nullVal = degDist(:,iNet);
        pval_deg(pp,iNet) =mean(abs(nullVal) >= abs(realVal));

        realVal = real_betw(pp,iNet);
        nullVal = betwDist(:,iNet);
        pval_betw(pp,iNet) = mean(abs(nullVal) >= abs(realVal));

        realVal = real_seg(pp,iNet);
        nullVal = segDist(:,iNet);
        pval_seg(pp,iNet) = mean(abs(nullVal) >= abs(realVal));


    end
end

% Save results if you wish
save(fullfile(OutputPath,'permutation_ttest_results_allMeasures.mat'), ...
    'real_ge','pval_ge','real_mod','pval_mod','real_deg','pval_deg','real_betw','pval_betw','real_seg','pval_seg',...
    'timepointPairs');

disp('Done! Permutation test results saved.');

%% ===================== HELPER FUNCTION =====================
% computeAllMeasures: given adjacency (or raw data) and networks, 
% returns the average measure per network
%
% If your input is raw NxROI data, we first do corrcoef, threshold negatives, zero diagonal
% If your input is already NxN adjacency, adjust below as needed.
function M = computeAllMeasures(data, networks)

    % If data is NxROI (raw), convert to NxN adjacency:
    % (We assume if size(data,1)==size(data,2), it might already be adjacency.)
    [r,c] = size(data);
    if r ~= c
        A = corrcoef(data);
    else
        A = ensure_fisher_z(data); 
    end
    
    % Threshold: remove negative edges, set diagonal=0
    A(eye(size(A))>0) = 0;
    A(A<0) = 0;

    
    
    nNets = numel(networks);
    
    % Preallocate
    M.ge   = zeros(1,nNets);
    M.deg  = zeros(1,nNets);
    M.seg  = zeros(1,nNets);

    node_deg_full = strengths_und(A);

    for iNet = 1:nNets
        cols = networks(iNet).Cols;
        subConn = A(cols, cols);

        % 1) Global Efficiency in that subgraph
        M.ge(iNet) = efficiency_wei(subConn);

        % 2) Degree => strengths_und
        M.deg(iNet) = mean(node_deg_full);

        % 3) System Segregation (Fisher-z)
        M.seg(iNet)  = system_segregation_z(A, cols);

    end
end


function seg = system_segregation_z(R_or_Z, cols)
    if numel(cols) < 2, seg = NaN; return; end
    others = setdiff(1:size(R_or_Z,1), cols);
    if isempty(others), seg = NaN; return; end

    Zw = ensure_fisher_z(R_or_Z(cols, cols));
    Zb = ensure_fisher_z(R_or_Z(cols, others));

    n  = numel(cols);
    Zw(1:n+1:end) = 0;
    zw = mean(Zw(triu(true(n),1)), 'all');
    zb = mean(Zb(:));

    if abs(zw) <= eps
        seg = NaN;
    else
        seg = (zw - zb) / zw;
    end
end

function Z = ensure_fisher_z(M)
    if any(abs(M(:)) > 1.0001)   % already z
        Z = M;
    else
        M = max(min(M,0.999999),-0.999999);
        Z = atanh(M);
    end
end