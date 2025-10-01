% Program for constructing single-subject maturational coupling matrices based on
% longitudinal change in cortical thickness
% Sperately for two age ranges: 0-2 and 2-10
% Follow the R script '0_clean_sbMCN_CT.R' under the same folder
% Need to run for train and test set seprately
% Budhachandra Khundrakpam, MNI, 25-11-2017
%% Yuyao Zhao, adapted for EBDS data, Sep 4 2024

clear;
clc;
% BCT path (update path as needed)
addpath(genpath('/Applications/MATLAB_R2025a.app/toolbox/BCT'));
% Path to the folder where your CSV files actually reside (train or test)
folderPath = '...';
% Specify the folder where you want to save the outputs and figures
outputFolder = '.../Outputs_Y2to10/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y2to10';
figureFolder = '.../Figures_Y2to10/reg_gender_ROIave_no_MVM_thresholdZ2_WB_AAL_Y2to10';

% Create the output folder if it doesn't already exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end
%% main loop
% ages 0-2
% Columns = [subjectID, age, gender, ROI1, ..., ROI30];
fileToRead = fullfile(folderPath, 'sbMCN_avg_z_sort_WB_no_MVM_Y2to10.csv');
Y = csvread(fileToRead, 1, 1);

% Get unique subject IDs and their counts
[uniqueSubjects, ~, idx] = unique(Y(:,1));

% Define the conditions for valid subjects
minStartPoints = [2, 4];
maxEndPoints = [8, 10];

% Group time points for each subject
subjectTimePoints = accumarray(idx, Y(:,2), [], @(x) {sort(x')});

% Filter subjects who have at least three time points,
% with the first time point at 2 or 4 years and the final time point at 8 or 10 years
validSubjects = uniqueSubjects(cellfun(@(x) ...
    numel(x) >= 3 && ismember(x(1), minStartPoints) && ismember(x(end), maxEndPoints), subjectTimePoints));

% Filter data for valid subjects
filteredData = Y(ismember(Y(:,1), validSubjects), :);
N_total = size(filteredData, 1);
N_ROIs = size(filteredData, 2) - 3;  % Adjust for ROI columns

% Initialize cell array to store averaged results
AveragedResults = {};  % Each entry will store {subjectID, averaged_cos_theta}
% Initialize an empty cell array to store results
Results_ind = [];
% Loop over valid subjects to calculate and average coupling
for subjID = validSubjects'
    subjData = filteredData(filteredData(:,1) == subjID, :);
    timePoints = subjData(:, 2);
    
    % Initialize a matrix to accumulate cosine matrices
    sum_cos_theta_001_12 = zeros(N_ROIs, N_ROIs);
    pairCount = 0;
    
    for i = 1:length(timePoints)-1
        j = i+1;
            T1 = timePoints(i);
            T2 = timePoints(j);
            CT = subjData([i, j], :);
            disp(['Constructing sbMCN for subject', num2str(subjID), ' between Year ',num2str(T1), ' and Year ', num2str(T2)]);
            % Calculating slopes
            Slope_001_12 = zeros(1, N_ROIs);
            for k = 1:N_ROIs
                Slope_001_12(k) = (CT(2, k+3) - CT(1, k+3)) / (T2 - T1);
            end
            
            % Calculating angles and cosines
            cos_theta_001_12 = zeros(N_ROIs, N_ROIs);  % Cosine matrix
            
            for p = 1:N_ROIs
                for q = 1:N_ROIs
                    m1 = Slope_001_12(p);
                    m2 = Slope_001_12(q);
                    theta_12 = atan(abs((m1 - m2) / (1 + (m1 * m2))));
                    cos_theta_001_12(p, q) = cos(theta_12);
                end
            end
            cos_theta_001_12 = cos_theta_001_12 - eye(N_ROIs, N_ROIs);
            
            % Calculate the row mean values for cos_theta_001_12
            row_mean_values = mean(cos_theta_001_12, 2)';
            mean_values = mean(cos_theta_001_12(:));
            
            % Accumulate the cosine matrix
            sum_cos_theta_001_12 = sum_cos_theta_001_12 + cos_theta_001_12;
            pairCount = pairCount + 1;
            
            % Normalize and save single MCN for each pair of time points
            MCM_ind = cos_theta_001_12;
%            MCM_z = zeros(size(MCM_ind));
%            M = nonzeros(MCM_ind(:));
%            ind = find(MCM_ind ~= 0);
%            MCM_z(ind) = (MCM_ind(ind)-mean(M))/std(M) ;
            
            Matrix_cdf = zeros(N_ROIs,N_ROIs) ;
            for ii = 1:N_ROIs
                for jj = 1:N_ROIs
                    [f,x]   = ecdf([MCM_ind(ii,setdiff(1:N_ROIs,ii))'; MCM_ind(jj,setdiff(1:N_ROIs,ii))']) ;
                    [~,ind] = min(abs(x-MCM_ind(ii,jj))) ;
                    Matrix_cdf(ii,jj) = f(ind) ;
                end
            end
            cos_theta_001_12 = Matrix_cdf;

            % Calculate efficiency and modularity
            ge = efficiency_wei(cos_theta_001_12);
            [community_structure,mod] = modularity_und(cos_theta_001_12);
                        % Node metrics on subW
            node_deg  = strengths_und(cos_theta_001_12);
            node_betw = betweenness_wei(cos_theta_001_12);
            node_eff  = efficiency_wei(cos_theta_001_12, 1); % nodal efficiency
            
            % Average across nodes in this network
            deg_lobes  = mean(node_deg);
            betw_lobes = mean(node_betw);
            eff_lobes  = mean(node_eff);
            % Append result for this pair of time points
            Results_ind = [Results_ind; subjID, T1, T2, ge, mod, mean_values, deg_lobes, betw_lobes, eff_lobes];
            
        
    end
    
    % Calculate the averaged cosine matrix for this subject
    averaged_cos_theta_001_12 = sum_cos_theta_001_12 / pairCount;
    
    % Store the averaged result with the subject ID
    AveragedResults{end+1} = {subjID, averaged_cos_theta_001_12};
end

%% Compute mean MCM for a group
MCM_temp = zeros(N_ROIs,N_ROIs) ;
MCM_sum  = zeros(N_ROIs,N_ROIs) ;
N_subject = length(validSubjects);
for i=1:N_subject
    MCM_temp(:,:) = AveragedResults{1,i}{1,2};
    MCM_sum = MCM_temp + MCM_sum ;
end
MCM_av = MCM_sum/N_subject ;

%% Normalization step for group
N = size(MCM_av,1) ;
Matrix_cdf = zeros(N,N) ;
for i = 1:N    
    for j = 1:N
        [f,x]   = ecdf([MCM_av(i,setdiff(1:N,i))'; MCM_av(j,setdiff(1:N,i))']) ;
        [~,ind] = min(abs(x-MCM_av(i,j))) ;
        Matrix_cdf(i,j) = f(ind) ;
    end
end
MCM_av_normalized = Matrix_cdf ;





dataTable = readtable(fileToRead);

% Extract the headers (variable names) from the 4th to 33rd columns
ROI_labels = dataTable.Properties.VariableNames(5:N_ROIs+4);

% Clean the labels and standardize their case
cleaned_labels = regexprep(ROI_labels, '^X_\d+_(.*?)_\d+$', '$1');


%% Figure of group netowrk
figure;
newOrder = 1:78;
% Reorder both rows and columns
MCM_av_normalized_plot=MCM_av_normalized(newOrder, newOrder);
imagesc(MCM_av_normalized_plot); 
colorbar;

% Set the ticks to correspond to the ROIs
%xticks(1:18);  % Set x-axis ticks at positions 1 to 30
%yticks(1:18);  % Set y-axis ticks at positions 1 to 30

% Assign the labels to the ticks
%xticklabels(cleaned_labels); 
%yticklabels(cleaned_labels);

% Optionally, rotate x-axis labels for better readability
%xtickangle(45);  % Rotate x-axis labels by 45 degrees


% Ensure labels are not interpreted as LaTeX
set(gca, 'TickLabelInterpreter', 'none');

% Save the current figure as a PNG file
saveas(gcf, fullfile(figureFolder,'group_averaged_sbMCN_z_sort_no_MVM.png'));



%% Normalization step for individual
MCM_subjects_norm = cell(1,N_subject);
for sublist = 1:N_subject
    disp(['Individual normalization for subject ',num2str(AveragedResults{1,sublist}{1,1})])
    MCM_ind = AveragedResults{1,sublist}{1,2};
    Matrix_cdf = zeros(N_ROIs,N_ROIs) ;
    for i = 1:N_ROIs
        for j = 1:N_ROIs
            [f,x]   = ecdf([MCM_ind(i,setdiff(1:N_ROIs,i))'; MCM_ind(j,setdiff(1:N_ROIs,i))']) ;
            [~,ind] = min(abs(x-MCM_ind(i,j))) ;
            Matrix_cdf(i,j) = f(ind) ;
        end
    end
    MCM_subjects_norm{1,sublist} = Matrix_cdf ;
end



%% save results
% Convert results to a table for better readability
ResultsTable_ind = array2table(Results_ind, 'VariableNames', {'SubjectID', 'T1', 'T2', 'ge', 'mod','mean','degree','betweenness','loce'});
% Save the results to a CSV file
writetable(ResultsTable_ind, fullfile(outputFolder,'sub_SCN_results_z_sort_no_MVM.csv'));
% Save all variables in the workspace to a file named 'workspace.mat'
save(fullfile(outputFolder,'workspace_z_sorted_no_MVM.mat'));
disp('Results have been saved csvs');

%% violin plot

% Extract the T1 and ge values
T1_values = ResultsTable_ind.T1;
ge_values = ResultsTable_ind.ge;

% Find unique time points
unique_timepoints = unique(T1_values);

% Initialize the figure
figure;
hold on;

% Define a scaling factor to control the width of each 'violin'
scaling_factor = 0.2;  

% Loop over each timepoint and plot the KDE (violin shape)
for i = 1:length(unique_timepoints)
    current_timepoint = unique_timepoints(i);
    
    % Extract ge values for the current timepoint
    ge_current = ge_values(T1_values == current_timepoint);
    
    % Perform kernel density estimation
    [density, xi] = ksdensity(ge_current);
    
    % Normalize the density for plotting
    density = density / max(density) * scaling_factor;
    
    % Plot the violin shape (left and right parts)
    fill(current_timepoint + [density -fliplr(density)], [xi fliplr(xi)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Scatter the actual ge values for this timepoint
    scatter(current_timepoint * ones(size(ge_current)), ge_current, 'k', 'filled');
end

% Calculate the mean ge values for each unique timepoint
mean_ge = arrayfun(@(t) mean(ge_values(T1_values == t)), unique_timepoints);

% Plot the mean values as black circles and connect them with a line
plot(unique_timepoints, mean_ge, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');

% Label the axes
xlabel('Timepoint (T1)');
ylabel('mean MCI');
title('Violin Plot of mean sbMCN by Timepoint (Manual)');

% Optional: Customize plot appearance
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Save the violin plot as an image
saveas(gcf, fullfile(figureFolder,'manual_violin_plot_ge_by_timepoint.png'));
disp('Manual violin plot with mean line saved as manual_violin_plot_ge_by_timepoint.png');
