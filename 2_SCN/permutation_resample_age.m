%% Permutation tests
% Zhao.Yuyao
clc; clear;
% Path to the folder where your CSV files actually reside
folderPath = '.../3_WB_matprep_regoverallave_globZ2_no_MVM_AAL';
% Specify the folder where you want to save the outputs
outputFolder = '.../Outputs/Y0to10_reg_gender_overallave_no_MVM_globZ2_AAL/Permutation_1000';
% Parameters
num_permutation = 1000;
years = [0, 1, 2, 4, 6, 8, 10];
num_years = length(years);

permtest_age = [];
weight_index = {};

% Create the output folder if it doesn't already exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Loop through each year and Calculate Connectivity
for i = 1:num_years-1
    year1 = years(i);
    year2 = years(i+1);
    filename1 = sprintf('T%d_DK_ID_resid_age.csv', i);
    filename2 = sprintf('T%d_DK_ID_resid_age.csv', i+1);
    T1 = readtable(fullfile(folderPath, filename1));
    T2 = readtable(fullfile(folderPath, filename2));
    subjectIDs1 = T1{:,1}; % Extract the Subject ID column
    subjectIDs2 = T2{:,1};
    data1 = table2array(T1(:, 2:end)); % Extract residuals
    data2 = table2array(T2(:, 2:end));
    only_T1 = setdiff(subjectIDs1, subjectIDs2);   % Only T1
    only_T2 = setdiff(subjectIDs2, subjectIDs1);   % Only T2
    both = intersect(subjectIDs1, subjectIDs2);    % Both
    % Number of subjects
    num_T1 = size(data1, 1);
    num_T2 = size(data2, 1);
    num_both = length(both);
    %combine the two sample
    PermAll = [data1; data2];
    num_all = size(PermAll,1);
    Perm1Num =1:num_all;
    % Create weights for permutation
    Weight1_1 = zeros(num_T1, 1);
    Weight1_2 = zeros(num_T2, 1);
    
    % Assign weights:
    % - 0.5 for subjects with both time points (assigned to each group)
    % - 0.25 for subjects with only one time point
    % for the first half
     % Weights for T1 section
    sample1 = datasample(both,round(num_both/2),'Replace',false);
    Weight1_1(ismember(subjectIDs1, sample1)) = 0.5;
    Weight1_1(ismember(subjectIDs1, only_T1)) = 0.25;
    % Weights for T2 section
    sample2 = both(~ismember(both,sample1));
    Weight1_2(ismember(subjectIDs2, sample2)) = 0.5;
    Weight1_2(ismember(subjectIDs2, only_T2)) = 0.25;
    Weight1 = [Weight1_1',Weight1_2'];
    Weight2 = 0.5 * ones(num_all, 1)' - Weight1;
    %PermTest1 for T1
    for x = 1:num_permutation
        if rem(x,2)==0
            Perm1Index = datasample(Perm1Num,num_T1,'Replace',false,'Weights',Weight1);
            for y = 1:num_T1
                resampCT(y,:)=PermAll(Perm1Index(1,y),:);
            end
            Perm1{1,x} = corrcoef(resampCT);
        else
            Perm1Index = datasample(Perm1Num,num_T1,'Replace',false,'Weights',Weight2);
            for y = 1:num_T1
                resampCT(y,:)=PermAll(Perm1Index(1,y),:);
            end
            Perm1{1,x} = corrcoef(resampCT);
        end
    end
    %PermTest2 for T2
    for x = 1:num_permutation
        if rem(x,2)==0
            Perm2Index = datasample(Perm1Num,num_T2,'Replace',false,'Weights',Weight1);
            for y = 1:num_T2
                resampCT(y,:)=PermAll(Perm2Index(1,y),:);
            end
            Perm2{1,x} = corrcoef(resampCT);
        else
            Perm2Index = datasample(Perm1Num,num_T2,'Replace',false,'Weights',Weight2);
            for y = 1:num_T2
                resampCT(y,:)=PermAll(Perm2Index(1,y),:);
            end
            Perm2{1,x} = corrcoef(resampCT);
        end
    end
    weight_index{end+1, 1} = Weight1;
    weight_index{end+1, 1} = Weight2;
    permtest_age = [permtest_age; Perm1; Perm2];
end 
save(fullfile(outputFolder,'permtest_age.mat'),'permtest_age'); 
Weight_table = table(weight_index);
writetable(Weight_table, fullfile(outputFolder,'weight_index.csv'));