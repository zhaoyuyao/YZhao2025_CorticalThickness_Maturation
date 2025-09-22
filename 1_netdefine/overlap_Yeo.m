%%  Mapping AAL atlas to Yeo 7-network
% Yuyao Zhao Jul 9 2025
clc;clear;

% --- Load Yeo (Yeo-7) atlas and labels ---
Yeo_file = '.../atlases_nii/2mm/Yeo-7_space-MNI152NLin6_res-2x2x2.nii';
Yeo_vol = niftiread(Yeo_file);

Yeo_labels_file = '.../labels/Yeo-7.csv';
Yeo_labels = readtable(Yeo_labels_file);  % Must have columns: index, name



%% --- Helper function ---
function assign_to_yeo(atlas_vol, atlas_labels, atlas_name, Yeo_vol, networkList)
    output_file = sprintf('.../overlap/%s_to_Yeo7_assignment.txt', atlas_name);
    fid_output  = fopen(output_file, 'w');
    fprintf(fid_output, '%s_Region\tBest_YeoNetwork\tOverlapPercent\tCorticalVoxel\n', atlas_name);

    for r = 1:height(atlas_labels)
        region_idx = atlas_labels{r,1};
        region_name = atlas_labels{r,2}{1};
        this_mask = (atlas_vol == region_idx);
        cortical_mask = (Yeo_vol ~= 0); % restrict to cortex
        total_vox = sum(this_mask(:) & cortical_mask(:));

        if total_vox <= 0
            fprintf('%s region: %s --> Best network: no overlapping\n', atlas_name, region_name);
            fprintf(fid_output, '%s\tno overlapping\t0\n', region_name);
            continue
        end

        bestNetwork = '';
        bestOverlapFraction = 0;

        for iNet = 1:size(networkList,1)
            netIdx = networkList{iNet,1};
            netName = networkList{iNet,2}{1};
            netMask = (Yeo_vol == netIdx);

            overlap_vox = sum(this_mask(:) & netMask(:));
            overlap_fraction = overlap_vox / total_vox;

            if overlap_fraction > bestOverlapFraction
                bestOverlapFraction = overlap_fraction;
                bestNetwork = netName;
            end
        end
        
        fprintf('%s region: %s --> Best network: %s (%.2f%% overlap)\n', atlas_name, region_name, bestNetwork, bestOverlapFraction*100);
        fprintf(fid_output, '%s\t%s\t%.2f\t%.2f\n', region_name, bestNetwork, bestOverlapFraction*100, total_vox);
    end
    fclose(fid_output);
end

%% --- Map AAL atlas ---
AAL_file = '.../atlases/2mm/AAL_space-MNI152NLin6_res-2x2x2.nii';
AAL_vol = niftiread(AAL_file);
AAL_labels_file = '.../labels/AAL.csv';
AAL_labels = readtable(AAL_labels_file); % columns: index, name

assign_to_yeo(AAL_vol, AAL_labels, 'AAL', Yeo_vol, Yeo_labels);
