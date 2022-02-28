% Image processing up to this point was either just background removal in
% ImageJ or co-localization thresholding leading to a merged image

%% Threshold and CC-based
clear all
parentfolder = '/Users/margaretschroeder/Dropbox (MIT)/2020.7.31_antigen_retrieval_crop/';
% pcontents = dir(parentfolder);
folders = {
    'AR_cav';
    'AR_homer';
    'AR_psd95';
    'noAR_cav';
    'noAR_homer';
    'noAR_psd95'
    };

%% Post-expansion staining
%use same connectivity for all right now
for fidx = 1:length(folders)
    data_post(fidx).data = segment_cc_coloc_fixedsize_vCrop_AR(parentfolder,folders{fidx},'post');
end

%% Pre-expansion staining

for fidx = 1:length(folders)
    data_pre(fidx).data = segment_cc_coloc_fixedsize_vCrop_AR(parentfolder,folders{fidx},'pre');
end

%% Compile results for copy/paste

clear data
proteins = {'cav','homer','psd95'};
vars = fieldnames(data_pre(1).data);

for vv = 1:length(vars)-1
    array = zeros(7,12);
    for ii = 1:length(data_post)
        var = vars(vv+1);
        protein = data_pre(ii).data.protein;
        struct_ind = find(contains(proteins,protein));
        
        if sum(array(struct_ind,1:3)) == 0 && sum(array(struct_ind,7:9)) == 0
            array(struct_ind,1:3) = getfield(data_pre(ii).data,vars{vv+1})';
            array(struct_ind,7:9) = getfield(data_post(ii).data,vars{vv+1})';
        else
            array(struct_ind,4:6) = getfield(data_pre(ii).data,vars{vv+1})';
            array(struct_ind,10:12) = getfield(data_post(ii).data,vars{vv+1})';
        end
    end
    data.(vars{vv+1}) = array;
%     if vv == 1
%         data = struct(vars{vv+1},array);
%     else
%         data = [data;struct(vars{vv+1},array)];
%     end
end

%% Compile results, protein-wise

subset_vars = {'signal_out_clean','signal_in_clean'};%{'total_vol_out','total_vol_in'};%%{'mean_vol_out','mean_vol_in'};%%

array = zeros(6,4);
clear data_protein
for ii = 1:length(data_post)
    var1 = subset_vars{1};
    var2 = subset_vars{2};
    
    protein = data_pre(ii).data.protein;
    disp(protein)
    try
        array = data_protein.(protein);
        array(4:6,1) = getfield(data_pre(ii).data,var1);
        array(4:6,2) = getfield(data_post(ii).data,var1);
        
        array(4:6,3) = getfield(data_pre(ii).data,var2);
        array(4:6,4) = getfield(data_post(ii).data,var2);
    catch
        array = zeros(6,4);
        array(1:3,1) = getfield(data_pre(ii).data,var1);
        array(1:3,2) = getfield(data_post(ii).data,var1);
        
        array(1:3,3) = getfield(data_pre(ii).data,var2);
        array(1:3,4) = getfield(data_post(ii).data,var2);
    end
    data_protein.(protein) = array;
    clear array
end

