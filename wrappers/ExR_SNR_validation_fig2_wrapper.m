% FILL IN THE INFORMATION UP HERE.
% Image processing up to this point was  background removal in ImageJ

%% Populate parameters
params.xystep = 0.1716; %um/voxel in x and y
params.zstep = 0.040; % um/voxel in z
params.conn = 26; %connectivity
params.layerfolders = {'L1','L23','L4'};
params.lowerlim_um = 0.100; %lower limit on size of a synaptic puncta in one linear dimension, in units of uM

%% Create list of folders for analysis

parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\idExmDecrowdingSegs\Decrowding crop image\images\';
pcontents = dir(parentfolder);
for ff = 3:length(pcontents)
    folders{ff-2,1} = pcontents(ff).name;
end

%% Run puncta extraction and signal/volume calculations on post-expansion stained images

for fidx = 1:length(folders)
    data_post(fidx).data = segment_cc_ExR_synapses_cropped(parentfolder,folders{fidx},'post',params);
end

%% Run puncta extraction and signal/volume calculations on pre-expansion stained images

for fidx = 1:length(folders)
    data_pre(fidx).data = segment_cc_ExR_synapses_cropped(parentfolder,folders{fidx},'pre',params);
end

%% Compile results for copy/paste into Graphpad Prism

clear data
proteins = {'bassoon','cacha','homer','psd95','rim','shank','syngap'};
vars = fieldnames(data_pre(1).data);
nreps = 3; %number of biological replicates
nlayers = length(params.layerfolders);
array = zeros(7,nreps*nlayers*2);

for vv = 1:length(vars)-1 %loop through all measured variables, except the protein name
    for ii = 1:length(data_post)
        var = vars(vv+1);
        protein = data_pre(ii).data.protein;
        struct_ind = find(contains(proteins,protein));
        if vv == 1
            disp(protein);
            disp(struct_ind);
        end
        
        if ii<8
            array(struct_ind,1:3) = getfield(data_pre(ii).data,vars{vv+1})';
            array(struct_ind,10:12) = getfield(data_post(ii).data,vars{vv+1})';
        elseif (8 <= ii) && (ii < 15)
            array(struct_ind,4:6) = getfield(data_pre(ii).data,vars{vv+1})';
            array(struct_ind,13:15) = getfield(data_post(ii).data,vars{vv+1})';
        else
            array(struct_ind,7:9) = getfield(data_pre(ii).data,vars{vv+1})';
            array(struct_ind,16:18) = getfield(data_post(ii).data,vars{vv+1})';
        end
    end
    data.(vars{vv+1}) = array;
end

%% Compile results for copy/past into Graphpad Prism, each protein individually

subset_vars = {'total_vol_out','total_vol_in'};%{'signal_out_clean','signal_in_clean'};%%%{'mean_vol_out','mean_vol_in'};%

array = zeros(9,4);
clear data_protein
sampleinds = [1:3;
    1:3;
    1:3;
    1:3];
for ii = 1:length(data_post)
    var1 = subset_vars{1};
    var2 = subset_vars{2};
    
    protein = data_pre(ii).data.protein;
    disp(protein)
    if ii < 8
        array = zeros(9,4);
        array(1:3,1) = getfield(data_pre(ii).data,var1);
        array(1:3,2) = getfield(data_post(ii).data,var1);
        
        array(1:3,3) = getfield(data_pre(ii).data,var2);
        array(1:3,4) = getfield(data_post(ii).data,var2);
    elseif (8 <= ii) && (ii < 15)
        array = data_protein.(protein);
        array(4:6,1) = getfield(data_pre(ii).data,var1);
        array(4:6,2) = getfield(data_post(ii).data,var1);
        
        array(4:6,3) = getfield(data_pre(ii).data,var2);
        array(4:6,4) = getfield(data_post(ii).data,var2);
    else
        array = data_protein.(protein);
        array(7:9,1) = getfield(data_pre(ii).data,var1);
        array(7:9,2) = getfield(data_post(ii).data,var1);
        
        array(7:9,3) = getfield(data_pre(ii).data,var2);
        array(7:9,4) = getfield(data_post(ii).data,var2);
    end
            
    data_protein.(protein) = array;
    clear array
end

