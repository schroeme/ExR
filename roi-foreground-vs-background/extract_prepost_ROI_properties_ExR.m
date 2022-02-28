function data = extract_prepost_ROI_properties_ExR(folder,params)
% UPDATE FUNCTION DESCRIPTION HERE
%

%load parameters
parentfolder = params.parentfolder;
savefolder = params.savefolder;
xystep = params.xystep;
zstep = params.zstep;
conn = params.conn;
layerfolders = params.layerfolders;
lowerlim_um = params.lowerlim_um;
xshifts = params.xshifts;
yshifts = params.yshifts;
zshifts = params.zshifts;

for lidx = 1:length(layerfolders) %loop through layers
    %Load the subfolder
    subfolder = [folder filesep layerfolders{lidx} filesep 'masks'];
    disp(subfolder)

    %% get filenames for all images in the folder
    prefnames_bin = dir([parentfolder subfolder '/Bin-pre*.tif']); %binary from test stain
    T = struct2table(prefnames_bin);
    sortedT = sortrows(T,'bytes');
    prefnames_bin = table2struct(sortedT);
    
    postfnames_bin = dir([parentfolder subfolder '/Bin-post*.tif']); %binary from test stain
    T = struct2table(postfnames_bin);
    sortedT = sortrows(T,'bytes');
    postfnames_bin = table2struct(sortedT);
    
    prefnames_raw = dir([parentfolder subfolder '/Raw-pre*_noBG.tif']);
    T = struct2table(prefnames_raw);
    sortedT = sortrows(T,'bytes');
    prefnames_raw = table2struct(sortedT);
    
    postfnames_raw = dir([parentfolder subfolder '/Raw-post*_noBG.tif']);
    T = struct2table(postfnames_raw);
    sortedT = sortrows(T,'bytes');
    postfnames_raw = table2struct(sortedT);
    
    %% loop through each cropped image in the folder
    for fidx = 1:length(prefnames_bin) 
        
        prefname_bin = prefnames_bin(fidx).name;
        postfname_bin = postfnames_bin(fidx).name;
        
        splitsname = split(prefname_bin,"-");
        layer = splitsname{3};
        samplesplit = split(splitsname{4},".");
        sample = samplesplit{1};
        
        %% load the images
        clear info I_pre_bin
        info=imfinfo([parentfolder subfolder filesep prefname_bin]);
        z_length = numel(info);

        for k = 1:z_length
            I_pre_bin(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        dim = size(I_pre_bin);

        clear info I_post_bin
        info=imfinfo([parentfolder subfolder filesep postfname_bin]);
        for k = 1:z_length
            I_post_bin(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end

        clear info I_pre_raw
        prefname_raw = ['Raw-preCh-' layer '-' sample '_noBG.tif'];
        info=imfinfo([parentfolder subfolder filesep prefname_raw]);
        for k = 1:dim(3)
            I_pre_raw(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        
        clear info I_post_raw
        prefname_raw = ['Raw-postCh-' layer '-' sample '_noBG.tif'];
        info=imfinfo([parentfolder subfolder filesep prefname_raw]);
        for k = 1:dim(3)
            I_post_raw(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        
        %% Binarize the images and count the number of puncta
        I_pre_filt = medfilt3(I_pre_bin,[5 5 3]);
        I_post_filt = medfilt3(I_post_bin,[5 5 3]);
        
        %apply size filtration
        lowerlim = ceil((lowerlim_um/xystep)*(lowerlim_um/xystep)*(lowerlim_um/zstep));
        I_pre_filt = bwareaopen(I_pre_filt, lowerlim);
        I_post_filt = bwareaopen(I_post_filt, lowerlim);
        
        %find connected components
        CC_pre = bwconncomp(I_pre_filt,conn);
        CC_post = bwconncomp(I_post_filt,conn);
        
        %calculate the number of puncta
        npre = CC_pre.NumObjects;
        n_pre(fidx,lidx) = npre;
        npost = CC_post.NumObjects;
        n_post(fidx,lidx) = npost;
        delta_npuncta(fidx,lidx) = npost-npre;
        delta_npuncta_normalized(fidx,lidx) = (npost-npre)/npre;
        
        if npre > 0 && npost > 0 %proceed only if we have at least one synaptic punctum

            %convert to double
            I_pre_raw=double(I_pre_raw);
            I_post_raw=double(I_post_raw);

            %first try min-max normalization
            [correlations_minmax,autocorr_pre_minmax,autocorr_post_minmax,~,~] = pairwise_corr_autocorr(I_pre_raw,...
                I_post_raw,'minmax',xshifts,yshifts,zshifts);
            
            %then Blanpied method normalization: mean-normalize and mask
            %also pull out the volume overlaps at the same time
            [correlations_blanpied,vol_overlap,vol_overlap_norm] = pairwise_corr_autocorr(I_pre_raw,...
                I_post_raw,'blanpied',xshifts,yshifts,zshifts,I_pre_filt,I_post_filt);
            
            vol_overlaps{fidx,lidx} = vol_overlap;
            vol_overlaps_norm{fidx,lidx}=vol_overlap_norm;

            correlations_norm_minmax{fidx,lidx}=correlations_minmax;
            correlations_norm_blanpied{fidx,lidx}=correlations_blanpied;

            autocorrs_norm_minmax_pre{fidx,lidx}=autocorr_pre_minmax;
            autocorrs_norm_minmax_post{fidx,lidx}=autocorr_post_minmax;

            filename{fidx,lidx} = prefname_bin;
            
        else
            vol_overlaps{fidx,lidx} = NaN;
            vol_overlaps_norm{fidx,lidx}= NaN;

            correlations_norm_minmax{fidx,lidx}= NaN;
            correlations_norm_blanpied{fidx,lidx}= NaN;

            autocorrs_norm_minmax_pre{fidx,lidx}= NaN;
            autocorrs_norm_minmax_post{fidx,lidx}= NaN;
            filename{fidx,lidx} = prefname_bin;
        end
    end
end

splits = split(subfolder,"_");
data.protein = splits{2};
data.n_pre = n_pre;
data.n_post = n_post;
data.delta_npuncta = delta_npuncta;
data.delta_npuncta_norm = delta_npuncta_normalized;

data.vol_overlap = vol_overlaps;
data.frac_vol_overlap = vol_overlaps_norm;

data.correlations_norm_minmax = correlations_norm_minmax;
data.correlation_blanpied = correlations_norm_blanpied;
data.autocorr_pre=autocorrs_norm_minmax_pre;
data.autocorr_post=autocorrs_norm_minmax_post;

data.fname = filename;


end
        
