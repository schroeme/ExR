function data_struct = segment_cc_ExR_synapses_cropped(parentfolder,folder,test_key,params)
% UPDATE THIS SECTION
% Creates connected-components based binary segmentation of a
% background-subtracted image and measurse volume and intensity of puncta
% located within and outside of a reference stain.
%
% INPUTS:
% - parentfolder
% - folder
% - test_key: keyword in the filename of the "test" channel to measure with
% respect to the reference channel (e.g., either pre or post expansion
% staining)
% - params: structure containing required parameters, including xystep,
% zstep, conn, and layerfolders
%
% OUTPUTS:
% - data_struct, structure containing measurements for each synapse


%voxel size
xystep = params.xystep;
zstep = params.zstep;
conn = params.conn;
layerfolders=params.layerfolders;
lowerlim_um=params.lowerlim_um;

signal_in_clean = zeros(length(layerfolders),1); %signal inside the synapses, background subtracted version
signal_out_clean = zeros(length(layerfolders),1); %signal inside the synapses, background subtracted version

true_pos = zeros(length(layerfolders),1);
false_pos = zeros(length(layerfolders),1);
nobjects_test_final = zeros(length(layerfolders),1);
nobjects_ref_final = zeros(length(layerfolders),1);

mean_vol = zeros(length(layerfolders),1);

mean_vol_in = zeros(length(layerfolders),1);
total_vol_in = zeros(length(layerfolders),1);

mean_vol_out = zeros(length(layerfolders),1);
total_vol_out = zeros(length(layerfolders),1);

for lidx = 1:length(layerfolders) %loop through layers
    %Load the subfolder
    subfolder = [folder filesep layerfolders{lidx} filesep 'masks'];

    %look at the images in that folder
    fnames = dir([parentfolder subfolder '/Bin-' test_key '*.tif']); %binary from test stain
    T = struct2table(fnames);
    sortedT = sortrows(T,'bytes');
    fnames = table2struct(sortedT);
    
    reffnames = dir([parentfolder subfolder '/Bin-ref*.tif']); %binary from ref stain
    T = struct2table(reffnames);
    sortedT = sortrows(T,'bytes');
    reffnames = table2struct(sortedT);
    
    cleanfnames = dir([parentfolder subfolder '/Raw-' test_key '*_noBG.tif']);
    T = struct2table(cleanfnames);
    sortedT = sortrows(T,'bytes');
    cleanfnames = table2struct(sortedT);
    
    %initialize empty variables
    signal_in_clean_syn = zeros(length(fnames),1);
    signal_out_clean_syn = zeros(length(fnames),1);
    true_pos_syn = zeros(length(fnames),1);
    total_vol_in_syn = zeros(length(fnames),1);
    total_vol_out_syn = zeros(length(fnames),1);
    false_pos_syn = zeros(length(fnames),1);
    mean_vol_syn = zeros(length(fnames),1);
    mean_vol_in_syn = zeros(length(fnames),1);
    mean_vol_out_syn = zeros(length(fnames),1);
    nobjects_test_syn = zeros(length(fnames),1);
    nobjects_ref_syn = zeros(length(fnames),1);
    
    for fidx = 1:length(fnames) %loop through each cropped synapse
        
        fname = fnames(fidx).name;
        splitsname = split(fname,"-");
        prepost = splitsname{2};
        layer = splitsname{3};
        samplesplit = split(splitsname{4},".");
        sample = samplesplit{1};
        
        clear info I
        info=imfinfo([parentfolder subfolder filesep fname]);
        z_length = numel(info);
        
        %load the binary mask file for the test channel
        for k = 1:z_length
            I(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        dim = size(I);

        %load the binary mask file for the reference channel
        clear refbininfo Iref
        refbinfname = reffnames(fidx).name;
        refbininfo=imfinfo([parentfolder subfolder filesep refbinfname]);
        for k = 1:dim(3)
            Iref(:,:,k) = imread(refbininfo(1).Filename, k, 'Info', refbininfo);
        end

        %load the raw (background-subtracted) intensity image file for the
        %test channel
        clear info_clean I_clean
        cleanfname = ['Raw-' prepost '-' layer '-' sample '_noBG.tif'];
        info_clean=imfinfo([parentfolder subfolder filesep cleanfname]);
        for k = 1:dim(3)
            I_clean(:,:,k) = imread(info_clean(1).Filename, k, 'Info', info_clean);
        end

        %apply a median filter to reference and test channels to get rid of
        %high intensity noise (likely from nonspecific staining)
        test_filt = medfilt3(I,[5 5 3]);
        ref_filt = medfilt3(Iref,[5 5 3]);  

        %apply lower filter before morphological closing
        %set lower limit on ROI size to minimize noise
        lowerlim = ceil((lowerlim_um/xystep)*(lowerlim_um/xystep)*(lowerlim_um/zstep));
        If_ref = bwareaopen(ref_filt, lowerlim); %filtration, removes binary objects LESS than this size
        If_test = bwareaopen(test_filt, lowerlim);

        %create segmentation image
        ref_seg = If_ref;
        test_seg = If_test;

        %test_seg = test_seg.*If_test;
        CCtest = bwconncomp(test_seg);
        nobjects_test = CCtest.NumObjects;
        testvols_final = regionprops3(CCtest,'Volume');
        testvols_final = testvols_final.Volume;
        
        se = strel('disk',6); %structuring element for dilating the reference seg
        Iseg = double(imdilate(ref_seg,se)); %dilate the reference segmentation
        Iseg_inv = 1-Iseg;

        CCfinal = bwconncomp(Iseg,conn); %connected components on filtered merged image
        nobjects_final = CCfinal.NumObjects;

        %signal_inside_raw = Iseg.*double(I_raw);
        %signal_outside_raw = Iseg_inv.*double(I_raw);
        nobjects_test_syn(fidx,1) = nobjects_test;
        nobjects_ref_syn(fidx,1) = nobjects_final;
        if nobjects_test > 0 && nobjects_final > 0
            bin_outside = Iseg_inv.*test_seg;
            bin_out_CC = bwconncomp(bin_outside,conn);
            outvols = regionprops3(bin_out_CC, 'Volume');
            outvols = outvols.Volume;

            bin_inside = Iseg.*test_seg;
            bin_in_CC = bwconncomp(bin_inside,conn);
            invols = regionprops3(bin_in_CC, 'Volume');
            invols = invols.Volume;

            mean_in_clean = regionprops3(logical(Iseg),I_clean,'MeanIntensity');
            mean_out_clean = regionprops3(logical(Iseg_inv),I_clean,'MeanIntensity');

            signal_in_clean_syn(fidx,1) = mean(mean_in_clean.MeanIntensity);
            signal_out_clean_syn(fidx,1) = mean(mean_out_clean.MeanIntensity);   

            true_pos_syn(fidx,1) = sum(bin_inside(:))/sum(test_seg(:));
            total_vol_in_syn(fidx,1) = sum(bin_inside(:));
            total_vol_out_syn(fidx,1) = sum(bin_outside(:));
            false_pos_syn(fidx,1) = sum(bin_outside(:))/sum(test_seg(:));
            mean_vol_syn(fidx,1) = mean(testvols_final);

            mean_vol_in_syn(fidx,1) = mean(invols);
            mean_vol_out_syn(fidx,1) = mean(outvols);
        elseif nobjects_test ==0 && nobjects_final > 0 %no test objects, but reference channel ok
            disp(['too few detected objects in test channel for ' folder fname])
            signal_in_clean_syn(fidx,1) = 0;
            signal_out_clean_syn(fidx,1) = 0;   

            true_pos_syn(fidx,1) = 0;
            total_vol_in_syn(fidx,1) = 0;
            total_vol_out_syn(fidx,1) = 0;
            false_pos_syn(fidx,1) = 0;
            mean_vol_syn(fidx,1) = 0;

            mean_vol_in_syn(fidx,1) = 0;
            mean_vol_out_syn(fidx,1) = 0;
        elseif nobjects_test > 0 && nobjects_final == 0
            disp(['too few detected objects in reference channel for ' folder fname])
            signal_in_clean_syn(fidx,1) = nan;
            signal_out_clean_syn(fidx,1) = nan;   

            true_pos_syn(fidx,1) = nan;
            total_vol_in_syn(fidx,1) = nan;
            total_vol_out_syn(fidx,1) = nan;
            false_pos_syn(fidx,1) = nan;
            mean_vol_syn(fidx,1) = nan;

            mean_vol_in_syn(fidx,1) = nan;
            mean_vol_out_syn(fidx,1) = nan;
        else 
            disp(['too few detected objects in reference or test channels for ' folder fname])
            signal_in_clean_syn(fidx,1) = nan;
            signal_out_clean_syn(fidx,1) = nan;   

            true_pos_syn(fidx,1) = nan;
            total_vol_in_syn(fidx,1) = nan;
            total_vol_out_syn(fidx,1) = nan;
            false_pos_syn(fidx,1) = nan;
            mean_vol_syn(fidx,1) = nan;

            mean_vol_in_syn(fidx,1) = nan;
            mean_vol_out_syn(fidx,1) = nan;
        end

    end
    signal_in_clean(lidx,1) = nanmean(signal_in_clean_syn);
    signal_out_clean(lidx,1) = nanmean(signal_out_clean_syn);
    true_pos(lidx,1) = nanmean(true_pos_syn);
    total_vol_in(lidx,1) = nanmean(total_vol_in_syn);
    total_vol_out(lidx,1) = nanmean(total_vol_out_syn);
    false_pos(lidx,1) = nanmean(false_pos_syn);
    mean_vol(lidx,1) = nanmean(mean_vol_syn);
    mean_vol_in(lidx,1) = nanmean(mean_vol_in_syn);
    mean_vol_out(lidx,1) = nanmean(mean_vol_out_syn);
    nobjects_test_final(lidx,1) = mean(nobjects_test_syn);
    nobjects_ref_final(lidx,1) = mean(nobjects_ref_syn);
end

splits = split(subfolder,"_");
data_struct.protein = splits{2};
data_struct.true_pos = true_pos;
data_struct.false_pos = false_pos;
data_struct.signal_in_clean = signal_in_clean;
data_struct.signal_out_clean = signal_out_clean;
data_struct.mean_vol = mean_vol;
data_struct.mean_vol_in = mean_vol_in;
data_struct.mean_vol_out = mean_vol_out;
data_struct.total_vol_in = total_vol_in;
data_struct.total_vol_out = total_vol_out;
data_struct.nobjects_test = nobjects_test_final;
data_struct.nobjects_ref = nobjects_ref_final;

