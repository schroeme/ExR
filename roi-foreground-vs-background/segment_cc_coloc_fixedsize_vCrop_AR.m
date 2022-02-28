function data_struct = segment_cc_coloc_fixedsize_vCrop_AR(parentfolder,folder,test_key)
% Creates connected-components based binary segmentation of a
% background-subtracted, registered image

%pre-expansion voxel size
xystep = 0.01144; %um/voxel in x and y
zstep = 0.026667; % um/voxel in z
conn = 26;

%layerfolders = {'L1','L23','L4'};

signal_in_clean = zeros(1,1); %signal inside the synapses, background subtracted version
signal_out_clean = zeros(1,1); %signal inside the synapses, background subtracted version

%Load the subfolder
subfolder = [folder filesep 'masks'];

%look at the images in that folder
fnames = dir([parentfolder subfolder '/Bin-' test_key '*.tif']); %binary from test stain
T = struct2table(fnames);
sortedT = sortrows(T,'bytes');
fnames = table2struct(sortedT);

reffnames = dir([parentfolder subfolder '/Bin-ref*.tif']); %binary from ref stain
T = struct2table(reffnames);
sortedT = sortrows(T,'bytes');
reffnames = table2struct(sortedT);

%rawfnames = dir([parentfolder subfolder '/Raw-' test_key '*_withBG.tif']); %binary from test stain
cleanfnames = dir([parentfolder subfolder '/Raw-' test_key '*.tif']);
T = struct2table(cleanfnames);
sortedT = sortrows(T,'bytes');
cleanfnames = table2struct(sortedT);

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
    otherinfo = splitsname{3};
    split2 = split(otherinfo,"_");
    protein = split2{1};
    AR = split2{2};
    numbersplit = split(split2{3},".");
    number = numbersplit{1};
    
    clear info I
    info=imfinfo([parentfolder subfolder filesep fname]);
    z_length = numel(info);

    for k = 1:z_length
        I(:,:,k) = imread(info(1).Filename, k, 'Info', info);
    end
    dim = size(I);

    clear refbininfo Iref
    refbinfname = reffnames(fidx).name;
    refbininfo=imfinfo([parentfolder subfolder filesep refbinfname]);
    z_length = numel(refbininfo);
    for k = 1:dim(3)
        Iref(:,:,k) = imread(refbininfo(1).Filename, k, 'Info', refbininfo);
    end

    clear info_clean I_clean
    cleanfname = ['Raw-' prepost '-' protein '_' AR '_' number '.tif'];
    info_clean=imfinfo([parentfolder subfolder filesep cleanfname]);
    for k = 1:dim(3)
        I_clean(:,:,k) = imread(info_clean(1).Filename, k, 'Info', info_clean);
    end

    %I(I>0)=1;
    %Iref(Iref>0)=1;
    test_filt = medfilt3(I,[5 5 3]);
    ref_filt = medfilt3(Iref,[5 5 3]);  

    %apply lower filter on both before closing, and upper filter after
    %closing
    %set lower limit on ROI size to minimize noise: .1 um in each dimension
    lowerlim = ceil((.12/xystep)*(.12/xystep)*(.12/zstep));
    If_ref = bwareaopen(ref_filt, lowerlim); %filtration, removes binary objects LESS than this size
    If_test = bwareaopen(test_filt, lowerlim);

%         CCref = bwconncomp(If_ref,conn);
%         %upperlim = 1^3*(1/xystep)*(1/xystep)*(1/zstep); %upper limit of 1.2 um^3 or 1um in each dimension
%         volref = regionprops3(CCref, 'Volume');
%         volref = volref.Volume;

%         CCtest = bwconncomp(If_test,conn);
%         voltest = regionprops3(CCtest, 'Volume');
%         voltest = voltest.Volume;

%         nobjects_test = CCtest.NumObjects;
%         indexes_test = regionprops3(CCtest,'SubarrayIdx');
%         indexes_test = indexes_test.SubarrayIdx;

    %create segmentation image
    ref_seg = If_ref;%zeros(size(I)); %black background
    test_seg = If_test;%zeros(size(I));

%         indexes = regionprops3(CCref,'SubarrayIdx');
%         indexes = indexes.SubarrayIdx;
%         nobjects_ref = CCref.NumObjects;

%         for cc_idx = 1:nobjects_test %loop through each putative synapse and add a bounding box for segmentation
%             if voltest(cc_idx) < upperlim
%                 test_seg(indexes_test{cc_idx,:})=1;
%             end
%         end

    %test_seg = test_seg.*If_test;
    CCtest = bwconncomp(test_seg);
    nobjects_test = CCtest.NumObjects;
    testvols_final = regionprops3(CCtest,'Volume');
    testvols_final = testvols_final.Volume;

%         for cc_idx = 1:nobjects_ref %loop through each putative synapse and add a bounding box for segmentation
%             if volref(cc_idx) < upperlim      
%                 ref_seg(indexes{cc_idx,:})=1;
%             end
%         end

    %test_seg = double(test_filt);
%         CCtest = bwconncomp(test_seg,conn);
%         testvols_final = regionprops3(CCtest, 'Volume');
%         testvols_final = testvols_final.Volume;

%         nobjects_test = CCtest.NumObjects;

    se = strel('disk',6); %structuring element for dilating the reference seg
    Iseg = double(imdilate(ref_seg.*If_ref,se)); %dilate the reference segmentation
    Iseg_inv = 1-Iseg;

    CCfinal = bwconncomp(Iseg,conn); %connected components on filtered merged image
    nobjects_final = CCfinal.NumObjects;
    nCC_final = nobjects_final;

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

        nCC_test = nobjects_test;

        %thresh = prctile(refvols_final,1);

%         combined = Iseg+test_seg;
%         figure();
%         imagesc(combined(:,:,1));

        %mean_in_raw = regionprops3(Iseg,I_raw,'MeanIntensity');
        %signal_in_raw(iidx,1) = mean_in_raw.MeanIntensity;
        %mean_out_raw = regionprops3(Iseg_inv,I_raw,'MeanIntensity');
        %signal_out_raw(iidx,1) = mean_out_raw.MeanIntensity;

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
% signal_in_clean = nanmean(signal_in_clean_syn);
% signal_out_clean = nanmean(signal_out_clean_syn);
% true_pos = nanmean(true_pos_syn);
% total_vol_in = nanmean(total_vol_in_syn);
% total_vol_out = nanmean(total_vol_out_syn);
% false_pos = nanmean(false_pos_syn);
% mean_vol = nanmean(mean_vol_syn);
% mean_vol_in = nanmean(mean_vol_in_syn);
% mean_vol_out = nanmean(mean_vol_out_syn);
% nobjects_test_final = mean(nobjects_test_syn);
% nobjects_ref_final = mean(nobjects_ref_syn);

splits = split(subfolder,"_");
data_struct.protein = splits{2};
data_struct.true_pos = true_pos_syn;
data_struct.false_pos = false_pos_syn;
data_struct.signal_in_clean = signal_in_clean_syn;
data_struct.signal_out_clean = signal_out_clean_syn;
data_struct.mean_vol = mean_vol_syn;
data_struct.mean_vol_in = mean_vol_in_syn;
data_struct.mean_vol_out = mean_vol_out_syn;
data_struct.total_vol_in = total_vol_in_syn;
data_struct.total_vol_out = total_vol_out_syn;
data_struct.nobjects_test = nobjects_test_syn;
data_struct.nobjects_ref = nobjects_ref_syn;

