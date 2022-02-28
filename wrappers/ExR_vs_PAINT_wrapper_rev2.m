%Wrapper for calling functions that calculate ExR vs. PAINT distortions,
%puncta properties, etc.
%Margaret Schroeder 09/07/2021

clear all;

params.parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\ExR_PAINT\unregistered\';

img2_names = {'well1_roi1_PAINT.tif';
    'well1_roi2_PAINT.tif';
    'well2_roi1_PAINT.tif';
    'well2_roi2_PAINT.tif';
    'well4_roi1_PAINT.tif';
    };

img1_names = {'well1_roi1_ExR_reg.tif';
    'well1_roi2_ExR_reg.tif';
    'well2_roi1_ExR_reg.tif';
    'well2_roi2_ExR_reg.tif';
    'well4_roi1_ExR_reg.tif';
    };

img_prefixes = {'well1_roi1';
    'well1_roi2';
    'well2_roi1';
    'well2_roi2';
    'well4_roi1';
    };

%Set the resolutions so we can map back to pixel indices:
XY_res = [0.1625/7.9 0.1625/7.9 0.1625/7.7 0.1625/7.7 0.1625/9.0]; %here divided by measured expansion factors
params.Z_res=1;

%% Sliding autocorrelation method for distortion analysis

params.parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\ExR_PAINT\cropped_rois_singles_unreg\';
%params.xystep = 0.01625; %um/voxel in x and y
params.savefolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\ExR_PAINT\sliding_autocorrelation\';
params.normalization = 'minmax';
params.error_channel = '1'; %channel index, channel on which to calculate error
params.xshifts = 0:19;
params.yshifts = 0:19;
params.zshifts=0;
error = [];

for ff=1:length(img_prefixes)
    params.xystep = XY_res(ff); %need to change this in every iteration now that all the images don't have same resolution
    params.fprefix = img_prefixes{ff};
    error_temp = measure_round_alignment_ExR_PAINT(params);
    error = [error; error_temp];
end


%% Loop through fovs and calculate number of puncta and mean puncta area
npunctas=[];
areas = [];
circs=[];
eccs=[];
eqdis=[];
perims=[];

for ff = 1:length(img1_names)
    params.XY_res = XY_res(ff);
    img1_name = img1_names{ff};
    img2_name = img2_names{ff};
    [npuncta,area,circ,ecc,eqdi,perim] = calculate_punctaprops_ExR_vs_PAINT(img1_name,img2_name,params);
    npunctas = [npunctas; npuncta];
    areas = [areas; area];
    circs = [circs; circ];
    eccs = [eccs; ecc];
    eqdis = [eqdis; eqdi];
    perims = [perims; perim];
end

%% Loop through cropped ROIs of PAIRS 
params.parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\ExR_PAINT\cropped_rois_doubles_unreg\';
params.lowerlim_um = 0.1;

delta_pw_distances = [];
delta_puncta_areas = [];
delta_npunctas = [];
dists_paint = [];
dists_exr = [];
for ff = 1:length(img_prefixes)
    img_prefix = img_prefixes{ff};
    params.XY_res = XY_res(ff);
    [delta_pw_distance,delta_puncta_area,delta_npuncta,dist_paint,dist_exr] = calculate_punctadistance_ExR_vs_PAINT(img_prefix,params);
    delta_pw_distances = [delta_pw_distances; delta_pw_distance];
    delta_puncta_areas = [delta_puncta_areas; delta_puncta_area];
    delta_npunctas = [delta_npunctas;delta_npuncta];
    dists_paint = [dists_paint; dist_paint];
    dists_exr = [dists_exr; dist_exr];
end

%% NOT USED IN SECOND REVISION- Loop through cropped ROIs of SINGLE puncta
% params.parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\ExR_PAINT\cropped_rois_singles_unreg\';
% 
% npunctas=[];
% areas = [];
% circs=[];
% eccs=[];
% eqdis=[];
% perims=[];
% 
% for ff = 1:length(img_prefixes)
%     img_prefix = img_prefixes{ff};
%     [npuncta,area,circ,ecc,eqdi,perim] = calculate_cropped_punctaprops_ExR_vs_PAINT(img_prefix,params);
%     npunctas = [npunctas; npuncta];
%     areas = [areas; area];
%     circs = [circs; circ];
%     eccs = [eccs; ecc];
%     eqdis = [eqdis; eqdi];
%     perims = [perims; perim];
% end

