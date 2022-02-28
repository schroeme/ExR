% WRAPPER FOR CALLING SCRIPTS FOR MEASURING DISTORTION/ERROR BETWEEN PRE-
% AND POST EXPANSION STAINING CHANNELS

%% Populate parameters. Edit this section for your dataset
params.parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\idExmDecrowdingSegs\Decrowding crop image\images\';
params.xystep = 0.01716; %um/voxel in x and y
params.zstep = 0.040; % um/voxel in z
params.conn = 26; %connectivity, used for finding connected components
params.layerfolders = {'L1','L23','L4'};
params.lowerlim_um = 0.100; %lower limit on size of a synaptic puncta in one linear dimension, in units of uM
%directory for saving .mat data files and .fig files containing images
params.savefolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\idExmDecrowdingSegs\Decrowding crop image\distortion\';
params.xshifts = 0:15;
params.yshifts = 0:15;
params.zshifts = 0:8;
params.proteins = {'homer','shank'};
params.proteins_plotnames = {'Homer1','Shank3'};
params.error_range_vals = 1:100; %range of reasonable error values (in nm) to use in error calculation
params.nreplicates = 3; %number of biological replicates

%% Create list of folders for analysis

%Select only protein combinations with good pre-expansion staining (or at
%least where volume/SNR between pre- and post-expansion staining are
%similar for this analysis

%put these folders in order of protein to make concatenation easier later
%on
folders = {
    'A4_homer_shank';
    'B4_homer_shank';
    'C4_homer_shank';
    'A7_shank_homer';
    'B7_shank_homer';
    'C7_shank_homer';    
    };

%% Run the alignment script, can be parallelized if desires

for fidx = 1:length(folders)
    data(fidx).data = extract_prepost_ROI_properties_ExR(folders{fidx},params);
end

% save the data
save([params.savefolder 'prepost_ROI_properties.mat'],'data')

%% Concatenate results from same protein

concat_data = prepost_distortion_concatenate_data_ExR(data,params);

%% Plot volume overlap as a function of shift distance

%Create 2D heatmaps for average and 7 randomly selected synapses of each
%protein
plot_volume_overlap_heatmaps(concat_data,params);

%% Correlation as a function of shift distance

%Create 2D heatmaps for average and 7 randomly selected synapses of each
%protein
plot_distortion_correlation_heatmaps_ExR(concat_data,params);

%% Compare the autocorrelation of pre- and post-expansion channels in x, y, and z

%Create 1D plots for average and 3 randomly selected synapses of each
%protein
plot_correlation_1D_ExR(concat_data,params);

%% Calculate distortion between pre- and post-expansion staining rounds
% This is calculated as the difference between the half-maximum shift
% values of pre-post correlation and pre-pre autocorrelation (or post-post
% autocorrelation) functions. See methods for more details.

%The error cell array will be organized as follows:
% - Rows correspond to proteins, in the order of the proteins in your
% params.proteins and the folders you listed
% - The first column is error with respect to pre-expansion's autocorrelation function
% - The second column is error with respect to post-expansion's
% autocorrelation function
% - The number of values in each matrix is the number of synapses (or ROIs)

error = measure_error_autocorr(concat_data,params);

%% LEGACY CODE BELOW - NOT USED IN MANUSCRIPT, BUT MAY BE OF INTEREST
% %find the distance at which normalized correlation falls to 50% of max
% 
% 
% %Create 2D heatmaps for average and 7 randomly selected synapses of each
% %protein
% selects_index = randi(105,[7,1]);
% xshifts = round([0:19]*17.16);
% yshifts = round([0:19]*17.16);
% zshifts = round([0:8]*40);
% 
% for pp = 1:length(concat_data)
%     mean_vol = nanmean(cat(4, concat_data(pp).correlations_norm_int{:}),4);
%     std_vol = std(cat(4, concat_data(pp).correlations_norm_int{:}),0,4);
%     
%     pointsx = mean_vol(:,1,1);
%     targetx = max(pointsx)*.99;
%     pointsy = mean_vol(1,:,1);
%     targety = max(pointsy)*.99;
%     temp = mean_vol(1,1,:);
%     pointsz = temp(:);
%     targetz = max(pointsz)*.99;
%     
%     dy = @(x,p1,p2,p3) 2*p1*x.^2 + 2*p2*x + p3;
%     
%     fitobjectx = fit(xshifts',pointsx,'poly3');
%     fitx = fitobjectx(range_vals);
%     [~,ind] = min(abs(fitx-targetx));
%     indx(pp,1) = range_vals(ind);
%     
%     fitobjecty = fit(yshifts',pointsy','poly3');
%     fity = fitobjecty(range_vals);
%     [~,ind] = min(abs(fity-targety));
%     indy(pp,1) = range_vals(ind);%range(inds(1));
%     
%     fitobjectz = fit(zshifts',pointsz,'poly3');
%     fitz = fitobjectz(range_vals);
%     [~,ind] = min(abs(fitz-targetz));
%     indz(pp,1) = range_vals(ind);%range(inds(1));
% 
% end