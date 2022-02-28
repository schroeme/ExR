clear all

%% Set folderpaths and parameters for each experiment
% 
% parentpath = 'C:\Users\BoydenLabber\Dropbox (MIT)\Fall 2019\SNB_Rotation_Margaret\2020.01_multiplexing\01.19\';
% roundpaths = {'R01\';'R02\';'R03\';'R04\'};
% regionfolders = {'R\';'L\'};

parentpath ='C:/Users/BoydenLabber/Dropbox (MIT)/BoydenLab/ExR_PAINT/';

targetpath=parentpath;

roundpaths = {
    'PAINT/'
    'ExR/0.05x PBS/MIPs/';
    };

%SET PARAMS
params.xystep =  0.1625/10; %um/pixel, adjusted for expansion factor
params.maxshiftumsmall = 2; %maximum shift between points between rounds, in um
params.pFlag=0; %plotting flag - change to 1 to see the plots
params.nchannels = 1;%5; %number of channels per round
params.nrounds = length(roundpaths);
params.DAPIind = 0;%5; %index of DAPI channel - we exclude this
params.regtype = 'similarity'; %type of registration for alignment AFTER SURF z-crop
params.normint = 1;%1; %whether or not to normalize image histogram during initial feature detection
params.gaussfiltsmall = 0;
params.nOctaves = 6; %number of octaves to use in SURF registration feature detection

fovdirs_last_round = dir([parentpath roundpaths{1} '*.tif*']); %only want fovs that made it all the way to the last round
nfovs = length(fovdirs_last_round);
% fovs = [1 3];

%% Run the registration code
for fovind=1:nfovs
    registration_ExR_2D(parentpath, targetpath, fovind, roundpaths, params);
end
