%% Count the number of cropped synapses in each folder 

clear all
%parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\idExmDecrowdingSegs\Decrowding crop image\';
parentfolder = '/Users/margaretschroeder/Dropbox (MIT)/BoydenLab/idExMDecrowdingSegs/Decrowding crop image/images/';
%% Create list of folders for analysis

pcontents = dir(parentfolder);
for ff = 4:length(pcontents)
    folders{ff-3,1} = pcontents(ff).name;
end

%% get the counts
for fidx = 1:length(folders)
    nsynapses{fidx,1} = count_nsynapses(parentfolder,folders{fidx},'pre');
end
