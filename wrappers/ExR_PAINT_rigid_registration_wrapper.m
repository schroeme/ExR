
%specify parentfolder and basename below
parentfolder = 'C:\Users\BoydenLabber\Dropbox (MIT)\BoydenLab\ExR_PAINT\unregistered\';
basenames = {'well1_roi1';
    'well1_roi2';
    'well2_roi1';
    'well2_roi2';
   'well4_roi1'
   };

%set type of registration - either SURF or intensity
regtype = 'intensity';

for ii = 1:length(basenames)
    ExR_rigid_registration_2D(parentfolder,basenames{ii},regtype,0)
end



