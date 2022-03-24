%%Open tiffs of separate channel stacks to show mult-color max projection.
%%This script assumes that each channel is separated into it's own stack
%%and is in .tiff format. This script is formatted for 3 channels
%%currently, and can be adjusted to fewer or more channels if desired.

clear;
file0 = ''; %Enter file name here, include proper path if current folder does not include file to be analyzed
mkdir(file0)

%Load stacks for each channel
d0g = double(loadtiff([file0 ''])); %Define which channel is green
d0r = double(loadtiff([file0 ''])); %Define which channel is red
d0b = double(loadtiff([file0 ''])); %Define which channel is blue

%maximal projection along z
d0gmax = double(max(d0g,[],3));   
d0rmax = double(max(d0r,[],3));
d0bmax = double(max(d0b,[],3));
imsizerow = size(d0g,1); imsizecol = size(d0g,2);
d0m = zeros(imsizerow, imsizecol, 3);

%Adjust contrast for each channel
t = max(0,d0rmax-100); t = t/1200; d0m(:,:,1) = t(:,:); %Red channel
t = max(0,d0gmax-150); t = t/3200; d0m(:,:,2) = t(:,:); %Green channel
t = max(0,d0bmax-100); t = t/1000; d0m(:,:,3) = t(:,:); %Blue channel

%Display max projection of all channels
figure; imshow(d0m,'border','tight')

%% Part2: Draw an ROI to pick up a synapse and save it as a new tiff stack

rect = getrect; % draw a box to select a ROI 
rectangle('Position',rect,'EdgeColor','white')
x1 = floor(rect(1)); x2 = floor(rect(1)+rect(3)+1); 
y1 = floor(rect(2)); y2 = floor(rect(2)+rect(4)+1); 
A0 = d0b(y1:y2,x1:x2,:);    %blue channel
B0 = d0r(y1:y2,x1:x2,:);    %red channel
C0 = d0g(y1:y2,x1:x2,:);    %green channel
x0 = x1; y0 = y1;

Amaxx = max(permute(A0,[3 2 1]),[],3);
Bmaxx = max(permute(B0,[3 2 1]),[],3);
Cmaxx = max(permute(C0,[3 2 1]),[],3);
sideview = zeros(size(Amaxx,1),size(Amaxx,2),3);
sideview(:,:,3) = Amaxx(:,:)/max(Amaxx(:));
sideview(:,:,1) = Bmaxx(:,:)/max(Bmaxx(:));
sideview(:,:,2) = Cmaxx(:,:)/max(Cmaxx(:));
figure; imshow(sideview,[100,800])    %show the side view of the synapse
% 
% Draw an ROI to pick up the synapse from sideview
rect = getrect; % draw a box to select a ROI 
rectangle('Position',rect,'EdgeColor','white')
x1 = floor(rect(1)); x2 = floor(rect(1)+rect(3)+1); 
y1 = floor(rect(2)); y2 = floor(rect(2)+rect(4)+1); 
x1 = max(1,x1); y1 = max(1,y1);
x2 = min(x2,size(sideview,2)); y2 = min(y2,size(sideview,1));
A = A0(:,x1:x2,y1:y2); B = B0(:,x1:x2,y1:y2); C = C0(:,x1:x2,y1:y2);
close

outputFileName = fullfile([file0 '\synapse_' num2str(x0) '_' num2str(y0) '.tif'])
A = A - min(A(:)); A = A/max(A(:));
B = B - min(B(:)); B = B/max(B(:));
C = C - min(C(:)); C = C/max(C(:));
for K=1:length(A(1, 1, :))
   imwrite(B(:, :, K), outputFileName, 'WriteMode', 'append', 'Compression','none');
   imwrite(C(:, :, K), outputFileName, 'WriteMode', 'append', 'Compression','none');
   imwrite(A(:, :, K), outputFileName, 'WriteMode', 'append', 'Compression','none');
end  %write the RGB stack of this synapse in the folder, two numbers represent ROI corner 

%% Part3: Read in the saved tiffs and do the cross correlation
% expand the matrix to make voxel cubic
Exp = 15;   %expansion factor, adjust as needed
list0 = dir();
res = [];       %array to keep the results
rmax = 20;      %Maximum shift distance, in pixels (defined below)
for kkk=3:length(list0)
    if list0(kkk).isdir == 0, continue; end
    d0 = list0(kkk).name;
    list = dir([d0 '\*.tif']);
for ii = 1:length(list)
    file1 = fullfile([d0 '\' list(ii).name]);
    test = loadtiff(file1);
    numberofimgs = size(test,3)/3;
    t = 1:numberofimgs;
    A = test(:,:,t*3);  %blue channel
    B = test(:,:,t*3-2); %red channel
    C = test(:,:,t*3-1); %green channel
    
%Subdivide pixel size (physical pixel size of camera chip divided by
%magnification of objective) and step size of stack in order to achieve
%isometric sub-voxels. 
    A1 = expand(A,[2,2,3]); 
    B1 = expand(B,[2,2,3]); 
    C1 = expand(C,[2,2,3]); 
    pixel = 83.3/Exp; %Pixel size defined here, in nm

%Shift clusters to be compared so that they are overlapped. If clusters
%are expected to be colocalized and occupy the same space (and thus require
%no shift), set xyzshift = [0,0,0]
xyzshift = []; 
distance = [20,180]; %Min and max distance allowed for xyz shift, in nm
flag = 0;
cc = get_corr_3dMatrix_final(A1, B1, pixel, rmax, xyzshift, distance, flag);
temp = [kkk, ii, 0, cc];
xyzshift = cc(1:3);
if mean(cc(10:15))~=0
    res = [res;temp];
end
end
end
save('CrossCorrelation.txt', 'res', '-ascii', '-tabs');

%Column 1: Folder number
%Column 2: Synapse number in that particular folder
%Column 3: Blank column
%Columns 4-6: xyz shift
%Column 7: Average intensity of channel 1
%Column 8: Average intensity of channel 2
%Column 9-end: Cross-correlation results
%% Part4: Read in the saved tiffs and do the autocorrelation
% expand the matrix to make voxel cubic
%Autocorrelation is calculated just as the cross-correlation except there
%is no xyzshift
Exp = 15;   %expansion factor, adjust as needed
list0 = dir();
res = [];       %array to keep the results
rmax = 20;      %Maximum shift distance, in pixels (defined below)
for kkk=3:length(list0)
    if list0(kkk).isdir == 0, continue; end
    d0 = list0(kkk).name;
    list = dir([d0 '\*.tif']);
for ii = 1:length(list)
    file1 = fullfile([d0 '\' list(ii).name]);
    test = loadtiff(file1);
    numberofimgs = size(test,3)/3;
    t = 1:numberofimgs;
    A = test(:,:,t*3);  %blue channel, 
    B = test(:,:,t*3-2); %red channel, 
    C = test(:,:,t*3-1); %green channel, 

%Subdivide pixel size (physical pixel size of camera chip divided by
%magnification of objective) and step size of stack in order to achieve
%isometric sub-voxels. 
    A1 = expand(A,[2,2,3]); 
    B1 = expand(B,[2,2,3]); 
    C1 = expand(C,[2,2,3]); 
    pixel = 83.3/Exp;   %Pixel size defined here, in nm

xyzshift = [0,0,0];     %No shift for autocorrelation
flag = 0;
cc = get_corr_3dMatrix_final(A1, A1, pixel, rmax, xyzshift, flag);
temp = [kkk, ii, 0, cc];
res = [res;temp];
end
end
save('AutoCorrelation.txt', 'res', '-ascii', '-tabs');

%Column 1: Folder number
%Column 2: Synapse number in that particular folder
%Column 3: Blank column
%Columns 4: Average intensity
%Column 5: Blank column
%Column 6-end: Auto-correlation results
%% Part5: Read in the saved tiffs and do protein enrichment analysis
% expand the matrix to make voxel cubic
Exp = 15;   %expansion factor, adjust as needed
list0 = dir();
resaa = [];       %array to keep the results
resab = []; resba = []; resbb = []; 
rmax = 180;     %Max distance over which to run the analysis, in nm
for kkk=3:length(list0)
    if list0(kkk).isdir == 0, continue; end
    d0 = list0(kkk).name;
    list = dir([d0 '\*.tif']);
for ii = 1:length(list)
    file1 = fullfile([d0 '\' list(ii).name]);
    test = loadtiff(file1);
    numberofimgs = size(test,3)/3;
    t = 1:numberofimgs;
    A = test(:,:,t*3);  %blue channel, 
    B = test(:,:,t*3-2); %red channel, 
    C = test(:,:,t*3-1); %green channel, 

%Subdivide pixel size (physical pixel size of camera chip divided by
%magnification of objective) and step size of stack in order to achieve
%isometric sub-voxels.     
    A1 = expand(A,[2,2,3]); 
    B1 = expand(B,[2,2,3]); 
    C1 = expand(C,[2,2,3]);
    pixel = 83.3/Exp;   %Pixel size defined here, in nm

%Shift clusters to be compared so that they are overlapped. If clusters
%are expected to be colocalized and occupy the same space (and thus require
%no shift), set xyzshift = [0,0,0]
xyzshift = []; 
distance = [20,180];    %Min and max distance allowed for xyz shift, in nm
flag = 0;
[Raa,Rab,Rba,Rbb] = get_enrichment_3dMatrix_final(A1, B1, pixel, rmax, xyzshift, distance,flag);
aatemp = [kkk,ii,0,Raa'];
abtemp = [kkk,ii,0,Rab'];
batemp = [kkk,ii,0,Rba'];
bbtemp = [kkk,ii,0,Rbb'];
resaa = [resaa;aatemp];    %A density along distance to peak of A
resab = [resab;abtemp];    %A density along distance to peak of B
resba = [resba;batemp];    %B density along distance to peak of A
resbb = [resbb;bbtemp];    %B density along distance to peak of B
end
end
save('AtoApeak.txt', 'resaa', '-ascii', '-tabs');
save('AtoBpeak.txt', 'resab', '-ascii', '-tabs');
save('BtoApeak.txt', 'resba', '-ascii', '-tabs');
save('BtoBpeak.txt', 'resbb', '-ascii', '-tabs');

%Column 1: Folder number
%Column 2: Synapse number in that particular folder
%Column 3: Blank column
%Column 4-end: Enrichment results
