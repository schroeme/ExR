function [Images,Imagesflat] = load_flatten_filter2D(parentfolder,folders,normalization,maxes,savename,subtractmean)

for fidx = 1:length(folders)%loop through all folders
    fnames = dir([parentfolder folders{fidx} filesep '*.tif']);
    Iflat = [];
    for imidx = 1:length(fnames) %loop through all files
        fname = fnames(imidx).name;
        clear info
        info=imfinfo([parentfolder folders{fidx} filesep fname]);
        z_length = numel(info);
        %numchunks = z_length/maxes(3);

        Im = [];
        for k = 1:z_length
            A = imread(info(1).Filename, k, 'Info', info);
            for color = 1:size(A,3)
                A(:,:,color) = medfilt2(A(:,:,color),'symmetric');
            end
            Im = cat(4,A,Im);
        end
        
        %take the max projection to create a 2D image
        Amip = max(Im, [], 4);
        I = max(Amip, [], 3); %combine the colors for now
        
        clear A Amip
        Iflat(imidx,:) = I(:);
        Images{imidx,fidx}=I;
    end
    Iflatall{fidx}=Iflat;
    clear Iflat
end

Imagesflat = Iflatall{1};
for fidx = 2:length(folders) %loop through all folders
    Imagesflat = [Imagesflat; Iflatall{fidx}];
end

% remove all zero rows and all zero columns
Imagesflat = Imagesflat(any(Imagesflat,2),:);

%subtract off the mean of each pixel
if subtractmean
    Imagesflat = Imagesflat - mean(Imagesflat,1);
end

% Save them
save(savename,'Images','-v7.3')
