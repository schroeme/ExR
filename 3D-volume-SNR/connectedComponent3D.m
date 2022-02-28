function newIm = connectedComponent3D(im, th)

% convert to binary image
binary_im = imbinarize(im, th);

% find connected components
cc = bwconncomp(binary_im, 26); %% connectivity
stats = regionprops3(cc, 'Volume', 'VoxelList');

newIm = zeros(size(im), 'logical');

if size(stats,1)
    % select LARGEST volume
    maxVolIdx = 0;
    maxVol = 0;
    for volume_idx = 1: size(stats,1)
        if maxVol < stats(volume_idx, 1).Volume
            maxVol = stats(volume_idx, 1).Volume;
            maxVolIdx = volume_idx;
        end
    end
    
    voxelList = stats(maxVolIdx, 2).VoxelList;
    if iscell(voxelList)
        voxelList = voxelList{:,:};
    end
    
    % visualize selected LARGEST area
    for pix_idx = 1: size(voxelList, 1)
        newIm(voxelList(pix_idx,2), voxelList(pix_idx,1), voxelList(pix_idx,3)) = true;
    end
end
% figure; labelvolshow(newIm);
% figure; sliceViewer(newIm);


end