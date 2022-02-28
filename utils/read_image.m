function [img,nz,QC] = read_image(path,chind,nchannels,se_rad)

    info = imfinfo(path);
    C = strsplit(info(1).ImageDescription,char(10));
    ch_str = strsplit(C{3},"=");
    z_str = strsplit(C{4},"=");
    nc = str2double(ch_str{2});
    nz = str2double(z_str{2});
    
    QC = 0;
    if nc ~= nchannels
        QC = 1;
        img = [];
        nz = [];
    else
        %nz = numel(info)/nchannels;
        img = zeros(info(1).Height,info(1).Width,nz);
        slices = chind:nchannels:numel(info);
        for j = 1:length(slices)
            img(:,:,j) = mat2gray(imread(path, slices(j)));
            if nargin > 3
                se = strel('disk',se_rad);
                img(:,:,j) = imtophat(img(:,:,j),se);
            end
        end

        img = medfilt3(img,[3 3 1]);
    end