close all; clear all; clc;

dirname = 'A5_homer(rb) pre 546 post 488_shank(gp)647/L1';

for im_idx = 1:1
    clear pre post ref;

    imname = fullfile(dirname, sprintf('%d.tif', im_idx));
    iminfo = imfinfo(imname);
    
    % process each channel
    for channel_idx = 1: size(iminfo, 1)
        % read image
        im = imread(imname, 'tiff', channel_idx);
%         figure; imshow(100*im); % show original image * 100 (intensity)
        
        % determine threshold % ------- IMPORTANT ------- *
        th = 50; %3 * median(im, 'all');
        
        % convert to binary image
        binary_im = im > th;
%         figure; imshow(binary_im); % show binary image
        
        % find connected components
        cc = bwconncomp(binary_im);
        stats = regionprops(cc, 'Area');
        
        % select LARGEST area
        maxAreaIdx = 1;
        maxArea = 0;
        for area_idx = 1: size(cc.PixelIdxList, 2)
            if maxArea < stats(area_idx).Area
                maxArea = stats(area_idx).Area;
                maxAreaIdx = area_idx;
            end
        end
        pixelIdx = cc.PixelIdxList(maxAreaIdx);
        pixelIdx = pixelIdx{1,1};
        
        % convert 1-D index to 2-D
        [y,x] = ind2sub(size(im), pixelIdx);
        
        % visualize selected LARGEST area
        newIm = zeros(size(im,1), size(im,2));
        for pix_idx = 1: size(pixelIdx, 1)
            newIm(y(pix_idx), x(pix_idx)) = 1;
        end
%         figure; imshow(newIm);
        if isfolder(fullfile('result', dirname)) == 0
            mkdir(fullfile('result', dirname));
        end
        imwrite(newIm, fullfile('result', dirname, ...
                            sprintf('%d_%d.tif', im_idx, channel_idx)));
        
        % store each result
        if channel_idx == 1
            pre = newIm;
        elseif channel_idx == 2
            post = newIm;
        elseif channel_idx == 3
            ref = newIm;
        end
    end
    
    % compute score
    denominator = sum(ref);
    score_pre = sum(pre .* ref) / denominator;
    score_post = sum(post .* ref) / denominator;
    
    % display result
    disp('score_pre: '); disp(score_pre * 100);
    disp('score_post: '); disp(score_post * 100);
    
    % write result to excel
    filename = 'result.xlsx';
    writematrix([score_pre, score_post], filename, ...
        'Sheet', 1, 'Range', sprintf('A%d:B%d', im_idx, im_idx));
    
end


