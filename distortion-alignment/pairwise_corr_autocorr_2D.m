function [correlations,autocorr_img1,autocorr_img2] = pairwise_corr_autocorr_2D(img1,img2,norm,xshifts,yshifts,zshifts)

correlations=zeros(length(xshifts),length(yshifts),length(zshifts));
autocorr_img1 = zeros(size(correlations));
autocorr_img2 = zeros(size(correlations));

img1 = double(img1);
img2 = double(img2);

dim1 = size(img1);
dim2 = size(img2);
    
if strcmp(norm,'minmax')
    img1 = (img1 - min(img1(:))) ./ (max(img1(:)) - min(img1(:)));
    img2 = (img2 - min(img2(:))) ./ (max(img2(:)) - min(img2(:)));
end
    
for xx = 1:length(xshifts)
    disp('Shifting in x...')
    xshift = xshifts(xx);
    for yy = 1:length(yshifts)
        yshift = yshifts(yy);
        
        img2_shifted = imtranslate(img2,[xshift yshift]);

        if xshift == 0 && yshift > 0
            edge = imtranslate(img2,[0 yshift-dim2(1)]);
        elseif xshift > 0 && yshift == 0
            edge = imtranslate(img2,[xshift-dim2(2) 0]);
        elseif xshift == 0 && yshift == 0
            edge = zeros(dim2);
        end

        img2_shifted = img2_shifted + edge;

        %shift img1
        img1_shifted = imtranslate(img1,[xshift yshift]);

        if xshift == 0 && yshift > 0
            edge = imtranslate(img1,[0 yshift-dim1(1)]);
        elseif xshift > 0 && yshift == 0
            edge = imtranslate(img1,[xshift-dim1(2) 0]);
        elseif xshift == 0 && yshift == 0
            edge = zeros(dim1);
        end

        img1_shifted = img1_shifted + edge;

        %correlation between images as a function of shift
        %distance in x,y, and z directions
        raw_corr = corr(double(img2_shifted(:)),double(img1(:)));
        correlations(xx,yy) = raw_corr;

        %normalized autocorrelations (want to compare these)
        autocorr_img1(xx,yy) = corr(img1(:),img1_shifted(:));
        autocorr_img2(xx,yy) = corr(img2(:),img2_shifted(:));

    end
end
end