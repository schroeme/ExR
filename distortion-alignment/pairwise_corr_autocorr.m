function [correlations,autocorr_img1,autocorr_img2,vol_overlap,vol_overlap_norm] = pairwise_corr_autocorr(img1,img2,norm,xshifts,yshifts,zshifts,img1_bin,img2_bin)

correlations=zeros(length(xshifts),length(yshifts),length(zshifts));
autocorr_img1 = zeros(size(correlations));
autocorr_img2 = zeros(size(correlations));
vol_overlap = zeros(size(correlations));
vol_overlap_norm = zeros(size(correlations));

img1 = double(img1);
img2 = double(img2);
    
if strcmp(norm,'minmax')
    img1 = (img1 - min(img1(:))) ./ (max(img1(:)) - min(img1(:)));
    img2 = (img2 - min(img2(:))) ./ (max(img2(:)) - min(img2(:)));
end

if strcmp(norm,'blanpied')
    masked_img1 = double(img1) .* double(img1_bin);
    mean_int_inmask = mean(mean(mean(masked_img1(masked_img1>0))));
    masked_img1_mean = masked_img1;
    masked_img1_mean(masked_img1_mean>0) = mean_int_inmask;

    masked_img2 = double(img2) .* double(img2_bin);
    mean_int_inmask = mean(mean(mean(masked_img2(masked_img2>0))));
    masked_img2_mean = masked_img2;
    masked_img2_mean(masked_img2_mean>0) = mean_int_inmask;

    img1 = double(masked_img1./masked_img1_mean);
    img2 = double(masked_img2./masked_img2_mean);
    img1(isnan(img1))=0;
    img2(isnan(img2))=0;
end
    
for xx = 1:length(xshifts)
    %disp('Shifting in x...')
    xshift = xshifts(xx);
    for yy = 1:length(yshifts)
        yshift = yshifts(yy);
        for zz = 1:length(zshifts)
            zshift = zshifts(zz);

            %shift img2
            img2_shifted = imtranslate(img2,[xshift yshift zshift]);
            edge = find_shift_edge(img2,xshift,yshift,zshift);
            img2_shifted = img2_shifted + edge;

            %shift img1
            img1_shifted = imtranslate(img1,[xshift yshift zshift]);
            edge = find_shift_edge(img1,xshift,yshift,zshift);
            img1_shifted = img1_shifted + edge;
            
            %correlation between images as a function of shift
            %distance in x,y, and z directions
            raw_corr = corr(double(img1(:)),double(img2_shifted(:)));
            correlations(xx,yy,zz) = raw_corr;

            %normalized autocorrelations (want to compare these)
            autocorr_img1(xx,yy,zz) = corr(img1(:),img1_shifted(:));
            autocorr_img2(xx,yy,zz) = corr(img2(:),img2_shifted(:));
            
            %calculate volume overlap and normalized volume overlap
            if (strcmp(norm,'volume')) && (nargin)>6
                img1_shifted_bin = imtranslate(img2_bin,[xshift yshift zshift]);
                edge = find_shift_edge(img2_bin,xshift,yshift,zshift);
                img1_shifted_bin = img1_shifted_bin + edge;
                vol_img1 = nnz(img1_bin);
                
                intersection = img1_shifted_bin & img2_bin;
                vol_overlap(xx,yy,zz) = nnz(intersection);
                frac_vol_overlap = nnz(intersection)/vol_img1;
                vol_overlap_norm(xx,yy,zz) = frac_vol_overlap;
            end

        end
    end
end
end