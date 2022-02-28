function [correlations,autocorr_img1,autocorr_img2] = pairwise_vol_overlap(img1,img2,xshifts,yshifts,zshifts)

correlations=zeros(length(xshifts),length(yshifts),length(zshifts));
autocorr_img1 = zeros(size(correlations));
autocorr_img2 = zeros(size(correlations));



dim = size(img1);
    
    
for xx = 1:length(xshifts)
    disp('Shifting in x...')
    xshift = xshifts(xx);
    for yy = 1:length(yshifts)
        yshift = yshifts(yy);
        for zz = 1:length(zshifts)
            zshift = zshifts(zz);

            %shift img2
            img2_shifted = imtranslate(img2,[xshift yshift zshift]);

            if xshift == 0 && yshift == 0 && zshift > 0
                edge = imtranslate(img2,[0 0 zshift-dim(3)]);
            elseif xshift == 0 && yshift > 0 && zshift > 0
                edge1 = imtranslate(img2,[0 yshift-dim(1) 0]);
                edge2 = imtranslate(img2,[0 0 zshift-dim(3)]);
                edge = edge1 + edge2;
            elseif xshift == 0 && yshift > 0 && zshift == 0
                edge = imtranslate(img2,[0 yshift-dim(1) 0]);
            elseif xshift > 0 && yshift == 0 && zshift == 0
                edge = imtranslate(img2,[xshift-dim(2) 0 0]);
            elseif xshift > 0 && yshift > 0 && zshift == 0
                edge1 = imtranslate(img2,[xshift-dim(2) 0 0]);
                edge2 = imtranslate(img2,[0 yshift-dim(1) 0]);
                edge = edge1 + edge2;
            elseif xshift > 0 && yshift > 0 && zshift > 0
                edge1 = imtranslate(img2,[0 0 zshift-dim(3)]);
                edge2 = imtranslate(img2,[0 yshift-dim(1) 0]);
                edge3 = imtranslate(img2,[xshift-dim(2) 0 0]);
                edge = edge1 + edge2 + edge3;
            elseif xshift == 0 && yshift == 0 && zshift == 0
                edge = zeros(dim);
            end

            img2_shifted = img2_shifted + edge;

            %shift img1
            img1_shifted = imtranslate(img1,[xshift yshift zshift]);

            if xshift == 0 && yshift == 0 && zshift > 0
                edge = imtranslate(img1,[0 0 zshift-dim(3)]);
            elseif xshift == 0 && yshift > 0 && zshift > 0
                edge1 = imtranslate(img1,[0 yshift-dim(1) 0]);
                edge2 = imtranslate(img1,[0 0 zshift-dim(3)]);
                edge = edge1 + edge2;
            elseif zshift == 0 && yshift > 0 && zshift == 0
                edge = imtranslate(img1,[0 yshift-dim(1) 0]);
            elseif xshift > 0 && yshift == 0 && zshift == 0
                edge = imtranslate(img1,[xshift-dim(2) 0 0]);
            elseif xshift > 0 && yshift > 0 && zshift == 0
                edge1 = imtranslate(img1,[xshift-dim(2) 0 0]);
                edge2 = imtranslate(img1,[0 yshift-dim(1)  0]);
                edge = edge1 + edge2;
            elseif xshift > 0 && yshift > 0 && zshift > 0
                edge1 = imtranslate(img1,[0 0 zshift-dim(3)]);
                edge2 = imtranslate(img1,[0 yshift-dim(1) 0]);
                edge3 = imtranslate(img1,[xshift-dim(2) 0 0]);
                edge = edge1 + edge2 + edge3;
            elseif xshift == 0 && yshift == 0 && zshift == 0
                edge = zeros(dim);
            end

            img1_shifted = img1_shifted + edge;
            
            intersection = pre_shifted_bin & I_post_bin;
            vol_overlaps(xx,yy,zz) = nnz(intersection);
            frac_vol_overlap = nnz(intersection)/vol_pre;
            frac_vol_overlaps(xx,yy,zz) = frac_vol_overlap;

        end
    end
end
end