function [diff1,diff2] = find_offset_from_correlations_2D(xshifts,yshifts,autocorr_img1,autocorr_img2,correlations)
%find the error/distance between two images by comparing their correlation
%to the autocorrelation functions

%Find difference between half-maximum shifts of img1-img1, img2-img2, and
%img1-img2
range_vals=0:.1:1000; %in nm

%img1-img1 - find target correlation value in each shift direction
pointsx_11 = autocorr_img1(:,1);
targetx_11 = max(pointsx_11)*.5;
pointsy_11 = autocorr_img1(1,:);
targety_11 = max(pointsy_11)*.5;

%2-2 - find target correlation value in each shift direction
pointsx_22 = autocorr_img2(:,1);
targetx_22 = max(pointsx_22)*.5;
pointsy_22 = autocorr_img2(1,:);
targety_22 = max(pointsy_22)*.5;

%1-2 - find target correlation value 
pointsx_12 = correlations(:,1);
targetx_12 = max(pointsx_12)*.5;
pointsy_12 = correlations(1,:);
targety_12 = max(pointsy_12)*.5;

%1-1
fitobjectx = fit(xshifts',pointsx_11,'poly3');
fitx = fitobjectx(range_vals);
[~,ind] = min(abs(fitx-targetx_11));
inds_11(1,1) = range_vals(ind);

fitobjecty = fit(yshifts',pointsy_11','poly3');
fity = fitobjecty(range_vals);
[~,ind] = min(abs(fity-targety_11));
inds_11(1,2) = range_vals(ind);%range(inds(1));

%2-2
fitobjectx = fit(xshifts',pointsx_22,'poly3');
fitx = fitobjectx(range_vals);
[~,ind] = min(abs(fitx-targetx_22));
inds_22(1,1) = range_vals(ind);

fitobjecty = fit(yshifts',pointsy_22','poly3');
fity = fitobjecty(range_vals);
[~,ind] = min(abs(fity-targety_22));
inds_22(1,2) = range_vals(ind);%range(inds(1));
    

%1-2
fitobjectx = fit(xshifts',pointsx_12,'poly3');
fitx = fitobjectx(range_vals);
[~,ind] = min(abs(fitx-targetx_12));
inds_12(1,1) = range_vals(ind);

fitobjecty = fit(yshifts',pointsy_12','poly3');
fity = fitobjecty(range_vals);
[~,ind] = min(abs(fity-targety_12));
inds_12(1,2) = range_vals(ind);%range(inds(1)); 

diffs_12_11 = inds_12(1,:) - inds_11(1,:);
diffs_12_22 = inds_12(1,:) - inds_22(1,:);

diff1 = mean([diffs_12_11(:)]);
diff2 = mean([diffs_12_22(:)]);

end

