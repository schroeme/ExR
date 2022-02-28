function diff = find_offset_from_correlations(xshifts,yshifts,zshifts,autocorr_img1,autocorr_img2,correlations)
%find the error/distance between two images by comparing their correlation
%to the autocorrelation functions

%Find difference between half-maximum shifts of img1-img1, img2-img2, and
%img1-img2
range_vals=0:.1:1000; %in nm

%img1-img1 - find target correlation value in each shift direction
pointsx_11 = autocorr_img1(:,1,1);
targetx_11 = max(pointsx_11)*.5;
pointsy_11 = autocorr_img1(1,:,1);
targety_11 = max(pointsy_11)*.5;
temp = autocorr_img1(1,1,:);
pointsz_11 = temp(:);
targetz_11 = max(pointsz_11)*0.5;

%2-2 - find target correlation value in each shift direction
pointsx_22 = autocorr_img2(:,1,1);
targetx_22 = max(pointsx_22)*.5;
pointsy_22 = autocorr_img2(1,:,1);
targety_22 = max(pointsy_22)*.5;
temp = autocorr_img2(1,1,:);
pointsz_22 = temp(:);
targetz_22 = max(pointsz_22)*.5;

%1-2 - find target correlation value 
pointsx_12 = correlations(:,1,1);
targetx_12 = max(pointsx_12)*.5;
pointsy_12 = correlations(1,:,1);
targety_12 = max(pointsy_12)*.5;
temp = correlations(1,1,:);
pointsz_12 = temp(:);
targetz_12 = max(pointsz_12)*.5;

%1-1
fitobjectx = fit(xshifts',pointsx_11,'poly3');
fitx = fitobjectx(range_vals);
[~,ind] = min(abs(fitx-targetx_11));
inds_11(1,1) = range_vals(ind);

fitobjecty = fit(yshifts',pointsy_11','poly3');
fity = fitobjecty(range_vals);
[~,ind] = min(abs(fity-targety_11));
inds_11(1,2) = range_vals(ind);%range(inds(1));

fitobjectz = fit(zshifts',pointsz_11,'poly3');
fitz = fitobjectz(range_vals);
[~,ind] = min(abs(fitz-targetz_11));
inds_11(1,3) = range_vals(ind);%range(inds(1));

%2-2
fitobjectx = fit(xshifts',pointsx_22,'poly3');
fitx = fitobjectx(range_vals);
[~,ind] = min(abs(fitx-targetx_22));
inds_22(1,1) = range_vals(ind);

fitobjecty = fit(yshifts',pointsy_22','poly3');
fity = fitobjecty(range_vals);
[~,ind] = min(abs(fity-targety_22));
inds_22(1,2) = range_vals(ind);%range(inds(1));

fitobjectz = fit(zshifts',pointsz_22,'poly3');
fitz = fitobjectz(range_vals);
[~,ind] = min(abs(fitz-targetz_22));
inds_22(1,3) = range_vals(ind);%range(inds(1));        

%1-2
fitobjectx = fit(xshifts',pointsx_12,'poly3');
fitx = fitobjectx(range_vals);
[~,ind] = min(abs(fitx-targetx_12));
inds_12(1,1) = range_vals(ind);

fitobjecty = fit(yshifts',pointsy_12','poly3');
fity = fitobjecty(range_vals);
[~,ind] = min(abs(fity-targety_12));
inds_12(1,2) = range_vals(ind);%range(inds(1));

fitobjectz = fit(zshifts',pointsz_12,'poly3');
fitz = fitobjectz(range_vals);
[~,ind] = min(abs(fitz-targetz_12));
inds_12(1,3) = range_vals(ind);%range(inds(1));   

diffs_12_11 = inds_12(1,:) - inds_11(1,:);
% diffs_12_22 = inds_12(1,:) - inds_22(1,:);
% 
% diff = mean([diffs_12_11(:); diffs_12_22(:)]);
diff = mean([diffs_12_11(:)]);
%diff = mean(diff);

end

