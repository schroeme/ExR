function [corrvals,tformsout,tformbest] = slicewise_corr_RB_intensity(img1,img2,gaussfilt,normint)

nslices = size(img1,3);
corrvals = zeros(nslices,1);
tformsmat = zeros(3,3,nslices);
scales=zeros(nslices,1);

for ss = 1:nslices %using reference channel, get transform that produces largest correlation
    slice1 = img1(:,:,ss);
    slice2 = img2(:,:,ss);
    
    if normint
        slice1 = (slice1 - mean(slice1(:)))/(std(slice1(:)));
        slice2 = (slice2 - mean(slice2(:)))/(std(slice2(:)));
    end
    
    [corrvals(ss,1),tformsout{ss}] = corr_transverse_slices_int(slice1,slice2,'rigid',0,'multimodal');
    tformsmat(:,:,ss) = tformsout{ss}.T;
    scales(ss,1) = det(tformsmat(:,:,ss))^(1/3);
end
    
    corrvals=corrvals(0.99 <= scales <= 1.01);
    [sorted,ind] = sort(corrvals,'descend');
    tformbest = tformsout{1};
    if any(sorted>=.80) %if any correlation is above 0.8, take the mean of all of those
        above80 = corrvals>=0.80;
        A = mean(tformsmat(:,:,above80),3);
    else %otherwise, take the maximum correlation transformation
        A = tformsmat(:,:,ind(1));
    end
    
    tformbest.T = A;
    
end