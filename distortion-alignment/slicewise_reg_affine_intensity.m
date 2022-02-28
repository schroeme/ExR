function [recovered,tformsout,tformbest] = slicewise_reg_affine_intensity(img1,img2,gaussfilt,tformin)

nslices = size(img1,3);
recovered = zeros(size(img1,1),size(img1,2),nslices);
tformsmat = zeros(3,3,nslices);
corrvals = zeros(nslices,1);
scales = zeros(nslices,1);

if nargin > 3 %we already have the transformations
    for ss = 1:nslices
        slice1 = img1(:,:,ss);
        slice2 = img2(:,:,ss);

        outputView = imref2d(size(slice1));
        recovered(:,:,ss) = imwarp(slice2,tformin,'OutputView',outputView,'interp','cubic');
    end
    tformsout={};
    recovered = uint16(recovered);
else
    %img2_filt = mat2gray(img2);
    %img2_filt = medfilt3(img2_filt,[5 5 3]);
    
    for ss = 1:nslices %registering reference channel, get transformation matrices
        slice1 = img1(:,:,ss);
        %slice2 = img2_filt(:,:,ss);
        slice2 = double(img2(:,:,ss));
        
        %try registering with initial guess first to save time
        [registered_temp,tform] = register_slices_INT(slice1,slice2,'affine',gaussfilt,'multimodal',0);
        
        tformsout{ss,1} = tform;
        scales(ss,1) = det(tform.T)^(1/3);
        corr_temp = corr(slice1(:),registered_temp(:));
        
        if 0.97 > scales(ss,1) || scales(ss,1) > 1.03 || corr_temp < 0.4 %if the scaling is too drastic or correlation is too low
            [registered_temp,tform] = register_slices_INT(slice1,registered_temp,'affine',gaussfilt,'multimodal',0);
            tformsout{ss,1} = tform;
            scales(ss,1) = det(tform.T)^(1/3);
        end
        
        corr_temp = corr(slice1(:),registered_temp(:));
        corrvals(ss,1) = corr_temp;
        tformsmat(:,:,ss) = tform.T;
        
    end
    
    [sorted,ind] = sort(corrvals,'descend');
    if all(sorted>=.90) %if correlation is super high for all, take mean
        A = mean(tformsmat,3);
%     elseif any(sorted>=.80) %if any correlation is above 0.8, take the mean of all of those
%         above80 = corrvals>=0.80;
%         A = mean(tformsmat(:,:,above80),3);
    else %otherwise, take the mean of the top n
        %A = tformsmat(:,:,ind(1:3));
        A = tformsmat(:,:,ind(1));
        %A = tformsmat(:,:,ind(1:5));
        %A = mean(A,3);
    end
    
    tformbest = tform;
    tformbest.T = A;
    
    for ss = 1:nslices %register the registration channel
        recovered(:,:,ss) = imwarp(img2(:,:,ss),tformbest,'OutputView',imref2d(size(slice1)),'interp','cubic');
    end
    recovered = uint16(recovered);
end



