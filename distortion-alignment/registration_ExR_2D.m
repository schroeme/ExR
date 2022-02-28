function registration_ExR_2D(parentpath, targetpath, fovind, roundpaths, params)

    %extract params
    xystep = params.xystep;
    
    maxshiftumsmall = params.maxshiftumsmall;
    maxshiftsmall = round(maxshiftumsmall*(1/xystep));
    pFlag=params.pFlag; %plotting flag - change to 1 to see the plots
    nrounds = params.nrounds; 
    regtype = params.regtype;
    normint = params.normint;    
    gaussfiltsmall = params.gaussfiltsmall;
    nOctaves = params.nOctaves;
    options.overwrite = true;
    options.message = false;
    warning('off','all');

    %% Running section - for each round

    fovdirs_R1 = dir([parentpath roundpaths{1} '*.tif*']);
    r1path = [fovdirs_R1(fovind).folder filesep fovdirs_R1(fovind).name]; %extract round 1 path for this fov
    disp(r1path);

    for roundind = 2:nrounds %loop through each round after the first
        ref1 = imread(r1path);
        ref1 = double(medfilt2(ref1,[3 3]));
        
        % Read in the images and rescale them
        fovdirs_Rn = dir([parentpath roundpaths{roundind} '*.tif*']);
        rnpath = [fovdirs_Rn(fovind).folder filesep fovdirs_Rn(fovind).name];
        disp(rnpath)
        
        for roundind = 2:nrounds %loop through each round after the first
            targetdir = [targetpath roundpaths{roundind} 'registered'];
            mkdir(targetdir)
            savepath = [targetdir filesep fovdirs_Rn(fovind).name(1:end-4) '_reg.tif'];
            refn_unfilt = imread(rnpath);
            refn = double(medfilt2(refn_unfilt,[3 3]));

            %register reference channel to the first round
            if strcmp(regtype,'translation') %translation registration using intensity
                %NEED TO UPDATE THIS PATYHWAY
                [registered,tforms] = slicewise_reg_TR_intensity(ref1,refn,normint,gaussfiltsmall);
            elseif strcmp(regtype,'rigid') %rigid body registration using intensity
                [corrvals,tformsout,tformbest] = slicewise_corr_RB_intensity(ref1,refn,gaussfiltsmall,normint);
                [registered,~] = slicewise_reg_RB_intensity(ref1,refn_unfilt,gaussfiltsmall,tformbest);
            elseif strcmp(regtype,'feature-affine') %feature-based registration using SURF
                %NEED TO UPDATE THIS PATHWAY
                [corrvals,tformsout,tformbest_feat] = slicewise_corr_SURF_intensity_transverse(ref1,refn,maxshiftsmall,nOctaves,pFlag,0,normint,gaussfiltsmall);
                [registered,~] = slicewise_reg_SURF_intensity(ref1,refn_unfilt,gaussfiltsmall,tformbest_feat);
                %[corrvals,tformsout,tformbest] = slicewise_corr_affine_intensity(ref1,refn,0,0);
                [registered,tformsout_affine,tformbest_affine] = slicewise_reg_affine_intensity(ref1,registered,gaussfiltsmall);
            elseif strcmp(regtype,'similarity')
                [corrvals,tformsout,tformbest] = slicewise_corr_similarity_intensity(ref1,refn,gaussfiltsmall);
                [registered,~] = slicewise_reg_similarity_intensity(ref1,refn_unfilt,gaussfiltsmall,tformbest);
            end

            saveastiff(registered, savepath, options);

        end
    end

end