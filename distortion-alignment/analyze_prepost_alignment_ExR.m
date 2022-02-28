function data_struct = analyze_prepost_alignment_ExR(folder,params)
% UPDATE FUNCTION DESCRIPTION HERE
%

%load parameters
parentfolder = params.parentfolder;
savefolder = params.savefolder;
xystep = params.xystep;
zstep = params.zstep;
conn = params.conn;
layerfolders = params.layerfolders;
lowerlim_um = params.lowerlim_um;

for lidx = 1:length(layerfolders) %loop through layers
    %Load the subfolder
    subfolder = [folder filesep layerfolders{lidx} filesep 'masks'];
    disp(subfolder)

    %get filenames for all images in the folder
    prefnames_bin = dir([parentfolder subfolder '/Bin-pre*.tif']); %binary from test stain
    T = struct2table(prefnames_bin);
    sortedT = sortrows(T,'bytes');
    prefnames_bin = table2struct(sortedT);
    
    postfnames_bin = dir([parentfolder subfolder '/Bin-post*.tif']); %binary from test stain
    T = struct2table(postfnames_bin);
    sortedT = sortrows(T,'bytes');
    postfnames_bin = table2struct(sortedT);
    
    prefnames_raw = dir([parentfolder subfolder '/Raw-pre*_noBG.tif']);
    T = struct2table(prefnames_raw);
    sortedT = sortrows(T,'bytes');
    prefnames_raw = table2struct(sortedT);
    
    postfnames_raw = dir([parentfolder subfolder '/Raw-post*_noBG.tif']);
    T = struct2table(postfnames_raw);
    sortedT = sortrows(T,'bytes');
    postfnames_raw = table2struct(sortedT);
    
    %loop through each cropped image in the folder
    for fidx = 1:length(prefnames_bin) 
        
        prefname_bin = prefnames_bin(fidx).name;
        postfname_bin = postfnames_bin(fidx).name;
        
        splitsname = split(prefname_bin,"-");
        layer = splitsname{3};
        samplesplit = split(splitsname{4},".");
        sample = samplesplit{1};
        
        %load all images
        clear info I_pre_bin
        info=imfinfo([parentfolder subfolder filesep prefname_bin]);
        z_length = numel(info);

        for k = 1:z_length
            I_pre_bin(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        dim = size(I_pre_bin);

        clear info I_post_bin
        info=imfinfo([parentfolder subfolder filesep postfname_bin]);
        for k = 1:z_length
            I_post_bin(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end

        clear info I_pre_raw
        prefname_raw = ['Raw-preCh-' layer '-' sample '_noBG.tif'];
        info=imfinfo([parentfolder subfolder filesep prefname_raw]);
        for k = 1:dim(3)
            I_pre_raw(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        
        clear info I_post_raw
        prefname_raw = ['Raw-postCh-' layer '-' sample '_noBG.tif'];
        info=imfinfo([parentfolder subfolder filesep prefname_raw]);
        for k = 1:dim(3)
            I_post_raw(:,:,k) = imread(info(1).Filename, k, 'Info', info);
        end
        
        
        vol_pre = nnz(I_pre_bin);
        vol_post = nnz(I_pre_bin);
        vol_pres(fidx,lidx) = vol_pre;
        vol_posts(fidx,lidx) = vol_post;
        total_vol = nnz(I_pre_bin) + nnz(I_post_bin);
        
        I_pre_filt = medfilt3(I_pre_bin,[5 5 3]);
        I_post_filt = medfilt3(I_post_bin,[5 5 3]);
        
        I_pre_filt = bwareaopen(I_pre_filt, 30);
        I_post_filt = bwareaopen(I_post_filt, 30);
        
        CC_pre = bwconncomp(I_pre_filt,conn);
        CC_post = bwconncomp(I_post_filt,conn);
        
        npre = CC_pre.NumObjects;
        n_pre(fidx,lidx) = npre;
        npost = CC_post.NumObjects;
        n_post(fidx,lidx) = npost;
        %disp(CC_post.NumObjects - CC_pre.NumObjects);
        
        if npre > 0 && npost > 0
            
            pre_bin_one = zeros(dim); %just keep largest puncta for some analysis - the largest
            if npre >1 
                sizes=zeros(npre,1);
                for oo = 1:npre
                    sizes(oo,1)= size(CC_pre.PixelIdxList{oo},1);
                end
                [~,ind_pre]=max(sizes);
                pre_bin_one(CC_pre.PixelIdxList{ind_pre})=1;
            else
                ind_pre =1;
                pre_bin_one = I_pre_bin;
            end
            
            post_bin_one = zeros(dim);%just keep largest puncta for some analysis - the largest
            if npost >1
                sizes=zeros(npost,1);
                for oo = 1:npost
                    sizes(oo,1)= size(CC_post.PixelIdxList{oo},1);
                end
                [~,ind_post]=max(sizes);
                post_bin_one(CC_post.PixelIdxList{ind_post})=1;
            else
                ind_post = 1;
                post_bin_one = I_post_bin;
            end
            
            centroids_pre = regionprops3(CC_pre,'Centroid');
            centroid_pre = centroids_pre.Centroid(ind_pre,:);
            centroids_post = regionprops3(CC_post,'Centroid');
            centroid_post = centroids_post.Centroid(ind_post,:);

            %linearization method one - maximum, mean, and minimum geodesic distance between boundaries
            if fidx == 10 && lidx == 2
                plotting = 1;
            else
                plotting = 0;
            end
            [mean_dist,max_dist,min_dist] = geodesic_btwn_edges(pre_bin_one,post_bin_one,plotting);
            mean_dist_geo(fidx,lidx) = mean_dist*xystep;
            max_dist_geo(fidx,lidx) = max_dist*xystep;
            min_dist_geo(fidx,lidx) = min_dist*xystep;
            %disp(mean_dist*xystep);

            %linearization method two - fit curve and take shortest
            %distance between lines tangent to the curves at different
            %points
            [mean_dist,max_dist,min_dist] = distance_between_curves(pre_bin_one,post_bin_one,plotting);
            mean_dist_curve(fidx,lidx) = mean_dist*xystep;
            max_dist_curve(fidx,lidx) = max_dist*xystep;
            min_dist_curve(fidx,lidx) = min_dist*xystep;
            %disp(mean_dist*xystep);
            
            %linearization method three - radial distances between puncta
            %centroids. convert to um using a spherical assumption
            centroid_dist = ((centroid_post(1)-centroid_pre(1))^2 + ...
                (centroid_post(2) - centroid_pre(2))^2 + ...
                (centroid_post(3) - centroid_pre(3))^2)^(1/2);
            centroid_dist_voxels(fidx,lidx) = centroid_dist;
            centroid_dist_um(fidx,lidx) = centroid_dist*mean([xystep xystep zstep]);
            %disp(centroid_dist*mean([xystep xystep zstep]));

            %allocate space for shift data
            vol_overlaps=zeros(length(xshifts),length(yshifts),length(zshifts));
            correlations=zeros(size(vol_overlaps));
            frac_vol_overlaps=zeros(size(vol_overlaps));
            correlations_norm_int=zeros(size(vol_overlaps));
            correlations_norm_vol=zeros(size(vol_overlaps));
            correlations_norm_both=zeros(size(vol_overlaps));
            correlations_tyler=zeros(size(vol_overlaps));
            intnorm_autocorr_pre=zeros(size(vol_overlaps));
            intnorm_autocorr_post=zeros(size(vol_overlaps));
            volnorm_autocorr_pre=zeros(size(vol_overlaps));
            volnorm_autocorr_post=zeros(size(vol_overlaps));
            area_overlap=zeros(size(vol_overlaps));
            area_overlap_norm=zeros(size(vol_overlaps));

            %normalize to min-max intensity
            I_pre_raw=double(I_pre_raw);
            I_post_raw=double(I_post_raw);
            pre_norm = (I_pre_raw - min(I_pre_raw(:))) ./ (max(I_pre_raw(:)) - min(I_pre_raw(:)));
            post_norm = (I_post_raw - min(I_post_raw(:))) ./ (max(I_post_raw(:)) - min(I_post_raw(:)));

            masked_post = double(I_post_raw) .* double(I_post_bin);
            mean_int_inmask = mean(mean(mean(masked_post(masked_post>0))));
            masked_post_mean = masked_post;
            masked_post_mean(masked_post_mean>0) = mean_int_inmask;

            for xx = 1:length(xshifts)
                xshift = xshifts(xx);
                for yy = 1:length(yshifts)
                    yshift = yshifts(yy);
                    for zz = 1:length(zshifts)
                        zshift = zshifts(zz);

                        %shift the pre-expansion channel
                        pre_shifted = imtranslate(I_pre_raw,[xshift yshift zshift]);
                        pre_shifted_bin = imtranslate(I_pre_bin,[xshift yshift zshift]);

                        if xshift == 0 && yshift == 0 && zshift > 0
                            edge = imtranslate(I_pre_raw,[0 0 zshift-dim(3)]);
                            edge_bin = imtranslate(I_pre_bin,[0 0 zshift-dim(3)]);
                        elseif xshift == 0 && yshift > 0 && zshift > 0
                            edge1 = imtranslate(I_pre_raw,[0 yshift-dim(1) 0]);
                            edge1_bin = imtranslate(I_pre_bin,[0 yshift-dim(1) 0]);

                            edge2 = imtranslate(I_pre_raw,[0 0 zshift-dim(3)]);
                            edge2_bin = imtranslate(I_pre_bin,[0 0 zshift-dim(3)]);

                            edge = edge1 + edge2;
                            edge_bin = edge1_bin + edge2_bin;
                        elseif zshift == 0 && yshift > 0 && zshift == 0
                            edge = imtranslate(I_pre_raw,[0 yshift-dim(1) 0]);
                            edge_bin = imtranslate(I_pre_bin,[0 yshift-dim(1) 0]);
                        elseif xshift > 0 && yshift == 0 && zshift == 0
                            edge = imtranslate(I_pre_raw,[xshift-dim(2) 0 0]);
                            edge_bin = imtranslate(I_pre_bin,[xshift-dim(2) 0 0]);
                        elseif xshift > 0 && yshift > 0 && zshift == 0
                            edge1 = imtranslate(I_pre_raw,[xshift-dim(2) 0 0]);
                            edge1_bin = imtranslate(I_pre_bin,[xshift-dim(2) 0 0]);

                            edge2 = imtranslate(I_pre_raw,[0 yshift-dim(1)  0]);
                            edge2_bin = imtranslate(I_pre_bin,[0 yshift-dim(1) 0]);

                            edge = edge1 + edge2;
                            edge_bin = edge1_bin + edge2_bin;
                        elseif xshift > 0 && yshift > 0 && zshift > 0
                            edge1 = imtranslate(I_pre_raw,[0 0 zshift-dim(3)]);
                            edge1_bin = imtranslate(I_pre_bin,[0 0 zshift-dim(3)]);

                            edge2 = imtranslate(I_pre_raw,[0 yshift-dim(1) 0]);
                            edge2_bin = imtranslate(I_pre_bin,[0 yshift-dim(1) 0]);

                            edge3 = imtranslate(I_pre_raw,[xshift-dim(2) 0 0]);
                            edge3_bin = imtranslate(I_pre_bin,[xshift-dim(2) 0 0]);

                            edge = edge1 + edge2 + edge3;
                            edge_bin = edge1_bin + edge2_bin + edge3_bin;
                        elseif xshift == 0 && yshift == 0 && zshift == 0
                            edge = zeros(dim);
                            edge_bin = uint8(zeros(dim));
                        end

                        pre_shifted = pre_shifted + edge;
                        pre_shifted_bin = pre_shifted_bin + edge_bin;

                        %shift the post-expansion channel
                        post_shifted = imtranslate(I_post_raw,[xshift yshift zshift]);
                        %post_shifted_bin = imtranslate(I_post_bin,[xshift yshift zshift]);

                        if xshift == 0 && yshift == 0 && zshift > 0
                            edge = imtranslate(I_post_raw,[0 0 zshift-dim(3)]);
                            %edge_bin = imtranslate(I_post_bin,[0 0 zshift-dim(3)]);
                        elseif xshift == 0 && yshift > 0 && zshift > 0
                            edge1 = imtranslate(I_post_raw,[0 yshift-dim(1) 0]);
                            %edge1_bin = imtranslate(I_post_bin,[0 yshift-dim(1) 0]);

                            edge2 = imtranslate(I_post_raw,[0 0 zshift-dim(3)]);
                            %edge2_bin = imtranslate(I_post_bin,[0 0 zshift-dim(3)]);

                            edge = edge1 + edge2;
                            %edge_bin = edge1_bin + edge2_bin;
                        elseif zshift == 0 && yshift > 0 && zshift == 0
                            edge = imtranslate(I_post_raw,[0 yshift-dim(1) 0]);
                            %edge_bin = imtranslate(I_post_bin,[0 yshift-dim(1) 0]);
                        elseif xshift > 0 && yshift == 0 && zshift == 0
                            edge = imtranslate(I_post_raw,[xshift-dim(2) 0 0]);
                            %edge_bin = imtranslate(I_post_bin,[xshift-dim(2) 0 0]);
                        elseif xshift > 0 && yshift > 0 && zshift == 0
                            edge1 = imtranslate(I_post_raw,[xshift-dim(2) 0 0]);
                            %edge1_bin = imtranslate(I_post_bin,[xshift-dim(2) 0 0]);

                            edge2 = imtranslate(I_post_raw,[0 yshift-dim(1)  0]);
                            %edge2_bin = imtranslate(I_post_bin,[0 yshift-dim(1) 0]);

                            edge = edge1 + edge2;
                            %edge_bin = edge1_bin + edge2_bin;
                        elseif xshift > 0 && yshift > 0 && zshift > 0
                            edge1 = imtranslate(I_post_raw,[0 0 zshift-dim(3)]);
                            %edge1_bin = imtranslate(I_post_bin,[0 0 zshift-dim(3)]);

                            edge2 = imtranslate(I_post_raw,[0 yshift-dim(1) 0]);
                            %edge2_bin = imtranslate(I_post_bin,[0 yshift-dim(1) 0]);

                            edge3 = imtranslate(I_post_raw,[xshift-dim(2) 0 0]);
                            %edge3_bin = imtranslate(I_post_bin,[xshift-dim(2) 0 0]);

                            edge = edge1 + edge2 + edge3;
                            %edge_bin = edge1_bin + edge2_bin + edge3_bin;
                        elseif xshift == 0 && yshift == 0 && zshift == 0
                            edge = zeros(dim);
                            %edge_bin = uint8(zeros(dim));
                        end

                        post_shifted = post_shifted + edge;
                        %post_shifted_bin = post_shifted_bin + edge_bin;

                        %volume overlap (and normalized to total volume) as a function of
                        %shift distance in x,y, and z directions
                        intersection = pre_shifted_bin & I_post_bin;
                        vol_overlaps(xx,yy,zz) = nnz(intersection);
                        frac_vol_overlap = nnz(intersection)/vol_pre;
                        frac_vol_overlaps(xx,yy,zz) = frac_vol_overlap;

                        %correlation between pre and post images as a function of shift
                        %distance in x,y, and z directions
                        raw_corr = corr(double(pre_shifted(:)),double(I_post_raw(:)));
                        correlations(xx,yy,zz) = raw_corr;

                        %normalize shifted to min/max intensity
                        pre_shifted_norm = (pre_shifted - min(pre_shifted(:))) / (max(pre_shifted(:)) - min(pre_shifted(:)));
                        post_shifted_norm = (post_shifted - min(post_shifted(:))) / (max(post_shifted(:)) - min(post_shifted(:)));

                        %normalize correlations to intensity and volume
                        norm_corr = corr(double(pre_shifted_norm(:)),double(post_norm(:)));
                        correlations_norm_int(xx,yy,zz) = norm_corr;
                        correlations_norm_vol(xx,yy,zz) = raw_corr ./ vol_pre;
                        correlations_norm_both(xx,yy,zz) = norm_corr ./ vol_pre;

                        masked_pre = double(pre_shifted) .* double(pre_shifted_bin);
                        mean_int_inmask = mean(mean(mean(masked_pre(masked_pre>0))));
                        masked_pre_mean = masked_pre;
                        masked_pre_mean(masked_pre_mean>0) = mean_int_inmask;
                        
                        masked_post = double(I_post_raw) .* double(I_post_bin);
                        mean_int_inmask = mean(mean(mean(masked_post(masked_post>0))));
                        masked_post_mean = masked_post;
                        masked_post_mean(masked_post_mean>0) = mean_int_inmask;
                        
                        norm_tyler1 = double(masked_post./masked_post_mean);
                        norm_tyler2 = double(masked_pre./masked_pre_mean);
                        norm_tyler1(isnan(norm_tyler1))=0;
                        norm_tyler2(isnan(norm_tyler2))=0;

                        correlation_tyler(xx,yy,zz) = corr(norm_tyler1(:),norm_tyler2(:));

                        %normalized autocorrelations (want to compare these)
                        autocorr_pre = corr(double(pre_shifted_norm(:)),double(pre_norm(:)));
                        autocorr_post = corr(double(post_shifted_norm(:)),double(post_norm(:)));

                        intnorm_autocorr_pre(xx,yy,zz) = autocorr_pre;
                        intnorm_autocorr_post(xx,yy,zz) = autocorr_post;

                        volnorm_autocorr_pre(xx,yy,zz) = autocorr_pre./vol_pre;
                        volnorm_autocorr_post(xx,yy,zz) = autocorr_post./vol_post;

                    end
                end
            end
        
            vol_overlaps_all{fidx,lidx} = vol_overlaps;
            correlations_all{fidx,lidx} = correlations;

            frac_vol_overlaps_all{fidx,lidx}=frac_vol_overlaps;

            correlations_norm_int_all{fidx,lidx}=correlations_norm_int;
            correlations_norm_vol_all{fidx,lidx}=correlations_norm_vol;
            correlations_norm_both_all{fidx,lidx}=correlations_norm_both;
            correlations_tyler_all{fidx,lidx}=correlation_tyler;

            intnorm_autocorr_pre_all{fidx,lidx}=intnorm_autocorr_pre;
            intnorm_autocorr_post_all{fidx,lidx}=intnorm_autocorr_post;
            volnorm_autocorr_pre_all{fidx,lidx}=volnorm_autocorr_pre;
            volnorm_autocorr_post_all{fidx,lidx}=volnorm_autocorr_post;

            area_overlap_all{fidx,lidx}=area_overlap;
            frac_area_overlap_all{fidx,lidx}=area_overlap_norm;
            filename{fidx,lidx} = prefname_bin;
            
        else
            vol_overlaps_all{fidx,lidx} = NaN;
            correlations_all{fidx,lidx} = NaN;

            frac_vol_overlaps_all{fidx,lidx}=NaN;

            correlations_norm_int_all{fidx,lidx}=NaN;
            correlations_norm_vol_all{fidx,lidx}=NaN;
            correlations_norm_both_all{fidx,lidx}=NaN;
            correlations_tyler_all{fidx,lidx}=NaN;

            intnorm_autocorr_pre_all{fidx,lidx}=NaN;
            intnorm_autocorr_post_all{fidx,lidx}=NaN;
            volnorm_autocorr_pre_all{fidx,lidx}=NaN;
            volnorm_autocorr_post_all{fidx,lidx}=NaN;

            area_overlap_all{fidx,lidx}=NaN;
            frac_area_overlap_all{fidx,lidx}=NaN;
            filename{fidx,lidx} = prefname_bin;
        end
    end
end

splits = split(subfolder,"_");
data_struct.protein = splits{2};
data_struct.vol_pre = vol_pres;
data_struct.vol_post = vol_posts;
data_struct.n_pre = n_pre;
data_struct.n_post = n_post;

data_struct.centroid_dist_voxels = centroid_dist_voxels;
data_struct.centroid_dist_um = centroid_dist_um;

data_struct.vol_overlap = vol_overlaps_all;
data_struct.frac_vol_overlap = frac_vol_overlaps_all;

data_struct.mean_dist_geo = mean_dist_geo;
data_struct.min_dist_geo = min_dist_geo;
data_struct.max_dist_geo = max_dist_geo;

data_struct.mean_dist_curve = mean_dist_curve;
data_struct.min_dist_curve = min_dist_curve;
data_struct.max_dist_curve = max_dist_curve;

data_struct.correlations = correlations_all;
data_struct.correlations_norm_int = correlations_norm_int_all;
data_struct.correlations_norm_vol = correlations_norm_vol_all;
data_struct.correlations_norm_both = correlations_norm_both_all;
data_struct.correlation_tyler = correlations_tyler_all;
data_struct.intnorm_autocorr_pre=intnorm_autocorr_pre_all;
data_struct.intnorm_autocorr_post=intnorm_autocorr_post_all;
data_struct.volnorm_autocorr_pre=volnorm_autocorr_pre_all;
data_struct.volnorm_autocorr_post=volnorm_autocorr_post_all;

data_struct.fname = filename;

        
