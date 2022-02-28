%FOR ABETA/KV7.2 COLOCALIZATION

clear all
folder = '/Users/margaretschroeder/Dropbox (MIT)/iterative direct ExM/20200716 AD Cluster Crops/ROIs/bin';

%Extract parameters
xystep = 0.01144;%71639; %um/voxel in x and y
zstep = 0.02667;%4; %um/voxel in z

fnames_all = dir([folder filesep '*.tiff']);
nchannels = 2;
nclusters = length(fnames_all)/nchannels;

npuncta = zeros(nclusters,nchannels);
frac_vol_occ = zeros(nclusters,nchannels);

max_puncta_vol = zeros(nclusters,nchannels);
mean_puncta_vol = zeros(nclusters,nchannels);
max_puncta_SA = zeros(nclusters,nchannels);
mean_puncta_SA = zeros(nclusters,nchannels);
max_eqdi = zeros(nclusters,nchannels);
max_paxlength = zeros(nclusters,nchannels*3);
frac_kv_in = zeros(nclusters,1);
kv_in = zeros(nclusters,1);
kv_out = zeros(nclusters,1);
cluster_size = zeros(nclusters,1);

count = 0;
for imidx = [1:2:length(fnames_all)] %loop through each cluster
    count = count+1;
    %load both channels for that cluster
    splits = strsplit(fnames_all(imidx).name,"_");
    fname_ch1 = [splits{1} '_' splits{2} '_' splits{3} '_bin.tiff'];
    fname_ch2 = [splits{1} '_' splits{2} '_' splits{3} '_bin.tiff'];
    fnames_tocopy{count,1} = fname_ch1;
    
    clear info1 I1 info2 I2
    info1=imfinfo([folder filesep fname_ch1]);
    info2=imfinfo([folder filesep fname_ch2]);
    h = info1.Width;
    w = info1.Height;
    z = length(info1);
    
    I1 = zeros(w,h,z);
    for k = 1:z
        I1(:,:,k) = imread(info1(1).Filename, k, 'Info', info1);
    end
    
    I2 = zeros(w,h,z);
    for k = 1:z
        I2(:,:,k) = imread(info2(1).Filename, k, 'Info', info2);
    end
  
    CC1 = bwconncomp(I1,26);%find all connected components in the image
    CC2 = bwconncomp(I2,26);%find all connected components in the image
    
    npuncta(count,1) = CC1.NumObjects;
    npuncta(count,2) = CC2.NumObjects;
    
   %if there is an object in either channel
    if npuncta(count,1) > 0 || npuncta(count,2) > 0
        
        cluster_size(count,1) = numel(I1);
        cluster_size(count,2) = numel(I2);

        puncta_vols1 = regionprops3(CC1,'Volume');
        puncta_vols1 = puncta_vols1.Volume;
        puncta_SAs1 = regionprops3(CC1,'SurfaceArea');
        puncta_SAs1 = puncta_SAs1.SurfaceArea;

        puncta_vols2 = regionprops3(CC2,'Volume');
        puncta_vols2 = puncta_vols2.Volume;
        puncta_SAs2 = regionprops3(CC2,'SurfaceArea');
        puncta_SAs2 = puncta_SAs2.SurfaceArea;

        eqdi1 = regionprops3(CC1,'EquivDiameter');
        eqdi1 = eqdi1.EquivDiameter;
        eqdi2 = regionprops3(CC2,'EquivDiameter');
        eqdi2 = eqdi2.EquivDiameter;

        paxlength1 = regionprops3(CC1,'PrincipalAxisLength');
        paxlength1 = paxlength1.PrincipalAxisLength;
        paxlength2 = regionprops3(CC2,'PrincipalAxisLength');
        paxlength2 = paxlength2.PrincipalAxisLength;

        %some interesting shape parameters
        max_puncta_vol(count,1) = max(puncta_vols1);
        mean_puncta_vol(count,1) = mean(puncta_vols1);
        max_puncta_vol(count,2) = max(puncta_vols2);
        mean_puncta_vol(count,2) = mean(puncta_vols2);
        max_puncta_SA(count,1) = max(puncta_SAs1);
        mean_puncta_SA(count,1) = mean(puncta_SAs1);
        max_puncta_SA(count,2) = max(puncta_SAs2);
        mean_puncta_SA(count,2) = mean(puncta_SAs2);
        max_eqdi(count,1) = max(eqdi1);
        max_eqdi(count,2) = max(eqdi2);
        if CC1.NumObjects > 1
            max_paxlength(count,1:3) = max(paxlength1);
        else
            max_paxlength(count,1:3) = paxlength1;
        end
        max_paxlength(count,4:6) = max(paxlength2);

        %volume overlap - what fraction of Kv is inside AB?
        intersection = I1 & I2;
        n_intersection = nnz(intersection);
        frac_kv_in(count,1) = n_intersection/nnz(I2);
        kv_in(count,1) = n_intersection;
        AB_in(count,1) = n_intersection;
        AB_frac(count,1) = n_intersection/nnz(I1);
        Kv_frac(count,1) = n_intersection/nnz(I2);
        kv_out(count,1) = nnz(I2) - n_intersection;
        AB_out(count,1) = nnz(I1) - n_intersection;
        
    end
           
end


%%
figure()
scatter3(max_puncta_vol(:,1),kv_out,kv_in,300,'k.')
xlabel('Size of largest AB puncta')
ylabel('Volume of Kv7.2 overlapped with AB')
zlabel('Volume of Kv7.2 outside of AB')

%% FOR ABETA/NAV COLOCALIZATION

clear all
folder = '/Users/margaretschroeder/Dropbox (MIT)/iterative direct ExM/20210506 AD Cluster Crops/ROIs/masks';

%Extract parameters
xystep = 0.01144;%71639; %um/voxel in x and y
zstep = 0.02667;%4; %um/voxel in z

fnames_all = dir([folder filesep '*.tiff']);
nchannels = 2;
nclusters = length(fnames_all)/nchannels;

npuncta = zeros(nclusters,nchannels);
frac_vol_occ = zeros(nclusters,nchannels);

max_puncta_vol = zeros(nclusters,nchannels);
mean_puncta_vol = zeros(nclusters,nchannels);
max_puncta_SA = zeros(nclusters,nchannels);
mean_puncta_SA = zeros(nclusters,nchannels);
max_eqdi = zeros(nclusters,nchannels);
max_paxlength = zeros(nclusters,nchannels*3);
frac_nav_in = zeros(nclusters,1);
nav_in = zeros(nclusters,1);
nav_out = zeros(nclusters,1);
cluster_size = zeros(nclusters,1);

count = 0;
for imidx = [1:2:length(fnames_all)] %loop through each cluster
    count = count+1;
    %load both channels for that cluster
    splits = strsplit(fnames_all(imidx).name,"_");
    fname_ch1 = [splits{1} '_' splits{2} '_ch1_bin.tiff'];
    fname_ch2 = [splits{1} '_' splits{2} '_ch2_bin.tiff'];
    fnames_tocopy{count,1} = fname_ch1;
    
    clear info1 I1 info2 I2
    info1=imfinfo([folder filesep fname_ch1]);
    info2=imfinfo([folder filesep fname_ch2]);
    h = info1.Width;
    w = info1.Height;
    z = length(info1);
    
    I1 = zeros(w,h,z);
    for k = 1:z
        I1(:,:,k) = imread(info1(1).Filename, k, 'Info', info1);
    end
    
    I2 = zeros(w,h,z);
    for k = 1:z
        I2(:,:,k) = imread(info2(1).Filename, k, 'Info', info2);
    end
  
    CC1 = bwconncomp(I1,26);%find all connected components in the image
    CC2 = bwconncomp(I2,26);%find all connected components in the image
    
    npuncta(count,1) = CC1.NumObjects;
    npuncta(count,2) = CC2.NumObjects;
    
   %if there is an object in either channel
    if npuncta(count,1) > 0 || npuncta(count,2) > 0
        
        cluster_size(count,1) = numel(I1);
        cluster_size(count,2) = numel(I2);

        puncta_vols1 = regionprops3(CC1,'Volume');
        puncta_vols1 = puncta_vols1.Volume;
        puncta_SAs1 = regionprops3(CC1,'SurfaceArea');
        puncta_SAs1 = puncta_SAs1.SurfaceArea;

        puncta_vols2 = regionprops3(CC2,'Volume');
        puncta_vols2 = puncta_vols2.Volume;
        puncta_SAs2 = regionprops3(CC2,'SurfaceArea');
        puncta_SAs2 = puncta_SAs2.SurfaceArea;

        eqdi1 = regionprops3(CC1,'EquivDiameter');
        eqdi1 = eqdi1.EquivDiameter;
        eqdi2 = regionprops3(CC2,'EquivDiameter');
        eqdi2 = eqdi2.EquivDiameter;

        paxlength1 = regionprops3(CC1,'PrincipalAxisLength');
        paxlength1 = paxlength1.PrincipalAxisLength;
        paxlength2 = regionprops3(CC2,'PrincipalAxisLength');
        paxlength2 = paxlength2.PrincipalAxisLength;

        %some interesting shape parameters
        if npuncta(count,1) > 0
            max_puncta_vol(count,1) = nanmax(puncta_vols1);
            mean_puncta_vol(count,1) = nanmean(puncta_vols1);
            max_puncta_SA(count,1) = nanmax(puncta_SAs1);
            mean_puncta_SA(count,1) = nanmean(puncta_SAs1);
            max_eqdi(count,1) = nanmax(eqdi1);
        end
        
        if npuncta(count,2) > 0
            max_puncta_vol(count,2) = nanmax(puncta_vols2);
            mean_puncta_vol(count,2) = nanmean(puncta_vols2);
            max_puncta_SA(count,2) = nanmax(puncta_SAs2);
            mean_puncta_SA(count,2) = nanmean(puncta_SAs2);
            max_eqdi(count,2) = nanmax(eqdi2);
        end
        
        if CC1.NumObjects > 1
            max_paxlength(count,1:3) = max(paxlength1);
        elseif CC1.NumObjects > 0
            max_paxlength(count,1:3) = paxlength1;
        else
            max_paxlength(count,1:3) = [nan nan nan];
        end
        
        if CC2.NumObjects > 1
            max_paxlength(count,4:6) = max(paxlength2);
        elseif CC2.NumObjects > 0
            max_paxlength(count,4:6) = paxlength2;
        else
            max_paxlength(count,4:6) = [nan nan nan];
        end

        %volume overlap - what fraction of Kv is inside AB?
        intersection = I1 & I2;
        n_intersection = nnz(intersection);
        frac_nav_in(count,1) = n_intersection/nnz(I1);
        nav_in(count,1) = n_intersection;
        AB_in(count,1) = n_intersection;
        AB_frac(count,1) = n_intersection/nnz(I2);
        nav_frac(count,1) = n_intersection/nnz(I1);
        nav_out(count,1) = nnz(I1) - n_intersection;
        AB_out(count,1) = nnz(I2) - n_intersection;
        
    end
           
end
