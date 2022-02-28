function concat_data = prepost_distortion_concatenate_data_ExR(data,params)
% If something breaks in here, it's probably because you don't have the
% same number of biological replicates and proteins as in the original ExR
% manuscript. You will have to manually edit this code.
%

% Manually edit this based on the number of biological replicates you have
% Apologies for this being clunky
%find indices of each protein of interest in data
nproteins = length(params.proteins);

if nproteins == 2
    p1_inds = 1:params.nreplicates; %indices of protein 1
    p2_inds = params.nreplicates+1:(nproteins*params.nreplicates);
    pinds = [p1_inds p2_inds];
else
    disp('More than 2 proteins - please manually adjust protein indices')
end

for pp = 1:length(params.proteins)
    r1 = pinds(pp);
    r2 = pinds(pp+1);
    r3 = pinds(pp+2);
    
    concat_data(pp).correlations_blanpied = [data(r1).data.correlation_blanpied(:);
        data(r2).data.correlation_blanpied(:);
        data(r3).data.correlation_blanpied(:)];
    
    concat_data(pp).correlations_norm_minmax = [data(r1).data.correlations_norm_minmax(:);
        data(r2).data.correlations_norm_minmax(:);
        data(r3).data.correlations_norm_minmax(:)];
    
    concat_data(pp).volo = [data(r1).data.frac_vol_overlap(:);
        data(r2).data.frac_vol_overlap(:);
        data(r3).data.frac_vol_overlap(:)];
    
    concat_data(pp).protein = params.proteins{pp};
    concat_data(pp).protein_plotnames = params.proteins_plotnames{pp};
    
    concat_data(pp).autocorr_pre = [data(r1).data.autocorr_pre(:);
        data(r2).data.autocorr_pre(:);
        data(r3).data.autocorr_pre(:)]; 
    
    concat_data(pp).autocorr_post = [data(r1).data.autocorr_post(:);
        data(r2).data.autocorr_post(:);
        data(r3).data.autocorr_post(:)]; 
    
    concat_data(pp).npuncta_pre = [data(r1).data.n_pre(:);
        data(r2).data.n_pre(:);
        data(r3).data.n_pre(:)]; 
    
    concat_data(pp).npuncta_post = [data(r1).data.n_post(:);
        data(r2).data.n_post(:);
        data(r3).data.n_post(:)]; 
    
    concat_data(pp).delta_npuncta = [data(r1).data.delta_npuncta(:);
        data(r2).data.delta_npuncta(:);
        data(r3).data.delta_npuncta(:)]; 
    
    concat_data(pp).delta_npuncta_norm = [data(r1).data.delta_npuncta_norm(:);
        data(r2).data.delta_npuncta_norm(:);
        data(r3).data.delta_npuncta_norm(:)]; 
    
    %go through and replace nans with empty cells
    for cc = 1:length(concat_data(pp).correlations_blanpied)
        if isnan(concat_data(pp).volo{cc})
            concat_data(pp).volo{cc} = [];
        end
        if isnan(concat_data(pp).autocorr_pre{cc})
            concat_data(pp).autocorr_pre{cc} = [];
        end
        if isnan(concat_data(pp).autocorr_post{cc})
            concat_data(pp).autocorr_post{cc} = [];
        end
        if isnan(concat_data(pp).correlations_blanpied{cc})
            concat_data(pp).correlations_blanpied{cc} = [];
        end
        if isnan(concat_data(pp).correlations_norm_minmax{cc})
            concat_data(pp).correlations_norm_int{cc} = [];
        end
    end
end

end

