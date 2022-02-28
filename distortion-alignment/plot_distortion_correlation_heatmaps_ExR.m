function plot_distortion_correlation_heatmaps_ExR(concat_data,params)
% Plot correlation as a function of shift distance

xshifts=params.xshifts*params.xystep*1000; %convert to nm
yshifts=params.yshifts*params.xystep*1000;
selects_index = randi(length(concat_data(1).autocorr_pre),[7,1]);

for pp = 1:length(concat_data)
    mean_vol = nanmean(cat(4, concat_data(pp).correlations_blanpied{:}),4);
    std_vol = std(cat(4, concat_data(pp).correlations_blanpied{:}),0,4);
    
    h = figure();
    subplot(2,4,1)
    imagesc(xshifts,yshifts,mean_vol(:,:,1));
    xlabel('Shift in x (nm)')
    ylabel('Shift in y (nm)')
    title('Mean');
    caxis([0 .8]) 
    
    subplot(2,4,2)
    imagesc(xshifts,yshifts,std_vol(:,:,1));
    xlabel('Shift in x (nm)')
    ylabel('Shift in y (nm)')
    title('Std Dev');
    caxis([0 .8]) 
    
    for ss = 1:length(selects_index)-1
        subplot(2,4,2+ss)
        imagesc(xshifts,yshifts,concat_data(pp).correlations_blanpied{selects_index(ss)}(:,:,1))
        xlabel('Shift in x (nm)')
        ylabel('Shift in y (nm)')
        title(['Synapse ' num2str(selects_index(ss))]); 
        caxis([0 .8])
        if ss == 2 || ss == 6
            colorbar()
        end
    end
    
    suptitle([concat_data(pp).protein_plotnames ' : Pixel-wise Correlation between Mean-Normalized, Masked Pre- and Post-expansion Channels'])
    %left bottom width height
    set(gcf,'Position',[100 100 800 400])
    savefig(h,[params.savefolder concat_data(pp).protein_plotnames 'correlations_blanpied.fig'])
    close(h)
end


