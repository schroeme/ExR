function plot_correlation_heatmaps_2D(correlations,img1name,img2name,xshifts,yshifts,targetfolder)
% Plot correlation as a function of shift distance
    
    h = figure();
    imagesc(xshifts,yshifts,correlations);
    xlabel('Shift in x (nm)')
    ylabel('Shift in y (nm)')
    title('Mean');
    %caxis([0 .8]) 
    
    suptitle([img1name ' / ' img2name ' : Pixel-wise Correlation between Min-Max Normalized Images'])
    %left bottom width height
    savefig(h,[targetfolder filesep img1name '_' img2name 'correlations.fig'])
    close(h)
end

