function plot_correlation_1D(correlations,autocorr_img1,autocorr_img2,...
    img1name,img2name,xshifts,yshifts,roiname,zshifts,targetfolder)
%Plot correlation vs. autocorrelation in 1D

yellow = [0.9290, 0.6940, 0.1250];
roi1name = [roiname '-' img1name];
roi2name = [roiname '-' img2name];

% x direction
if zshifts
    h = figure();
    subplot(3,1,1)
    plot(xshifts,autocorr_img1(:,1,1),'Color',yellow,'LineWidth',3)
    hold on
    plot(xshifts,autocorr_img2(:,1,1),'m','LineWidth',3)
    plot(xshifts,correlations(:,1,1),'k','LineWidth',3)
    hold off
    xlabel('Shift in x (nm)')
    ylabel('Pixel-wise Correlation')
    legend([roi1name '-' roi1name],[roi2name '-' roi2name],[roi1name '-' roi2name])
    %axis([0 300 -.1 1])

    subplot(3,1,2)
    plot(yshifts,autocorr_img1(1,:,1),'Color',yellow,'LineWidth',3)
    hold on
    plot(yshifts,autocorr_img2(1,:,1),'m','LineWidth',3)
    plot(yshifts,correlations(1,:,1),'k','LineWidth',3)
    hold off
    xlabel('Shift in y (nm)')
    ylabel('Pixel-wise Correlation')
    legend([roi1name '-' roi1name],[roi2name '-' roi2name],[roi1name '-' roi2name])
    %axis([0 300 -.1 1])

    subplot(3,1,3)
    temp = autocorr_img1(1,1,:);
    plot(zshifts,temp(:),'Color',yellow,'LineWidth',3)
    hold on
    temp = autocorr_img2(1,1,:);
    plot(zshifts,temp(:),'m','LineWidth',3)
    temp = correlations(1,1,:);
    plot(zshifts,temp(:),'k','LineWidth',3)
    hold off
    xlabel('Shift in z (nm)')
    ylabel('Pixel-wise Correlation')
    legend([roi1name '-' roi1name],[roi2name '-' roi2name],[roi1name '-' roi2name])
    %axis([0 300 -.1 1])

    %set(gcf,'Position',[10 10 800 1000])
    suptitle([img1name ' / ' img2name ' : Pixel-wise Correlations/Autocorrelations between Min-Max Normalized Images'])
    savefig(h,[targetfolder filesep roiname '_' img1name '_' img2name 'correlations_autocorrs_1D.fig'])
    close(h)
else
    h = figure();
    subplot(2,1,1)
    plot(xshifts,autocorr_img1(:,1),'Color',yellow,'LineWidth',3)
    hold on
    plot(xshifts,autocorr_img2(:,1),'m','LineWidth',3)
    plot(xshifts,correlations(:,1),'k','LineWidth',3)
    hold off
    xlabel('Shift in x (nm)')
    ylabel('Pixel-wise Correlation')
    legend([roi1name '-' roi1name],[roi2name '-' roi2name],[roi1name '-' roi2name])
    %axis([0 300 -.1 1])

    subplot(2,1,2)
    plot(yshifts,autocorr_img1(1,:),'Color',yellow,'LineWidth',3)
    hold on
    plot(yshifts,autocorr_img2(1,:),'m','LineWidth',3)
    plot(yshifts,correlations(1,:),'k','LineWidth',3)
    hold off
    xlabel('Shift in y (nm)')
    ylabel('Pixel-wise Correlation')
    legend([roi1name '-' roi1name],[roi2name '-' roi2name],[roi1name '-' roi2name])
    %axis([0 300 -.1 1])

    %set(gcf,'Position',[10 10 800 1000])
    suptitle([img1name ' / ' img2name ' : Pixel-wise Correlations/Autocorrelations between Min-Max Normalized Images'])
    savefig(h,[targetfolder filesep roiname '_' img1name '_' img2name 'correlations_autocorrs_1D.fig'])
    close(h)
end

    

end

