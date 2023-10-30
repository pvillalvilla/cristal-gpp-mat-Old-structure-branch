figure;
[~,outputPath_size] = size(filesBulk.outputPath);
plots_folder = [filesBulk.outputPath filesBulk.filename_L1B(outputPath_size+1:end-3) '\bursts\'];
if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
end
max_val=zeros(1,length(L1BS_buffer(i_surf_stacked).burst_index(1):L1BS_buffer(i_surf_stacked).burst_index(end)));
max_pos=zeros(1,length(L1BS_buffer(i_surf_stacked).burst_index(1):L1BS_buffer(i_surf_stacked).burst_index(end)));
for i_burst_plots = L1BS_buffer(i_surf_stacked).burst_index(1):L1BS_buffer(i_surf_stacked).burst_index(end)
    subplot(2,3,[1:2 4:5]);
	mesh(abs(fft((fftshift(squeeze(L1A_buffer(i_burst_plots).beams_focused_shifted(:,:)),2)).', chd.N_samples_sar ).').^2); hold all;
    beam_selected=(abs(fft((fftshift(squeeze(L1A_buffer(i_burst_plots).beams_focused_shifted(L1BS_buffer(i_surf_stacked).beam_index(find(L1BS_buffer(i_surf_stacked).burst_index==i_burst_plots)),:)),2)).', chd.N_samples_sar).').^2);
    plot3(1:256,zeros(1,256)+ L1BS_buffer(i_surf_stacked).beam_index(find(L1BS_buffer(i_surf_stacked).burst_index==i_burst_plots)),reshape(beam_selected,[size(beam_selected) 1]),'or');
    set(gca,'XLim',[1 256],'FontSize',12);
    set(gca,'YLim',[1 64],'FontSize',12);
    set(gca,'ZLim',[0 5e-10],'FontSize',12);
    figlabels('Samples','Beams','Power',['Burst ' num2str(i_burst_plots) ', Selected Beam: ' num2str(L1BS_buffer(i_surf_stacked).beam_index(find(L1BS_buffer(i_surf_stacked).burst_index==i_burst_plots)))],12)
    view(60,50)
    hold off;
    [max_val(i_burst_plots),max_pos(i_burst_plots)]=max(max(abs(fft((ifftshift(squeeze(L1A_buffer(i_burst_plots).beams_focused_shifted(:,:)),2)).', chd.N_samples_sar ).').^2));
    subplot(2,3,3);
    plot(90:i_burst_plots,max_val(90:i_burst_plots));hold all; plot(i_burst_plots,max_val(i_burst_plots),'or');hold off;
    set(gca,'XLim',[89 321],'FontSize',12);
    set(gca,'YLim',[0 5e-10],'FontSize',12);
    figlabels('Bursts','Power','',['Burst: ' num2str(i_burst_plots)  ],12);
    subplot(2,3,6);
    plot(90:i_burst_plots,max_pos(90:i_burst_plots));hold all; plot(i_burst_plots,max_pos(i_burst_plots),'or');hold off;
    set(gca,'XLim',[89 321],'FontSize',12);
    set(gca,'YLim',[1 256],'FontSize',12);
    figlabels('Bursts','Sample','',['Burst ' num2str(i_burst_plots)],12)
    figName = [plots_folder  num2str(i_burst_plots, '%03d') '.png'];
    print(gcf, figName,'-dpng');
%     saveas (gcf,[figName,'.png']);
%     saveas (h,['./results/plots/Burst_' num2str(i_burst_plots, '%03d') '.png']);

    
    
    
end
