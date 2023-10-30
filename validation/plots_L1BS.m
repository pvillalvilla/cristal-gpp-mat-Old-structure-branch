SAVE      = 1;
VISIBLE   = 1;
font_size = 20;
plotting_stacks = 1;
set_default_plot;

[~,outputPath_size] = size(filesBulk.outputPath);

plots_folder = [filesBulk.outputPath filesBulk.filename_L1B(outputPath_size+1:end-3) '/'];
if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
end

waveforms(i_surf_stacked,:) = L1B.wfm_cor_i2q2;

if(strcmp(cnf.processing_mode,'SIN'))
    phase_diff(i_surf_stacked,:) = L1B.phase_difference;
    coherence(i_surf_stacked,:) = L1B.coherence;

end

if(plotting_stacks)
    if(i_surf_stacked==1)
        h=figure;
    end
    
    AIR = sum(abs(L1BS_buffer(i_surf_stacked).beams_rng_cmpr).');
    subplot(2,2,2)
    plot(AIR,1:size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1));
%     h1 = axes;
    set(gca, 'Ydir', 'reverse');
%     set(gca,'XLim',[0  4e-4],'FontSize',font_size);
    set(gca,'YLim',[1 size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1)],'FontSize',font_size); 
    figlabels('Power','Beam','',['AIR Mean: ' num2str(L1B.stack_centre_index) '  sigma: ' num2str(L1B.stack_std_index) ' [Doppler beams]' ],font_size);
    power_fitted = max(AIR).* L1B.stack_agauss_index .* exp (-((1:size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1)) - L1B.stack_centre_index).^2 /(2*L1B.stack_std_index.^2));
    hold all; plot(power_fitted,1:size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1),'.-')
    legend('Power','Gauss fit');
    hold off;
    subplot(2,2,4)
    plot(sum(abs(L1BS_buffer(i_surf_stacked).beams_rng_cmpr)));
    set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size);
    figlabels('Range bin','Power','',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' RIR'],font_size);
%     set(gca,'YLim',[0  1e-3],'FontSize',font_size);
    hold off;
    
    [central, valid_beams] = TRP_beam_analysis(L1BS_buffer(i_surf_stacked));
    [~,alignment]=max(abs(L1BS_buffer(i_surf_stacked).beams_rng_cmpr).');
    [Slope_coef_A1,SS,MUMU] = polyfit(valid_beams,alignment(valid_beams),1);
    Slope_A1 = polyval(Slope_coef_A1, valid_beams,SS,MUMU); %Slope_coef_A1(1) related with datation error
    Stack_noise_A1=std(alignment(valid_beams)-Slope_A1)*chd.T0_nom*cst.c/2; % units [m]
    Stack_alignment_A1=(Slope_A1(end)-Slope_A1(1))/length(valid_beams)*chd.T0_nom*cst.c/2; % units [m/beam]

    subplot(2,2,3)
    plot(alignment,1:size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1));
%     h1 = axes;
    set(gca, 'Ydir', 'reverse');
    set(gca,'YLim',[1 size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1)],'FontSize',font_size);
    set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',font_size);
    figlabels('Max pos','Beam','',['alignment: ' num2str(Stack_alignment_A1,'%1f') ' [m/beam] and noise: ' num2str(Stack_noise_A1,'%1f') ' [m]'],font_size);
   
    hold all; 
    plot(Slope_A1,(valid_beams));
    legend('Max Pos','Slope');
    hold off;
    subplot(2,2,1)
    imagesc(1: chd.N_samples_sar, 1:size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1),abs(L1BS_buffer(i_surf_stacked).beams_rng_cmpr));
    set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);
    set(gca,'YLim',[1 size(L1BS_buffer(i_surf_stacked).beams_rng_cmpr,1)],'FontSize',font_size);
    colormap('jet'); 
    figlabels('Samples','Beam','Power',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' after alignment'],font_size);
    hold all;
    plot(Slope_A1,(valid_beams),'r');
    plot(L1B.stack_mask_vector,1:length(L1B.stack_mask_vector),'g');
    hold off;
    
    
    if SAVE == 1
        figName = [plots_folder 'Stack_after_alignment_' num2str(L1BS_buffer(i_surf_stacked).surf_counter,'%0.3d')];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
%         saveas (gcf,[figName,'.fig']);
    end 
    if(i_surf_stacked==N_total_surf_loc)
        close(h);
    end
    
    
end

if i_surf_stacked == N_total_surf_loc
    surface_locations = [[L1BS_buffer(:).lon_surf] ;[L1BS_buffer(:).lat_surf]]; 
    sat_locations = [[L1A_buffer(:).lon_sar_surf] ;[L1A_buffer(:).lat_sar_surf]]; 

    figure;plot(surface_locations(1,:),surface_locations(2,:),'bo'); hold all; plot(sat_locations(1,:),sat_locations(2,:),'.r');

    figure; mesh(10.*log10(waveforms.')); 
    figlabels('Record index L1B','Samples','Power [dB]','L1B Waveforms' ,font_size);
    set(gca,'YLim',[1  chd.N_samples_sar*cnf.zp_fact_range ],'FontSize',font_size);
    set(gca,'XLim',[1 N_total_surf_loc],'FontSize',font_size);
    colormap('jet'); colorbar; view(50,20);
    if SAVE == 1
        figName = [plots_folder '1_Scaled_waveform_L1B_internal'];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
    end
    
    
    if(strcmp(cnf.processing_mode,'SIN'))

        figure; mesh(phase_diff.'); 
        figlabels('Record index L1B','Samples','Phase diff [rad]','L1B Phase Diff' ,font_size);
        set(gca,'YLim',[1  chd.N_samples_sar*cnf.zp_fact_range ],'FontSize',font_size);
        set(gca,'XLim',[1 N_total_surf_loc],'FontSize',font_size);
        colormap('jet'); colorbar; view(50,20);
        if SAVE == 1
            figName = [plots_folder '2_Phase diff_L1B_internal'];
            print(gcf, figName,'-dpng');
            saveas (gcf,[figName,'.png']);
        end
        figure; mesh(coherence.'); 
        figlabels('Record index L1B','Samples','Coherence','L1B Coherence' ,font_size);
        set(gca,'YLim',[1  chd.N_samples_sar*cnf.zp_fact_range ],'FontSize',font_size);
        set(gca,'XLim',[1 N_total_surf_loc],'FontSize',font_size);
        colormap('jet'); colorbar; view(50,20);
        if SAVE == 1
            figName = [plots_folder '3_Coherence_L1B_internal'];
            print(gcf, figName,'-dpng');
            saveas (gcf,[figName,'.png']);
        end
    end
    
    
    
    
end