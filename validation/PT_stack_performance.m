%% PLOTTING
mission='PICE';
SAVE      = 1;
VISIBLE   = 1;
font_size = 20;
[~,outputPath_size] = size(filesBulk.outputPath);
plots_folder = [filesBulk.outputPath filesBulk.filename_L1B(outputPath_size+1:end-3) '/'];
tag_pt ='_1';
pos=strfind(filesBulk.CHD_file,'description');
tag_pt=filesBulk.CHD_file(pos+11:end-4);
if(isempty(tag_pt))
    tag_pt ='_1';
end
if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);

end


set(0,'defaultFigureVisible','on');
zp_fact_trp=1;


non_centered_spectra_1=fftshift(squeeze(L1BS_buffer(i_surf_stacked).beams_surf(:,:)),2);
beams_zp_fft_plot = fft([non_centered_spectra_1(:,1:chd.N_samples_sar/2),...
                       zeros(size(non_centered_spectra_1,1),(zp_fact_trp-1)*chd.N_samples_sar),...
                       non_centered_spectra_1(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                       chd.N_samples_sar * zp_fact_trp,2);
non_centered_spectra_2=fftshift(squeeze(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:)),2);
beams_zp_fft_aligned_plot = fft([non_centered_spectra_2(:,1:chd.N_samples_sar/2),...
                       zeros(size(non_centered_spectra_2,1),(zp_fact_trp-1)*chd.N_samples_sar),...
                       non_centered_spectra_2(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                       chd.N_samples_sar * zp_fact_trp,2);



h=figure;
k=surf(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:size(beams_zp_fft_plot,1),abs(beams_zp_fft_plot));
set(k, 'edgecolor','none');view(62,30);
set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
colormap('jet'); colorbar;
figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' before alignment' ' PT #' tag_pt(end:end) ],font_size)
if SAVE == 1
    figName = [plots_folder 'Stack_before_alignment' ' PT #' tag_pt(end:end) ];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
    saveas (gcf,[figName,'.fig']);
end
plot(L1BS_buffer(i_surf_stacked).shift-1,1:L1BS_buffer(i_surf_stacked).N_beams_stack);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
set(gca,'XLim',[- chd.N_samples_sar   chd.N_samples_sar ],'FontSize',font_size);
figlabels('Shift [samples]','Beams','',['Range corrections' ' PT #' tag_pt(end:end) ],font_size);
if SAVE == 1
    figName = [plots_folder 'Alignment' ' PT #' tag_pt(end:end) ];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
    saveas (gcf,[figName,'.fig']);
end
k=surf(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS_buffer(i_surf_stacked).N_beams_stack,abs(beams_zp_fft_aligned_plot));
set(k, 'edgecolor','none');
colormap('jet'); colorbar;view(62,30);
set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' aligned' ' PT #' tag_pt(end:end) ],font_size)
if SAVE == 1
    figName = [plots_folder 'Stack_after_alignment' ' PT #' tag_pt(end:end) ];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
    saveas (gcf,[figName,'.fig']);
end




zp_fact_trp=512;

non_centered_spectra=fftshift(squeeze(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:)),2);
    beams_zp_fft = fft([non_centered_spectra(:,1:chd.N_samples_sar/2),...
                       zeros(max(L1BS_buffer(i_surf_stacked).N_beams_stack),(zp_fact_trp-1)*chd.N_samples_sar),...
                       non_centered_spectra(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                       chd.N_samples_sar * zp_fact_trp,2);
                   
% samples_to_check=1:size(beams_zp_fft,2);
% i_sample_corr=0;
% if(strcmp(tag_pt(end:end),'2') || strcmp(tag_pt(end:end),'5') || strcmp(tag_pt(end:end),'8')) % central column PT
    i_sample_corr=(chd.N_samples_sar/2-28)*zp_fact_trp-1;
    samples_to_check=(chd.N_samples_sar/2-28)*zp_fact_trp:(chd.N_samples_sar/2+28)*zp_fact_trp;    
% elseif(strcmp(tag_pt(end:end),'3') || strcmp(tag_pt(end:end),'6') || strcmp(tag_pt(end:end),'9'))
%     i_sample_corr=220*zp_fact_trp-1;
%     samples_to_check=220*zp_fact_trp:250*zp_fact_trp;

% end

[max_val,max_pos]=max(abs(beams_zp_fft(:,samples_to_check)).');
max_pos = max_pos/zp_fact_trp-1/zp_fact_trp+i_sample_corr/zp_fact_trp; %reference from [0 to N] IMPORTANT!!
plot(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,sum(abs(beams_zp_fft)));
set(gca,'XLim',[1  chd.N_samples_sar],'FontSize',font_size);
% set(gca,'XLim',[mean(max_pos)-5  mean(max_pos)+5 ],'FontSize',font_size);


% %% NC reading testing
% look_i_samples_rx1_V = real(L1BS_buffer(i_surf_stacked).beam_geo_corr);
% look_q_samples_rx1_V = imag(L1BS_buffer(i_surf_stacked).beam_geo_corr);
% look_i_samples_rx2_V = real(L1BS_buffer(i_surf_stacked).beam_geo_corr_2);
% look_q_samples_rx2_V = imag(L1BS_buffer(i_surf_stacked).beam_geo_corr_2);
% 
% clear i_scale_factor q_scale_factor iq_scale_factor_rx1 look_i_samples_rx1 look_q_samples_rx1
% for i_look=1:size(look_i_samples_rx1_V,1)
%         i_scale_factor(i_look) = double(max(abs(look_i_samples_rx1_V(i_look,:))) / (2^7-1));
%         q_scale_factor(i_look) = double(max(abs(look_q_samples_rx1_V(i_look,:))) / (2^7-1));
%         iq_scale_factor_rx1(i_look)= max([i_scale_factor(i_look) q_scale_factor(i_look)])';
%         look_i_samples_rx1(i_look,:) = double((look_i_samples_rx1_V(i_look,:) ./ iq_scale_factor_rx1(i_look)));
%         look_q_samples_rx1(i_look,:) = double((look_q_samples_rx1_V(i_look,:) ./ iq_scale_factor_rx1(i_look)));
% end
% 
% % Data read
% look_i_samples_nc=ncread('CRA_IR_1S_HR_SIO_20280101T193201_20280101T193206_20230525T154740______________ISRD_SIR____TST.nc', 'data/ku/look_i_samples_rx1');
% look_q_samples_nc=ncread('CRA_IR_1S_HR_SIO_20280101T193201_20280101T193206_20230525T154740______________ISRD_SIR____TST.nc', 'data/ku/look_q_samples_rx1');
% iq_scale_factor_nc=ncread('CRA_IR_1S_HR_SIO_20280101T193201_20280101T193206_20230525T154740______________ISRD_SIR____TST.nc', 'data/ku/iq_scale_factor_rx1');
% % Waveform building
% wfm=(double(look_i_samples_nc(:,1:555,62)) + 1i*double(look_q_samples_nc(:,1:555,62))).*iq_scale_factor_nc(1:555,62).';
% wfm=wfm.';
% % Zero padding stack
% non_centered_spectra_nc=fftshift(squeeze(wfm),2);
% beams_zp_fft_nc = fft([non_centered_spectra_nc(:,1:chd.N_samples_sar/2),...
%     zeros(max(L1BS_buffer(i_surf_stacked).N_beams_stack),(zp_fact_trp-1)*chd.N_samples_sar),...
%     non_centered_spectra_nc(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
%     chd.N_samples_sar * zp_fact_trp,2);
% % Samples to check
% i_sample_corr=(chd.N_samples_sar/2-28)*zp_fact_trp-1;
% samples_to_check=(chd.N_samples_sar/2-28)*zp_fact_trp:(chd.N_samples_sar/2+28)*zp_fact_trp;
% % Maximum calculation in the zero padded stack
% [max_val_nc,max_pos_nc]=max(abs(beams_zp_fft_nc(:,samples_to_check)).');
% max_pos_nc = max_pos_nc/zp_fact_trp-1/zp_fact_trp+i_sample_corr/zp_fact_trp; %reference from [0 to N] IMPORTANT!!
% 
% [Slope_coef_A1_nc,SS,MUMU] = polyfit(valid_beams,max_pos_nc(valid_beams),1);
% Slope_A1_nc = polyval(Slope_coef_A1_nc, valid_beams,SS,MUMU); %Slope_coef_A1(1) related with datation error
% Stack_noise_A1=std(max_pos_nc(valid_beams)-Slope_A1_nc)*chd.T0_nom*cst.c/2; % units [m]
% Stack_alignment_A1=(Slope_A1_nc(end)-Slope_A1_nc(1))*chd.T0_nom*cst.c/2; % units [m]
% 
% 
% %look time testing
% look_time_nc=ncread('CRA_IR_1S_HR_SIO_20280101T193201_20280101T193206_20230525T154740______________ISRD_SIR____TST.nc', 'data/ku/look_time');


%%
figlabels('Samples','Power','',['L1B waveform #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) 'zero padding: ' num2str(zp_fact_trp) ' PT #' tag_pt(end:end) ],font_size)
    if SAVE == 1
        figName = [plots_folder 'L1B waveform' ' PT #' tag_pt(end:end) ];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
    end

if(strcmp(cnf.processing_mode,'SIN'))
   
    
    non_centered_spectra_2=fftshift(squeeze(L1BS_buffer(i_surf_stacked).beam_geo_corr_2(:,:)),2);
    beams_zp_fft_2 = fft([non_centered_spectra(:,1:chd.N_samples_sar/2),...
                       zeros(max(L1BS_buffer(i_surf_stacked).N_beams_stack),(zp_fact_trp-1)*chd.N_samples_sar),...
                       non_centered_spectra_2(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                       chd.N_samples_sar * zp_fact_trp,2);


    [max_val_2,max_pos_2]=max(abs(beams_zp_fft_2(:,samples_to_check)).');
    max_pos_2 = max_pos_2/zp_fact_trp-1/zp_fact_trp+i_sample_corr/zp_fact_trp; %reference from [0 to N] IMPORTANT!!
    zp_fact_trp=1;
    
    stack_phase_diff = compute_phase_difference (L1BS_buffer(i_surf_stacked),beams_zp_fft,beams_zp_fft_2, max_pos, max_pos_2);
%     phase_diff = mean ()
    
    B = abs(chd.y_cog_ant_2- chd.y_cog_ant);
    angle_meas = asin (chd.wv_length * (L1B.phase_difference)./(2* pi* B));
    plot(angle_meas/pi*180); hold all;
    plot(round(mean(max_pos)),angle_meas(round(mean(max_pos)))/pi*180,'rO');
    set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);
    figlabels('Samples','AoA [deg]','',['Measured Angle of Arrival: ' num2str(angle_meas(round(mean(max_pos)))/pi*180) ' PT #' tag_pt(end:end) ],font_size);
%     set(gca,'YLim',[angle_meas(round(mean(max_pos)))/pi*180-0.0001  angle_meas(round(mean(max_pos)))/pi*180+0.0001 ],'FontSize',font_size);
    set(gca,'YLim',[-0.6  +0.6],'FontSize',font_size);
    if SAVE == 1
        figName = [plots_folder 'AoA_L1B' ' PT #' tag_pt(end:end) ];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
    end 
    k=surf(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS_buffer(i_surf_stacked).N_beams_stack,abs(L1BS_buffer(i_surf_stacked).phase_diff));
    set(k, 'edgecolor','none');
    colormap('jet'); colorbar;view(62,30);
    set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);
    set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
    figlabels('Samples','Beams','Phase diff [rad]',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' aligned' ' PT #' tag_pt(end:end) ],font_size)
    if SAVE == 1
        figName = [plots_folder 'Phase difference stack' ' PT #' tag_pt(end:end) ];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
        saveas (gcf,[figName,'.fig']);
    end  
    k=surf(abs(L1BS_buffer(i_surf_stacked).coherence));
    set(k, 'edgecolor','none');
    colormap('jet'); colorbar;view(62,30);
    set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);
    set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
    figlabels('Samples','Beams','Coherence',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' aligned' ' PT #' tag_pt(end:end) ],font_size)
    if SAVE == 1
        figName = [plots_folder 'Coherence stack' ' PT #' tag_pt(end:end) ];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
        saveas (gcf,[figName,'.fig']);
    end   
    zp_fact_trp=512;   
end


%% GEOPHYSICAL CORRECTIONS

alt_rate_sat=L1BS_buffer(i_surf_stacked).alt_rate_sat;


%% IN CODE %%
wet=0; %to be readed from the data.
geophysical_correction=wet; %09/04/2016


%%
[central, valid_beams] = TRP_beam_analysis(L1BS_buffer(i_surf_stacked));

% init_beam=25; end_beam=217; % %07/11/2015 CR2
% init_beam=29; end_beam=222; % %30/04/2012 CR2
% init_beam=21; end_beam=216; % %29/05/2012 CR2
% init_beam=37; end_beam=217; % %10/05/2013 CR2
% valid_beams=init_beam:end_beam;
% % valid_beams = valid_beams(find(valid_beams~=150));
% valid_beams = valid_beams(find(valid_beams~=115));

[Slope_coef_A1,SS,MUMU] = polyfit(valid_beams,max_pos(valid_beams),1);
Slope_A1 = polyval(Slope_coef_A1, valid_beams,SS,MUMU); %Slope_coef_A1(1) related with datation error
Stack_noise_A1=std(max_pos(valid_beams)-Slope_A1)*chd.T0_nom*cst.c/2; % units [m]
Stack_alignment_A1=(Slope_A1(end)-Slope_A1(1))/length(valid_beams)*chd.T0_nom*cst.c/2; % units [m/beam]

if(strcmp(cnf.processing_mode,'SIN'))
    [Slope_coef_A2,SS,MUMU] = polyfit(valid_beams,max_pos_2(valid_beams),1);
    Slope_A2 = polyval(Slope_coef_A2, valid_beams,SS,MUMU); %Slope_coef_A1(1) related with datation error
    Stack_noise_A2=std(max_pos_2(valid_beams)-Slope_A2)*chd.T0_nom*cst.c/2; % units [m]
    Stack_alignment_A2=(Slope_A2(end)-Slope_A2(1))/length(valid_beams)*chd.T0_nom*cst.c/2; % units [m/beam]
    
end


subplot(2,1,1);
plot(1:L1BS_buffer(i_surf_stacked).N_beams_stack,max_val);
set(gca,'XLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
figlabels('Beams','Power','',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' zero padding: ' num2str(zp_fact_trp) ' PT #' tag_pt(end:end) ],font_size)
if(strcmp(cnf.processing_mode,'SIN'))
    hold all;plot(1:L1BS_buffer(i_surf_stacked).N_beams_stack,max_val_2);
end
subplot(2,1,2);
plot(1:L1BS_buffer(i_surf_stacked).N_beams_stack,max_pos);
set(gca,'XLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',font_size);
set(gca,'YLim',[max_pos(floor(L1BS_buffer(i_surf_stacked).N_beams_stack/2))-5 max_pos(floor(L1BS_buffer(i_surf_stacked).N_beams_stack/2))+5],'FontSize',font_size);
% set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',font_size);

hold all;
plot((valid_beams),Slope_A1);
figlabels('Beams','Samples','',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' alignment: ' num2str(Stack_alignment_A1) ' [m/beam] and noise: ' num2str(Stack_noise_A1) ' [m]; zero padding: ' num2str(zp_fact_trp) ' ' ' PT #' tag_pt(end:end) ],font_size)
legend('Max Pos','Slope');
if(strcmp(cnf.processing_mode,'SIN'))
    plot(1:L1BS_buffer(i_surf_stacked).N_beams_stack,max_pos_2);
    plot((valid_beams),Slope_A2);
    legend('Max Pos Rx1','Slope Rx1','Max Pos Rx2','Slope Rx2');
    figlabels('Beams','Samples','',[' Alignment Rx1: ' num2str(abs(Stack_alignment_A1)*1e3,'%4f') ' [mm/beam] and Noise: ' num2str(Stack_noise_A1*1e3,'%4f') ' [mm] ' ' PT #' tag_pt(end:end) ; ' Alignment Rx2: ' num2str(abs(Stack_alignment_A2)*1e3,'%4f') ' [mm/beam] and Noise: ' num2str(Stack_noise_A2*1e3,'%4f') ' [mm] ' ' PT #' tag_pt(end:end) ],font_size)
%' Alignment Rx1: ' num2str(Stack_alignment_A1) ' [m/beam] and Noise: ' num2str(Stack_noise_A1) ' [m] ' ' PT #' tag_pt(end:end) 
%' Alignment Rx2: ' num2str(Stack_alignment_A2) ' [m/beam] and Noise: ' num2str(Stack_noise_A2) ' [m] ' ' PT #' tag_pt(end:end) 
end

hold off

if SAVE == 1
    figName = [plots_folder 'Stack_alignment' ' PT #' tag_pt(end:end) ];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
    saveas (gcf,[figName,'.fig']);
end


%% COMPUTING RESULTS
TRP_int_delay = 0; % Svalbard


p = lla2ecef([chd.lat_trp,chd.lon_trp,chd.alt_trp],cst.flat_coeff,cst.semi_major_axis);
    x_TRP = p(:,1).';
    y_TRP = p(:,2).';
    z_TRP = p(:,3).';
TRP_coord = [x_TRP,y_TRP,z_TRP];
% TRP_coord = [L1BS_buffer(i_surf_stacked).x_surf,L1BS_buffer(i_surf_stacked).y_surf,L1BS_buffer(i_surf_stacked).z_surf];
window_corr = L1BS_buffer(i_surf_stacked).win_delay_surf;
time_stack = window_corr +((max_pos-(chd.N_samples_sar/2))*chd.T0_nom);

range_stack = time_stack *cst.c/2  + geophysical_correction-TRP_int_delay/2;
range_without_slant =range_stack-(L1BS_buffer(i_surf_stacked).slant_range_corr(:).')*chd.T0_nom*cst.c/2;

for i_beam=1:length(range_without_slant)
    time_beam(i_beam)   = L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(i_beam)).time;
%     time_burst (i_beam) = L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(i_beam)).time_rx_1st;
    SAT_coord           = [L1BS_buffer(i_surf_stacked).x_sar_sat_beam(i_beam),L1BS_buffer(i_surf_stacked).y_sar_sat_beam(i_beam),L1BS_buffer(i_surf_stacked).z_sar_sat_beam(i_beam)];
    dif                 = (SAT_coord(1)-TRP_coord(1))^2 + (SAT_coord(2)-TRP_coord(2))^2 + (SAT_coord(3)-TRP_coord(3))^2;
    theo_range(i_beam) =  sqrt(dif);
    if(strcmp(cnf.processing_mode,'SIN'))
    
        angle_theo(i_beam) = asin (L1BS_buffer(i_surf_stacked).dist/theo_range(i_beam))*180/cst.pi;
    end
end
if(strcmp(cnf.processing_mode,'SIN'))
    angle_theo_mean=-mean(angle_theo);
else
    angle_theo_mean=0;
end
theo_range_aligned=theo_range+(L1BS_buffer(i_surf_stacked).slant_range_corr(:).')*chd.T0_nom*cst.c/2;

% time interpolation:
interp_fact         = 0.001;
delta_time          = mean(diff(time_beam));
delta_time_interp   = delta_time*interp_fact;
times_STACK         = min(time_beam):delta_time_interp:max(time_beam);
% times_BURSTS        = min(time_burst):delta_time_interp:max(time_burst);

[range_theo_coef,SS,MUMU] = polyfit(time_beam,theo_range,2);
range_theo_interp = polyval(range_theo_coef, times_STACK,SS,MUMU); 
range_theo_interp_2 = interp1(time_beam, theo_range,times_STACK,'splines');
[range_meas_coef,SS,MUMU] = polyfit(time_beam(valid_beams),range_without_slant(valid_beams),2);
range_meas_interp = polyval(range_meas_coef, times_STACK,SS,MUMU); 
range_meas_interp_2 = interp1(time_beam(valid_beams), range_without_slant(valid_beams),times_STACK);

deltaR = fitfunctions_yaxis (range_meas_interp, range_theo_interp);
deltaT = -1*fitfunctions_xaxis ((range_meas_interp-deltaR), range_theo_interp); %--20110711-- reduce Y diff
if deltaT==0
  deltaT = fitfunctions_xaxis (range_theo_interp, range_meas_interp-deltaR);
end   
range_error_A1_v2 = deltaR; %[m] 
datation_error_A1_v2 = deltaT*delta_time_interp *1e6;%[micros]
[range_theo,pos_theo]=min(range_theo_interp);
[range_meas,pos_meas]=min(range_meas_interp);
datation_error_A1_v1 = (times_STACK(pos_meas)-times_STACK(pos_theo))*1e6; %[micros]
range_error_A1_v1=range_meas-range_theo;
elevation = L1BS_buffer(i_surf_stacked).alt_sat-L1BS_buffer(i_surf_stacked).win_delay_surf*cst.c/2-((max_pos-(chd.N_samples_sar/2))*chd.T0_nom*cst.c/2) + TRP_int_delay/2- geophysical_correction;


[l1b_max_val,l1b_max_pos]=max(L1B.wfm_cor_i2q2);
geoid_correction=geoidheight(L1BS_buffer(i_surf_stacked).lat_sat,L1BS_buffer(i_surf_stacked).lon_sat); 
geoid_correction_trp=geoidheight(chd.lat_trp,chd.lon_trp); 

sat_elevation_1b=L1BS_buffer(i_surf_stacked).alt_sat-L1BS_buffer(i_surf_stacked).win_delay_surf*cst.c/2 + geoid_correction;

 elevation_try=L1BS_buffer(i_surf_stacked).alt_sat-L1BS_buffer(i_surf_stacked).win_delay_surf*cst.c/2-((l1b_max_pos-(chd.N_samples_sar/2))*chd.T0_nom*cst.c/2) + TRP_int_delay/2- geophysical_correction+ geoid_correction;


disp(['Range error fitting: ' num2str(range_error_A1_v2*1e3) ' mm & minimum value: ' num2str(range_error_A1_v1*1e3) ' mm & aligned ranges: ' num2str(mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3) ' mm']);
disp(['Datation error fitting: ' num2str(datation_error_A1_v2) ' microseconds & minimum value: ' num2str(datation_error_A1_v1) ' microseconds']);
disp(['Datation accuracy: ' num2str(delta_time_interp*1e6) ' microseconds']);
disp(['Alignment: ' num2str(Stack_alignment_A1*1000) ' [mm/beam]']);
disp([' Noise: ' num2str(Stack_noise_A1*1000) ' [mm]' ]);
disp([' zero padding: ' num2str(zp_fact_trp)]);
disp(['Range bias from the multilooked stack (L1B-HR): ' num2str(elevation(floor(length(elevation)/2))-chd.alt_trp) ' cm']); %JPL:: to check if this is the correct way to compute it
angle_meas_mean = 0;
if(strcmp(cnf.processing_mode,'SIN'))
    angle_meas_mean =    angle_meas(round(mean(max_pos)))/pi*180;
    B = abs(chd.y_cog_ant_2- chd.y_cog_ant);
    disp(['AoA measured: ' num2str(angle_meas(round(mean(max_pos)))/pi*180) ' [deg] ']);
    disp(['AoA theoreical: ' num2str(angle_theo_mean) ' [deg] ']);
    disp(['AoA bias ' num2str((angle_meas_mean-angle_theo_mean)*360) 'arcseconds' ]);
    jump_AoA     = asin ((chd.wv_length * (2*pi)./(2*pi*B)))/pi*180;
    
end
%% View side parameter%%
switch mission
    case 'PICE'
        if ( sum([L1BS_buffer(1:5).lat_surf]) > sum([L1BS_buffer(6:10).lat_surf]) ) %Descending
            if (L1BS_buffer(i_surf_stacked).lon_surf > chd.lon_trp)
                sat_view = 'Right';
            else
                sat_view = 'Left';
            end
        else %Ascending
            if (L1BS_buffer(i_surf_stacked).lon_surf > chd.lon_trp)
                sat_view = 'Left';
            else
                sat_view = 'Right';
            end
        end
    case 'CR2'
        CR2_mode = filename_L1A(13:15);
        %{
        if strcmp(CR2_mode,'SIN')
            if ( sum([L1A_buffer(1:5).lat_surf]) > sum([L1A_buffer(6:10).lat_surf]) ) %Descending
                if (L1A_buffer(pos).lon_surf > chd.lon_trp)
                    sat_view = 'Right';
                else
                    sat_view = 'Left';
                end
            else %Ascending
                if (L1A_buffer(pos).lon_surf > chd.lon_trp)
                    sat_view = 'Left';
                else
                    sat_view = 'Right';
                end
            end
        elseif strcmp(CR2_mode,'SAR')
            if ( sum([L1A_buffer(1:5).lat_sar_surf]) > sum([L1A_buffer(6:10).lat_sar_surf]) ) %Descending
                if (L1A_buffer(pos).lon_sar_surf > chd.lon_trp)
                    sat_view = 'Right';
                else
                    sat_view = 'Left';
                end
            else %Ascending
                if (L1A_buffer(pos).lon_sar_surf > chd.lon_trp)
                    sat_view = 'Left';
                else
                    sat_view = 'Right';
                end
            end    
        end
        %}
end



%% Writing 
fidResults=fopen([plots_folder 'results_PT' tag_pt '.csv'],'a+');
cycle='';
%{
if chd.lat_trp>78
    results_file=dir('./../*Svalbard.csv');%Svalbard
else
    results_file=dir('./../*Crete.csv');%Crete
end

fidResults=fopen(['../'  results_file.name],'a+');
%}
fprintf(fidResults,'Cycle; Range error fitting [mm]; Range error minimum value [mm]; Range error aligned ranges [mm]; Datation error fitting [microseconds]; Datation error minimum value[microseconds]; alignment [mm/beam]; noise [mm]; AoA meas [deg]; AoA theo [deg]; wet[mm];;geophysical_correction [mm]; zero padding; n_beams; dist;\n');
%fprintf(fidResults,'%d; %s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; %s; %d; %d; %d;  \n',...
switch mission
    case 'PICE'
        if(strcmp(cnf.processing_mode,'SIN'))
        fprintf(fidResults,'%d; %s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; %s; %d; %d; %d;\n',...
            cycle,...
            (range_error_A1_v2*1e3),...
            (range_error_A1_v1*1e3),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos),
            (angle_meas),(angle_theo_mean),(wet),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),L1BS_buffer(i_surf_stacked).dist);%,...
            %(mode_L1),(mode_L2),(TRP_int_delay*1000));%,L1BS_buffer(i_surf_stacked).dist); %,instrument_range_correction_tx_rx,USO_correction)
        %lla2kml(['./'  currenfolder '_L1BS'], [L1BS_buffer.lat_surf], [L1BS_buffer.lon_surf]);
        %lla2kml_tour(['./'  currenfolder '_tour'],L1BS_buffer);
        %lla2kml(['./' currenfolder '_wet_min'], L1A_buffer(wet_min_pos).lat_sar_sat, L1A_buffer(wet_min_pos).lon_sar_sat, L1A_buffer(wet_min_pos).alt_sar_surf,'.','yellow');
        %lla2kml(['./' currenfolder '_dry_min'], L1A_buffer(dry_min_pos).lat_sar_sat, L1A_buffer(dry_min_pos).lon_sar_sat, L1A_buffer(dry_min_pos).alt_sar_surf,'.','green');
        %lla2kml(['./' currenfolder '_pos'], L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(pos)).lat_sar_sat, L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(pos)).lon_sar_sat, L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(pos)).alt_sar_surf,'.','red');
        elseif (strcmp(cnf.processing_mode,'SAR'))
            angle_meas=0;
            fprintf(fidResults,'%d; %s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; %s; %d; %d; %d;\n',...
            cycle,...
            (range_error_A1_v2*1e3),...
            (range_error_A1_v1*1e3),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos),
            (angle_meas),(angle_theo_mean),(wet),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),L1BS_buffer(i_surf_stacked).dist);        
        end       
    case 'CR2'
        fprintf(fidResults,'%s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; \n',...
            currenfolder,...
            (range_error_A1_v2),...
            (range_error_A1_v1),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos), (wet), (wet_min_pos), (dry), (dry_min_pos), (iono),...
            (pos), (wet), (dry), (iono),...
            (solid_earth),(geocentric_tide),(ocean_loading),(geophysical_correction),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),...
            (TRP_int_delay*1000),L1BS_buffer(i_surf_stacked).dist);%,sat_view); %,instrument_range_correction_tx_rx,USO_correction)
        
    case 'LRM'
        fprintf(fidResults,'%s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; \n',...
            currenfolder,...
            (range_error_A1_v2),...
            (range_error_A1_v1),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos), (wet), (wet_min_pos), (dry), (dry_min_pos), (iono),...
            (pos), (wet), (dry), (iono),...
            (solid_earth),(geocentric_tide),(ocean_loading),(geophysical_correction),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),...
            (TRP_int_delay*1000),L1BS_buffer(i_surf_stacked).dist);

end
fclose(fidResults);


h=figure;
subplot(1,2,1);
plot(range_without_slant'); hold all; plot(theo_range');
figlabels('Beams','Range [m]','','Range' ,14);
legend('Range meas', 'Range theo');
subplot(1,2,2);plot(range_without_slant-theo_range);
figlabels('Beams','Range [m]','',['Error = Range meas - Range theo, zp:' num2str(zp_fact_trp)] ,14);
                                               %L1BS_buffer(i_surf_stacked).win_delay_surf-TRP_int_delay/cst.c*2;
elevation = L1BS_buffer(i_surf_stacked).alt_sat-L1BS_buffer(i_surf_stacked).win_delay_surf*cst.c/2-((max_pos-(chd.N_samples_sar/2))*chd.T0_nom*cst.c/2) + TRP_int_delay/2- geophysical_correction;
   
disp(elevation(floor(length(elevation)/2))-chd.alt_trp);
saveas (h,[plots_folder 'Stack_range_bias_' num2str(L1BS_buffer(i_surf_stacked).surf_counter) '.png']);
close(h);

% 
