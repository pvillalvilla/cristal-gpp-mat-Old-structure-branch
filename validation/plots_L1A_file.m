function plots_L1A_file(filesBulk, cnf, chd, cst)

SAVE      = 1;
VISIBLE   = 0;
font_size = 20;
plotting_bursts = 0;
set_default_plot;
N_bursts = get_num_record(filesBulk.filename_L1A,'nb');

[~,outputPath_size] = size(filesBulk.outputPath);

plots_folder = [filesBulk.outputPath filesBulk.filename_L1A(outputPath_size+1:end-3) '/'];
if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
    mkdir([plots_folder '\bursts']);
end

burst_counter = double(ncread(filesBulk.filename_L1A, 'burst_count_prod_l1a_echo'));
data_record_time   = double(ncread(filesBulk.filename_L1A, 'time_l1a_echo'));
com_position_vector = double([ncread(filesBulk.filename_L1A, 'x_pos_l1a_echo') ncread(filesBulk.filename_L1A, 'y_pos_l1a_echo') ncread(filesBulk.filename_L1A, 'z_pos_l1a_echo')].');
com_velocity_vector = double([ncread(filesBulk.filename_L1A, 'x_vel_l1a_echo') ncread(filesBulk.filename_L1A, 'y_vel_l1a_echo') ncread(filesBulk.filename_L1A, 'z_vel_l1a_echo')].');
com_altitude_rate = double(ncread(filesBulk.filename_L1A, 'alt_l1a_echo'));
% interferomter_baseline_direction_vector = ncread(filesBulk.filename_L1A, 'interferomter_baseline_direction_vector');
off_nadir_pitch_angle = ncread(filesBulk.filename_L1A, 'roll_mispointing_l1a_echo');
off_nadir_roll_angle = ncread(filesBulk.filename_L1A, 'pitch_mispointing_l1a_echo');
off_nadir_yaw_angle = ncread(filesBulk.filename_L1A, 'yaw_mispointing_l1a_echo');
% isp_type = ncread(filesBulk.filename_L1A, 'isp_type');
% instrument_mode_id = ncread(filesBulk.filename_L1A, 'instrument_mode_id');
% instrument_configuration_flags = ncread(filesBulk.filename_L1A, 'instrument_configuration_flags');
% instrument_status_flags = ncread(filesBulk.filename_L1A, 'instrument_status_flags');
% power_scaling_to_antenna = ncread(filesBulk.filename_L1A, 'power_scaling_to_antenna');
h0_initial_altitude_instruction = ncread(filesBulk.filename_L1A, 'h0_applied_l1a_echo');
cor2_altitude_rate_estimation = ncread(filesBulk.filename_L1A, 'cor2_applied_l1a_echo');
% coarse_altitude_instruction = ncread(filesBulk.filename_L1A, 'coarse_altitude_instruction');
% fine_altitude_instruction = ncread(filesBulk.filename_L1A, 'fine_altitude_instruction');
tracker_range_calibrated = ncread(filesBulk.filename_L1A, 'range_l1a_echo');
burst_repetition_interval = ncread(filesBulk.filename_L1A, 'bri_l1a_echo');
pulse_repetition_interval = ncread(filesBulk.filename_L1A, 'pri_l1a_echo');
ambiguity_rank = ncread(filesBulk.filename_L1A, 'ambiguity_order_l1a_echo');
% number_of_pulses_in_1_radar_cycle = ncread(filesBulk.filename_L1A, 'number_of_pulses_in_1_radar_cycle');
% altimeter_clock = ncread(filesBulk.filename_L1A, 'altimeter_clock');#
% scale_i_samples = double(ncread(filesBulk.filename_L1A, 'rx1_i_scale_factor'))/10^0.3;
% scale_q_samples = double(ncread(filesBulk.filename_L1A, 'rx1_q_scale_factor'))/10^0.3;
rx1_complex_waveforms_i_samples = double(ncread(filesBulk.filename_L1A, 'rx1_complex_waveforms_i_samples'));
rx1_complex_waveforms_q_samples = double(ncread(filesBulk.filename_L1A, 'rx1_complex_waveforms_q_samples'));
if(strcmp(cnf.processing_mode,'SIN'))
    rx2_complex_waveforms_i_samples = double(ncread(filesBulk.filename_L1A, 'rx2_complex_waveforms_i_samples'));
    rx2_complex_waveforms_q_samples = double(ncread(filesBulk.filename_L1A, 'rx2_complex_waveforms_q_samples'));
end


%% -----------------------------
lla = ecef2lla(com_position_vector.',cst.flat_coeff,cst.semi_major_axis);
figure; plot3(lla(:,2),lla(:,1),lla(:,3));
figlabels('Longitude [degrees]','Latitude [degrees]','Altitude [m]','Orbit',font_size);
if SAVE == 1
    figName = [plots_folder '1_Orbit_L1A'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; subplot(2,1,1);plot(tracker_range_calibrated);
figlabels('Record index L1A','Tracker Range [m]','','',font_size);
subplot(2,1,2);plot(lla(:,3)-tracker_range_calibrated);
figlabels('Record index L1A','Tracker Elevation [m]','','',font_size);
if SAVE == 1
    figName = [plots_folder '2_Tracker_range_L1A'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; plot(diff(lla(:,3)));
figlabels('Record index L1A','Altitude variation [m]','','',font_size);
if SAVE == 1
    figName = [plots_folder '3_Altitude_difference_L1A'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; plot(sqrt(com_velocity_vector(1,:).^2+com_velocity_vector(2,:).^2+com_velocity_vector(3,:).^2));
figlabels('Record index L1A','Velocity [m/s]','','',font_size);
if SAVE == 1
    figName = [plots_folder '4_Velocity_vector_L1A'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
% figure; subplot(2,1,1);plot(scale_i_samples);
% figlabels('Record index L1A','I Scale Factor ','','',font_size);
% subplot(2,1,2);plot(scale_q_samples);
% figlabels('Record index L1A','Q Scale Factor ','','',font_size);
% if SAVE == 1
%     figName = [plots_folder '9_Scaling_factors_L1A'];
%     print(gcf, figName,'-dpng');
%     saveas (gcf,[figName,'.png']);
% end



%% -----------------------------
zp = 1;
wfm_isp_av              = zeros(N_bursts,chd.N_samples_sar);
wfm_isp_av_2              = zeros(N_bursts,chd.N_samples_sar);
wd_shift                = zeros(1,N_bursts);
wfm_shift               = zeros(N_bursts,chd.N_samples_sar*zp);
figure;
max_val=zeros(1,N_bursts);
max_pos=zeros(1,N_bursts);
for i_burst_plot = 1:N_bursts
    Waveforms_rx1(:,:) =          squeeze(rx1_complex_waveforms_i_samples(:,:,i_burst_plot)) ...
                            + 1i.*squeeze(rx1_complex_waveforms_q_samples(:,:,i_burst_plot));
                        
     if(strcmp(cnf.processing_mode,'SIN'))
             Waveforms_rx2(:,:) =          squeeze(rx2_complex_waveforms_i_samples(:,:,i_burst_plot)) ...
                            + 1i.*squeeze(rx2_complex_waveforms_q_samples(:,:,i_burst_plot));
     end
    % FFT & power
           
        wfmAUX2 = abs((fft(ifftshift(squeeze(Waveforms_rx1(:,:))),chd.N_samples_sar*zp)).').^2/(chd.N_samples_sar);
        if(strcmp(cnf.processing_mode,'SIN'))
            wfmAUX2_2 = abs((fft(ifftshift(squeeze(Waveforms_rx2(:,:))),chd.N_samples_sar*zp)).').^2/(chd.N_samples_sar);
        end
    if(plotting_bursts) 
        if(i_burst_plot>89 && i_burst_plot<321)
    %     wfmAUX2 = ((abs(Waveforms_rx1)).').^2/(chd.N_samples_sar);
            subplot(2,3,[1:2 4:5]);
%             mesh(abs(fft(((fftshift(fft((fftshift(squeeze(Waveforms_rx1(:,:)),1)),[],2),2)))).').^2/(chd.N_samples_sar));
            mesh(abs(fft(((fftshift(fft(((squeeze(Waveforms_rx1(:,:)))),[],2),2)))).').^2/(chd.N_samples_sar));
            colorbar;
            
            set(gca,'XLim',[1 256],'FontSize',12);
            set(gca,'XLim',[1 1024],'FontSize',12); %JPLZ: redefined limits to match data
            set(gca,'YLim',[1 64],'FontSize',12);
%             set(gca,'ZLim',[0 1.5e-10],'FontSize',12);
            figlabels('Samples','Beams','Power',['Burst ' num2str(i_burst_plot)],12)
            view(60,50);
            hold off;
            [max_val(i_burst_plot),max_pos(i_burst_plot)]=max(max(abs(fft(((fftshift(fft(((squeeze(Waveforms_rx1(:,:)))),[],2),2)))).').^2/(chd.N_samples_sar)));
            subplot(2,3,3);
            plot(90:i_burst_plot,max_val(90:i_burst_plot));hold all; plot(i_burst_plot,max_val(i_burst_plot),'or');hold off;
            set(gca,'XLim',[89 321],'FontSize',12);
%             set(gca,'YLim',[0 1.5e-10],'FontSize',12);
            figlabels('Bursts','Power','',['Burst ' num2str(i_burst_plot)],12);
            subplot(2,3,6);
            plot(90:i_burst_plot,max_pos(90:i_burst_plot));hold all; plot(i_burst_plot,max_pos(i_burst_plot),'or');hold off;
            set(gca,'XLim',[89 321],'FontSize',12);
            set(gca,'YLim',[1 256],'FontSize',12);
            set(gca,'YLim',[1 512],'FontSize',12); %JPLZ: redefined limits to match data
            figlabels('Bursts','Sample','',['Burst ' num2str(i_burst_plot)],12)
            figName = [plots_folder '/bursts/' num2str(i_burst_plot, '%03d') '_burst_focursed.png'];
            print(gcf, figName,'-dpng');
            % Average
        end
    end
    for i_sample = 1:chd.N_samples_sar*zp
        %             wfm_isp_av(i_burst,i_sample) = sum(wfmAUX(:,i_sample))/chd.N_pulses_burst;
        wfm_isp_av(i_burst_plot,i_sample) = sum(wfmAUX2(:,i_sample))./chd.N_pulses_burst;
        if(strcmp(cnf.processing_mode,'SIN'))
             wfm_isp_av_2(i_burst_plot,i_sample) = sum(wfmAUX2_2(:,i_sample))./chd.N_pulses_burst;
        end
    end
end

%% -----------------------------
figure; imagesc(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_isp_av.')
set(gca,'XLim',[1 N_bursts])
set(gca,'YLim',[1 chd.N_samples_sar*zp])
xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
if SAVE == 1
    figName = [plots_folder '5_Bursts_L1A'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end

if(strcmp(cnf.processing_mode,'SIN'))
    figure; imagesc(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_isp_av_2.')
    set(gca,'XLim',[1 N_bursts])
    set(gca,'YLim',[1 chd.N_samples_sar*zp])
    xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
    ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
    if SAVE == 1
        figName = [plots_folder '5_Bursts_L1A_rx2'];
        print(gcf, figName,'-dpng')
        saveas (gcf,[figName,'.png'])
    end
    
end
%% -----------------------------
reference =  200;
for i_burst_plot = 1:N_bursts
    wd_shift(i_burst_plot) = (tracker_range_calibrated(i_burst_plot))/cst.c*2* chd.T0_nom + 20;
    wfm_cal_gain_wdcorr(i_burst_plot,:) = circshift(wfm_isp_av(i_burst_plot,:),[0,round(wd_shift(i_burst_plot)*zp)]);
end
% figure; imagesc(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_cal_gain_wdcorr.');
% set(gca,'XLim',[1 N_bursts]) % --> eixos
% set(gca,'YLim',[1 chd.N_samples_sar*zp]) % --> eixos
% xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
% ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
% colorbar;
% if SAVE == 1
%     figName = [plots_folder '6_Bursts_L1A_tracker_range_aligned'];
%     print(gcf, figName,'-dpng')
%     saveas (gcf,[figName,'.png'])
% end
%% -----------------------------
sum_burst = sum(wfm_isp_av.');
figure; plot(sum_burst);
set(gca,'XLim',[1 N_bursts]) % --> eixos
xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Accumulated energy')
if SAVE == 1
    figName = [plots_folder '7_Burst_Energy_L1A'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end
%% -----------------------------
sum_range = sum(wfm_isp_av);
figure; plot(sum_range);
set(gca,'XLim',[1 chd.N_samples_sar*zp]) % --> eixos
xlab = get(gca,'XLabel'); set(xlab,'String','Range samples')
ylab = get(gca,'YLabel'); set(ylab,'String','Accumulated energy')
if SAVE == 1
    figName = [plots_folder '8_Range_Energy_L1A'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end

close all;

