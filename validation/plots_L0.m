function plots_L0(filesBulk, cnf, chd, cst)

SAVE      = 1;
VISIBLE   = 1;
font_size=20;
set_default_plot;
plotting_bursts = 0;

plots_folder = [filesBulk.outputPath filesBulk.filename_L0(end-56:end-3) '/'];

N_bursts = get_num_record(filesBulk.filename_L0,'nb');

if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
    mkdir([plots_folder '\bursts']);
end
burst_counter = ncread(filesBulk.filename_L0, 'burst_counter');
data_record_time   = ncread(filesBulk.filename_L0, 'data_record_time');
com_position_vector = ncread(filesBulk.filename_L0, 'com_position_vector');
com_velocity_vector = ncread(filesBulk.filename_L0, 'com_velocity_vector');
com_altitude_rate = ncread(filesBulk.filename_L0, 'com_altitude_rate');
interferomter_baseline_direction_vector = ncread(filesBulk.filename_L0, 'interferomter_baseline_direction_vector');
off_nadir_pitch_angle = ncread(filesBulk.filename_L0, 'off_nadir_pitch_angle');
off_nadir_roll_angle = ncread(filesBulk.filename_L0, 'off_nadir_roll_angle');
off_nadir_yaw_angle = ncread(filesBulk.filename_L0, 'off_nadir_yaw_angle');
isp_type = ncread(filesBulk.filename_L0, 'isp_type');
instrument_mode_id = ncread(filesBulk.filename_L0, 'instrument_mode_id');
instrument_configuration_flags = ncread(filesBulk.filename_L0, 'instrument_configuration_flags');
instrument_status_flags = ncread(filesBulk.filename_L0, 'instrument_status_flags');
power_scaling_to_antenna = ncread(filesBulk.filename_L0, 'power_scaling_to_antenna');
h0_initial_altitude_instruction = ncread(filesBulk.filename_L0, 'h0_initial_altitude_instruction');
cor2_altitude_rate_estimation = ncread(filesBulk.filename_L0, 'cor2_altitude_rate_estimation');
coarse_altitude_instruction = ncread(filesBulk.filename_L0, 'coarse_altitude_instruction');
fine_altitude_instruction = ncread(filesBulk.filename_L0, 'fine_altitude_instruction');
tracker_range_calibrated = ncread(filesBulk.filename_L0, 'tracker_range_calibrated');
burst_repetition_interval = ncread(filesBulk.filename_L0, 'burst_repetition_interval');
pulse_repetition_interval = ncread(filesBulk.filename_L0, 'pulse_repetition_interval');
ambiguity_rank = ncread(filesBulk.filename_L0, 'ambiguity_rank');
number_of_pulses_in_1_radar_cycle = ncread(filesBulk.filename_L0, 'number_of_pulses_in_1_radar_cycle');
altimeter_clock = ncread(filesBulk.filename_L0, 'altimeter_clock');
rx1_complex_waveforms_i_samples = ncread(filesBulk.filename_L0, 'rx1_complex_waveforms_i_samples');
rx1_complex_waveforms_q_samples = ncread(filesBulk.filename_L0, 'rx1_complex_waveforms_q_samples');
rx2_complex_waveforms_i_samples = ncread(filesBulk.filename_L0, 'rx2_complex_waveforms_i_samples');
rx2_complex_waveforms_q_samples = ncread(filesBulk.filename_L0, 'rx2_complex_waveforms_q_samples');
%% -----------------------------
lla = ecef2lla(com_position_vector.',cst.flat_coeff,cst.semi_major_axis);
figure; plot3(lla(:,2),lla(:,1),lla(:,3));
figlabels('Longitude [degrees]','Latitude [degrees]','Altitude [m]','Orbit',font_size);
if SAVE == 1
    figName = [plots_folder '1_Orbit_L0'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; subplot(2,1,1);plot(tracker_range_calibrated);
figlabels('Record index L0','Tracker Range [m]','','',font_size);
subplot(2,1,2);plot(lla(:,3)-tracker_range_calibrated);
figlabels('Record index L0','Tracker Elevation [m]','','',font_size);
if SAVE == 1
    figName = [plots_folder '2_Tracker_range_L0'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; plot(diff(lla(:,3)));
figlabels('Record index L0','Altitude variation [m]','','',font_size);
if SAVE == 1
    figName = [plots_folder '3_Altitude_difference_L0'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; plot(sqrt(com_velocity_vector(1,:).^2+com_velocity_vector(2,:).^2+com_velocity_vector(3,:).^2));
figlabels('Record index L0','Velocity [m/s]','','',font_size);
if SAVE == 1
    figName = [plots_folder '4_Velocity_vector_L0'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
Waveforms_rx1 = rx1_complex_waveforms_i_samples + 1i.*rx1_complex_waveforms_q_samples;
zp = 1;
wfm_isp_av              = zeros(N_bursts,chd.N_samples_sar);
wd_shift                = zeros(1,N_bursts);
wfm_shift               = zeros(N_bursts,chd.N_samples_sar*zp);
N_bursts=size(Waveforms_rx1,3);
figure;
max_val=zeros(1,N_bursts);
max_pos=zeros(1,N_bursts);




for i_burst_plot = 1:N_bursts
    % FFT & power
    wfmAUX2 = abs((fft(ifftshift(squeeze(Waveforms_rx1(:,:,i_burst_plot))),chd.N_samples_sar*zp)).').^2/(chd.N_samples_sar);
%       wfmAUX2 = abs(squeeze(Waveforms_rx1(:,:,i_burst_plot)).').^2/(chd.N_samples_sar);
   if(plotting_bursts)
       if(i_burst_plot>89 && i_burst_plot<321)
       subplot(2,3,[1:2 4:5]);
            mesh(abs(fft(((fftshift(fftshift(fft((squeeze(Waveforms_rx1(:,:,i_burst_plot))),[],2),2),1)))).').^2/(chd.N_samples_sar));
%             mesh(abs(fft((((fft(((squeeze(Waveforms_rx1(:,:,i_burst_plot)))),[],2))))).').^2/(chd.N_samples_sar));
            colorbar;
            
            set(gca,'XLim',[1 256],'FontSize',12);
            %set(gca,'XLim',[1 1024],'FontSize',12); %JPLZ: added new limits ti natch the data correctly
            set(gca,'YLim',[1 64],'FontSize',12);
            %set(gca,'ZLim',[0 1.2e-10],'FontSize',12);
            %set(gca,'ZLim',[0 1e-5],'FontSize',12); %JPLZ: added new limits ti natch the data correctly
            figlabels('Samples','Beams','Power',['Burst ' num2str(i_burst_plot)],12)
            view(60,50);
            hold off;
            [max_val(i_burst_plot),max_pos(i_burst_plot)]=max(max(abs(fft(((fftshift(fftshift(fft((squeeze(Waveforms_rx1(:,:,i_burst_plot))),[],2),2),1)))).').^2/(chd.N_samples_sar)));
            subplot(2,3,3);
            plot(90:i_burst_plot,max_val(90:i_burst_plot));hold all; plot(i_burst_plot,max_val(i_burst_plot),'or');hold off;
            set(gca,'XLim',[89 321],'FontSize',12);
            %set(gca,'YLim',[0 1.2e-10],'FontSize',12);
            %set(gca,'YLim',[0 1e-5],'FontSize',12); %JPLZ: added new limits ti natch the data correctly
            figlabels('Bursts','Power','',['Burst ' num2str(i_burst_plot)],12);
            subplot(2,3,6);
            plot(90:i_burst_plot,max_pos(90:i_burst_plot));hold all; plot(i_burst_plot,max_pos(i_burst_plot),'or');hold off;
            set(gca,'XLim',[89 321],'FontSize',12);
            set(gca,'YLim',[1 256],'FontSize',12);
            set(gca,'YLim',[1 512],'FontSize',12); %JPLZ: added new limits ti natch the data correctly
            figlabels('Bursts','Sample','',['Burst ' num2str(i_burst_plot)],12)
            figName = [plots_folder '/bursts/' num2str(i_burst_plot, '%03d') '_burst_focursed.png'];
            print(gcf, figName,'-dpng');
       end% Average
   end
    for i_sample = 1:chd.N_samples_sar*zp
        %             wfm_isp_av(i_burst,i_sample) = sum(wfmAUX(:,i_sample))/chd.N_pulses_burst;
        wfm_isp_av(i_burst_plot,i_sample) = sum(wfmAUX2(:,i_sample))/chd.N_pulses_burst;
    end
end

%% -----------------------------
figure; imagesc(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_isp_av.')
set(gca,'XLim',[1 N_bursts])
set(gca,'YLim',[1 chd.N_samples_sar*zp])
xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
if SAVE == 1
    figName = [plots_folder '5_Bursts_L0'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end
%% -----------------------------
reference =  200;
for i_burst_plot = 1:N_bursts
    wd_shift(i_burst_plot) = (tracker_range_calibrated(i_burst_plot))/cst.c*2* chd.bw + 20;
    wfm_cal_gain_wdcorr(i_burst_plot,:) = circshift(wfm_isp_av(i_burst_plot,:),[0,round(wd_shift(i_burst_plot)*zp)]);
end
figure; imagesc(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_cal_gain_wdcorr.');
set(gca,'XLim',[1 N_bursts]) % --> eixos
set(gca,'YLim',[1 chd.N_samples_sar*zp]) % --> eixos
xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
colorbar;
if SAVE == 1
    figName = [plots_folder '6_Bursts_L0_tracker_range_aligned'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end
%% -----------------------------
sum_burst = sum(wfm_isp_av.');
figure; plot(sum_burst);
set(gca,'XLim',[1 N_bursts]) % --> eixos

xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Accumulated energy')
if SAVE == 1
    figName = [plots_folder '7_Burst_Energy_L0'];
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
    figName = [plots_folder '8_Range_Energy_L0'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end



close all