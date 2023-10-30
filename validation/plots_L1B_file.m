function plots_L1B_file(filesBulk, cnf, chd, cst)

SAVE      = 1;
VISIBLE   = 1;
font_size = 20;
plotting_bursts = 0;
set_default_plot;
N_bursts = get_num_record(filesBulk.filename_L1B,'nb');

[~,outputPath_size] = size(filesBulk.outputPath);

plots_folder = [filesBulk.outputPath filesBulk.filename_L1B(outputPath_size+1:end-3) '/'];
if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
end

com_position_vector = double([ncread(filesBulk.filename_L1A, 'x_pos_l1a_echo') ncread(filesBulk.filename_L1A, 'y_pos_l1a_echo') ncread(filesBulk.filename_L1A, 'z_pos_l1a_echo')].');
com_velocity_vector = double([ncread(filesBulk.filename_L1A, 'x_vel_l1a_echo') ncread(filesBulk.filename_L1A, 'y_vel_l1a_echo') ncread(filesBulk.filename_L1A, 'z_vel_l1a_echo')].');
com_altitude_rate = double(ncread(filesBulk.filename_L1A, 'alt_l1a_echo'));

lat = double(ncread(filesBulk.filename_L1B, 'lat_l1b_echo'));
lon = double(ncread(filesBulk.filename_L1B, 'lon_l1b_echo'));
alt = double(ncread(filesBulk.filename_L1B, 'alt_l1b_echo'));
x_vel = double(ncread(filesBulk.filename_L1B, 'x_vel_l1b_echo'));
y_vel = double(ncread(filesBulk.filename_L1B, 'y_vel_l1b_echo'));
z_vel =  double(ncread(filesBulk.filename_L1B, 'z_vel_l1b_echo'));
alt_rate = double(ncread(filesBulk.filename_L1B, 'orb_alt_rate_l1b_echo'));
i2q2_meas_ku_l1b_echo = double(ncread(filesBulk.filename_L1B, 'i2q2_meas_ku_l1b_echo'));

waveform_scale_factor = double(ncread(filesBulk.filename_L1B, 'waveform_scale_factor_l1b_echo'));
if(strcmp(cnf.processing_mode,'SIN'))
    phase_diff_meas_ku_l1b_echo= double(ncread(filesBulk.filename_L1B, 'phase_diff_meas_ku_l1b_echo'));
    coherence_meas_ku_l1b_echo = double(ncread(filesBulk.filename_L1B, 'coherence_meas_ku_l1b_echo'));
    i2q2_meas_ku_l1b_echo_2 = double(ncread(filesBulk.filename_L1B, 'i2q2_meas_ku_l1b_echo_2'));
    waveform_scale_factor_2 = double(ncread(filesBulk.filename_L1B, 'waveform_scale_factor_l1b_echo_2'));
end
range_ku_l1b_echo = ncread(filesBulk.filename_L1B, 'range_ku_l1b_echo');
L1B_written = readanyNETCDF_V1(filesBulk.filename_L1B);
range_scale_factor = L1B_written.attributes.range_ku_l1b_echo.scale_factor;
range_offset = L1B_written.attributes.range_ku_l1b_echo.add_offset;

% number_of_pulses_in_1_radar_cycle = ncread(filesBulk.filename_L1A, 'number_of_pulses_in_1_radar_cycle');
% altimeter_clock = ncread(filesBulk.filename_L1A, 'altimeter_clock');#
% scale_i_samples = double(ncread(filesBulk.filename_L1A, 'rx1_i_scale_factor'))/10^0.3;
% scale_q_samples = double(ncread(filesBulk.filename_L1A, 'rx1_q_scale_factor'))/10^0.3;
%rx1_complex_waveforms_i_samples = double(ncread(filesBulk.filename_L1A, 'rx1_complex_waveforms_i_samples'));
%rx1_complex_waveforms_q_samples = double(ncread(filesBulk.filename_L1A, 'rx1_complex_waveforms_q_samples'));

%%%
%{
L1B_written = readanyNETCDF_V1(filesBulk.filename_L1B);
figure;
plot(((L1B_written.data.range_ku_l1b_echo)*L1B_written.attributes.range_ku_l1b_echo.scale_factor)+L1B_written.attributes.range_ku_l1b_echo.add_offset)
figure;
mesh( (double(L1B_written.data.waveform_scale_factor_l1b_echo) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B_written.data.i2q2_meas_ku_l1b_echo))
%}
%%

%% -----------------------------
figure; plot3(lat,lon,alt);
figlabels('Longitude [degrees]','Latitude [degrees]','Altitude [m]','Orbit',font_size);
if SAVE == 1
    figName = [plots_folder '1_Orbit_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; plot(diff(lat));
figlabels('Record index L1B','Altitude variation [m]','','',font_size);
if SAVE == 1
    figName = [plots_folder '3_Altitude_difference_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; plot(sqrt(x_vel.^2+y_vel.^2+z_vel.^2));
figlabels('Record index L1B','Velocity [m/s]','','',font_size);
if SAVE == 1
    figName = [plots_folder '4_Velocity_vector_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
figure; mesh(((waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(i2q2_meas_ku_l1b_echo));
%figure; mesh(10.*log10(((waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(i2q2_meas_ku_l1b_echo)));
figlabels('Record index L1B','Samples','Power','L1B Waveforms' ,font_size);
set(gca,'YLim',[1  chd.N_samples_sar*cnf.zp_fact_range ],'FontSize',font_size);
set(gca,'XLim',[1 size(i2q2_meas_ku_l1b_echo,2)],'FontSize',font_size);
colormap('jet'); colorbar; view(128,34);
if SAVE == 1
    figName = [plots_folder 'Scaled_waveform_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end

if(strcmp(cnf.processing_mode,'SIN'))
    figure; mesh(double(coherence_meas_ku_l1b_echo));
    figlabels('Record index L1B','Samples','Power','L1B Coherence' ,font_size);
    set(gca,'YLim',[1  chd.N_samples_sar*cnf.zp_fact_range ],'FontSize',font_size);
    set(gca,'XLim',[1 size(i2q2_meas_ku_l1b_echo,2)],'FontSize',font_size);
    colormap('jet'); colorbar; view(128,34);
    if SAVE == 1
        figName = [plots_folder 'Coherence_L1B'];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
    end
    figure; mesh(double(phase_diff_meas_ku_l1b_echo));
    figlabels('Record index L1B','Samples','Phase','L1B Phase Diff' ,font_size);
    set(gca,'YLim',[1  chd.N_samples_sar*cnf.zp_fact_range ],'FontSize',font_size);
    set(gca,'XLim',[1 size(i2q2_meas_ku_l1b_echo,2)],'FontSize',font_size);
    colormap('jet'); colorbar; view(128,34);
    if SAVE == 1
        figName = [plots_folder 'Phase_diff_L1B'];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
    end

end

%% -----------------------------
figure; plot(((range_ku_l1b_echo)*range_scale_factor)+range_offset);
figlabels('Tracker Range L1B','Range [m]','','',font_size);
if SAVE == 1
    figName = [plots_folder '5_Tracker_range_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
% 87% of peak, first value

scaled_waveforms = ((waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(i2q2_meas_ku_l1b_echo);
[max_value,max_pos] = peak_retracker(scaled_waveforms);
max_pos=max_pos/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;

parab     = fitgaus(scaled_waveforms.');
parab_max = parab(:,1)/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
maxpos_gaus = parab_max';

PP = peakiness(scaled_waveforms,cnf);
delta_h = compute_snow_depth_correction (PP,chd);
delta_h =0;
i_th= 0.95; %Sea Ice Ku
  
    [tfmr_val,tfmr_pos] = threshold_first_maximum_restracker(scaled_waveforms,i_th,cnf);
    tfmr_pos=tfmr_pos/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
    height_tfmr        = ((alt - (range_ku_l1b_echo-(chd.N_samples_sar/2 - tfmr_pos.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');
    height_peak     = ((alt - (range_ku_l1b_echo-(chd.N_samples_sar/2 - max_pos.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');
    height_gaussian = ((alt - (range_ku_l1b_echo-(chd.N_samples_sar/2 - maxpos_gaus.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');

    
i_th= 0.87; % Ocean
% 87% of peak, first value
waveform_shape = ((waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(i2q2_meas_ku_l1b_echo);
[max_value,max_pos]=(max(((waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(i2q2_meas_ku_l1b_echo)));

[~,records] = size(max_value);
pc87_value=zeros(1,records);
pc87_pos=zeros(1,records);
for i=1:records
    [pc87_value(i),pc87_pos(i)] = min( abs(waveform_shape(1:max_pos,i)-(i_th*max_value(i))) );
end

    figure; plot(height_peak);
    hold all; plot(height_tfmr); plot(height_gaussian);
    figlabels('L2 Records','Elevation [m]','','',font_size);
    legend(['Peak > Mean: ' num2str(mean(height_peak)) ' m , Std: ' num2str(std(height_peak)) ' m'],...
        ['TFMR ' num2str(i_th) ' > Mean: ' num2str(mean(height_tfmr)) ' m , Std: ' num2str(std(height_tfmr)) ' m'],...
        ['Gaussian > Mean: ' num2str(mean(height_gaussian)) ' m , Std: ' num2str(std(height_gaussian)) ' m']);

    if SAVE == 1
        figName = [plots_folder '7_Retrackers_TFMR_' num2str(i_th) '_L1B'];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
    end

    


close all
%% -----------------------------
figure; plot(max_value);
figlabels('Peak Power','Power','','',font_size);
if SAVE == 1
    figName = [plots_folder '8_Peak_power_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end
%% -----------------------------
close all;