function [L2] = wfm_characteristics (L1B, L2, cnf_L2, chd)
% ** flag_waveform ??
if (cnf_L2.flag_waveform == 0) % Rx1 (always used for Ka band)
    L2.power_waveform=L1B.power_waveform_rx1;
    L2.waveform_scale_factor= L1B.waveform_scale_factor_rx1;
elseif (cnf_L2.flag_waveform == 1) % Rx2 (never used for Ka band)
    L2.power_waveform=L1B.power_waveform_rx2;
    L2.waveform_scale_factor= L1B.waveform_scale_factor_rx2;
elseif (cnf_L2.flag_waveform == 2) % Combined (never used for Ka band)
    L2.power_waveform=L1B.power_waveform_comb;
    L2.waveform_scale_factor= L1B.waveform_scale_factor_comb;
end
if strcmp(chd.band , 'Ka')% No selection among antenna to be performed for Ka band
    L2.off_nadir_pitch_angle= L1B.off_nadir_pitch_angle';
    L2.off_nadir_roll_angle= L1B.off_nadir_roll_angle';
    L2.off_nadir_yaw_angle= L1B.off_nadir_yaw_angle';
elseif (cnf_L2.flag_waveform == 0) % Rx1
    L2.off_nadir_pitch_angle= L1B.off_nadir_pitch_angle_ant1';
    L2.off_nadir_roll_angle= L1B.off_nadir_roll_angle_ant1';
    L2.off_nadir_yaw_angle= L1B.off_nadir_yaw_angle_ant1';
elseif (cnf_L2.flag_waveform == 1) % Rx2
    L2.off_nadir_pitch_angle= L1B.off_nadir_pitch_angle_ant2';
    L2.off_nadir_roll_angle= L1B.off_nadir_roll_angle_ant2';
    L2.off_nadir_yaw_angle= L1B.off_nadir_yaw_angle_ant2';
elseif (cnf_L2.flag_waveform == 2) % Combined
    L2.off_nadir_pitch_angle = (L1B.off_nadir_pitch_angle_ant1'+L1B.off_nadir_pitch_angle_ant2')/2;
    L2.off_nadir_roll_angle = (L1B.off_nadir_roll_angle_ant1'+L1B.off_nadir_roll_angle_ant2')/2;
    L2.off_nadir_yaw_angle = (L1B.off_nadir_yaw_angle_ant1'+L1B.off_nadir_yaw_angle_ant2')/2;
end
L2.scaled_waveform = L2.power_waveform.*L2.waveform_scale_factor;
L2.scaled_waveform = L2.scaled_waveform';
abs_waveform = 1:1:chd.N_samples_sar; %samples_ov in HR mode

%1.------ Noise Floor
L2.mean_noise = mean(L2.scaled_waveform(cnf_L2.noise_first_sample_physical_rtk:cnf_L2.noise_last_sample_physical_rtk,:),1);
L2.std_floor  = std(L2.scaled_waveform(cnf_L2.noise_first_sample_physical_rtk:cnf_L2.noise_last_sample_physical_rtk,:),1);
%2.------ Peakiness
L2.noise_floor = mean(L2.scaled_waveform(cnf_L2.noise_first_sample_physical_rtk: cnf_L2.noise_last_sample_physical_rtk),1);
[L2.peakiness,wvf_max] = peakiness(L2.scaled_waveform, L2.noise_floor, cnf_L2,chd); % peakiness

%3.------ Prominence
% waveform_db = 10*log10(L2.scaled_waveform);
% [~,~,~,L2.prominence] = findpeaks(waveform_db); %**
% %4.------Num Peaks
% L2.number_of_peaks = length(L2.prominence);
% %5.------ Slope of the trailing edge
% fist_te = min([ind_max, cnf_L2.trailing_edge_first_sample]);
% L2.slope_trailing_edge =polyfit(abs_waveform(fist_te:cnf_L2.trailing_edge_last_sample),...
%     L2.scaled_waveform(fist_te:cnf_L2.trailing_edge_last_sample), 1);
%6.------ PSNR
L2.peak_noise_ratio = wvf_max./L2.noise_floor;

end