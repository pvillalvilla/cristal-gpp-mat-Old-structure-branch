% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: - Correct every echo of a burst from power and phase variations
% within a burst using instrumental calibration corrections (CAL1).
%  -Correct waveforms for the CAL2 provided by instrumental calibration
%   corrections (CAL4).
%
% INPUTs:
%
%  -N_ku_pulses_burst_chd                         uc	CHD file (�2.2.1)
%  -N_samples_sar_chd                           us	CHD file (�2.2.1)
%  -fai_fine_shift_number_chd                   us	CHD file (�2.2.1)
%  -gain_corr_instr_sar                     dB	do	Instrument Gain (�7.4)
%  -burst_phase_array_cor_cal1_sar          rad	do	CAL1 selection (�5.1)
%  -burst_power_array_cor_cal1_sar_rep      dB	do	CAL1 selection (�5.1)
%  -wfm_cal2_science_sar                    FFT ss	CAL2 selection (�5.2)
%  -fai_sar                                 s	do	Preliminary Window Delay (�7.3)
%  -wfm_sar_reversed                        FFT do  OnBoard reversion (�7.43.1)
%
% OUTPUTs:
%  -burst_phase_array_cor_cal1_sar          rad do
%  -burst_power_array_cor_cal1_sar          dB  do
%  -wfm_cal2_science_sar                    FFT do
%  -wfm_cal_corrected                       FFT do
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WFM_COR_OUT] = waveforms_correction (CH1, CH2, CAL1, CAL4, cnf, chd)

% Band selection is selected outside as the input depends on the band
% ONBRD_REV_OUT and INSTR_GAIN_OUT

%% Apply gain correction
wfm_gain_cor_ch1 = CH1.wfm_reversed ./ 10^( CH1.instr_gain / 20 );

%% CAL1 Application
if cnf.flag_cal1_intraburst_corrections
    wfm_cal1_cor_ch1 = wfm_gain_cor_ch1 .* 10.^( CAL1.burst_power_array_cor / 20 )...
        .* exp( 1i*CAL1.burst_phase_array_cor );
    
else
    wfm_cal1_cor_ch1 = wfm_gain_cor_ch1;
end

%% SARIn
if strcmp(CH2.band, 'Ku2')
    % Apply Gain correction
    wfm_gain_cor_ch2 = CH2.wfm_reversed ./ 10^( CH2.instr_gain / 20 );
    
    % CAL1 Application
    if cnf.flag_cal1_intraburst_corrections
        wfm_cal1_cor_ch2 = wfm_gain_cor_ch2 .* 10.^( CAL1.burst_power_array_cor / 20 )...
            .* exp( 1i*CAL1.burst_phase_array_cor );
        
    else
        wfm_cal1_cor_ch2 = wfm_gain_cor_ch2;
    end
    
    % CAL4 Application
    phase_difference_corr = CAL4.phase_diff_cor + chd.ext_phase_difference;
    wfm_cal4_cor = wfm_cal1_cor_ch2 .* exp ( 1i*phase_difference_corr);
    
    
    %% Waveforms alignment
    samples_pb_chd = zeros(chd.N_pulses_burst,1);
    samples_pb_chd(:) = 1:chd.N_pulses_burst;
    if cnf.flag_height_rate_application
        samples_sar = 1:chd.N_samples_sar;
    else
        samples_sar = 1:chd.N_samples_rmc_onboard;
    end
    
    % Computation path delay
    delta_win_delay = CH2.win_delay - CH1.win_delay;
    % Range shift frequency channel 2
    wfm_shifted_ch2 = wfm_cal4_cor .* exp ( 1i*2*cst.pi * 2 *...
        samples_pb_chd * delta_win_delay * CH2.pri .*...
        ( samples_sar - length(samples_sar)/2) ./ length(samples_sar) );
    % Update window delay channel 2 to match channel 1
    CH2.win_delay = CH1.win_delay;
    
    % Output
    WFM_COR_OUT.wfm_shifted_ch2  = wfm_shifted_ch2;
    WFM_COR_OUT.win_delay_ch2 = CH2.win_delay;
    WFM_COR_OUT.band_ch2 = CH1.band;

    
    
else
    % Output
    WFM_COR_OUT.wfm_cor = wfm_cal1_cor_ch1;
    WFM_COR_OUT.band_ch1 = CH1.band;

    
    
end

end
