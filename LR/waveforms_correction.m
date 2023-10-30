% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% ---------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WFM_COR_OUT] = waveforms_correction (CH1, CH2, cnf, chd)

% Band selection is selected outside as the input depends on the band
instr_gain_lrm_ch1 = CH1.instr_gain_lrm;

%% Instrument gain application
wfm_cor_lrm_ch1 = CH1.wfm_i2q2_lrm ./ ( 10^(instr_gain_lrm_ch1/10) *...
    chd.N_pulses_burst * chd.N_bursts_rc );

%% PNR computation
thermal_noise_ch1 = sum(wfm_cor_lrm_ch1(cnf.thermal_noise_first_range_bin:...
    cnf.thermal_noise_first_range_bin+cnf.thermal_noise_width-1) );
PNR_ch1 = max(wfm_cor_lrm_ch1) ./ thermal_noise_ch1;

%% SARIn
if strcmp(CH2.band, 'Ku2')
    % Instrument gain application
    instr_gain_lrm_ch2 = CH2.instr_gain_lrm;
    wfm_cor_lrm_ch2 = CH2.wfm_i2q2_lrm ./ ( 10^(instr_gain_lrm_ch2/10 ) *...
        chd.N_pulses_burst * chd.N_bursts_rc );
    
    % PNR computation
    thermal_noise_ch2 = sum(wfm_cor_lrm_ch2(cnf.thermal_noise_first_range_bin:...
        cnf.thermal_noise_first_range_bin+cnf.thermal_noise_width-1) );
    PNR_ch2 = max(wfm_cor_lrm_ch2) ./ thermal_noise_ch2;
    
    % Waveforms alignment
    wfm_combined = ( wfm_cor_lrm_ch1 + wfm_cor_lrm_ch2 ) / 2;
        
    WFM_COR_OUT.wfm_cor_lrm_ch2 = wfm_cor_lrm_ch2;
    WFM_COR_OUT.thermal_noise_ch2 = thermal_noise_ch2;
    WFM_COR_OUT.PNR_ch2 = PNR_ch2;
    WFM_COR_OUT.chd.band_ch2 = CH2.band;
    WFM_COR_OUT.wfm_combined = wfm_combined;
      
end

%% Output
WFM_COR_OUT.wfm_cor_lrm_ch1 = wfm_cor_lrm_ch1;
WFM_COR_OUT.thermal_noise_ch1 = thermal_noise_ch1;
WFM_COR_OUT.PNR_ch1 = PNR_ch1;
WFM_COR_OUT.band_ch1 = CH1.band;



end
