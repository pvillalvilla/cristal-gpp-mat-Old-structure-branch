% Instrument Gain correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% ---------------------------------------------------------
% Objective: compute gain/attenuation applied to SAR wfms
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [INSTR_GAIN_OUT] = instrument_gain_corr(ISP,Process_ID,chd)

if( strcmp(chd.band, 'Ku') )
    att_conv = attenuator_table_att_ku1(ISP.attcode);
    gain_isp = ISP.gain_ku1;
    shift_isp = ISP.shift_ku1;
    g_scaling = ISP.g_scaling_ku1;
    power_cal1 = ISP.power_cal1_ku1;
    antenna_gain_two_way = 2 * chd.antenna_gain_ku1;
    
elseif( strcmp(chd.band, 'Ku2') )
    att_conv = attenuator_table_att_ku2(ISP.attcode);
    gain_isp = ISP.gain_ku2;
    shift_isp = ISP.shift_ku2;
    g_scaling = ISP.g_scaling_ku2;
    power_cal1 = ISP.power_cal1_ku2;
    antenna_gain_two_way = chd.antenna_gain_ku1 + chd.antenna_gain_ku2;
    
elseif( strcmp(chd.band, 'Ka') )
    att_conv = attenuator_table_att_ka(ISP.attcode);
    gain_isp = ISP.gain_ka;
    shift_isp = ISP.shift_ka;
    g_scaling = ISP.g_scaling_ka;
    power_cal1 = ISP.power_cal1_ka;
    antenna_gain_two_way = 2 * chd.antenna_gain_ka;
    
end

if strcmp(Process_ID, 'RAW')
    var_digital_gain = 20*log10( gain_isp / (2^11 * 2^shift_isp) );
    residual_onboard_proc = chd.residual_onboard_proc_sar_raw;
    
elseif strcmp(Process_ID, 'RMC')
    var_digital_gain = 20*log10( 1 / 2^shift_isp );
    residual_onboard_proc = chd.residual_onboard_proc_sar_rmc;

end

% in dB
instr_gain = -att_conv + var_digital_gain - g_scaling + ...
    power_cal1 + residual_onboard_proc - antenna_gain_two_way;

%% Output
INSTR_GAIN_OUT.att_conv = att_conv;
INSTR_GAIN_OUT.instr_gain = instr_gain;
INSTR_GAIN_OUT.var_digital_gain = var_digital_gain;
INSTR_GAIN_OUT.band = ISP.band;

end