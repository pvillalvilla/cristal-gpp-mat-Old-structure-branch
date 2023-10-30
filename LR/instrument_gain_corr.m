% Instrument Gain correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% ---------------------------------------------------------
% Objective: compute gain/attenuation applied to SAR wfms
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [INSTR_GAIN_OUT] = instrument_gain_corr(ISP,chd)

if( strcmp(ISP.band, 'Ku') )
    att_conv_lrm = ISP.attenuator_table_att_ku1(ISP.attcode_lrm);
    gain_isp = ISP.gain_ku1;
    shift_isp = ISP.shift_ku1;
    g_scaling = ISP.g_scaling_ku1;
    power_cal1 = ISP.power_cal1_ku1;
    antenna_gain_two_way = 2*chd.antenna_gain_ku1;
    
elseif( strcmp(ISP.band, 'Ku2') )
    att_conv_lrm = ISP.attenuator_table_att_ku2(ISP.attcode_lrm);
    gain_isp = ISP.gain_ku2;
    shift_isp = ISP.shift_ku2;
    g_scaling = ISP.g_scaling_ku2;
    power_cal1 = ISP.power_cal1_ku2;
    antenna_gain_two_way = chd.antenna_gain_ku1 + chd.antenna_gain_ku2;
    
elseif( strcmp(ISP.band, 'Ka') )
    att_conv_lrm = ISP.attenuator_table_att_ka(ISP.attcode_lrm);
    gain_isp = ISP.gain_ka;
    shift_isp = ISP.shift_ka;
    g_scaling = ISP.g_scaling_ka;
    power_cal1 = ISP.power_cal1_ka;
    antenna_gain_two_way = 2*chd.antenna_gain_ka;
    
end

var_digital_gain_lrm = 10*log10( gain_isp / (2^11 * 2^shift_isp) );

instr_gain_lrm = -att_conv_lrm + var_digital_gain_lrm - g_scaling + ...
    power_cal1 + chd.residual_onboard_proc_lrm - antenna_gain_two_way;

%% Output
INSTR_GAIN_OUT.att_conv_lrm = att_conv_lrm;
INSTR_GAIN_OUT.instr_gain_lrm = instr_gain_lrm;
INSTR_GAIN_OUT.var_digital_gain_lrm = var_digital_gain_lrm;
INSTR_GAIN_OUT.band = ISP.band;



end