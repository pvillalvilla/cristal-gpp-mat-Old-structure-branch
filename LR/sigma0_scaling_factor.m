% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: - The sigma-0 scaling factor is used by level 2 
% processing for computing the backscatter coefficient of
% the surface where the echo is reflected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SCL_FACT_OUT] = sigma0_scaling_factor(chd, cst, alt_lrm)

% Band selection
if( strcmp(ISP.band, 'Ku') )
    antenna_gain_tx = chd.antenna_gain_tx_ku1;
    antenna_gain_rx = chd.antenna_gain_rx_ku1;
    wv_length = cst.c / chd.freq_ku;
    BW = chd.bw_ku;
    
elseif( strcmp(ISP.band, 'Ku2') )
    antenna_gain_tx = chd.antenna_gain_tx_ku1;
    antenna_gain_rx = chd.antenna_gain_rx_ku2;
    wv_length = cst.c / chd.freq_ku;
    BW = chd.bw_ku;
    
elseif( strcmp(ISP.band, 'Ka') )
    antenna_gain_tx = chd.antenna_gain_tx_ka;
    antenna_gain_rx = chd.antenna_gain_rx_ka;
    wv_length = cst.c / chd.freq_ka;
    BW = chd.bw_ku;
    
end

s0_scaling_factor_i = 4*cst.pi / ( antenna_gain_tx * ...
    antenna_gain_rx * wv_length^2 );

s0_scaling_factor_lrm = s0_scaling_factor_i * 16*cst.pi * BW *...
    alt_lrm^3 * ( (cst.earth_radius + alt_lrm) / cst.earth_radius ) / cst.c;

s0_scaling_factor_lrm_dB = 10 * log10( s0_scaling_factor_lrm );

%% Ouputs
SCL_FACT_OUT.s0_scaling_factor_lrm = s0_scaling_factor_lrm;
SCL_FACT_OUT.s0_scaling_factor_lrm_dB = s0_scaling_factor_lrm_dB;

end