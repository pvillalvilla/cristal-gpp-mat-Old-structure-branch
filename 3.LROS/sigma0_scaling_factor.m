% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: - The sigma-0 scaling factor is used by level 2 
% processing for computing the backscatter coefficient of
% the surface where the echo is reflected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SCL_FACT_OUT] = sigma0_scaling_factor(L1A, ORBIT,chd, cst)

% Band selection
if( strcmp(L1A.band, 'Ku') )
    antenna_gain_tx = chd.antenna_gain_tx_ku1;
    antenna_gain_rx = chd.antenna_gain_rx_ku1;
    wv_length = cst.c / chd.freq_ku;
    BW = chd.bw_ku;
    
elseif( strcmp(L1A.band, 'Ku2') )
    antenna_gain_tx = chd.antenna_gain_tx_ku1;
    antenna_gain_rx = chd.antenna_gain_rx_ku2;
    wv_length = cst.c / chd.freq_ku;
    BW = chd.bw_ku;
    
elseif( strcmp(L1A.band, 'Ka') )
    antenna_gain_tx = chd.antenna_gain_tx_ka;
    antenna_gain_rx = chd.antenna_gain_rx_ka;
    wv_length = cst.c / chd.freq_ka;
    BW = chd.bw_ku;
    
end

s0_scaling_factor_i = 4*cst.pi / ( antenna_gain_tx * ...
    antenna_gain_rx * wv_length^2 );

s0_scaling_factor_lros = s0_scaling_factor_i * 16*cst.pi * BW *...
    ORBIT.alt_lros^3 * ( (cst.earth_radius + ORBIT.alt_lros) / cst.earth_radius ) / cst.c;

s0_scaling_factor_lros_dB = 10 * log10( s0_scaling_factor_lros );

%% Ouputs
SCL_FACT_OUT.s0_scaling_factor_lros = s0_scaling_factor_lros;
SCL_FACT_OUT.s0_scaling_factor_lros_dB = s0_scaling_factor_lros_dB;

end