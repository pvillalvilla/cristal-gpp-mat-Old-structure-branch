% PRELIMINARY WINDOW DELAY ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: The computation of the waveforms preliminary window delay.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WIN_DELAY_OUT] = window_delay (T0_lrm_pre_dat,ATTITUDE,ORBIT,CAL1, ISP, cnf, chd, cst)

%% 0.Selection of OB/CB
if( strcmp(chd.meas_mode, 'OPEN_BURST') )
    Npb = chd.N_pulses_burst;
    
elseif( strcmp(chd.meas_mode, 'CLOSED_BURST') )
    Npb = chd.prf/chd.brf; % BRI / PRI
    
end 

%% Applied or GNSS
if( ISP.gnss )
    cor2 = ISP.cor2_gnss_lrm;
    h0 = ISP.h0_lrm_gnss;
    
else
    cor2 = ISP.cor2_lrm;
    h0 = ISP.h0_lrm;
    
end 

%% Attitude corrections
z_cog_ant1_corr = - chd.x_cog_ant * sin(ATTITUDE.pitch_pre_dat) + chd.y_cog_ant * sin(ATTITUDE.roll_pre_dat) + ...
    chd.z_cog_ant * ( cos(ATTITUDE.pitch_pre_dat) + cos(ATTITUDE.roll_pre_dat) );

z_cog_ant2_corr = -chd.x_cog_ant_2 * sin(ATTITUDE.pitch_pre_dat) - chd.y_cog_ant_2 * sin(ATTITUDE.roll_pre_dat) + ...
    chd.z_cog_ant * ( cos(ATTITUDE.pitch_pre_dat) + cos(ATTITUDE.roll_pre_dat) );

delta_CoGant_1_ku1_ka = 2 * z_cog_ant1_corr / cst.c;
delta_CoGant_2_ku2 = ( (z_cog_ant1_corr + z_cog_ant2_corr )/2 ) / cst.c;

% Band selection
if( strcmp(ISP.band, 'Ku') )
    delta_CoGant = delta_CoGant_1_ku1_ka;
    attenuator_table_delay = ISP.attenuator_table_delay_ku1;
    attcode_lrm = ISP.attcode_lrm_ku1;
    ext_delay_ground = chd.ext_delay_ground_ku1;
    wv_length = cst.c / chd.freq_ku;
    bw = chd.bw_ku;

elseif( strcmp(ISP.band, 'Ku2') )
    delta_CoGant = delta_CoGant_2_ku2;
    attenuator_table_delay = ISP.attenuator_table_delay_ku2;
    attcode_lrm = ISP.attcode_lrm_ku2;
    ext_delay_ground = chd.ext_delay_ground_ku2;
    wv_length = cst.c / chd.freq_ku;
    bw = chd.bw_ku;
    
elseif( strcmp(ISP.band, 'Ka') )
    delta_CoGant = delta_CoGant_1_ku1_ka;
    attenuator_table_delay = ISP.attenuator_table_delay_ka;
    attcode_lrm = ISP.attcode_lrm_ka;
    ext_delay_ground = chd.ext_delay_ground_ka;
    wv_length = cst.c / chd.freq_ka;
    bw = chd.bw_ka;
    
end

%% Computation of H0

for i_burst = 0:chd.N_bursts_rc-1
    for i_pulse = 0:chd.N_pulses_burst-1
        
        cor2n(i_burst+1,i_pulse+1) = (i_pulse + i_burst*Npb) * cor2;
        
        if cor2n(i_burst+1, i_pulse+1) > 0
            cor2_inc(i_burst+1, i_pulse+1) = floor(cor2n(i_burst+1, i_pulse+1)/ (Npb*chd.N_bursts_rc)  );
        else
            cor2_inc(i_burst+1, i_pulse+1) = ceil(cor2n(i_burst+1, i_pulse+1)/ (Npb*chd.N_bursts_rc) );
        end        
    end
end

h_lrm = h0* chd.h0_cor2_unit_conv + cor2_inc;
hn_mean_lrm = mean( h_lrm(:) );

%% Computation of RX delay and CoG to antenna correction
% T0_lrm_pre_dat: real altimeter period calling the USO_Selection
rx_delay_lrm = hn_mean_lrm * T0_lrm_pre_dat / (chd.T0_h0_unit_conv * chd.h0_cor2_unit_conv)...
    + delta_CoGant;

rx_delay_lrm_nom = hn_mean_lrm * T0_nom / (chd.T0_h0_unit_conv * chd.h0_cor2_unit_conv)...
    + delta_CoGant;

%% Doppler correction (not applied to window delay, only for monitoring)
if (cnf.flag_doppler_correction)
    % alt_rate: computed in Orbit_selection using real altimeter period
    doppler_corr_lrm = cnf.doppler_range_correction_sign * (-2 * ORBIT.alt_rate /...
        (wv_length * bw / chd.pulse_length) );
else
    doppler_corr_lrm = 0;
end

%% Instrument corrections
int_delay_att_lrm= attenuator_table_delay(1 + attcode_lrm);
instr_delay = CAL1.int_delay_cor + ext_delay_ground;

%% Final window delay
ref_sample = ISP.fs - chd.central_sample_matched_filter + ...
    cnf.tracker_range_L1B_reference_sample;

ref_sample_shift = ref_sample * T0_lrm_pre_dat;
ref_sample_shift_nom = ref_sample * chd.T0_nom;

win_delay_lrm = rx_delay_lrm + instr_delay + doppler_corr_lrm +...
    ref_sample_shift - int_delay_att_lrm + chd.residual_delay;
win_delay_lrm_nom = rx_delay_lrm_nom + instr_delay + doppler_corr_lrm +...
    ref_sample_shift_nom - int_delay_att_lrm + chd.residual_delay;

% The reference pulse is the first Ku/Ka pulse, which is the first pulse in the burst.

%% Outputs
WIN_DELAY_OUT.win_delay_lrm = win_delay_lrm;
WIN_DELAY_OUT.win_delay_lrm_nom = win_delay_lrm_nom;
WIN_DELAY_OUT.band = ISP.band;
WIN_DELAY_OUT.chd.meas_mode = chd.meas_mode;
WIN_DELAY_OUT.ISP.gnss = ISP.gnss;

end
