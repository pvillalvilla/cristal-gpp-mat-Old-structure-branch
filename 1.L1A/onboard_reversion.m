% ONBOARD REVERSION ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
%
% ---------------------------------------------------------
% Objective: - reverse the onboard RMC for the SAR RMC data.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ONBRD_REV_OUT] = onboard_reversion(T0_pre_dat, WIN_DELAY_OUT, ORBIT,ISP,chd,cnf,cst)


%% Band selection
if( strcmp(ISP.band, 'Ku') )
    wfm_i = ISP.wfm_i_ku1;
    wfm_q = ISP.wfm_q_ku1;
      
    freq = chd.freq_ku;
    
elseif( strcmp(ISP.band, 'Ku2') )
    wfm_i = ISP.wfm_i_ku2;
    wfm_q = ISP.wfm_q_ku2;

    freq = chd.freq_ku;
    
elseif( strcmp(ISP.band, 'Ka') )
    wfm_i = ISP.wfm_i_ka;
    wfm_q = ISP.wfm_q_ka;
    
    freq = chd.freq_ka;   
    
end

if( strcmp(ISP.pid, 'RMC') )
    
    if( strcmp(ISP.band, 'Ku') )
        rmc_matrix_i_samples = chd.rmc_matrix_i_samples_ku;
        rmc_matrix_q_samples = chd.rmc_matrix_q_samples_ku;
        rmc_weights_i_samples = chd.rmc_weights_i_samples_ku;
        rmc_weights_q_samples = chd.rmc_weights_q_samples_ku;

    elseif( strcmp(ISP.band, 'Ka') )
        rmc_matrix_i_samples = chd.rmc_matrix_i_samples_ka;
        rmc_matrix_q_samples = chd.rmc_matrix_q_samples_ka;
        rmc_weights_i_samples = chd.rmc_weights_i_samples_ka;
        rmc_weights_q_samples = chd.rmc_weights_q_samples_ka;

    end
    %%
    wfm_1 = zeros(chd.N_pulses_burst,chd.N_samples_rmc_onboard);
    
    samples_pb_chd = zeros(chd.N_pulses_burst,1);
    samples_pb_chd(:) = 1:chd.N_pulses_burst;
    samples_rmc_onboard = 1:chd.N_samples_rmc_onboard;
    samples_sar = 1:chd.N_samples_sar;
    
    % 0 <= i_sample_start <= 127
    % N_samples_rmc = 128
    % N_samples_sar = 256
    % N_samples_rmc_onboard = 512
    
    %% 0. Group real and imag of I&Q SAR RMC waveforms from ISP
    wfm_rmc_iq = wfm_i + 1i*wfm_q;
    wfm_rmc_iq = reshape(wfm_rmc_iq,chd.N_pulses_burst,chd.N_samples_rmc);
    
    %% 1. Truncation reversion
    % assumed that first sample is i_sample_start = 0
    wfm_1(:,chd.i_sample_start+1:chd.i_sample_start + chd.N_samples_rmc) = wfm_rmc_iq;
    
    %% 2. FFT range dimension: frequency to time
    wfm_2 = fftshift(fft(wfm_1,chd.N_samples_rmc_onboard,2),2);
    
    %% 3. RMC correction
    rmc_matrix = rmc_matrix_i_samples + 1i*rmc_matrix_q_samples;
    rmc_matrix = reshape(rmc_matrix ,chd.N_pulses_burst,chd.N_samples_rmc_onboard);
    wfm_3 = wfm_2 .* conj(rmc_matrix);
    
    %% 4. IFFT azimuth dimension: doppler beams to pulses
    wfm_4 = ifft(fftshift(wfm_3,1),chd.N_pulses_burst,1);
    
    %% 5. Doppler centroid reversion
    delta_alt = ISP.vertical_speed.*(ISP.pri);
    delta_tau = (2.*delta_alt/cst.c);
    
    wfm_5 = wfm_4 .* exp( -1i*2*cst.pi * samples_pb_chd *...
        delta_tau * freq );
    
    %% 6. Radial speed correction reversion
    wfm_6 = wfm_5 .* exp( -1i*2*cst.pi * samples_pb_chd *...
        delta_tau / T0_pre_dat .*...
        ( samples_rmc_onboard - chd.N_samples_rmc_onboard/2 ) ./ chd.N_samples_rmc_onboard );
    
    %% 7. CAI/FAI reversion
    wfm_7 = wfm_6 .* exp( 1i*2*cst.pi*( ...
        ( (WIN_DELAY_OUT.cai_namb).' - WIN_DELAY_OUT.cai_namb(1) ) * chd.cai_cor2_unit_conv...
        + ( (WIN_DELAY_OUT.fai).' - WIN_DELAY_OUT.fai(1) ) / chd.h0_cor2_unit_conv/chd.T0_h0_unit_conv )...
        / T0_pre_dat.*...
        ( samples_rmc_onboard - chd.N_samples_rmc_onboard/2 ) ./ chd.N_samples_rmc_onboard );
    
    %% 8. Undo azimuth weighting
    rmc_weights = (rmc_weights_i_samples + 1i*rmc_weights_q_samples).';
    wfm_8 = wfm_7 ./ rmc_weights;
    
    %% 9. IFFT range dimension: time to frequency
    wfm_9 = circshift( fliplr(wfm_8) , 1 , 2 );
    wfm_9 = 1/chd.N_samples_rmc_onboard * wfm_9(:,1:2:chd.N_samples_rmc_onboard);

else % RAW
    wfm_9 = reshape( wfm_i + 1i*wfm_q, chd.N_pulses_burst , chd.N_samples_rmc_onboard );
    
end

%%
if cnf.flag_height_rate_application
    
    %% 10. CAI/FAI re-reversion
    wfm_10 = wfm_9 .* exp( 1i*2*cst.pi*( ...
        ( (WIN_DELAY_OUT.cai_namb).' - WIN_DELAY_OUT.cai_namb(1) ) * chd.cai_cor2_unit_conv...
        + ( (WIN_DELAY_OUT.fai).' - WIN_DELAY_OUT.fai(1)) / chd.h0_cor2_unit_conv/chd.T0_h0_unit_conv )...
        / T0_pre_dat...
        .* ( samples_sar - chd.N_samples_sar ) ./ chd.N_samples_sar );
    
    % Now pulses have same window delay
    
    %% 11. Altitude rate application
    wfm_reversed = wfm_10 .* exp( 1i*2*cst.pi * 2 * ORBIT.alt_rate_sar_sat * samples_pb_chd .*...
        ISP.pri / cst.c / T0_pre_dat .* ...
        ( samples_sar - chd.N_samples_sar/2) ./ chd.N_samples_sar );
    
    %% 12. Window delay update
    alt_rate_wd_corr = ORBIT.alt_rate_sar_sat * samples_pb_chd .* ...
        ISP.pri * 2 / cst.c;
    alt_rate_wd_corr_nom = ORBIT.alt_rate_sar_sat * samples_pb_chd .* ...
        ISP.pri * chd.T0_nom / T0_pre_dat * 2 / cst.c;
    
    win_delay_updt = WIN_DELAY_OUT.win_delay_ref + mean(alt_rate_wd_corr(:));
    win_delay_updt_nom = WIN_DELAY_OUT.win_delay_ref_nom + mean(alt_rate_wd_corr_nom(:));
    
    uso_drift = win_delay_updt - win_delay_updt_nom;
    
else
    wfm_reversed = wfm_9;
    uso_drift = WIN_DELAY_OUT.win_delay - WIN_DELAY_OUT.win_delay_nom;
    % window delays remain the same as output from window_delay.m
    
end



%% Output
ONBRD_REV_OUT.wfm_reversed = wfm_reversed;
ONBRD_REV_OUT.uso_drift = uso_drift;
ONBRD_REV_OUT.win_delay_updt = win_delay_updt;
ONBRD_REV_OUT.band = ISP.band;
ONBRD_REV_OUT.pid = ISP.pid;



end
