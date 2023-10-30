% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: The frequency domain waveforms are converted
% into range bins (time domain) by FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RC_OUT] = range_compression(ORBIT,L1A, chd, cnf)
% The range compression is performed at burst level
% as the input is an L1A record

% Doppler correction applied to WD
doppler_corr = cnf.doppler_range_correction_sign * (-2 * ORBIT.alt_rate /...
        (chd.wv_length * chd.bw / chd.pulse_length) );
    
win_delay = L1A.win_delay + doppler_corr;

% Range frequency to range time
wfm_RC = fft(fftshift( L1A.wfm , 2 ), chd.N_samples* cnf.zp_fact_range, 2 );

% Burst averaging
wfm_av_burst = mean( wfm_RC , 1); % abs or .^2 ?

% Output
RC_OUT.wfm_RC = wfm_RC;
RC_OUT.wfm_av_burst = wfm_av_burst;
RC_OUT.win_delay = win_delay;
RC_OUT.band = L1A.band;
end