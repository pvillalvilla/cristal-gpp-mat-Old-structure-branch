% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L.
% --------------------------------------------------------
% ---------------------------------------------------------
% Objective: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AV_OUT] = waveform_averaging(L1A, chd, cnf)
% L1A: array L1A for each record

% which variables to average?
% alt_rate, x_vel_sat_sar (y and z), xyz, sigma0
% lat, lon, WDelay, thr1 (temperature), doppler

% N = 10 ;
% L1A = struct;
% for i = 1:N
%     L1A(i).N_bursts_rc = N;
%     L1A(i).alt_rate = rand(1);
%     L1A(i).win_delay = 1e-7*rand(1);
%     L1A(i).x_vel_sat = 1e3*rand(1);
%     L1A(i).y_vel_sat = 1e3*rand(1);
%     L1A(i).z_vel_sat = 1e3*rand(1);
%     L1A(i).s0_scaling_factor_lros = 20*rand(1);
%     
% end

alt_rate_av_rc = 0;
win_delay_av_rc = 0;
x_vel_sat_av_rc = 0;
y_vel_sat_av_rc = 0;
z_vel_sat_av_rc = 0;
sigma0_av_rc = 0;

for i_burst = 1:L1A(1).N_bursts_rc
    wfm_av_rc = wfm_av_rc + L1A(i_burst).wfm_av_burst;
    alt_rate_av_rc = alt_rate_av_rc + L1A(i_burst).alt_rate;
    win_delay_av_rc = win_delay_av_rc + L1A(i_burst).win_delay;
    x_vel_sat_av_rc = x_vel_sat_av_rc + L1A(i_burst).x_vel_sat;
    y_vel_sat_av_rc = y_vel_sat_av_rc + L1A(i_burst).y_vel_sat;
    z_vel_sat_av_rc = z_vel_sat_av_rc + L1A(i_burst).z_vel_sat;
    sigma0_av_rc = sigma0_av_rc + 10^(L1A(i_burst).s0_scaling_factor_lros / 10);
    % doppler
    % x,y,z
    % lat,lon
    % thr1 (temp)
    
end

AV_OUT.wfm_av_rc = wfm_av_rc / L1A(1).N_bursts_rc;
AV_OUT.alt_rate_av_rc = alt_rate_av_rc / L1A(1).N_bursts_rc;
AV_OUT.win_delay_av_rc = win_delay_av_rc / L1A(1).N_bursts_rc;
AV_OUT.x_vel_sat_av_rc = x_vel_sat_av_rc / L1A(1).N_bursts_rc;
AV_OUT.y_vel_sat_av_rc = y_vel_sat_av_rc / L1A(1).N_bursts_rc;
AV_OUT.z_vel_sat_av_rc = z_vel_sat_av_rc / L1A(1).N_bursts_rc;
AV_OUT.sigma0_av_rc = 10*log10( sigma0_av_rc / L1A(1).N_bursts_rc );
% doppler
% x,y,z
% lat,lon
% thr1 (temp)

     
end