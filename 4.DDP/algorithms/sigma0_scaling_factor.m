%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% 9.11 WAVEFORMS SCALING FACTOR %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of the waveforms scaling factor is to provide the L2
% processing with the waveform scaling factor in order to compute the
% backscatter coefficient of the surface from which the echo is reflected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v1.1 added S3 case. Sigma0 computed from the one comming from L1A. See
% "adaptnetCDF2internal.m"


%____________________________For each stack________________________________

function [L1BS , L1B]   = sigma0_scaling_factor (L1A, L1BS, L1B, cnf, chd, cst)




 

L1BS.wfm_scaling_factor_beam = zeros (1,chd.N_max_beams_stack);
azimuth_distance = zeros (1,chd.N_max_beams_stack);
range_distance = zeros (1,chd.N_max_beams_stack);
surface_area = zeros (1,chd.N_max_beams_stack);
L1B.wfm_scaling_factor = zeros (1,1);
    
    % 1. Preliminary factors computation
    %*********************************************************
    %**** Moved to Instrument Gain Corrections algorithm *****
    %*********************************************************
%     %   > Instrument attenuation averaging
%     aux = 0;
%     for 1:L1BS.N_beams_stack = 1:L1BS.N_beams_stack
%         aux = aux + att_corr_instr_sar(L1BS.burst_index1:L1BS.N_beams_stack));
%     end
%     att_corr_instr_sar_mean = aux / L1BS.N_beams_stack;
    %*********************************************************
    %*********************************************************
    %*********************************************************
    
    %   > Surface area (of every beam)
    % Pre-allocating memory
        norm_vel_sat(1:L1BS.N_beams_stack) = sqrt(L1BS.x_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.y_vel_sat_beam(1:L1BS.N_beams_stack).^2+ L1BS.z_vel_sat_beam(1:L1BS.N_beams_stack).^2);
        azimuth_distance (1:L1BS.N_beams_stack) = (1+L1BS.range_sat_surf(1:L1BS.N_beams_stack)./cst.earth_radius)*chd.wv_length .* L1BS.range_sat_surf(1:L1BS.N_beams_stack) ./ L1BS.pri_sar_sat_beam(1:L1BS.N_beams_stack) ./ (2 * (norm_vel_sat) .* chd.N_pulses_burst);
%         range_distance (1:L1BS.N_beams_stack) = 2 * sqrt (cst.c * L1BS.range_sat_surf(1:L1BS.N_beams_stack) * (1/(chd.pulse_length * chd.chirp_slope_chd)).*...
%                                     cst.earth_radius./(cst.earth_radius + L1BS.range_sat_surf(1:L1BS.N_beams_stack)) );
        range_distance (1:L1BS.N_beams_stack) = 2 * sqrt (cst.c * L1BS.range_sat_surf(1:L1BS.N_beams_stack) * (chd.PTR_width).*...
                                    cst.earth_radius./(cst.earth_radius + L1BS.range_sat_surf(1:L1BS.N_beams_stack)) );
        if cnf.hamming_window
            wf=1.486*0.92;%widening factor as denoted by Dinardo
        else
            wf=1.0;
        end
        %area_0.886*wf
        surface_area(1:L1BS.N_beams_stack) = (wf.*azimuth_distance(1:L1BS.N_beams_stack) .* range_distance(1:L1BS.N_beams_stack)).*0.886;
    
    surface_area_LRM= cst.pi*(0.5* min(range_distance(find(range_distance(:)))))^2;
    
    % 3. Compute the waveform scaling factor with external contributions for each beam
    aux = 0;
%         L1BS.wfm_scaling_factor_beam (1:L1BS.N_beams_stack) = 10*log10(64) + 30*log10(cst.pi)...
%             + 40*log10(L1BS.range_sat_surf(1:L1BS.N_beams_stack)) - 10*log10(chd.power_tx_ant) - 2*(chd.antenna_gain)...
%             - 20*log10(chd.wv_length) - 10*log10(mean(surface_area(1:L1BS.N_beams_stack)))- 10*log10(64);
        
        %added by EM 30.11.2015: consider the surface and not the mean over
        %the different beams & compensate for norm. in fft range
        %compression & TBP or pulse compression gain
        L1BS.wfm_scaling_factor_beam (1:L1BS.N_beams_stack) = 10*log10(64) + 30*log10(cst.pi)...
            + 40*log10(L1BS.range_sat_surf(1:L1BS.N_beams_stack)) - 10*log10(chd.power_tx_ant) - 2*(chd.antenna_gain)...
            - 20*log10(chd.wv_length) - 10*log10(surface_area(1:L1BS.N_beams_stack));
    
    
%         L1B.wfm_scaling_factor = 10*log10(64) + 30*log10(cst.pi)...
%             + 40*log10(mean(L1BS.range_sat_surf(:))) - 10*log10(chd.power_tx_ant) - 2*(chd.antenna_gain)...
%             - 20*log10(chd.wv_length) - 10*log10(min(surface_area(:)));
%         
%         aux = aux + L1BS.wfm_scaling_factor_beam(1:L1BS.N_beams_stack);
    
    
    
    % 4. Compute the averaged waveform scaling factor with external contributions
    L1B.wfm_scaling_factor = mean(L1BS.wfm_scaling_factor_beam (1:L1BS.N_beams_stack));


    % 2. Compute scaling factor only with instrumental factors
%     L1B.wfm_scaling_factor_instr = att_corr_instr_sar_mean + 10*log10(4) + 10*log10(cst.pi)...
%         -(power_tx_ant_chd) - 2*(chd.antenna_gain) - 20*log10(chd.wv_length) - (tx_rx_gain_ground_chd);
    
    L1B.wfm_scaling_factor_instr = 10*log10(4) + 10*log10(cst.pi)...
        -10*log10(chd.power_tx_ant) - 2*(chd.antenna_gain) - 20*log10(chd.wv_length)-10*log10(64) ;
    

    
    
    
    
    
end


