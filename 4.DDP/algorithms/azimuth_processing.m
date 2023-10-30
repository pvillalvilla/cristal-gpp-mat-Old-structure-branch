%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: The purpose of the azimuth processing is to steer the beams to the
% different surface locations.
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
%
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% (JP López-Zaragoza) 1.1, 2023/03/15, Updated approximate and exact method functions to adapt to the new beam steering method
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


function [L1A] = azimuth_processing (L1A,i_burst_focussed,N_bursts, cnf, chd, cst)

% Pre-allocating memory
wfm_cal_corrected_azw   = zeros (chd.N_pulses_burst, chd.N_samples_sar);
wfm_beams_focused       = zeros (chd.N_pulses_burst*cnf.zp_fact_azimut, chd.N_samples_sar);
L1A.beams_focused_shifted   = zeros (chd.N_pulses_burst*cnf.zp_fact_azimut, chd.N_samples_sar);

if(strcmp(cnf.processing_mode,'SIN'))
    wfm_cal_corrected_azw_2   = zeros (chd.N_pulses_burst, chd.N_samples_sar);
    wfm_beams_focused_2       = zeros (chd.N_pulses_burst*cnf.zp_fact_azimut, chd.N_samples_sar);
    L1A.beams_focused_shifted_2   = zeros (chd.N_pulses_burst*cnf.zp_fact_azimut, chd.N_samples_sar);
end

k0      = 2*cst.pi / chd.wv_length;

vel_sat = [L1A.x_vel_sat_sar, L1A.y_vel_sat_sar, L1A.z_vel_sat_sar];
vel_sat_norm = norm(vel_sat); clear vel_sat

%% 1. AZIMUTH WEIGHTING
if(exist('azimuth_weighting_filename_chd','var'))
    [start_angle_azw, N_weights_azw, weights_azw, angles_azw] = read_weighting(azimuth_weighting_filename_chd);
else
    weights_azw = ones(1,chd.N_pulses_burst);
end
%
% Missing code to interpolate the weighting to apply it properly
%
wfm_cal_corrected_azw(:,:) = L1A.wfm_cal_gain_corrected(:,:);
if(strcmp(cnf.processing_mode,'SIN'))
    wfm_cal_corrected_azw_2(:,:) = L1A.wfm_cal_gain_corrected_2(:,:);
end

%% 2. METHOD SELECTION

if(cnf.force_exact_method)
    method = 1; %Exact method
else
    method = 0; %Approximate method
end

if method == 1
    %% -2.1. EXACT METHOD
    
    noninf_beams=isfinite(L1A.beam_ang);
    noninfbeamindexes=find(noninf_beams);
    
    for i_beam = 1:(length(noninfbeamindexes))
        
        % Pre-azimuth FFT phase shift
        wfm_pulses_received_shifted = L1A.wfm_cal_gain_corrected .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(noninfbeamindexes(i_beam)))* L1A.pri_sar * (0:chd.N_pulses_burst-1)).'*ones(1,chd.N_samples_sar)) ;
        
        if(cnf.hamming_window)
            % Hamming
            wfm_pulses_received_shifted(:,:) = wfm_pulses_received_shifted(:,:) .* (cnf.hamm_win * ones(1,chd.N_samples_sar));
            
        end
        
        %test images
        %figure; imagesc(abs(fft(fftshift(fft(L1A.wfm_cal_gain_corrected),1).')).');
        %figure; imagesc(abs(fft(fftshift(fft(wfm_pulses_received_shifted),1).')).');
        %figure; imagesc(abs(fft(fftshift(fft(wfm_beams_received_shifted),1).')).');
        %[Max,idx_max_stack]=max(abs(fft(fftshift(fft(wfm_pulses_received_shifted),1).')));
        %[Min,idx_center_stack]=min(idx_max_stack)
        %figure;plot(idx_max_stack)
        
        % FFT
        wfm_beams_received_shifted = fft(wfm_pulses_received_shifted,chd.N_pulses_burst*cnf.zp_fact_azimut)/sqrt(chd.N_pulses_burst*cnf.zp_fact_azimut);
        wfm_beams_received_shifted = fftshift(wfm_beams_received_shifted,1);
        
        % Select the focused beam (the central one) and store it
        wfm_beams_focused (noninfbeamindexes(i_beam), :) = wfm_beams_received_shifted(L1A.beam_ang_index,:);
        
        if(strcmp(cnf.processing_mode,'SIN'))
            
            wfm_pulses_received_shifted_2 = L1A.wfm_cal_gain_corrected_2 .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(noninfbeamindexes(i_beam)))* L1A.pri_sar * (0:chd.N_pulses_burst-1)).'*ones(1,chd.N_samples_sar)) ;
            
            if(cnf.hamming_window)
                % Hamming antenna 2
                wfm_pulses_received_shifted_2 = wfm_pulses_received_shifted_2 .* (cnf.hamm_win * ones(1,chd.N_samples_sar));
            end
            % FFT
            wfm_beams_received_shifted_2 = fft(wfm_pulses_received_shifted_2,chd.N_pulses_burst*cnf.zp_fact_azimut)/sqrt(chd.N_pulses_burst*cnf.zp_fact_azimut);
            wfm_beams_received_shifted_2 = fftshift(wfm_beams_received_shifted_2,1);
            % Select the focused beam (the central one) and store it
            wfm_beams_focused_2 (noninfbeamindexes(i_beam), :) = wfm_beams_received_shifted_2(L1A.beam_ang_index,:);
            
        end
        
    end
    L1A.beams_focused_shifted(:,:) = wfm_beams_focused(:,:);
    %L1A.beams_focused_shifted(:,:) = fftshift(wfm_beams_focused(:,:),1);
    if(strcmp(cnf.processing_mode,'SIN'))
        L1A.beams_focused_shifted_2(:,:) = wfm_beams_focused_2(:,:);
    end
    
elseif method == 0
    %% -2.2. APPROXIMATE METHOD
    
    while(L1A.beam_ang( L1A.beam_ang_index)==Inf) %JPLZ testing only inf condition
        if cnf.shift_focusing && i_burst_focussed<N_bursts/2
            L1A.beam_ang_index= L1A.beam_ang_index+1; %JPLZ: testing new beam focusing at first bursts, we now look for the nonINF values to the right
        else
            L1A.beam_ang_index= L1A.beam_ang_index-1;
        end
    end
    
    %     noninf_beams=isfinite(L1A.beam_ang);
    %     noninfbeamindexes=find(noninf_beams);
    %     middlebeamindex=ceil(mean(noninfbeamindexes));
    %     L1A.beam_ang_index=middlebeamindex;
    
    %L1A.beam_ang_index=33;
    wfm_pulses_received_shifted = L1A.wfm_cal_gain_corrected .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(L1A.beam_ang_index))* L1A.pri_sar * (0:chd.N_pulses_burst-1)).'*ones(1,chd.N_samples_sar)) ;
    %test images
    %figure; imagesc(abs(fft(fftshift(fft(L1A.wfm_cal_gain_corrected),1).')).');
    %figure; imagesc(abs(fft(fftshift(fft(wfm_pulses_received_shifted),1).')).');
    %[Max,idx_max_stack]=max(abs(fft(fftshift(fft(wfm_pulses_received_shifted),1).')));
    %[Min,idx_center_stack]=min(idx_max_stack)
    %figure;plot(idx_max_stack)
    
    % Pre-allocating memory
    if(cnf.hamming_window)
        % Hamming
        wfm_pulses_received_shifted(:,:) = wfm_pulses_received_shifted(:,:) .* (cnf.hamm_win * ones(1,chd.N_samples_sar));
    end
    % FFT process (!!!!!DIVIDE BY EFFECTIVE NUMBER OF PULSES!!!)
    % right now_burst 2 each RC skipped, burst1 only 60 effectNpb, burst3(here 2 cause B#2 is skipped) 63 effectNpb
    if (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst_focussed==1||any(i_burst_focussed==1+11*int16(1:59)))
       effect_Npb=60;  
        %      effect_Npb=chd.N_pulses_burst; 
    elseif (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst_focussed==2||any(i_burst_focussed==2+11*int16(1:59)))
       effect_Npb=63; 
       %       effect_Npb=chd.N_pulses_burst; 
    else
       effect_Npb=chd.N_pulses_burst;    
    end
    
    wfm_beams_focused = fft(wfm_pulses_received_shifted,chd.N_pulses_burst*cnf.zp_fact_azimut)/sqrt(effect_Npb*cnf.zp_fact_azimut);
    L1A.beams_focused_shifted(:,:) = fftshift(wfm_beams_focused,1);
    %don't need to do any swap
    
    if(strcmp(cnf.processing_mode,'SIN'))
        wfm_pulses_received_shifted_2  = L1A.wfm_cal_gain_corrected_2 .* (exp(-2i*k0*vel_sat_norm*cos(L1A.beam_ang(L1A.beam_ang_index))* L1A.pri_sar * (0:chd.N_pulses_burst-1)).'*ones(1,chd.N_samples_sar));
        
        if(cnf.hamming_window)
            % Hamming antenna 2
            wfm_pulses_received_shifted_2 = wfm_pulses_received_shifted_2 .* (cnf.hamm_win * ones(1,chd.N_samples_sar));
        end
        
        % FFT process
        wfm_beams_focused_2 = fft(wfm_pulses_received_shifted_2,chd.N_pulses_burst*cnf.zp_fact_azimut)/sqrt(effect_Npb*cnf.zp_fact_azimut);
        L1A.beams_focused_shifted_2(:,:) = (fftshift(wfm_beams_focused_2,1));
        
        %testing plots
        %figure;imagesc(abs(fft(fftshift(wfm_beams_focused.'))).');
        %figure;imagesc(abs(fft(wfm_pulses_received_shifted.')).');
    end
    
end

end












