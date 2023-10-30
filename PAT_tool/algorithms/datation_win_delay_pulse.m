%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Computation of the timing and window delay per pulse
% pri and T0 per pulse are also provided as well as the ID (RAW, RMC) of the process
% per pulse
%
% Calling: 
% INPUTs:
%   L1A_buffer: Buffer of L1A structures of all the bursts
%   cnf:        Configuration parameters structure
%   chd:        Characterization parameters structure
%   cst:        Constant parameters structure
%
% OUTPUTs:
%   time:       Datation time per pulse for all the bursts (single vector)
%   win_delay:  Window delay per pulse for all the bursts (single vector)
%   pri:        Pulse repetition interval per pulse for all the bursts 
%   T0:         Sampling period per pulse for all the bursts
%   Process_ID_pulse: Idenftifier of the pulse type for all bursts
%
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%            Albert Garcia / isardSAT
%            Ferran Gibert / isardSAT
%
% ----------------------------------------------------------
% Version record:
% v 1.1 2019/06/14 - Enable possibility of forcing processing over desired and 
%                  coordinates (cnf.FFt.force_POI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,win_delay,pri,T0,Process_ID_pulse,L1A_buffer] = datation_win_delay_pulse(L1A_buffer,cnf,chd,cst)

N_bursts       = length([L1A_buffer(:).time]);
N_total_pulses = N_bursts*chd.N_pulses_burst;
time_dumm      = zeros(N_bursts,chd.N_pulses_burst);
win_delay_dumm = zeros(N_bursts,chd.N_pulses_burst);
pri_dumm       = zeros(N_bursts,chd.N_pulses_burst);
T0_dumm        = zeros(N_bursts,chd.N_pulses_burst);
Process_ID_dumm= zeros(N_bursts,chd.N_pulses_burst);

for i_burst=1:N_bursts
    
    %% WINDOW DELAY COMPUTATION
    if ~cnf.trp_flag  && ~cnf.FFt.force_POI
      % H0 and COR2 shall be used to compute window delay assuming h0 and
      % need to know wether H0 is for each burst and at first pulse
      % need to know COR2 is changing
      % units of COR2 and H0
              
    else
        %for time being and transponder case assume replicated values
        %within the burst
        win_delay_dumm(i_burst,:) = L1A_buffer(i_burst).win_delay*ones(1,chd.N_pulses_burst);             
    end            
    
    %% PRI
    pri_dumm(i_burst,:) = L1A_buffer(i_burst).pri_sar*ones(1,chd.N_pulses_burst);        
    
    %% T0
    T0_dumm(i_burst,:) = L1A_buffer(i_burst).T0_sar*ones(1,chd.N_pulses_burst);  
    
    %% DATATION COMPUTATION
    % assuming the time RX of the first pulse
    L1A_buffer(i_burst).time = L1A_buffer(i_burst).time_rx_1st;
    % compute time the 1st pulse hits surface and propagate it for the
    % different pulses
    time_dumm(i_burst,:) = L1A_buffer(i_burst).time_rx_1st -(L1A_buffer(i_burst).win_delay/2)+((0:chd.N_pulses_burst-1)).*L1A_buffer(i_burst).pri_sar; %

% %     %ADRIAN test S3 approach, setting 1st time value to 0
%       time_first_burst = L1A_buffer(1).time_rx_1st;
% %                  time_dumm(i_burst,:) = 0*time_first_burst + ... 
% %                                          (-(chd.N_pulses_burst-1)/2:1:(chd.N_pulses_burst-1)/2)*L1A_buffer(i_burst).pri_sar + ...
% %                                          1*L1A_buffer(i_burst).pri_sar + ...
% %                                          1/chd.brf*(i_burst - 1);  %       
% %                  time_dumm(i_burst,:) = 0*time_first_burst + ... 
% %                                          (-(chd.N_pulses_burst-1)/2:1:(chd.N_pulses_burst-1)/2)*L1A_buffer(i_burst).pri_sar + ...
% %                                          1*L1A_buffer(i_burst).pri_sar + ...
% %                                          1/chd.brf*(i_burst - 1);  % 
%                  time_dumm(i_burst,:) = 0*time_first_burst + ...
%                                          (-(chd.N_pulses_burst-1)/2:1:(chd.N_pulses_burst-1)/2)*L1A_buffer(i_burst).pri_sar + ...
%                                          1/chd.brf*(i_burst - 1);  %

    
    %% PROCESS ID: RAW OR RMC
    Process_ID_dumm(i_burst,:) = double(L1A_buffer(i_burst).inst_id_sar_isp)*ones(1,chd.N_pulses_burst); %-1 fill value, 0 RAW and 1: RMC 
    
end
%reshape as a single vector of all pulses
time      = reshape(time_dumm.',[1,N_total_pulses]);
win_delay = reshape(win_delay_dumm.',[1,N_total_pulses]);
pri       = reshape(pri_dumm.',[1,N_total_pulses]);
T0        = reshape(T0_dumm.',[1,N_total_pulses]);
Process_ID_pulse = reshape(Process_ID_dumm.',[1,N_total_pulses]);

end