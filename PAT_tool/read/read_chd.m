%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Read CHD xml file
%
% Calling: GPPICE (input_folder, aux_folder, output_folder, options)
% INPUTs:
%
%
% OUTPUTs:
%
%
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
%
%v1.1 added default values, NamberofSamples,  MeasurementMode,  chd.meas_mode
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chd = read_chd(CHD_file, cst)

fidXML = fopen(CHD_file,'r');
%Default
chd.N_samples_sar           = 256;  % SAR Samples
chd.N_samples_sar_rmc       = 128;  % SAR RMC Samples
chd.uso_freq_nom = 10e6;
chd.alt_freq_multiplier = 39.5;
chd.antenna_gain        = 0;
chd.tx1_sar             = 0;
chd.ext_delay_ground    = 0;
chd.int_delay_ground    = 0;
chd.power_tx_ant        = 1;
chd.simulation          = '';
% ASSUMED FROM ARESYS PM4 (updated if exists in the chd file)
chd.antenna_beamwidth_along_track        = 1.06*pi/180; 
chd.antenna_beamwidth_across_track       = 1.1992*pi/180; 

chd.N_cal_pulses_burst  = 0;
chd.N_c_pulses_burst    = 0;


while( ~feof(fidXML) )
    myLine = fgetl(fidXML);
    
    if(strfind(myLine,'<PulseLength'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.pulse_length    = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [seconds]
    elseif(strfind(myLine,'<Bandwidth'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.bw              = str2double(myLine(separator_start(1)+1:separator_end(2)-1));% [Hz]
    elseif(strfind(myLine,'<NumberofSamples'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.N_samples_sar   = str2double(myLine(separator_start(1)+1:separator_end(2)-1));% [Sampless]
    elseif(strfind(myLine,'<PulseEnergy'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.power_tx_ant    = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [Joules]
    elseif(strfind(myLine,'<PulseSamplingRate'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.alt_freq        = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % % [Hz]
        chd.T0_nom          = 1/str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % % [seconds]
        
    elseif(strfind(myLine,'<PulseStartFrequency'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.start_freq      = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [rad] should be needed as the phase departure
    elseif(strfind(myLine,'<PulseStartPhase'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.phase_departure = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [rad] should be needed as the phase departure
    elseif(strfind(myLine,'<AntennaGain'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.antenna_gain = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [rad] should be needed as the phase departure
    elseif(strfind(myLine,'<MeasurementMode'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.meas_mode       = (myLine(separator_start(1)+1:separator_end(2)-1)); % [Hz]
        switch chd.meas_mode
            case 'OPEN_BURST'
                %chd.N_ku_pulses_burst   = 65;
                %chd.N_bursts_cycle_sar  = 10;    % Bursts in a cycle
                % JPLZ: adjust to new chronogram (no holes)
                chd.N_ku_pulses_burst   = 64;
                chd.N_bursts_cycle_sar  = 12;    % Bursts in a cycle
                chd.SWST                = 2.6037973915238429e-05; % from Lisa email 22/01/2019
            case 'CLOSED_BURST'
                chd.N_ku_pulses_burst   = 64;
                chd.N_bursts_cycle_sar  = 4;    % Bursts in a cycle 
                chd.SWST                = 2.4745193075612004e-05; % from Lisa email 22/01/2019
        end
        
        
        
    elseif(strfind(myLine,'<AlongTrack3dBApertureOneWay'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.antenna_beamwidth_along_track             = str2double(myLine(separator_start(1)+1:separator_end(2)-1))*pi/180; % [rad]
    elseif(strfind(myLine,'<AcrossTrack3dBApertureOneWay'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.antenna_beamwidth_across_track             = str2double(myLine(separator_start(1)+1:separator_end(2)-1))*pi/180;% [rad]
    elseif(strfind(myLine,'<AcquisitionPRF'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.prf             = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [Hz]
        chd.pri_nom         = 1/chd.prf;    %[s]
    

    elseif(strfind(myLine,'<BurstRepetitionFrequency'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.brf             = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [Hz]
		chd.burst_duration_sar      = 1/str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [Hz]
     
    elseif(strfind(myLine,'<fc_hz'))
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.freq        = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [Hz]
        chd.wv_length   = cst.c/chd.freq;
        if(chd.freq > 23575000000)
            chd.band='Ka';
        elseif(chd.freq < 23575000000)
            chd.band='Ku';
        end
    elseif(strfind(myLine,'<NumberOfRxChains'))
        fgetl(fidXML);
        fgetl(fidXML);
        elseif(strfind(myLine,'<Value'))
        separator_start     = strfind(myLine,'>');
        separator_end       = strfind(myLine,'<');
        chd.num_rx_chains   = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % % [Hz]
    
        
        
      
    elseif(strfind(myLine,'Location of TxRx antenna'))
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.x_cog_ant       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.y_cog_ant       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.z_cog_ant       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
   elseif(strfind(myLine,'Location of RxOnly antenna'))
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.x_cog_ant_2       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.y_cog_ant_2       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.z_cog_ant_2       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
   % ---------------------------- LAND ICE coordinates ----------------------------
    elseif(strfind(myLine,'<Target>SnowIce'))
        chd.simulation = 'LI';
    elseif(strfind(myLine,'<TargetType>Distributed'))
        chd.simulation = 'LI';
    elseif(strfind(myLine,'<Zone3>'))
        chd.simulation = 'SI';
    elseif(strfind(myLine,'<LayerDescription>Snow'))
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.snow_height = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.snow_depth  = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        
    elseif(strfind(myLine,'<LayerDescription>Ice'))
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.ice_freeboard = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        myLine = fgetl(fidXML);
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.ice_thickness  = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        
        
    % ---------------------------- PT coordinates ----------------------------
    elseif(strfind(myLine,'<ECEF_X'))
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.x_trp       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        chd.simulation = 'PT';
   elseif(strfind(myLine,'<ECEF_Y'))
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.y_trp       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
   elseif(strfind(myLine,'<ECEF_Z'))
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.z_trp       = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        chd.simulation  = 'PT';
   % ---------------------------- OCEAN properties ----------------------------
   elseif(strfind(myLine,'<SSH'))
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.ssh         = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]
        if(~strcmp(chd.simulation,'SI'))
            chd.simulation  = 'OC';
        end
   elseif(strfind(myLine,'<SWH'))
        separator_start = strfind(myLine,'>');
        separator_end   = strfind(myLine,'<');
        chd.swh         = str2double(myLine(separator_start(1)+1:separator_end(2)-1)); % [m]     
   end
    
    
end

chd.chirp_slope =  chd.bw/chd.pulse_length;


if(isfield(chd,'N_ku_pulses_burst'))
    disp(chd.meas_mode)
else
    chd.N_ku_pulses_burst = 64;
    chd.N_bursts_cycle_sar  = 4;    % Bursts in a cycle 
    chd.SWST                = 2.4745193075612004e-05; % from Lisa email 22/01/2019
    disp('N_ku_pulses_burst force to 64');
    chd.meas_mode = 'CLOSED_BURST';
end
    
chd.N_pulses_burst      = chd.N_ku_pulses_burst;   % HR pulses in a burst

chd.burst_duration_sar  = chd.N_pulses_burst / chd.prf;
chd.N_pri     = 0;

chd.alt_freq_multiplier = chd.alt_freq /chd.T0_nom;
chd.pri_T0_unit_conv = 4;
chd.h0_cor2_unit_conv = 16;
chd.T0_h0_unit_conv = 64;
chd.cai_cor2_unit_conv = 4096;

chd.i_sample_start = 1;

chd.fai_shift_number = 1024;

chd.N_max_beams_stack = (chd.N_bursts_cycle_sar+3) .* chd.N_pulses_burst;

chd.RMC_start_sample =1;

chd.PTR_width = 2.819e-9; % Using same as in CR2

fclose(fidXML);


end