function ISP = read_isp(filename)

fid = fopen(filename,'r','b'); % open the input file for reading
dirf=dir(filename); filesize=dirf.bytes;
frewind(fid); % go to beg of file

tic

i_rec=0;

while ftell(fid)<filesize
    
    i_rec=i_rec+1;
    
    if mod(i_rec,1000)==0
        disp(['Reading record ' num2str(i_rec) '  -------  Timing: ' datestr(now)]);
    end
    
    %% read ISP Packet Header
    % Packet ID
    ISP(i_rec).packet_header.version_number = fread(fid,1,'ubit3'); % byte 1
    ISP(i_rec).packet_header.type = fread(fid,1,'ubit1'); % byte 1
    ISP(i_rec).packet_header.data_field_header_flag = fread(fid,1,'ubit1'); % byte 1
    ISP(i_rec).packet_header.mbs_process_id = fread(fid,1,'ubit3'); % byte 1
    ISP(i_rec).packet_header.lsb_process_id = fread(fid,1,'uint8'); % byte 2
    % Packet Sequence Control
    ISP(i_rec).packet_header.grouping_flag = fread(fid,1,'ubit2'); % byte 3
    ISP(i_rec).packet_header.sequence_count = fread(fid,1,'ubit14'); % byte 3-4
    % Packet Length
    ISP(i_rec).packet_header.packet_length = fread(fid,1,'uint16'); % bytes 5-6
    
    %% read ISP Data Field Header
    % Data Field Header
    ISP(i_rec).data_field_header.pus_version_number = fread(fid,1,'ubit3'); % byte 7
    fread(fid,1,'ubit4'); % spare bits, byte 7
    ISP(i_rec).data_field_header.service_type = fread(fid,1,'uint8'); % byte 8
    ISP(i_rec).data_field_header.message_subtype = fread(fid,1,'uint8'); % byte 9
    ISP(i_rec).data_field_header.message_type_count = fread(fid,1,'uint16'); % byte 10-11 Message Type Counter
    ISP(i_rec).data_field_header.destination_identifier = fread(fid,1,'uint16'); % byte 12-13 Destination identifier
    ISP(i_rec).data_field_header.coarse_time = fread(fid,1,'uint32'); % bytes 14-17
    ISP(i_rec).data_field_header.fine_time = fread(fid,1,'ubit24'); % bytes 18-20
    
    switch ISP(i_rec).packet_header.lsb_process_id
        
        % ------------------------------ TM_CAL1_LRM -------------------------------------
        case hex2dec('B2')
            if i_rec == 1
                disp('Reading TM_CAL1_LRM');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');   % Byte 21 (CAL1 Configuration) -- bit B0..B3 reserved
            ISP(i_rec).radar_pattern = fread(fid, 1, 'ubit1');   %#ok<*AGROW> % Byte 21 (CAL1 Configuration) -- bit B04 Radar pattern 0: CB, 1: OB
            ISP(i_rec).instr_mode = fread(fid, 1, 'ubit3');   % Byte 21 (CAL1 Configuration) --  Instrument mode 0b000 : CAL1 LRM mode
            fread(fid, 1, 'ubit3');   % Byte 22 (CAL1 Configuration) -- bit B0..B2 reserved
            ISP(i_rec).cal_path = fread(fid, 1, 'ubit3');   % Byte 22 (CAL1 Configuration) -- bit B3..B5 Calibration path 0b000 : CalTx, 0b001 : CalRx, 0b010 : CalComp, 0b011 : CalDiff1, 0b100 : CalDiff2
            ISP(i_rec).band = fread(fid, 1, 'ubit2');        % Byte 22 -- CAL1 Band B6B7: Band
            ISP(i_rec).c_rdb_cal1_ncal = fread(fid, 1, 'uint16');  % Byte 23-24 -- Number of calibrations C_RDB_CAL1_NCAL
            ISP(i_rec).calibration_counter = fread(fid, 1, 'uint16');  % Byte 25-26 -- Calibration counter (from 1 to C_RDB_CAL1_NCAL)
            ISP(i_rec).ncycle = fread(fid, 1, 'uint16');               % Byte 27-28 -- NCYCLE
            ISP(i_rec).tracking_cyc_counter = fread(fid, 1, 'uint16'); % Byte 29-30 -- Tracking cycle counter (from 1 to NCYCLE)
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');            % Byte 31-32 -- NIMP (unsigned integer, LSB = 1)
            ISP(i_rec).pri = fread(fid, 1, 'uint16');            % Byte 33-34 -- PRI (LSB = 4 / 600 MHz ? 6.67 ns)
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 35-38 -- Distance command H0 applied (LSB = (1 / 600 MHz / 64) ?26.04 ps)
            ISP(i_rec).att_attenuation = fread(fid, 1, 'uint8');       % Byte 39 -- Applied ATT attenuation for the Ku or Ka band
            ISP(i_rec).shift_t = fread(fid, 1, 'uint8');       % Byte 40 -- ADERx band shift of the tracking channel : SHIFT_T
            ISP(i_rec).gain_t = fread(fid, 1, 'uint16');      % Byte 41-42 -- ADERx band gain of the tracking channel : GAIN_T
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_t = fread(fid, 1, 'uint16');          % Byte 73-74 -- First sample for the distance selection of the tracking channel : FS_T
            ISP(i_rec).ns_t = fread(fid, 1, 'uint16');          % Byte 75-76 -- Number of samples for the distance selection of the tracking channel : NS_T
            for samples=1:256                            % Byte 77-588 -- I2Q2 sample
                ISP(i_rec).i2q2(samples) = fread(fid, 1, 'uint16');
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 589-590 -- PEC
            %                     for count=1:2
            %                         disp(['Pointer in file is at byte ', num2str(ftell(fid))]);
            %                         ISP = read_isp_cal1_lrm_sarcb(fid, i_rec);
            %                     end
            
            % ------------------------------- TM_CAL1_RAW -------------------------------------
        case hex2dec('B3')
            if i_rec == 1
                disp('Reading TM_CAL1_RAW');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');   % Byte 21 (CAL1 Configuration) -- bit B0..B3 reserved
            ISP(i_rec).radar_pattern = fread(fid, 1, 'ubit1');   % Byte 21 (CAL1 Configuration) -- bit B04 Radar pattern 0: CB, 1: OB
            ISP(i_rec).instr_mode = fread(fid, 1, 'ubit3');   % Byte 21 (CAL1 Configuration) --  Instrument mode 0b000 : CAL1 LRM mode, 0b010 : CAL1 SAR ATT mode
            fread(fid, 1, 'ubit3');   % Byte 22 (CAL1 Configuration) -- bit B0..B2 reserved
            ISP(i_rec).cal_path = fread(fid, 1, 'ubit3');   % Byte 22 (CAL1 Configuration) -- bit B3..B5 Calibration path 0b000 : CalTx, 0b001 : CalRx, 0b010 : CalComp, 0b011 : CalDiff1, 0b100 : CalDiff2
            ISP(i_rec).band = fread(fid, 1, 'ubit2');        % Byte 22 -- CAL1 Band B6B7: Band
            ISP(i_rec).num_of_claibrations = fread(fid, 1, 'uint16');  % Byte 23-24 -- Number of calibrations C_RDB_CAL1_NCAL
            ISP(i_rec).calibration_counter = fread(fid, 1, 'uint16');  % Byte 25-26 -- Calibration counter (from 1 to C_RDB_CAL1_NCAL)
            ISP(i_rec).ncycle = fread(fid, 1, 'uint16');               % Byte 27-28 -- NCYCLE
            ISP(i_rec).tracking_cyc_counter = fread(fid, 1, 'uint16'); % Byte 29-30 -- Tracking cycle counter (from 1 to NCYCLE)
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');            % Byte 31-32 -- NIMP (unsigned integer, LSB = 1)
            ISP(i_rec).pri = fread(fid, 1, 'uint16');            % Byte 33-34 -- PRI (LSB = 4 / 600 MHz ? 6.67 ns)
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 35-38 -- Distance command H0 applied (LSB = (1 / 600 MHz / 64) ?26.04 ps)
            ISP(i_rec).att_attenuation = fread(fid, 1, 'uint8');       % Byte 39 -- Applied ATT attenuation for the Ku or Ka band
            ISP(i_rec).shift_t = fread(fid, 1, 'uint8');       % Byte 40 -- ADERx band shift of the tracking channel : SHIFT_T
            ISP(i_rec).gain_t = fread(fid, 1, 'uint16');      % Byte 41-42 -- ADERx band gain of the tracking channel : GAIN_T
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_s = fread(fid, 1, 'uint16');          % Byte 73-74 -- First sample for the distance selection of the science channel : FS_S (unsigned integer, LSB = 1)
            ISP(i_rec).ns_s = fread(fid, 1, 'uint16');          % Byte 75-76 -- Number of samples for the distance selection of the science channel : NS_S (unsigned integer, LSB = 1)
            fseek(fid,1,'cof');                                 % Byte 77 -- Reserved
            ISP(i_rec).burst_num = fread(fid, 1, 'uint8');      % Byte 78 -- Burst number
            for pulses=1:64                                     % Byte 79-32846 -- I2Q2 sample
                for samples=1:256
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(pulses, samples) = complex(i, q);
                end
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 32847 to 32848 -- PEC
            
            % ------------------------------- TM_CAL1_RMC -------------------------------------
        case hex2dec('B4')
            if i_rec == 1
                disp('Reading TM_CAL1_RMC');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit5');   % Byte 21 (CAL1 Configuration) -- bit B0..B4 reserved
            ISP(i_rec).instr_mode = fread(fid, 1, 'ubit3');   % Byte 21 (CAL1 Configuration) --  Instrument mode 0b011 : CAL1 RMC mode
            fread(fid, 1, 'ubit3');   % Byte 22 (CAL1 Configuration) -- bit B0..B2 reserved
            ISP(i_rec).cal_path = fread(fid, 1, 'ubit3');   % Byte 22 (CAL1 Configuration) -- bit B3..B5 Calibration path 0b000 : CalTx, 0b001 : CalRx, 0b010 : CalComp, 0b011 : CalDiff1, 0b100 : CalDiff2
            ISP(i_rec).band = fread(fid, 1, 'ubit2');        % Byte 22 -- CAL1 Band B6B7: Band = 0b00 : Ku1 band, 0b01 : Ku2 band, 0b10 : Ka band
            ISP(i_rec).ncycle = fread(fid, 1, 'uint16');               % Byte 23-24 -- NCYCLE
            ISP(i_rec).tracking_cyc_counter = fread(fid, 1, 'uint16'); % Byte 25-26 -- Counter of tracking cycle (from 1 to NCYCLE)
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');            % Byte 27-28 -- NIMP (unsigned integer, LSB = 1)
            ISP(i_rec).pri = fread(fid, 1, 'uint16');            % Byte 29-30 -- PRI (LSB = 4 / 600 MHz ? 6.67 ns)
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 applied (LSB = (1 / 600 MHz / 64) ?26.04 ps)
            ISP(i_rec).att_attenuation = fread(fid, 1, 'uint8'); % Byte 35 -- Applied ATT attenuation for the Ku or Ka band
            fseek(fid,1,'cof');                             % Byte 36 -- Reserved
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_r = fread(fid, 1, 'uint16');          % Byte 67-68 -- First sample for the distance selection of the RMC channel : FS_R (unsigned integer, LSB = 1)
            ISP(i_rec).ns_r = fread(fid, 1, 'uint16');          % Byte 69-70 -- Number of samples for the distance selection of the RMC channel : NS_R (unsigned integer, LSB = 1)
            fseek(fid,1,'cof');                          % Byte 71 --  Reserved
            ISP(i_rec).burst_num = fread(fid, 1, 'uint8');      % Byte 72 -- Burst number
            fseek(fid,1,'cof');                          % Byte 73-74 --  Reserved
            ISP(i_rec).shift_r = fread(fid, 1, 'uint8');       % Byte 73-74 -- SHIFT_R band shift computed by ADERx and transmitted directly from ADERx to TWIST FPGA
            for pulses=1:64                              % Byte 75-16458 -- I2Q2 sample
                for samples=1:128
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(pulses, samples) = complex(i, q);
                end
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 16458 to 16460 -- PEC
            
            %----------------------------------- TM_CAL1_INSTR ---------------------------------------
        case hex2dec('B5')
            if i_rec == 1
                disp('Reading TM_CAL1_INSTR');
            end
            for n = 1:2
                if n == 2
                    disp(['Pointer in file is at byte ', num2str(ftell(fid))]);
                    fseek(fid,6,'cof');    %primary packet
                    fseek(fid,14,'cof');   %secondary packet
                end
                
                % Record --------------------------------------------------------------
                ISP(i_rec).configuration = fread(fid, 1, 'uint8');   % Byte 21 -- CAL1 Configuration
                ISP(i_rec).band = fread(fid, 1, 'uint8');        % Byte 22 -- CAL1 INSTR Configuration - B6B7: Band for INSTR mode
                ISP(i_rec).ncycle = fread(fid, 1, 'uint16');               % Byte 23-24 -- NCYCLE
                ISP(i_rec).num_tracking_cyc = fread(fid, 1, 'uint16'); % Byte 25-26 -- Number of tracking cycles performed = N
                ISP(i_rec).nimp = fread(fid, 1, 'uint16');            % Byte 27-28 -- NIMP (unsigned integer, LSB = 1)
                ISP(i_rec).pri = fread(fid, 1, 'uint16');            % Byte 29-30 -- PRI (LSB = 4 / 600 MHz ? 6.67 ns)
                ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 applied (LSB = (1 / 600 MHz / 64) ?26.04 ps)
                ISP(i_rec).att_attenuation = fread(fid, 1, 'uint8'); % Byte 35 -- Applied ATT attenuation for the Ku or Ka band
                ISP(i_rec).shift_i = fread(fid, 1, 'uint8');      % Byte 36 -- ADERx shift of the CAL1-INSTR channel : SHIFT_I
                ISP(i_rec).gain_i = fread(fid, 1, 'uint16');      % Byte 37-38 -- ADERx gain of the CAL1-INSTR channel : GAIN_I
                ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
                ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
                ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
                ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
                ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
                ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
                ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
                ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
                ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
                ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
                ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
                ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
                ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
                ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
                ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
                fseek(fid,1,'cof');
                ISP(i_rec).nTm = fread(fid, 1, 'uint8');            % Byte 70 -- TM number
                for samples=1:16384                                 % Byte 71-32838 -- I2Q2 sample
                    if ISP(i_rec).nTm == 2
                        samples = samples + 16384; %#ok<FXSET>
                    end
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(samples) = complex(i, q);
                end
                ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 32839 to 32840 -- PEC
            end
            
            %----------------------------------- TM_CAL2_RAW ---------------------------------------
        case hex2dec('B6')
            if i_rec == 1
                disp('Reading TM_CAL2_RAW');
            end
            % Record --------------------------------------------------------------
            fseek(fid,1,'cof');       % Byte 21 -- Reserved
            fread(fid, 1, 'ubit5');   % Byte 22 (CAL2 Configuration) -- bit B0..B4 reserved
            ISP(i_rec).radar_pattern = fread(fid, 1, 'ubit1');   % Byte 22 (CAL2 Configuration) -- bit B04 Radar pattern 0: CB, 1: OB
            ISP(i_rec).band = fread(fid, 1, 'ubit2');   % Byte 22 (CAL2 Configuration) --  Band
            ISP(i_rec).ncycle = fread(fid, 1, 'uint16');               % Byte 23-24 -- NCYCLE
            ISP(i_rec).num_tracking_cyc = fread(fid, 1, 'uint16'); % Byte 25-26 -- Number of the tracking cycle (from 1 to NCYCLE)
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');            % Byte 27-28 -- NIMP (unsigned integer, LSB = 1)
            ISP(i_rec).pri = fread(fid, 1, 'uint16');            % Byte 29-30 -- PRI (LSB = 4 / 600 MHz ? 6.67 ns)
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 applied (LSB = (1 / 600 MHz / 64) ?26.04 ps)
            ISP(i_rec).att_attenuation = fread(fid, 1, 'uint8'); % Byte 35 -- Applied ATT attenuation for the Ku or Ka band
            ISP(i_rec).shift_s = fread(fid, 1, 'uint8');      % Byte 36 -- ADERx band shift of the science channel : SHIFT_S
            ISP(i_rec).gain_s = fread(fid, 1, 'uint16');      % Byte 37-38 -- ADERx band gain of the science channel : GAIN_S
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_s = fread(fid, 1, 'uint16');          % Byte 69-70 -- First sample for the distance selection of the science channel : FS_S (unsigned integer, LSB = 1)
            ISP(i_rec).ns_s = fread(fid, 1, 'uint16');          % Byte 70-71 -- Number of samples for the distance selection of the science channel : NS_S (unsigned integer, LSB = 1)
            fseek(fid,1,'cof');
            ISP(i_rec).nBurst = fread(fid, 1, 'uint8');             % Byte 122 -- Burst number
            for pulses=1:64                              % Byte 75-32842 -- I2Q2 sample
                for samples=1:256
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(pulses, samples) = complex(i, q);
                end
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 32843 to 32844 -- PEC
            
            %------------------------------------ TM_ACQ ------------------------------
        case hex2dec('B7')
            if i_rec == 1
                disp('Reading TM_ACQ');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit5');   % Byte 21 (Tracking configuration) -- bit B0..B4 reserved
            ISP(i_rec).closed_loop_gain = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- B5 Closed loop gain = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) --  B6 Acquisition = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) --  B7 DEM MRAM read access enabled/disabled = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            fread(fid, 1, 'ubit6');   % Byte 22 (Acquisition phase) -- bit B0..B5 reserved
            ISP(i_rec).acq_phase = fread(fid, 1, 'ubit2');   % Byte 22 (Acquisition phase) -- Acquisition phase B6B7 = 0b00 : Noise estimation, 0b01 : Search phase, 0b10 : Positioning phase, 0b11 : Locking phase
            fseek(fid,1,'cof');         % Byte 23 -- Reserved
            ISP(i_rec).tracking_cyc_search_phase = fread(fid, 1, 'uint8'); % Byte 24 -- Tracking cycle counter for the search phase
            ISP(i_rec).tracking_cyc_position_phase = fread(fid, 1, 'uint8');  % Byte 25 -- Tracking cycle counter for the positioning phase
            ISP(i_rec).tracking_cyc_locking_phase = fread(fid, 1, 'uint8');   % Byte 26 -- Tracking cycle counter for the locking phase
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');            % Byte 27-28 -- NIMP (unsigned integer, LSB = 1)
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 29-30 -- PRI (LSB = 4 / 600 MHz ? 6.67 ns)
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');      % Byte 31-32  -- Ambiguity rank (unsigned integer, LSB = 1)
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');      % Byte 33-36 -- Distance command H0 applied on current Tracking Cycle
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');      % Byte 37-40 -- Distance command COR2 applied on current Tracking Cycle
            ISP(i_rec).threshold = fread(fid, 1, 'uint32');      % Byte 41-44 -- Detection threshold computed on Tracking Cycle (N-2)
            ISP(i_rec).nav_oop = fread(fid, 1, 'uint32');      % Byte 45-48 -- On Orbit Position NAV_OOP computed by navigation processing (LSB = 10-6 deg)
            ISP(i_rec).h0_nav = fread(fid, 1, 'uint32');      % Byte 49-52 -- Distance command H0_NAV computed by navigation processing (LSB = (1 / 600 MHz / 64) ?26.04 ps)
            ISP(i_rec).cor2_nav = fread(fid, 1, 'int32');      % Byte 53-56 -- Distance command COR2_NAV computed by navigation processing (signed integer) (LSB = (1 / 600 MHz / 1024) ?1.63 ps)
            ISP(i_rec).snr = fread(fid, 1, 'uint16');      % Byte 57-58  -- Estimated Signal to Noise Ratio computed on current Tracking Cycle (LSB = 1)
            ISP(i_rec).gates = fread(fid, 1, 'uint16');      % Byte 59-60  -- Estimated number of equivalent gates computed on current Tracking Cycle (LSB = 1)
            ISP(i_rec).att_attenuation = fread(fid, 1, 'uint8');       % Byte 61 -- Applied ATT attenuation for the Ku or Ka band
            ISP(i_rec).shift_t_ku = fread(fid, 1, 'uint8');       % Byte 62 -- ADERx Ku band shift of the tracking channel : SHIFT_T_KU
            ISP(i_rec).gain_t_ku = fread(fid, 1, 'uint16');      % Byte 63-64 -- ADERx Ku band gain of the tracking channel : GAIN_T_KU
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 65-68 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 69-71 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 72 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 73-76 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 77-80 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 81-84 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 85-88 -- Orbit pos
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_t = fread(fid, 1, 'uint16');          % Byte 119-120 -- First sample for the distance selection of the tracking channel : FS_T (unsigned integer, LSB = 1)
            ISP(i_rec).ns_t = fread(fid, 1, 'uint16');          % Byte 121-122 -- Number of samples for the distance selection of the tracking channel : NS_T (unsigned integer, LSB = 1)
            ISP(i_rec).ad = fread(fid, 1, 'uint16');            % Byte 123-124 -- Number of samples for the distance accumulation of the tracking channel : AD (unsigned integer, LSB = 1)
            for samples=1:256                           % Byte 125-636 -- I2Q2 sample
                ISP(i_rec).i2q2(samples) = fread(fid, 1, 'uint16');
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 637 to 638 -- PEC
            
            %------------------------------------ TM_ECHO_LRM_SAR_CB ------------------
        case hex2dec('B8')
            if i_rec == 1
                disp('Reading TM_ECHO_LRM_SAR_CB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');   % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank (unsigned integer, LSB = 1)
            ISP(i_rec).band = fread(fid, 1, 'uint8');           % Byte 29 -- Band
            fread(fid, 1, 'ubit7');                      % Byte 30 Reserved
            ISP(i_rec).loss_track = fread(fid, 1, 'ubit1');     % Byte 30 Loss of track criterion computed on Tracking Cycle (N-2) in OL mode B7 = 0 : normal, 1 : Loss of track
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
            ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
            ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
            ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
            ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
            ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
            ISP(i_rec).distance_err = fread(fid, 1, 'int32');   % Byte 55-58 -- Distance error
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
            ISP(i_rec).shift_t = fread(fid, 1, 'uint8');        % Byte 60 -- ADERx Ku or Ka band shift of the tracking channel : SHIFT_T
            ISP(i_rec).gain_t = fread(fid, 1, 'uint16');        % Byte 61-62 -- ADERx Ku or Ka band gain of the tracking channel : GAIN_T
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 63-66 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-86 -- Orbit pos
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_t = fread(fid, 1, 'uint16');            % Byte 117-118 -- First sample for the distance selection of the tracking channel : FS_T
            ISP(i_rec).ns_t = fread(fid, 1, 'uint16');            % Byte 119-120 -- Number of samples for the distance selection of the tracking channel : NS_T
            for samples=1:256                            % Byte 121-632 -- I2Q2 sample
                ISP(i_rec).i2q2(samples) = fread(fid, 1, 'int16');
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 633-634 -- PEC
            
            %----------------------------------- TM_ECHO_CAL_SAR_CB ----------------------
        case hex2dec('B9')
            if i_rec == 1
                disp('Reading TM_ECHO_CAL_SAR_CB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');   % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 27-30 -- Distance command H0 applied
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 31 -- Att
            ISP(i_rec).shift_s = fread(fid, 1, 'uint8');        % Byte 32 -- ADERx shift of the science channel for the calibration pulse : SHIFT_S
            ISP(i_rec).gain_s_cal = fread(fid, 1, 'uint16');        % Byte 33-34 -- ADERx gain of the science channel for the calibration pulse : GAIN_S_CAL
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_s = fread(fid, 1, 'uint16');            % Byte 65-66 -- First sample for the distance selection of the science channel : FS_S
            ISP(i_rec).ns_s = fread(fid, 1, 'uint16');            % Byte 67-68 -- Number of samples for the distance selection of the tracking channel : NS_T
            for samples=1:256                            % Byte 69-580 -- I2Q2 sample CALTx Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calTx_ku1(pulses, samples) = complex(i, q);
            end
            for samples=1:256                            % Byte 581-1092 -- I2Q2 sample CALComp Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calComp_ku1(pulses, samples) = complex(i, q);
            end
            for samples=1:256                            % Byte 1093-1604 -- I2Q2 sample CALRx Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calRx_ku1(pulses, samples) = complex(i, q);
            end
            for samples=1:256                            % Byte 1605-2116 -- I2Q2 sample  CALTx Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calTx_ka(pulses, samples) = complex(i, q);
            end
            for samples=1:256                            % Byte 2117-2628 -- I2Q2 sample  CALComp Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calComp_ka(pulses, samples) = complex(i, q);
            end
            for samples=1:256                            % Byte 2629-3140 -- I2Q2 sample  CALRx Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calRx_ka(pulses, samples) = complex(i, q);
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 3141-3142 -- PEC
            
            %----------------------------------- TM_ECHO_RMC_SAR_CB -----------------
        case hex2dec('BA')
            %                  if i_rec == 1
            disp('Reading TM_ECHO_RMC_SAR_CB');
            %                   end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');   % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank (unsigned integer, LSB = 1)
            ISP(i_rec).band = fread(fid, 1, 'uint8');           % Byte 29 -- Band
            fread(fid, 1, 'ubit7');                      % Byte 30 Reserved
            ISP(i_rec).loss_track = fread(fid, 1, 'ubit1');     % Byte 30 Loss of track criterion computed on Tracking Cycle (N-2) in OL mode B7 = 0 : normal, 1 : Loss of track
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
            ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
            ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
            ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
            ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
            ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
            ISP(i_rec).distance_err = fread(fid, 1, 'int32');   % Byte 55-58 -- Distance error
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
            fseek(fid,1,'cof');                          % Byte 60 -- Reserved
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 61-66 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-84 -- Orbit pos
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_r = fread(fid, 1, 'uint16');            % Byte 115-116 -- First sample for the distance selection of the RMC channel : FS_R
            ISP(i_rec).ns_r = fread(fid, 1, 'uint16');            % Byte 117-118 -- Number of samples for the distance selection of the RMC channel : NS_R
            fseek(fid,1,'cof');                            % Byte 119 -- Reserved
            ISP(i_rec).nBurst = fread(fid, 1, 'uint8');           % Byte 120 -- Burst number (from 1 to 4)
            fseek(fid,1,'cof');                            % Byte 121-122 -- Reserved
            ISP(i_rec).shift_r = fread(fid, 1, 'uint8');          % Byte 121-122 -- Band shift computed by ADERx: SHIFT_R
            for pulses = 1:64
                for samples=1:128                          % Byte 123-16506 -- I2Q2 sample
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(pulses, samples) = complex(i, q);
                end
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 16507 to 16508 -- PEC
            
            %----------------------------------- TM_ECHO_RAW_SAR_CB ------------------
        case hex2dec('BB')
            if i_rec == 1 || mod(i_rec, 100) == 0
                disp('Reading TM_ECHO_RAW_SAR_CB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');   % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank
            ISP(i_rec).band = fread(fid, 1, 'uint8');           % Byte 29 -- Band
            fread(fid, 1, 'ubit7');                      % Byte 30 Reserved
            ISP(i_rec).loss_track = fread(fid, 1, 'ubit1');     % Byte 30 Loss of track criterion computed on Tracking Cycle (N-2) in OL mode B7 = 0 : normal, 1 : Loss of track
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
            ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
            ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
            ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
            ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
            ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
            ISP(i_rec).distance_err = fread(fid, 1, 'int32');  % Byte 55-58 -- Distance error
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
            ISP(i_rec).shift_s = fread(fid, 1, 'uint8');        % Byte 60 -- ADERx Ku or Ka band shift of the science channel : SHIFT_S
            ISP(i_rec).gain_s = fread(fid, 1, 'uint16');          % Byte 61-62 -- ADERx Ku or Ka band gain of the science channel : GAIN_S
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 63-66 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-86 -- Orbit pos
            ISP(i_rec).thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_s = fread(fid, 1, 'uint16');          % Byte 117cat-118 -- FS_S
            ISP(i_rec).ns_s = fread(fid, 1, 'uint16');          % Byte 119-120 -- NS_S
            fseek(fid,1,'cof');                          % Byte 121 -- Reserved
            ISP(i_rec).nBurst = fread(fid, 1, 'uint8');         % Byte 122 -- Burst number
            for pulses = 1:64
                for samples=1:256                          % Byte 123-32890 -- I2Q2 sample
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(pulses, samples) = complex(i, q);
                end
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 32891-32892 -- PEC
            
            %--------------------------------------------- SARin ------------------------
            %-------------------------------------- TM_ECHO_LRM_SARin_CB ---------------
        case hex2dec('BC')
            if i_rec == 1
                disp('Reading TM_ECHO_LRM_SARin_CB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');                      % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank
            ISP(i_rec).band = fread(fid, 1, 'ubit1');           % Byte 29 -- Band
            fread(fid, 1, 'ubit7');                      % Byte 30 -- bit B0..B6 reserved
            ISP(i_rec).loss_track = fread(fid, 1, 'uint8');     % Byte 30 -- Loss of track B7 =  0 : normal, 1 : Loss of track
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
            ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
            ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
            ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
            ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
            ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
            ISP(i_rec).distance_err = fread(fid, 1, 'int32');  % Byte 55-58 -- Distance error
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
            ISP(i_rec).shift_t = fread(fid, 1, 'uint8');          % Byte 60 -- ADERx Ku or Ka band shift of the tracking channel : SHIFT_T
            ISP(i_rec).gain_t = fread(fid, 1, 'uint16');          % Byte 61-62 -- ADERx Ku or Ka band gain of the tracking channel : GAIN_T
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 63-66 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-86 -- Orbit pos
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs = fread(fid, 1, 'uint16');            % Byte 117-118 -- FS_T
            ISP(i_rec).ns = fread(fid, 1, 'uint16');            % Byte 119-120 -- NS_T
            for samples=1:2048                           % Byte 121-4216 -- I2Q2 sample
                ISP(i_rec).i2q2(samples) = fread(fid, 1, 'uint16');
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 4217-4218 -- PEC
            
            %-------------------------------------- TM_ECHO_CAL_SARin_CB ----------------
        case hex2dec('BD')
            if i_rec == 1
                disp('Reading TM_ECHO_CAL_SARin_CB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');                      % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 27-30 -- Distance command H0 applied
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 31 -- Att
            ISP(i_rec).shift_s = fread(fid, 1, 'uint8');        % Byte 32 -- ADERx shift of the science channel for the calibration pulse : SHIFT_S
            ISP(i_rec).gain_s_cal = fread(fid, 1, 'uint16');    % Byte 33-34 -- ADERx gain of the science channel for the calibration pulse : GAIN_S_CAL
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs_s = fread(fid, 1, 'uint16');            % Byte 65-66 -- FS_S
            ISP(i_rec).ns_s = fread(fid, 1, 'uint16');            % Byte 67-68 -- NS_S
            for samples=1:256                              % Byte 69-580 -- I2Q2 sample CALTx Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calTx_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 581-1092 -- I2Q2 sample CALComp Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calComp_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 1093-1604 -- I2Q2 sample CALRx Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calRx_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 1605-2116 -- I2Q2 sample CALdiff Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calDiff_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 2117-2628 -- I2Q2 sample CALdiff Ku2
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calDiff_ku2(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 2629-3140 -- I2Q2 sample CALTx Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calTx_ka(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 3141-3652 -- I2Q2 sample CALComp Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calComp_ka(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 3653-4164 -- I2Q2 sample CALRx Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calRx_ka(samples) = complex(i, q);
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 4165-4166 -- PEC
            
            %-------------------------------------- TM_ECHO_RAW_SARin_CB ----------------
        case hex2dec('BE')
            if i_rec == 1 || mod(i_rec, 100) == 0
                disp('Reading TM_ECHO_RAW_SARin_CB');
            end
            for n = 1:4
                if n > 1
                    
                    fseek(fid,6,'cof');    %primary packet
                    fseek(fid,14,'cof');   %secondary packet
                end
                % Record --------------------------------------------------------------
                fread(fid, 1, 'ubit4');                      % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
                ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
                ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
                ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
                ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
                ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
                ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
                ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
                ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank
                ISP(i_rec).band = fread(fid, 1, 'uint8');           % Byte 29 -- Band
                fread(fid, 1, 'ubit7');                      % Byte 30 Reserved
                ISP(i_rec).loss_track = fread(fid, 1, 'ubit1');     % Byte 30 Loss of track criterion computed on Tracking Cycle (N-2) in OL mode B7 = 0 : normal, 1 : Loss of track
                ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
                ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
                ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
                ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
                ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
                ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
                ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
                ISP(i_rec).distance_err = fread(fid, 1, 'int32');  % Byte 55-58 -- Distance error
                ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
                ISP(i_rec).shift = fread(fid, 1, 'uint8');          % Byte 60 -- Shift
                ISP(i_rec).gain = fread(fid, 1, 'uint16');          % Byte 61-62 -- Gain
                ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 63-66 -- Coarse Time
                fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
                fineTime2 = fread(fid, 1, 'uint8');
                binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
                ISP(i_rec).fine_time = bin2dec(binaryFineTime);
                ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
                ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
                ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
                ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
                ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-86 -- Orbit pos
                ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
                ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
                ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
                ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
                ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
                ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
                ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
                ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
                ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
                ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
                ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
                ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
                ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
                ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
                ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
                ISP(i_rec).fs = fread(fid, 1, 'uint16');            % Byte 117-118 -- FS_S
                ISP(i_rec).ns = fread(fid, 1, 'uint16');            % Byte 119-120 -- NS_S
                ISP(i_rec).nTm = fread(fid, 1, 'uint8');            % Byte 121 -- TM number
                ISP(i_rec).cal4_occurrence = fread(fid, 1, 'ubit1'); % Byte 122 -- B0 CAL4 occurrence = 0: No CAL4
                ISP(i_rec).nBurst = fread(fid, 1, 'ubit7');         % Byte 122 -- Burst number
                for p=1:16                                     % Byte 123-32890 -- I2Q2 sample
                    pulses = p + (n-1)*16;
                    for s=1:1024
                        samples = s + (n-1)*1024;
                        i = fread(fid, 1, 'int8');
                        q = fread(fid, 1, 'int8');
                        ISP(i_rec).iq(pulses, samples) = complex(i, q);
                    end
                end
                ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 32891-32892 -- PEC
            end
            
            %-------------------------------------- TM_ECHO_LRM_SARin_OB ----------------
        case hex2dec('BF')
            if i_rec == 1
                disp('Reading TM_ECHO_LRM_SARin_OB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');                      % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank
            ISP(i_rec).band = fread(fid, 1, 'uint8');           % Byte 29 -- Band
            fread(fid, 1, 'ubit7');                      % Byte 30 Reserved
            ISP(i_rec).loss_track = fread(fid, 1, 'ubit1');     % Byte 30 Loss of track criterion computed on Tracking Cycle (N-2) in OL mode B7 = 0 : normal, 1 : Loss of track
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
            ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
            ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
            ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
            ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
            ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
            ISP(i_rec).distance_err = fread(fid, 1, 'int32');  % Byte 55-58 -- Distance error
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
            ISP(i_rec).shift = fread(fid, 1, 'uint8');          % Byte 60 -- Shift
            ISP(i_rec).gain = fread(fid, 1, 'uint16');          % Byte 61-62 -- Gain
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 63-66 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-86 -- Orbit pos
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs = fread(fid, 1, 'uint16');            % Byte 117-118 -- FS_S
            ISP(i_rec).ns = fread(fid, 1, 'uint16');            % Byte 119-120 -- NS_S
            for samples=1:256                                    % Byte 121-632 -- I2Q2 sample
                ISP(i_rec).i2q2(samples) = fread(fid, 1, 'uint16');
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 633-634 -- PEC
            
            %-------------------------------------- TM_ECHO_CAL_SARin_OB ----------------
        case hex2dec('C0')
            if i_rec == 1
                disp('Reading TM_ECHO_CAL_SARin_OB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');                      % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 27-30 -- Distance command H0 applied
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 31 -- Att
            ISP(i_rec).shift_s = fread(fid, 1, 'uint8');        % Byte 32 -- Shift
            ISP(i_rec).gain_s_cal = fread(fid, 1, 'uint16');    % Byte 33-34 -- Gain
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs = fread(fid, 1, 'uint16');            % Byte 65-66 -- FS_S
            ISP(i_rec).ns = fread(fid, 1, 'uint16');            % Byte 67-68 -- NS_S
            for samples=1:256                                   % Byte 69-580 -- I2Q2 sample CALTx Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calTx_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 581-1092 -- I2Q2 sample CALComp Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calComp_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 1093-1604 -- I2Q2 sample CALRx Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calRx_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 1605-2116 -- I2Q2 sample CALdiff Ku1
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calDiff_ku1(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 2117-2628 -- I2Q2 sample CALdiff Ku2
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calDiff_ku2(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 2629-3140 -- I2Q2 sample CALTx Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calTx_ka(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 3141-3652 -- I2Q2 sample CALComp Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calComp_ka(samples) = complex(i, q);
            end
            for samples=1:256                              % Byte 3653-4164 -- I2Q2 sample CALRx Ka
                i = fread(fid, 1, 'int8');
                q = fread(fid, 1, 'int8');
                ISP(i_rec).iq_calRx_ka(samples) = complex(i, q);
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');       % Byte 4165-4166 -- PEC
            
            %-------------------------------------- TM_ECHO_RAW_SARin_OB ----------------
        case hex2dec('C1')
            if i_rec == 1
                disp('Reading TM_ECHO_RAW_SARin_OB');
            end
            % Record --------------------------------------------------------------
            fread(fid, 1, 'ubit4');                      % Byte 21 (Tracking configuration) -- bit B0..B3 reserved
            ISP(i_rec).gain_open_loop = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Gain in Open Loop bit B4 = 0 : Closed Loop gain, 1 : fixed gain
            ISP(i_rec).loop_gain_ref = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Loop gain reference B5 = 0 : nominal value, 1 : nominal value with back-off
            ISP(i_rec).acq = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- Acquisition B6 = 0 : no acquisition, 1 : acquisition
            ISP(i_rec).dem_mram_read = fread(fid, 1, 'ubit1');   % Byte 21 (Tracking configuration) -- DEM MRAM read access enabled/disabled B7 = 0 : DEM MRAM read access enabled, 1 : DEM MRAM read access disabled
            ISP(i_rec).mode_id = fread(fid, 1, 'uint8');        % Byte 22 -- Mode identifier
            ISP(i_rec).nimp = fread(fid, 1, 'uint16');          % Byte 23-24 -- NIMP
            ISP(i_rec).pri = fread(fid, 1, 'uint16');           % Byte 25-26 -- PRI
            ISP(i_rec).ambiguity_rank = fread(fid, 1, 'uint16');% Byte 27-28 -- Ambiguity rank
            ISP(i_rec).band = fread(fid, 1, 'uint8');           % Byte 29 -- Band
            fread(fid, 1, 'ubit7');                      % Byte 30 Reserved
            ISP(i_rec).loss_track = fread(fid, 1, 'ubit1');     % Byte 30 Loss of track criterion computed on Tracking Cycle (N-2) in OL mode B7 = 0 : normal, 1 : Loss of track
            ISP(i_rec).h0 = fread(fid, 1, 'uint32');            % Byte 31-34 -- Distance command H0 DEM
            ISP(i_rec).cor2 = fread(fid, 1, 'int32');          % Byte 35-38 -- Distance command COR2 DEM
            ISP(i_rec).h0_gnss = fread(fid, 1, 'uint32');       % Byte 39-42 -- Distance command H0
            ISP(i_rec).cor2_gnss = fread(fid, 1, 'int32');     % Byte 43-46 -- Distance command COR2
            ISP(i_rec).treshold = fread(fid, 1, 'uint32');      % Byte 47-50 -- Treshold
            ISP(i_rec).snr = fread(fid, 1, 'uint16');           % Byte 51-52 -- SNR
            ISP(i_rec).gates = fread(fid, 1, 'uint16');         % Byte 53-54 -- Gates
            ISP(i_rec).distance_err = fread(fid, 1, 'int32');  % Byte 55-58 -- Distance error
            ISP(i_rec).att = fread(fid, 1, 'uint8');            % Byte 59 -- Att
            ISP(i_rec).shift = fread(fid, 1, 'uint8');          % Byte 60 -- Shift
            ISP(i_rec).gain = fread(fid, 1, 'uint16');          % Byte 61-62 -- Gain
            ISP(i_rec).coarse_time = fread(fid, 1, 'uint32');   % Byte 63-66 -- Coarse Time
            fineTime1 = fread(fid, 1, 'uint16');         % Byte 67-69 -- Fine Time
            fineTime2 = fread(fid, 1, 'uint8');
            binaryFineTime = [dec2bin(fineTime1) dec2bin(fineTime2)];
            ISP(i_rec).fine_time = bin2dec(binaryFineTime);
            ISP(i_rec).nav_bullet = fread(fid, 1, 'uint8');     % Byte 70 -- Nav bullet
            ISP(i_rec).alt_ell = fread(fid, 1, 'uint32');       % Byte 71-74 -- Alt ellipsoid
            ISP(i_rec).alt_rate_ell = fread(fid, 1, 'uint32');  % Byte 75-78 -- Alt rate ellipsoid
            ISP(i_rec).orbit_num = fread(fid, 1, 'uint32');     % Byte 79-82 -- Orbit num
            ISP(i_rec).orbit_pos = fread(fid, 1, 'uint32');     % Byte 83-86 -- Orbit pos
            ISP(i_rec).thermistances.thr0 = fread(fid, 1, 'uint16');          % Thermistances 0
            ISP(i_rec).thermistances.thr1 = fread(fid, 1, 'uint16');          % Thermistances 1
            ISP(i_rec).thermistances.thr2 = fread(fid, 1, 'uint16');          % Thermistances 2
            ISP(i_rec).thermistances.thr3 = fread(fid, 1, 'uint16');          % Thermistances 3
            ISP(i_rec).thermistances.thr4 = fread(fid, 1, 'uint16');          % Thermistances 4
            ISP(i_rec).thermistances.thr5 = fread(fid, 1, 'uint16');          % Thermistances 5
            ISP(i_rec).thermistances.thr6 = fread(fid, 1, 'uint16');          % Thermistances 6
            ISP(i_rec).thermistances.thr7 = fread(fid, 1, 'uint16');          % Thermistances 7
            ISP(i_rec).thermistances.thr8 = fread(fid, 1, 'uint16');          % Thermistances 8
            ISP(i_rec).thermistances.thr9 = fread(fid, 1, 'uint16');          % Thermistances 9
            ISP(i_rec).thermistances.thr10 = fread(fid, 1, 'uint16');         % Thermistances 10
            ISP(i_rec).thermistances.thr11 = fread(fid, 1, 'uint16');         % Thermistances 11
            ISP(i_rec).thermistances.thr12 = fread(fid, 1, 'uint16');         % Thermistances 12
            ISP(i_rec).thermistances.thr13 = fread(fid, 1, 'uint16');         % Thermistances 13
            ISP(i_rec).thermistances.thr14 = fread(fid, 1, 'uint16');         % Thermistances 14
            ISP(i_rec).fs = fread(fid, 1, 'uint16');            % Byte 117-118 -- FS_S
            ISP(i_rec).ns = fread(fid, 1, 'uint16');            % Byte 119-120 -- NS_S
            fseek(fid,1,'cof');                                 % Byte 121 -- Reserved
            ISP(i_rec).nBurst = fread(fid, 1, 'uint8');         % Byte 122 -- Burst number (from 1 to 12)
            for pulses = 1:64
                for samples=1:256                          % Byte 123-32890 -- I2Q2 sample
                    i = fread(fid, 1, 'int8');
                    q = fread(fid, 1, 'int8');
                    ISP(i_rec).iq(pulses, samples) = complex(i, q);
                end
            end
            ISP(i_rec).pec = fread(fid, 1, 'uint16');           % Byte 633-634 -- PEC
            
            
            
    end
    
    
end

disp(['The ISP binary file has been read. Filesize is ' num2str(filesize) '. Reading position: byte nº ' num2str(ftell(fid)) '.']);
fclose (fid);

toc
end
