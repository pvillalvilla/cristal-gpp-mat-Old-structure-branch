%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Process from ISP to L1A
%
% Calling: 
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
% v1.1 added orbit struct computation and onboard reversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filesBulk = L1A_chain (filesBulk, chd, cnf, cst, options)


%% init variables
i_burst=1;
N_bursts = get_num_record(filesBulk.filename_L0,'nb');
if(N_bursts < chd.N_bursts_cycle_sar*chd.N_pulses_burst/2) % Condition to avoid processing L1A files with less records than the ones needed o fill half stack
    disp(['Not enough records inside mask to create a complete stack for the file ' filesBulk.filename_L0]);
    fprintf(filesBulk.fid_log,'%s\n' ,['Not enough records inside mask to create a complete stack for the file ' filesBulk.filename_L0]);
    return;
end

if(cnf.writting_flag(1))
    if cnf.skip_burst2_xRC && (strcmp(chd.meas_mode,'OPEN_BURST'))
        [filesBulk] = create_NetCDF_L1A(filesBulk,N_bursts-60, cnf, chd, cst); %manually create a l1a with the bursts dedicated to CAL4 (in SARIN OB) skipped
    elseif cnf.skip_burst2_xRC && (strcmp(chd.meas_mode,'CLOSED_BURST'))       
        [filesBulk] = create_NetCDF_L1A(filesBulk,N_bursts-20, cnf, chd, cst); %manually create a l1a with the bursts dedicated to CAL4 (in SARIN CB) skipped
    else
        [filesBulk] = create_NetCDF_L1A(filesBulk,N_bursts, cnf, chd, cst);
    end
end
[L1A]           = create_L1A_struct(cnf, chd);
[L1A_buffer]    = create_L1A_struct(cnf, chd);
[L1AP]          = create_L1AP_struct(cnf, chd, N_bursts);
[ORBIT]         = read_orbit (filesBulk);


bursts_skipped=0;
while(i_burst<=N_bursts)
    %% READING and ADAPTING
    if (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'CLOSED_BURST')) && (i_burst==2||i_burst==42||i_burst==82||i_burst==122||i_burst==162||i_burst==202||i_burst==242||i_burst==282||i_burst==322||i_burst==362||i_burst==402||i_burst==442||i_burst==482||i_burst==522||i_burst==562||i_burst==602||i_burst==642||i_burst==682||i_burst==722||i_burst==762)
        %CAL4 SARIN CB pulse, SKIP
        i_burst=i_burst+1;
        bursts_skipped=bursts_skipped+1;
        
    % TESTING: skip all bursts (1,2,3 of each RC) with CAL pulses    
%     elseif (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst==1||any(i_burst==1+12*int16(1:59)))
%         i_burst=i_burst+1;
%         bursts_skipped=bursts_skipped+1;
    elseif cnf.skip_burst2_xRC && (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst==2||any(i_burst==2+12*int16(1:59)))
        i_burst=i_burst+1;
        bursts_skipped=bursts_skipped+1;
%     elseif (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst==3||any(i_burst==3+12*int16(1:59)))
%         i_burst=i_burst+1;
%         bursts_skipped=bursts_skipped+1;
    end
    
    [L1A,filesBulk]     = read_adapt_L0_record(filesBulk,L1A,i_burst,N_bursts, cnf, chd, cst);
    
    % Set to zero the waveforms of CAL4 in the OB pulses, but don't skip burst
    if (strcmp(cnf.processing_mode,'SIN')) && (strcmp(chd.meas_mode,'OPEN_BURST')) && (i_burst==2||any(i_burst==2+12*int16(1:59)))

        L1A.wfm_cal_gain_corrected(10,:)=0;
        L1A.wfm_cal_gain_corrected_2(10,:)=0;
        % trying to set it to just thermal noise (like previous burst) instead of 0
%         L1A.wfm_cal_gain_corrected(10,:)=L1A.wfm_cal_gain_corrected(9,:);
%         L1A.wfm_cal_gain_corrected_2(10,:)= L1A.wfm_cal_gain_corrected_2(9,:);
        %bursts_skipped=bursts_skipped+1;
        %(strcmp(chd.meas_mode,'CLOSED_BURST')) && (i_burst==2||14 26
    end
    
    if cnf.starting_from_netcdf
        [L1A,L1AP]          = preliminary_datation_old (L1A, L1AP, cnf, chd, cst);
        [L1A,L1AP]          = window_delay_old (L1A,L1AP, cnf, chd, cst);
        [L1A,L1AP]          = final_datation_old (L1A,L1AP,cnf, chd, cst, ORBIT);
    else
        [L1A,L1AP]          = preliminary_datation (L1A, L1AP, cnf, chd, cst);
        [L1A,L1AP]          = window_delay (L1A,L1AP, cnf, chd, cst);
        [L1A,L1AP]          = final_datation (L1A,L1AP,cnf, chd, cst, ORBIT);
    end
    
    if(cnf.onboard_reversion_flag) %RMC
        [L1A,chd]               = onboard_reversion_old(L1A,filesBulk, chd,cnf,cst);
    end
    
    %     % Flagging bursts with OB CAL pulses
    %     if any(1:3:N_bursts == i_burst)
%         L1A.OB_cal_flag = 1; % 1st burst of TC, it has 4 holes
%         L1A.OB_cal_mask=[0 0 0 0 ones(1,60)];
%     elseif any(2:3:N_bursts == i_burst) 
%         L1A.OB_cal_flag = 2; % 2nd burst of TC, it has 3 holes
%         L1A.OB_cal_mask=[ones(1,7) 0 0 0 ones(1,54)];
%     elseif any(3:3:N_bursts == i_burst)
%         L1A.OB_cal_flag = 3; % 3rd burst of TC, it has 1 holes
%         L1A.OB_cal_mask=[ones(1,15) 0 ones(1,48)];
%     else
%         L1A.OB_cal_flag = 4; % 4 to 12 burst of TC, it has no holes
%         L1A.OB_cal_mask=ones(1,64);
%     end
    
        % Flagging bursts with OB CAL pulses
%     if any(1:3:N_bursts == i_burst) 
%         L1A.OB_cal_flag = 1; % 1st burst of TC, it has 4 holes
%         L1A.OB_cal_mask=[0 0 0 0 ones(1,60)];
%     elseif any(2:3:N_bursts == i_burst) 
%         L1A.OB_cal_flag = 2; % 2nd burst of TC, it has 3 holes
%         L1A.OB_cal_mask=[ones(1,7) 0 0 0 ones(1,54)];
%     else any(3:3:N_bursts == i_burst)
%         L1A.OB_cal_flag = 3; % 3rd burst of TC, it has 1 holes
%         L1A.OB_cal_mask=[ones(1,15) 0 ones(1,48)];
% %     else
% %         L1A.OB_cal_flag = 4; % 4 to 12 burst of TC, it has no holes
% %         L1A.OB_cal_mask=ones(1,64);
%     end

    L1A_buffer(i_burst-bursts_skipped) = L1A;

    
    %% WRITING
    if(cnf.writting_flag(1)==1)
        if i_burst==1
            filesBulk.ncid_L1A = netcdf.open(filesBulk.filename_L1A,'WRITE');
        end

        write_NetCDF_L1A (filesBulk,L1A_buffer(i_burst-bursts_skipped),i_burst-bursts_skipped, cnf, chd, cst);
        if i_burst==N_bursts
            
%             ncwriteatt(filesBulk.filename_L1A,'/','creation_time',datestr(now));
%             ncwriteatt(filesBulk.filename_L1A,'/','data_info','data simulated by ARESYS and processed by isardSAT');
            netcdf.close(filesBulk.ncid_L1A);
        end
    end
    i_burst=i_burst+1;
end

end