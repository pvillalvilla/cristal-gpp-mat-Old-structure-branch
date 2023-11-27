%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Process from L1A to L1B
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
%            JP López-Zaragoza / isardSAT
%---------------------------
% Version  record
% (JP López-Zaragoza) 1.1, 2023/03/15,  Updated the condition to process a new surface based on the new beam angles surface focussing
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filesBulk = DDP_chain (filesBulk, chd, cnf, cst, options)
%% init variables
i_burst             = 1;
i_surf              = 1;
i_burst_focussed    = 1;
i_surf_stacked      = 1;
exit_gpp            = 0;
N_bursts = get_num_record(filesBulk.filename_L1A,'nb');
N_surfs_loc_estimated = floor(N_bursts/chd.N_bursts_cycle_sar*5.0);
burst_indexes_in_memory = 1: N_bursts;
if(N_bursts < chd.N_bursts_cycle_sar*chd.N_pulses_burst/2) % Condition to avoid processing L1A files with less records than the ones needed o fill half stack
    disp(['Not enough records inside mask to create a complete stack for the file ' filesBulk.filename_L1A]);
    fprintf(filesBulk.fid_log,'%s\n' ,['Not enough records inside mask to create a complete stack for the file ' filesBulk.filename_L1A]);
    return;
end
if(cnf.writting_flag(2))
    [filesBulk] = create_NetCDF_L1Bs_HR(filesBulk,N_surfs_loc_estimated, N_bursts,chd.N_max_beams_stack, cnf, chd, cst,0);   
end
if(cnf.writting_flag(3))
    [filesBulk] = create_NetCDF_L1B(filesBulk,N_surfs_loc_estimated, N_bursts, cnf, chd, cst,0);
end
[L1A]           = create_L1A_struct(cnf, chd);
[L1A_buffer]    = create_L1A_struct(cnf, chd);
[L1BS_buffer]   = create_L1BS_struct(cnf, chd);
L1B             = [];

while(i_burst<=N_bursts)
    %% READING   
    [L1A,filesBulk]  = read_L1A_record(filesBulk,L1A,i_burst,N_bursts, cnf, chd, cst);
    L1A_buffer(i_burst) = L1A;   
    %% PROCESSING
    [L1BS_buffer, out_surf, cnf] = surface_locations (L1A_buffer,L1BS_buffer,N_bursts, i_burst, i_surf, cnf, chd, cst);
    if(i_surf < out_surf)
        i_surf = out_surf;
        N_total_surf_loc = i_surf-1;
        new_surface = 1;
    end
    if(i_surf > 64)
        if(i_burst_focussed==1)% handle final bursts
            burst_margin = i_burst;
        end
        
        [L1A_buffer(i_burst_focussed)]      = beam_angles (L1A_buffer(i_burst_focussed),   L1BS_buffer, N_total_surf_loc,i_surf_stacked, i_burst_focussed,N_bursts, cnf, chd, cst);
        [L1A_buffer(i_burst_focussed)]      = azimuth_processing    (L1A_buffer(i_burst_focussed),i_burst_focussed,N_bursts, cnf, chd, cst);
        
        i_burst_focussed = i_burst_focussed+1;
               
        nnz_beams=find(L1A_buffer(i_burst_focussed-1).surf_loc_index(:));
        if( L1A_buffer(i_burst_focussed-1).surf_loc_index(nnz_beams(1))>i_surf_stacked )
            stack_processing;
            i_surf_stacked = i_surf_stacked +1;
            if(i_surf_stacked>N_surfs_loc_estimated)
                disp('exit process');
                exit_gpp = 1;
                break;
            end
        end                      
    end
    
    i_burst=i_burst+1;
end

%% PROCESSING Lasts Bursts
last_burst=i_burst_focussed;
if(~exit_gpp)
    for i_burst_focussed = last_burst:N_bursts
       
        [L1A_buffer(i_burst_focussed)]      = beam_angles (L1A_buffer(i_burst_focussed),   L1BS_buffer, N_total_surf_loc,i_surf_stacked,i_burst_focussed,N_bursts, cnf, chd, cst);
        [L1A_buffer(i_burst_focussed)]      = azimuth_processing    (L1A_buffer(i_burst_focussed),i_burst_focussed,N_bursts, cnf, chd, cst);
             
        clear nnz_beams
        nnz_beams=find(L1A_buffer(i_burst_focussed).surf_loc_index(:));
        if( L1A_buffer(i_burst_focussed).surf_loc_index(nnz_beams(1))>i_surf_stacked )
            stack_processing;
            i_surf_stacked = i_surf_stacked +1;
            if(i_surf_stacked>N_surfs_loc_estimated)
                disp('exit process');
                exit_gpp = 1;
                break;
            end
        end
    end
end

% TEST PLOTS
%indx_non_zeroBeams=[L1A_buffer(:).N_beams]>0;
%surf_loc_index_matrix    = reshape([L1A_buffer(indx_non_zeroBeams).surf_loc_index],[chd.N_pulses_burst.*cnf.zp_fact_azimut,nnz(indx_non_zeroBeams)]).';
%figure;imagesc(surf_loc_index_matrix)
%figure;plot(surf_loc_index_matrix(:,64))

 %% PROCESSING Lasts Stacks
last_stack = i_surf_stacked;
if(~exit_gpp)
    for i_surf_stacked = last_stack:N_total_surf_loc
       
        stack_processing;
        if(i_surf_stacked>N_surfs_loc_estimated)
            disp('exit process');
            exit_gpp = 1;
            break;
        end
    end
end

%% ----------- Resize the netCDF file ------------------------------

for i=1:N_total_surf_loc
    beams_x_stack(i)=max(L1BS_buffer(i).N_beams_stack);
end
max_beams_x_stack=max(beams_x_stack);
if(cnf.writting_flag(2))
    netcdf.close(filesBulk.ncid_L1Bs); % close the already open netcdf
    filesBulk.filename_netCDF_2remove_Bs=filesBulk.filename_L1Bs;
    [files] = create_NetCDF_L1Bs_HR(filesBulk,i_surf_stacked-1,N_total_surf_loc-1,max_beams_x_stack, cnf, chd, cst,1);
    resize_NetCDF_L1Bs_HR(files,i_surf_stacked-1,max_beams_x_stack,cnf,chd,cst);
    delete(filesBulk.filename_netCDF_2remove_Bs);
end
if(cnf.writting_flag(3))
    netcdf.close(filesBulk.ncid); % close the already open netcdf
%     filesBulk.filename_netCDF_2remove=filesBulk.filename_L1B;
%     %[files] = create_NetCDF_L1B(filesBulk,i_surf_stacked-1,N_total_surf_loc-1, cnf, chd, cst);
%     [filesBulk] = create_NetCDF_L1B(filesBulk,i_surf_stacked-1,N_total_surf_loc-1, cnf, chd, cst,1);
%     resize_NetCDF_L1B(filesBulk,i_surf_stacked-1,cnf,chd,cst);
%     delete(filesBulk.filename_netCDF_2remove);
end