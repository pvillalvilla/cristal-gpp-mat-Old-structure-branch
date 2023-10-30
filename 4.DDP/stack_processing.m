%% Algorithms related with Stack Processing
% Stacking
% Geometry Corrections
% Range Transformation
% Stack Extension
% Multilooking

% v1.1 REmoved commented lines and activated sigma0 scaling factor.

[L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked), cnf, chd, cst);
[L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked), cnf, chd, cst);
[L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked), cnf, chd, cst);
if cnf.extend_window
    [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked), cnf, chd, cst);
end
[L1BS_buffer(i_surf_stacked),L1B]   = multilooking          (L1BS_buffer(i_surf_stacked), cnf, chd, cst);

if(L1BS_buffer(i_surf_stacked).surface_type_flag == 4)
    PT_stack_performance;
end

if cnf.include_wfms_aligned
    
    if i_surf_stacked== 1
        
        win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
        alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
    end
    [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref, cnf, chd, cst);
end

if cnf.starting_from_netcdf
    [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor_old (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B, cnf, chd, cst);
else
    [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B, cnf, chd, cst);
end

%% Writting routines
if(cnf.writting_flag(2))
    if i_surf_stacked==1
        filesBulk.ncid_L1Bs = netcdf.open(filesBulk.filename_L1Bs,'WRITE');
    end
    write_NetCDF_L1Bs_HR(filesBulk,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked,cnf, chd, cst);
end
if(cnf.writting_flag(3))
    if i_surf_stacked==1
        filesBulk.ncid = netcdf.open(filesBulk.filename_L1B,'WRITE');
    end
    write_NetCDF_L1B (filesBulk,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked-1,cnf, chd, cst);
end

indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1BS_buffer(i_surf_stacked).burst_index));

if(cnf.verify_L1BS_internal)
    plots_L1BS;
end

if(~isempty(indexes_to_remove))
    [L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
    burst_indexes_in_memory(indexes_to_remove)=[];
    
end
[L1BS_buffer(i_surf_stacked)] = empty_L1BS_struct(L1BS_buffer(i_surf_stacked));



if(L1BS_buffer(i_surf_stacked).surface_type_flag == 4)
    %TRP location being processed
    disp(['Computing the results for the record ' num2str(i_surf_stacked)]);    
end