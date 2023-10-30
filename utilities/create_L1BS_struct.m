function [L1BS] = create_L1BS_struct (cnf, chd)




L1BS         = [];
L1BS = setfield_new(L1BS,'surf_counter',0);
L1BS = setfield_new(L1BS,'lat_surf',0);
L1BS = setfield_new(L1BS,'lon_surf',0);
L1BS = setfield_new(L1BS,'alt_surf',0);
L1BS = setfield_new(L1BS,'x_surf',0);
L1BS = setfield_new(L1BS,'y_surf',0);
L1BS = setfield_new(L1BS,'z_surf',0);
L1BS = setfield_new(L1BS,'alt_rate_sat',0);
L1BS = setfield_new(L1BS,'time_surf',0);
L1BS = setfield_new(L1BS,'lat_sat',0);
L1BS = setfield_new(L1BS,'lon_sat',0);
L1BS = setfield_new(L1BS,'alt_sat',0);
L1BS = setfield_new(L1BS,'x_sat',0);
L1BS = setfield_new(L1BS,'y_sat',0);
L1BS = setfield_new(L1BS,'z_sat',0);
L1BS = setfield_new(L1BS,'x_vel_sat',0);
L1BS = setfield_new(L1BS,'y_vel_sat',0);
L1BS = setfield_new(L1BS,'z_vel_sat',0);
L1BS = setfield_new(L1BS,'roll_surf',0);
L1BS = setfield_new(L1BS,'pitch_surf',0);
L1BS = setfield_new(L1BS,'yaw_surf',0);
L1BS = setfield_new(L1BS,'win_delay_surf',0);
L1BS = setfield_new(L1BS,'surface_type_flag',[]);
L1BS = setfield_new(L1BS,'win_delay_surf_aligned',0);
L1BS = setfield_new(L1BS,'burst_index',0);
L1BS = setfield_new(L1BS,'ProcessID_surf',0);
L1BS = setfield_new(L1BS,'doppler_ang_surf',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'beam_ang_surf',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'ProcessID_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'lat_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'lon_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'alt_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'x_vel_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'y_vel_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'z_vel_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'x_sar_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'y_sar_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'z_sar_sat_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'win_delay_beam',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'pri_sar_sat_beam',zeros (1,chd.N_max_beams_stack));

L1BS = setfield_new(L1BS,'Gap_flag',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'pointing_ang_surf',1./zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'T0_sar_surf',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'beams_surf',zeros(chd.N_max_beams_stack,chd.N_samples_sar));
L1BS = setfield_new(L1BS,'N_beams_stack',0);
L1BS = setfield_new(L1BS,'beam_index',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'look_ang_surf',1./zeros (1,chd.N_max_beams_stack));

if(strcmp(cnf.processing_mode,'SIN'))
    L1BS = setfield_new(L1BS,'beams_surf_2',zeros(chd.N_max_beams_stack,chd.N_samples_sar));
end

% geometry corrections
L1BS = setfield_new(L1BS,'doppler_corr',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'range_sat_surf',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'slant_range_corr_time',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'slant_range_corr',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'wd_corr',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'shift',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'shift_coarse',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'shift_fine',zeros (1,chd.N_max_beams_stack));
L1BS = setfield_new(L1BS,'N_windows',0);
L1BS = setfield_new(L1BS,'N_windows_back',0);
L1BS = setfield_new(L1BS,'N_windows_fore',0);
L1BS = setfield_new(L1BS,'good_samples',zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
L1BS = setfield_new(L1BS,'beam_geo_corr',zeros(chd.N_max_beams_stack,chd.N_samples_sar));
L1BS = setfield_new(L1BS,'beam_ref',zeros(chd.N_max_beams_stack,chd.N_samples_sar));
L1BS = setfield_new(L1BS,'beam_ref2',zeros(chd.N_max_beams_stack,chd.N_samples_sar));

if(strcmp(cnf.processing_mode,'SIN'))
    L1BS = setfield_new(L1BS,'beam_geo_corr_2',zeros(chd.N_max_beams_stack,chd.N_samples_sar));
end
% range transformation
L1BS = setfield_new(L1BS,'beams_rng_cmpr',zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
L1BS = setfield_new(L1BS,'beams_rng_cmprIQ',zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
if(strcmp(cnf.processing_mode,'SIN'))
    L1BS = setfield_new(L1BS,'beams_rng_cmpr_2' ,zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
    L1BS = setfield_new(L1BS,'beams_rng_cmprIQ_2',zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
    L1BS = setfield_new(L1BS,'phase_diff'       ,zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
    L1BS = setfield_new(L1BS,'coherence'        ,zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
    L1BS = setfield_new(L1BS,'coherence_mask'   ,zeros(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
	if(cnf.coherency_mask)
        L1BS = setfield_new(L1BS,'beams_rng_cmpr_Coh' ,[]);
        L1BS = setfield_new(L1BS,'beams_rng_cmpr_2_Coh' ,[]);
        L1BS = setfield_new(L1BS,'beams_rng_cmprIQ_Coh' ,[]);
        L1BS = setfield_new(L1BS,'beams_rng_cmprIQ_2_Coh' ,[]);
        L1BS = setfield_new(L1BS,'phase_diff_Coh' ,[]);

    end
end


% stack masking
L1BS = setfield_new(L1BS,'stack_mask',ones(chd.N_max_beams_stack,chd.N_samples_sar*cnf.zp_fact_range));
% multi-looking
% sigma0
L1BS = setfield_new(L1BS,'wfm_scaling_factor_beam',zeros (1,chd.N_max_beams_stack));



% L1BS = setfield_new(L1BS,'beams_focused_shifted',0);
% if(strcmp(cnf.processing_mode,'SIN'))
%     L1BS = setfield_new(L1BS,'beams_focused_shifted_2',zeros(N_ku_pulses_burst_chd,chd.N_samples_sar));
% end



end

