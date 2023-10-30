
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT Ltd. 
% --------------------------------------------------------
%
% Ground Processor Prototype for altimetry missions:
% CryoSat-2 SARIn
% CryoSat-2 SAR
% Sentinel 6 RAW
% Sentinel 6 RMC
% 
% ---------------------------------------------------------
% L0 (ISP)+ orbit file + attitude file + meteo file% Inputs: 

% L1A (FBR in case of CryoSat-2)
% L1B-S
% L1B
% 
% Output
%   NetCDF L1A, L1B-S, L1B and/or L2
%
% ----------------------------------------------------------
% 
% Authors:  Albert Garcia-Mondejar  / isardSAT
%           Eduard Makhoul          / isardSAT
%           Roger Escola Jane       / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1.0 2016/04/xx   First version based on S6 GPP
% v2.0 2016/04/xx   Changed reading routines for each record
% v2.1 2016/04/xx   Changed for bulk processing. Files are opened one time within read_L1A_record
% v2.2 2016/05/12   Changed FBR for L1A. Opening file moved inside read_L1A_record fuction
% v2.3 2016/05/17   Include the possibility to activate or deactivate the alignment of the L1B waveforms w.r.t first one
%					Stacking of the remaining last surfaces
% v2.4 2016/05/18 	Add writting L1BS
% v2.5 2016/05/20 	Inlucded return to skip processing when no records are found within the mask
% v2.6 2016/05/20   Included extend_window_cnf option

function L1B_processing_ACDC_bis(filesBulk,i_fileL1A_input,targz_option_active,varargin)
    
	if(nargin<3 || nargin>(3+2*1))
		error('Wrong number of input parameters');
	end
	p = inputParser;
	p.addParamValue('proc_bsln_id',{''},@(x)ischar(x));
	p.parse(varargin{:});
	proc_bsln_id=char(p.Results.proc_bsln_id);
	clear p;

	global PLRM_mode ref_burst N_files files2delete_indexs mode
    global N_bursts_cycle_chd N_ku_pulses_burst_chd
    global extend_window_cnf include_wfms_aligned fromSARin
    
    writing_L1B = 1;
	writing_L1BS = 0;
    writing_L1B_pLRM = 1;
    exit_gpp = 0;
    
    i_RC = 0;
    clear wfms wfms_wdcorr
    
    %% Sanity check
    if PLRM_mode == 1
        writing_L1B = 0;
    else
        writing_L1B_pLRM = 0;
    end

    %%
	%EM: added 03.10.2016
    global ACDC_application_cnf cnf_p_ACDC figures_visible_cnf first
    global bw_ku_chd zp_fact_range_cnf N_samples c_cst zp_fact_azimut_cnf
    
    if figures_visible_cnf
        set(0, 'DefaultFigureVisible', 'on');
    else
        set(0, 'DefaultFigureVisible', 'off');
    end

    %global attributes of the main processing options
    global hamming_window_cnf height_rate_application_cnf FAI_application_cnf CAL2_flag_cnf CAL1p2p_flag_cnf
    global apply_stack_mask_cnf avoid_beams_mask_allzeros_cnf avoid_noisy_beams_cnf method_noise_est_cnf param_std_noise_thres_cnf noise_estimation_window mask_look_angles_CR2_cnf mask_look_angles_CR2_max_chd mask_look_angles_CR2_min_chd use_zeros_cnf
    global force_exact_method_cnf  sign_dopp_corr sign_beamforming S3_outer_beams_mask_cnf S3_outer_beams_mask_in_cnf S3_outer_beams_mask_end_cnf
    t6 = tic;
    %% Define the files structure for the specific FBR to be processed
    % including all other files (cnf,cst,chd,....) not DBL nor HDR related to the specific
    % folder (in order no to change the functions)
    files.inputPath =filesBulk.inputPath;
    files.resultPath =filesBulk.resultPath;
    files.inputFiles =filesBulk.inputFiles(~(filesBulk.filterDATAFILES));
    if targz_option_active
        %untar the file
        untar([files.inputPath filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name],files.inputPath);
        filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name=strrep(filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name,'TGZ','DBL');
    end
    
    %add the current L1A to be processed
    files.inputFiles(end+1) = filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input));
    files.indexaDirs      =   find(([files.inputFiles.isdir]));
    files.indexFiles      =   find(not([files.inputFiles.isdir]));
    files.nFiles          =   length(files.indexFiles);
    
    filename_L1A=filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name;
    disp(strcat('Processing L1A: ',filename_L1A))
    disp('....')
    
    progressbar('Burst processing','Surfaces computation','Burst focussing','Stack processing');
    i_surf              = 1;
    i_burst_focused     = 1;
    i_surf_stacked      = 1;
    % previous_gap_flag   = 0; % indicates that the burst has a gap unable to be fullfilled
 
    %% Read inputs
    [steps,meteo,files,first_burst,final_burst] = read_inputs(files,'proc_bsln_id',proc_bsln_id);
    N_bursts    = final_burst-first_burst+1;
    
    if (N_bursts < N_bursts_cycle_chd*N_ku_pulses_burst_chd)
		disp(['Not enough records inside mask to create a complete stack for the file ' files.filename_L1A])
        files2delete_indexs(i_fileL1A_input) = 1;
		return;
    else
        N_files = N_files + 1;
    end
    
    
% % %     %% fix
% % %     final_burst = first_burst + 2569 - 1;
% % %     N_bursts    = final_burst-first_burst+1;
    
    
    %%
    original_burst_index    =   first_burst:final_burst;
    burst_indexes_in_memory = 1: N_bursts;
    N_surfs_loc_estimated = floor(N_bursts/N_bursts_cycle_chd*2.0);
	disp(strcat('# Estimated surfaces',{' '},num2str(N_surfs_loc_estimated)))    
	if(writing_L1B)
		[files] = create_NetCDF_L1B_S3(files,N_surfs_loc_estimated);
    end
    if(writing_L1BS)
		[files] = create_NetCDF_L1BS_S3(files,N_surfs_loc_estimated);   
    end
    if writing_L1B_pLRM
        [files] = create_NetCDF_L1B_pLRM(files,floor(N_bursts/N_bursts_cycle_chd));
    end
    
    [L1A]           = create_L1A_struct;
    [L1A_buffer]    = create_L1A_struct;
    [L1BS_buffer]   = create_L1BS_struct;
    [L1B]           = create_L1B_struct;
    
    first = 1; %% used for displaying when the gaps are present (used in the check_continuity funcion)
    i_burst=1;
    
    %fix
%     N_bursts = floor(N_bursts/2);
    
    N_bursts_original=N_bursts;
%     N_bursts_original=2569;
    while(i_burst<=N_bursts)
% % %         i_burst
        progressbar(i_burst/N_bursts,[],[],[]);
        %% BURST Processor
        if (steps(1))
            [ISP]  = read_ISP_record(files,first_burst,original_burst_index(i_burst),N_bursts); %To be build
            [L1A]  = pre_datation    (ISP);
            [L1A]  = pre_win_delay   (L1A,ISP);
            [L1A]  = final_burst_datation   (L1A,ISP);
            [L1A]  = onboard_reversion   (L1A,ISP);
            [L1A]  = instrument_gain_corr   (ISP,L1A);
            [L1A]  = waveforms_correction   (ISP,L1A);
        else
            
            [L1A,files,CAI,FAI]  = read_L1A_record(files,L1A,original_burst_index(i_burst),i_burst);
            if(i_burst > 1) %&& strcmp(mode,'SAR'))
                %if there is are missing burst records in the L1A file, jumps
                %in time, they are created to provide continuity to the stacks.
                [L1A,original_burst_index,N_bursts] = check_continuity (L1A_buffer(end),L1A,original_burst_index,i_burst,N_bursts);
            end
            if(L1A.confi_block_degraded==1&&L1A.days==0) % handle degraded bursts
                [L1A,end_of_file,i_burst]  = fill_empty_burst(files,L1A_buffer(end),L1A,original_burst_index(i_burst),N_bursts_original,i_burst,first_burst);
                if(end_of_file==1)
                    N_bursts=i_burst-1;
                    break;
                end
            end
        end
        L1A_buffer(i_burst) = L1A;

        
        %%
        if(steps(2))
            switch mode
                case 'SAR'
                    [L1BS_buffer, out_surf] = surface_locations_SAR (L1A_buffer,L1BS_buffer, i_burst, i_surf);
                case 'SIN'
                    [L1BS_buffer, out_surf] = surface_locations_SARin (L1A_buffer,L1BS_buffer, i_burst, i_surf);
            end
            if(i_surf < out_surf)
                progressbar([],i_surf/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd,[],[]);
                i_surf = out_surf;
                N_total_surf_loc = i_surf-1;
                new_surface = 1;
                %             disp([i_burst N_total_surf_loc]);
            end
            
			if(i_surf > 64)
                
                if(i_burst_focused==1)% handle final bursts
                    burst_margin = i_burst;
                end
                progressbar([],[],i_burst_focused/N_bursts,[]);
                
                [L1A_buffer(i_burst_focused)]      = beam_angles (L1A_buffer(i_burst_focused),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
                [L1A_buffer(i_burst_focused)]      = azimuth_processing    (L1A_buffer(i_burst_focused));
                i_burst_focused = i_burst_focused+1;
                
                
                if(L1A_buffer(i_burst_focused-1).surf_loc_index(1)~=i_surf_stacked)
                    %             if(i_burst_focused > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
                    progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                    
                    
                    [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));

                    
                    [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked),i_surf_stacked);


                    [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));


                    if extend_window_cnf
                        [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
                    end

                    
                    [L1BS_buffer(i_surf_stacked),L1B]   = multilooking_new          (L1BS_buffer(i_surf_stacked));
                    

                    
                    if fromSARin
                        %average_ch1_ch2

                        L1B.wfm_cor_i2q2_sar_ku = 0.5*(L1B.wfm_cor_i2q2_sar_ku + L1B.wfm_cor_i2q2_sar_ku_2);

                    end
                    
                    if include_wfms_aligned
                        %surface alignment using window delay reference of first
                        %surface
                        if i_surf_stacked== 1
                            %reference surface
                            win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                            alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                        end
%                         [L1BS_buffer(i_surf_stacked),L1B]   = extend_L1B(L1B,L1BS_buffer(i_surf_stacked));
                        [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
                        wfms_wdcorr(i_surf_stacked,:) = L1B.wfm_cor_i2q2_sar_ku_wdcorr;
                    end
                    wfms(i_surf_stacked,:) = L1B.wfm_cor_i2q2_sar_ku;
                                        
                    
                    [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
                    
                    
					if(writing_L1B)
                        if i_surf_stacked==1
                            files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                        end
                        prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
					end
					if(writing_L1BS)
                        if i_surf_stacked==1
                        	files.ncid_BS = netcdf.open(files.filename_netCDF_BS,'WRITE');
                        end
						prepare_NetCDF_L1BS_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
					end 

                    i_surf_stacked = i_surf_stacked +1;
                    if(i_surf_stacked>N_surfs_loc_estimated)
                        disp('exit process');
                        exit_gpp=1;
                        break;
                    end
                    
                end
            end
        
        elseif (((L1A_buffer(i_burst).burst_sar_ku_fbr == N_bursts_cycle_chd)&&strcmp(mode,'SAR')&&(i_burst>=N_bursts_cycle_chd)) || (strcmp(mode,'SIN')))

            %% pLRM chain  
            i_RC = i_RC + 1;
                L1B(i_RC) = plrm_chain(L1A_buffer(i_burst-N_bursts_cycle_chd+1:i_burst),i_RC);

                sigma0_scal(i_RC) = L1B.sigma0_scaling_factor_RC;

                if writing_L1B_pLRM
                    if i_RC == 1
                        files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
                    end
                    switch mode
                        case 'SAR'
                            prepare_NetCDF_L1B_pLRM (files,L1A_buffer(i_burst-N_bursts_cycle_chd+ref_burst),L1B(i_RC),i_RC);
                        case 'SIN'
                            prepare_NetCDF_L1B_pLRM (files,L1A_buffer(i_burst),L1B(i_RC),i_RC);
                    end
                end
                if i_burst >= N_bursts
                    disp('exit process');
                    exit_gpp = 1;
                    break;
                end
%             else
%                 i_RC = i_RC - 1;
%             end
        end
    i_burst=i_burst+1;
    end
    
% % %     if PLRM_mode
% % %         wfms            = circshift(wfms,[0,30]);
% % %         [wfms]          = extend_L1B_PLRM(wfms);
% % %         [wfms_wdcorr]   = surface_win_delay_alignment_PLRM (L1B,wfms);
% % %         wfms_wdcorr     = circshift(wfms_wdcorr,[0,30]);
% % %     end
    
%     h=figure; imagesc(wfms);
%     print(h,['wfms_',num2str(i_fileL1A_input)],'-dpng');
%     saveas (h,['wfms_',num2str(i_fileL1A_input),'.fig'])
%     h=figure; imagesc(wfms_wdcorr);
%     print(h,['wfms_',num2str(i_fileL1A_input),'_wdcorr'],'-dpng');
%     saveas (h,['wfms_wdcorr_',num2str(i_fileL1A_input),'.fig'])
    
    %%
    %close the input binary file
    fclose(files.fid);
    last_burst=i_burst_focused;
    % Last bursts and stacks
    if(~exit_gpp) & (PLRM_mode == 0)
        for i_burst_focused = last_burst:N_bursts            
            progressbar([],[],i_burst_focused/N_bursts,[]);
            [L1A_buffer(i_burst_focused)]      = beam_angles (L1A_buffer(i_burst_focused),   L1BS_buffer, N_total_surf_loc,i_surf_stacked);
            [L1A_buffer(i_burst_focused)]      = azimuth_processing    (L1A_buffer(i_burst_focused));
            
            if(L1A_buffer(i_burst_focused).surf_loc_index(1)~=i_surf_stacked)
                %             if(i_burst_focused > (N_bursts_cycle_chd*N_ku_pulses_burst_chd))
                progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
                [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
                [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
                
                
                
                %--------------
%                 h=figure; imagesc(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr.'),1).').^2);
%                 title(strcat('Stack #',num2str(i_surf_stacked)));
%                 print('-dpng',['./results/plots/dedop/beam_geo_corr_surf',num2str(i_surf_stacked),'.png']);
%                 close(h)
                %--------------
                
                
                
                [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
                if extend_window_cnf
                    [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
                end
                [L1BS_buffer(i_surf_stacked),L1B]   = multilooking_new          (L1BS_buffer(i_surf_stacked));
                if include_wfms_aligned
                    if i_surf_stacked== 1
                        %reference surface
                        win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                        alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                    end
                    [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
                end
                [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
                
                %added EM: 30.10.2016
                %% ---------------- ACDC APPLICATION ----------------------
                if ACDC_application_cnf
                    if i_surf_stacked==1
                        %% ----------------- FITTING PARAMS INITIALIZE ----------------------------
                        %--------------------------------------------------------------------------
                        if cnf_p_ACDC.rou_flag
                            if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou cnf_p_ACDC.ini_Pu];
                            else
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou];
                            end
                        else
                            if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4 cnf_p_ACDC.ini_Pu];
                            else
                                fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4];
                            end
                        end
                        fit_params_ini_ACDC=[];
                        if cnf_p_ACDC.ini_Hs_rou_sliding_win_opt==1
                            accumulated_sigmaz_conv=0;
                            accumulated_SWH_ACDC=0;
                            accumulated_SSH_ACDC=0;
                            accumulated_amp_ACDC=0;
                        end
                        i_surf_fitted=0;
                        %% ----------------- INITIALIZE RANGE_INDEXATION --------------------------
                        %--------------------------------------------------------------------------
                        switch cnf_p_ACDC.range_index_method
                            %N_samples corresponds to non-zero padded ones
                            case 'conventional'
                                range_index=1:(N_samples*zp_fact_range_cnf);
                                delta_range=c_cst/(2.0*bw_ku_chd*zp_fact_range_cnf);
                            case 'resampling'
                                range_index=interp((1:N_samples),zp_fact_range_cnf);
                                delta_range=c_cst/(2.0*bw_ku_chd);
                        end
                        %% ----------------- LOAD THE LUTs FUNC_F0/F1 ------------------------------
                        %--------------------------------------------------------------------------
                        if (cnf_p_ACDC.lut_flag)
                            %-------------- load func_f0 --------------------------------------
                            load('./inputs/LUT_f0.mat','func_f0')
                            switch cnf_p_ACDC.power_wfm_model
                                case 'complete'
                                    load('./inputs/LUT_f1.mat','func_f1')
                                case 'simple'
                                    %
                                otherwise
                                    error('No valid power waveform model')
                            end
                        end
                        geo_param_buffer_conv=[];
                        geo_param_buffer_ACDC=[];
                    end
                    [L1B,geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,accumulated_amp_ACDC,i_surf_fitted,flag_exit]  = ACDC_stack_bis (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B,...
                        range_index,delta_range,i_surf_stacked,...
                        geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,...
                        accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,accumulated_amp_ACDC,i_surf_fitted,...
                        func_f0,func_f1,files);
                    
%                     Hs_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.Hs_conv;
%                     Hs_ACDC(i_surf_stacked)=L1B.ACDC.Hs;
%                     SSH_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.SSH_conv;
%                     SSH_ACDC(i_surf_stacked)=L1B.ACDC.SSH;
                end
                
                %% Writting routines
					
				if(writing_L1B)
					if i_surf_stacked==1
						files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
					end
					prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
				end
				if(writing_L1BS)
					if i_surf_stacked==1
						files.ncid_BS = netcdf.open(files.filename_netCDF_BS,'WRITE');
					end
					prepare_NetCDF_L1BS_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
				end 
                
% % %                 indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1BS_buffer(i_surf_stacked).burst_index));
% % %                 
% % %                 if(~isempty(indexes_to_remove))
% % %                     [L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
% % %                     burst_indexes_in_memory(indexes_to_remove)=[];
% % %                 end
% % %                 [L1BS_buffer(i_surf_stacked)] = empty_L1BS_struct(L1BS_buffer(i_surf_stacked));
                i_surf_stacked = i_surf_stacked +1;
                if(i_surf_stacked>N_surfs_loc_estimated)
                    disp('exit process');
                    break;
                end
                
            end
        end
        %stack the last surfaces 
        last_stack = i_surf_stacked;
        % Last stacks
        for i_surf_stacked = last_stack:N_total_surf_loc
                        
            progressbar([],[],[],i_surf_stacked/N_bursts/zp_fact_azimut_cnf*N_bursts_cycle_chd);
            [L1BS_buffer(i_surf_stacked)]       = stacking              (L1A_buffer,L1BS_buffer(i_surf_stacked));
            [L1BS_buffer(i_surf_stacked)]       = geometry_corrections  (L1BS_buffer(i_surf_stacked));
            
            
            
            %--------------
%             h=figure; imagesc(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr.'),1).').^2);
%             title(strcat('Stack #',num2str(i_surf_stacked)));
%             print('-dpng',['./results/plots/dedop/beam_geo_corr_surf',num2str(i_surf_stacked),'.png']);
%             close(h)
            %--------------
            
            
            
            [L1BS_buffer(i_surf_stacked)]       = range_transformation  (L1BS_buffer(i_surf_stacked));
            if extend_window_cnf
                [L1BS_buffer(i_surf_stacked)] = extend_stack(L1BS_buffer(i_surf_stacked));
            end
            [L1BS_buffer(i_surf_stacked),L1B]   = multilooking_new      (L1BS_buffer(i_surf_stacked));
            if include_wfms_aligned
                if i_surf_stacked== 1
                    %reference surface
                    win_delay_surf_ref=L1BS_buffer(i_surf_stacked).win_delay_surf;
                    alt_sat_ref=L1BS_buffer(i_surf_stacked).alt_sat;
                end
                [L1BS_buffer(i_surf_stacked),L1B]   = surface_win_delay_alignment (L1BS_buffer(i_surf_stacked),L1B,win_delay_surf_ref,alt_sat_ref);
            end
            [L1BS_buffer(i_surf_stacked),L1B]   = sigma0_scaling_factor (L1BS_buffer(i_surf_stacked),L1B);
            
            %added EM: 30.10.2016
            %% ---------------- ACDC APPLICATION ----------------------
            if ACDC_application_cnf
                if i_surf_stacked==1
                    %% ----------------- FITTING PARAMS INITIALIZE ----------------------------
                    %--------------------------------------------------------------------------
                    if cnf_p_ACDC.rou_flag
                        if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou cnf_p_ACDC.ini_Pu];
                        else
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_rou];
                        end
                    else
                        if strcmp(cnf_p_ACDC.multilook_option,'Cris_old')
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4 cnf_p_ACDC.ini_Pu];
                        else
                            fit_params_ini_conv      =   [cnf_p_ACDC.ini_Epoch cnf_p_ACDC.ini_Hs/4];
                        end
                    end
                    fit_params_ini_ACDC=[];
                    if cnf_p_ACDC.ini_Hs_rou_sliding_win_opt==1
                        accumulated_sigmaz_conv=0;
                        accumulated_SWH_ACDC=0;
                        accumulated_SSH_ACDC=0;
                        accumulated_amp_ACDC=0;
                    end
                    i_surf_fitted=0;
                    %% ----------------- INITIALIZE RANGE_INDEXATION --------------------------
                    %--------------------------------------------------------------------------
                    switch cnf_p_ACDC.range_index_method
                        %N_samples corresponds to non-zero padded ones
                        case 'conventional'
                            range_index=1:(N_samples*zp_fact_range_cnf);
                            delta_range=c_cst/(2.0*bw_ku_chd*zp_fact_range_cnf);
                        case 'resampling'
                            range_index=interp((1:N_samples),zp_fact_range_cnf);
                            delta_range=c_cst/(2.0*bw_ku_chd);
                    end
                    %% ----------------- LOAD THE LUTs FUNC_F0/F1 ------------------------------
                    %--------------------------------------------------------------------------
                    if (cnf_p_ACDC.lut_flag)
                        %-------------- load func_f0 --------------------------------------
                        load('./inputs/LUT_f0.mat','func_f0')
                        switch cnf_p_ACDC.power_wfm_model
                            case 'complete'
                                load('./inputs/LUT_f1.mat','func_f1')
                            case 'simple'
                                %
                            otherwise
                                error('No valid power waveform model')
                        end
                    end
                    geo_param_buffer_conv=[];
                    geo_param_buffer_ACDC=[];
                end
                [L1B,geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,accumulated_amp_ACDC,i_surf_fitted,flag_exit]  = ACDC_stack_bis (L1A_buffer,L1BS_buffer(i_surf_stacked),L1B,...
                    range_index,delta_range,i_surf_stacked,...
                    geo_param_buffer_conv,geo_param_buffer_ACDC,fit_params_ini_conv,fit_params_ini_ACDC,...
                    accumulated_sigmaz_conv,accumulated_SWH_ACDC,accumulated_SSH_ACDC,accumulated_amp_ACDC,i_surf_fitted,...
                    func_f0,func_f1,files);
                
%                 Hs_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.Hs_conv;
%                 Hs_ACDC(i_surf_stacked)=L1B.ACDC.Hs;
%                 SSH_preliminary_ACDC(i_surf_stacked)=L1B.ACDC.SSH_conv;
%                 SSH_ACDC(i_surf_stacked)=L1B.ACDC.SSH;
            end
            
            %% Writting routines
			if(writing_L1B)
				if i_surf_stacked==1
					files.ncid = netcdf.open(files.filename_netCDF,'WRITE');
				end
				prepare_NetCDF_L1B_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
			end
			if(writing_L1BS)
				if i_surf_stacked==1
					files.ncid_BS = netcdf.open(files.filename_netCDF_BS,'WRITE');
				end
				prepare_NetCDF_L1BS_S3 (files,L1A_buffer,L1BS_buffer(i_surf_stacked), L1B,i_surf_stacked);
            end 
            
        
% % %             indexes_to_remove = find(burst_indexes_in_memory(1:end) < min(L1BS_buffer(i_surf_stacked).burst_index));
% % %             
% % %             if(~isempty(indexes_to_remove))
% % %                 [L1A_buffer(burst_indexes_in_memory(indexes_to_remove))] = empty_L1A_struct(L1A_buffer(burst_indexes_in_memory(indexes_to_remove)));
% % %                 burst_indexes_in_memory(indexes_to_remove)=[];
% % %             end
% % %             [L1BS_buffer(i_surf_stacked)] = empty_L1BS_struct(L1BS_buffer(i_surf_stacked));
%             i_surf_stacked = i_surf_stacked +1;
            if(i_surf_stacked>N_surfs_loc_estimated)
                disp('exit process');
                break;
            end
            
            
        end
        
    end
    
%     h=figure; imagesc(wfms);
    %     print(h,['wfms_',num2str(i_fileL1A_input)],'-dpng');
    %     saveas (h,['wfms_',num2str(i_fileL1A_input),'.fig'])
%         h=figure; imagesc(wfms_wdcorr);
    %     print(h,['wfms_',num2str(i_fileL1A_input),'_wdcorr'],'-dpng');
    %     saveas (h,['wfms_wdcorr_',num2str(i_fileL1A_input),'.fig'])
    
    %Checking
    global mission
%     file_ws=strcat(strcat(files.resultPath,'data/'),mission,'_SR_1_SRA____',files.sph.product_info.product_id(20:20+30),'_temp_ws.mat')
%     save(file_ws,'win_delay_bursts','win_delay_bursts_h0cor2','lat_bursts','lon_bursts','instrument_range_correction_tx_rx_bursts','instrument_range_correction_rx_bursts')
   
    
    %% ----------- Resize the netCDF file ------------------------------
    if(writing_L1B)
        netcdf.close(files.ncid); % close the already open netcdf
        files.filename_netCDF_2remove=files.filename_netCDF;
        [files] = create_NetCDF_L1B_S3(files,i_surf_stacked-1);
        resize_NetCDF_L1B_S3(files,i_surf_stacked-1);
        delete(files.filename_netCDF_2remove);
        
    end
    if(writing_L1BS)
        netcdf.close(files.ncid_BS); % close the already open netcdf
        files.filename_netCDF_2remove_BS=files.filename_netCDF_BS;
        [files] = create_NetCDF_L1BS_S3(files,i_surf_stacked-1);
        resize_NetCDF_L1BS_S3(files,i_surf_stacked-1);
        delete(files.filename_netCDF_2remove_BS);
    end
    if(writing_L1B_pLRM)
        netcdf.close(files.ncid); % close the already open netcdf
        files.filename_netCDF_2remove=files.filename_netCDF;
        [files] = create_NetCDF_L1B_pLRM(files,i_RC-1);
        resize_NetCDF_L1B_PLRM(files,i_RC-1);
        delete(files.filename_netCDF_2remove);
    end
    if targz_option_active
        %remove the files DBL and HDR
        delete([files.inputPath filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name]); %DBL
        delete([files.inputPath strrep(filesBulk.inputFiles(filesBulk.indexFilesL1A(i_fileL1A_input)).name,'DBL','HDR')]); %HDR        
    end
    clear files;
    fclose all;
    close all;
    time = toc(t6);
    minutes_processing = floor(time/60);
    secs_processing = time - minutes_processing*60;
    strcat('Processing L1A: ',filename_L1A)
    disp(['Processing time for ',filename_L1A,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
    
    
    
    %% Retracker Processor
    if(steps(4))
        [L2]   = retracker(L1B);
    end
end