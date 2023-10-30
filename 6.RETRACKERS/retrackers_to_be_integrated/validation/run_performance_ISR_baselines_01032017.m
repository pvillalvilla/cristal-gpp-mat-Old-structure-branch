function [SSH_std,SWH_std,sigma0_std,amp_fit_std,...
    max_wvfm_std,ESA_SSH_std,ESA_sigma0_std,...
    STL_SSH_std,STL_sigma0_std,STL_SWH_std,...
    SSH_std_mean,SWH_std_mean,sigma0_std_mean,COR_std_mean,...
    ESA_SSH_std_mean,ESA_sigma0_std_mean,...
    STL_SSH_std_mean,STL_sigma0_std_mean,STL_SWH_std_mean,...
    SSH_mean,SWH_mean,sigma0_mean,...
    amp_fit_mean,max_wvfm_mean,hr_mean,nb_mean,...
    COR_mean,LAT_median,LON_median,...
    ESA_SSH_mean,ESA_sigma0_mean,...
    STL_SSH_mean,STL_sigma0_mean,STL_SWH_mean,...
    SSH,SWH,sigma0,COR,...
    ESA_L2_SSH_r1,ESA_L2_sigma0_r1,...
    STL_L2_SSH,STL_L2_sigma0,STL_L2_SWH,...
    amp_fit,max_wvfm,hr,nb,...
    ESA_hr_mean,ESA_hr,ESA_nb_mean,ESA_nb,...
    SSH_smooth,SWH_smooth,sigma0_smooth,COR_smooth,...
    SSH_RMSE,SWH_RMSE,sigma0_RMSE,COR_RMSE,...
    ESA_L2_SSH_r1_smooth,ESA_L2_sigma0_r1_smooth,ESA_SSH_RMSE,ESA_sigma0_RMSE,...
    STL_L2_SSH_smooth,STL_L2_SWH_smooth,STL_L2_sigma0_smooth,...
    STL_SSH_RMSE,STL_sigma0_RMSE,STL_SWH_RMSE,...
    SSH_mean_error_ESA_ISR,sigma0_mean_error_ESA_ISR,...
    SSH_mean_error_STL_ISR,SWH_mean_error_STL_ISR,sigma0_mean_error_STL_ISR,...
    num_surfaces]=run_performance_ISR_baselines_01032017(N_baselines,filesBulk,input_path_L2_ISR_bs,input_path_L2_ESA,input_path_L1_ESA,input_path_L2_STL,input_path_L1_ISR_bs,...
                                                                                                            name_bs,filename_mask_KML,...
                                                                                                            flag_outliers_removal,type_outliers_removal,outlier_percentil_low,outlier_percentil_high,IQR_times,hampel_wind,hampel_sigma,...
                                                                                                            smooth_param,win_size_detrending,...
                                                                                                            i_files_valid,plot_fits_downsampling,marker_ESA,color_ESA,marker_STL,color_STL,marker_bs,color_bs,sh_name_nc,...
                                                                                                            path_comparison_results,generate_plots,print_file,res_fig,file_ext,legend_fontsize,textbox_fontsize,annotation_box_active,...
                                                                                                            geo_mask,generate_kml,filter_land_surf_type)

%define some output variables as empty if not available the information
if cellfun(@isempty,input_path_L1_ISR_bs)
    amp_fit_std=[];
    max_wvfm_std=[];
    amp_fit_mean=[];
    max_wvfm_mean=[];
    hr_mean=[];
    nb_mean=[];
        amp_fit=[];
    max_wvfm=[];
    hr=[];
    nb=[];
end

if isempty(input_path_L1_ESA)
    ESA_hr_mean=[];
    ESA_hr=[];    
    ESA_nb_mean=[];
    ESA_nb=[];
end

if isempty(input_path_L2_ESA)
    ESA_SSH_std=[];
    ESA_sigma0_std=[];
    ESA_SSH_std_mean=[];
    ESA_sigma0_std_mean=[];
    ESA_SSH_mean=[];
    ESA_sigma0_mean=[];
    ESA_L2_SSH_r1=[];
    ESA_L2_sigma0_r1=[];
    SSH_mean_error_ESA_ISR=[];
    sigma0_mean_error_ESA_ISR=[];
    ESA_L2_SSH_r1_smooth=[];
    ESA_L2_sigma0_r1_smooth=[];
    
    ESA_SSH_RMSE=[];
    ESA_sigma0_RMSE=[];
end

if isempty(input_path_L2_STL)
    STL_SSH_std=[];
    STL_sigma0_std=[];
    STL_SWH_std=[];
    STL_SSH_std_mean=[];
    STL_sigma0_std_mean=[];
    STL_SWH_std_mean=[];
    STL_SSH_mean=[];
    STL_sigma0_mean=[];
    STL_SWH_mean=[];
    STL_L2_SSH=[];
    STL_L2_sigma0=[];
    STL_L2_SWH=[];
    SSH_mean_error_STL_ISR=[];
    SWH_mean_error_STL_ISR=[];
    sigma0_mean_error_STL_ISR=[];
    STL_L2_SSH_smooth=[];
    STL_L2_SWH_smooth=[];
    STL_L2_sigma0_smooth=[];
    
    STL_SSH_RMSE=[];
    STL_sigma0_RMSE=[];
    STL_SWH_RMSE=[];
end

if ~smooth_param
    SSH_smooth=[];
    SWH_smooth=[];
    sigma0_smooth=[];
    COR_smooth=[];    
    SSH_RMSE=[];
    SWH_RMSE=[];
    sigma0_RMSE=[];
    COR_RMSE=[];
end

                                      
for b=1:N_baselines
    
    %% ----------- Taking the files names -----------------------------------
    if b==1
        filename_L2_ISR=char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name);
        data_string=filename_L2_ISR(17:17+30);
        disp(num2str(i_files_valid));
        disp(data_string);
        
        %% --------------- READ ESA product -----------------------------------
        %try
        if ~isempty(input_path_L2_ESA)
            input_L2_ESA_Files   = dir(fullfile(char(input_path_L2_ESA),['*' data_string(1:15) '*.DBL']));
            if isempty(input_L2_ESA_Files)
                %add one second to initial time acquisition
                init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
                if isempty(input_L2_ESA_Files)
                    %add one second to initial time acquisition
                    end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                    input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
                    if isempty(input_L2_ESA_Files)
                        input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
                    end
                end
            end
            filename_L2_ESA=strcat(char(input_path_L2_ESA),input_L2_ESA_Files.name);
            [~,CS2]=Cryo_L2_read(filename_L2_ESA);
            s=size(CS2.MEA.surf_height_r1_20Hz);
            records_db=s(2);
            num_bursts_db=s(1);
            ESA_num_surfaces=records_db*num_bursts_db;
            %-------------- Geometry parameters ---------------------------------------
            ESA_lat_surf=reshape(CS2.MEA.LAT_20Hz,[1,ESA_num_surfaces]);
            ESA_lon_surf=reshape(CS2.MEA.LON_20Hz,[1,ESA_num_surfaces]);
            
            
            %-------------- Geophysical parameters ------------------------------------
            ESA_L2_SSH_r1=reshape(CS2.MEA.surf_height_r1_20Hz,[1,ESA_num_surfaces]);
            ESA_L2_sigma0_r1=reshape(CS2.MEA.backsc_sig_r1_20Hz,[1,ESA_num_surfaces]);
            clear CS2;
        end
        %---------- read the radial velocity information L1------------
        if ~isempty(input_path_L1_ESA)
            input_L1_ESA_Files   = dir(fullfile(char(input_path_L1_ESA),['*' data_string(1:15) '*.DBL']));
            if isempty(input_L1_ESA_Files)
                %add one second to initial time acquisition
                init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                input_L1_ESA_Files   = dir(fullfile(input_path_L1_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
                if isempty(input_L1_ESA_Files)
                    %add one second to initial time acquisition
                    end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                    input_L1_ESA_Files   = dir(fullfile(input_path_L1_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
                    if isempty(input_L1_ESA_Files)
                        input_L1_ESA_Files   = dir(fullfile(input_path_L1_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
                    end
                end
            end
            filename_L1_ESA=strcat(char(input_path_L1_ESA),input_L1_ESA_Files.name);
            [~,CS2]=Cryo_L1b_read(filename_L1_ESA);
            ESA_hr=reshape(CS2.GEO.H_rate,[1,ESA_num_surfaces]);
            ESA_nb=reshape(CS2.SAR.N_averaged_echoes,[1,ESA_num_surfaces]);
            ESA_flag_type_surface=reshape(ones(num_bursts_db,1)*(CS2.COR.surf_type),[1,ESA_num_surfaces]);
            
            
            clear CS;
        end
        
        %------------- Read Starlab file ------------------------------
        if ~isempty(input_path_L2_STL)
            input_L2_STL_Files   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
            filename_L2_STL=char(input_L2_STL_Files.name);
            filename_L2_STL=strcat(char(input_path_L2_STL),filename_L2_STL);
            STL_L2_SSH=double(ncread(filename_L2_STL,'ssh_corr').');
            STL_L2_sigma0=double(ncread(filename_L2_STL,'sigma0_L2').');
            STL_L2_SWH=double(ncread(filename_L2_STL,'swh').');
            
        end
        
    else
        input_L2_ISR_Files   = dir(fullfile(char(input_path_L2_ISR_bs(b)),['*' data_string(1:15) '*.nc']));
        filename_L2_ISR=char(input_L2_ISR_Files.name);
        clear input_L2_ISR_b2_Files;
    end
    filename_L2_ISR=strcat(char(input_path_L2_ISR_bs(b)),filename_L2_ISR);
    
    if ~cellfun(@isempty,input_path_L1_ISR_bs)
        input_L1_ISR_Files   = dir(fullfile(char(input_path_L1_ISR_bs(b)),['*' data_string(1:15) '*.nc']));
        filename_L1_ISR=char(input_L1_ISR_Files.name);
        filename_L1_ISR=strcat(char(input_path_L1_ISR_bs(b)),filename_L1_ISR);
    end
    
    %% ------------- Reading the information for each file ----------------
    % ********************* BASELINE -1 ***************************************
    if ~isempty(strfind(lower(char(name_bs(b))),'acdc'))
        %------------------- ACDC product -------------------------------------
        %---------------- Geometry variables --------------------------------------
        %Assuming they should have the same surfaces b1 and b2
        if b==1
            ISR_lat_surf=double(ncread(filename_L2_ISR,'lat_l1b_echo_sar_ku')).';
            ISR_lon_surf=double(ncread(filename_L2_ISR,'lon_l1b_echo_sar_ku')).';
        end
        
        %------------- Geophysical parameters ---------------------------------
        SSH_dumm=double(ncread(filename_L2_ISR,'ssh_ACDC_20_ku')).';
        SWH_dumm=double(ncread(filename_L2_ISR,'swh_ACDC_20_ku')).';
        sigma0_dumm=double(ncread(filename_L2_ISR,'sig0_ACDC_20_ku')).';
        COR_dumm=double(ncread(filename_L2_ISR,'Pearson_corr_ACDC_20_ku')).';
        
        %------------ read waveforms information ----------------------
        if ~cellfun(@isempty,input_path_L1_ISR_bs)
            ISD_i2q2_meas=double(ncread(filename_L2_ISR,'i2q2_meas_ku_ACDC_echo_sar_ku').');
            ISD_waveform_scale=ncread(filename_L2_ISR,'waveform_scale_factor_ACDC_echo_sar_ku');
            ISD_N_samples=length(ISD_i2q2_meas(1,:));
            %apply scaling factor to waveforms
            ISD_i2q2_meas=ISD_i2q2_meas.*repmat(ISD_waveform_scale,1,ISD_N_samples);
            max_wvfm_dumm=(max(ISD_i2q2_meas,[],2)).';
            clear ISD_i2q2_meas ISD_waveform_scale;
            Pu_dumm=10.^(double(ncread(filename_L2_ISR,'Pu_analytical_ACDC_20_ku')).'./10);
            amp_fit_dumm=Pu_dumm./max_wvfm_dumm;
            if b==1
                hr=double(ncread(filename_L1_ISR,'orb_alt_rate_l1b_echo_sar_ku')).'; %radial velocity
                ISD_flag_type_surface=double(ncread(filename_L1_ISR,'surf_type_l1b_echo_sar_ku').');
            end
            nb_dumm=double(ncread(filename_L1_ISR,'nb_stack_l1b_echo_sar_ku')).'; %number of beams
        end
        
    else
        %--------------------- L2 product -------------------------------------
        %---------------- Geometry variables ----------------------------------
        %Assuming they should have the same surfaces b1 and b2
        if b==1
            ISR_lat_surf=double(ncread(filename_L2_ISR,'lat_20_ku')).';
            ISR_lon_surf=double(ncread(filename_L2_ISR,'lon_20_ku')).';
        end
        
        %------------- Geophysical parameters ---------------------------------
        if isempty(strfind(lower(char(name_bs(b))),'thresh'))
            %analytical retracker SAR ocean model
            SSH_dumm=double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_analytical_SWH_MSSfixed_20_ku'))).';
            SWH_dumm=double(ncread(filename_L2_ISR,'swh_analytical_SWH_MSSfixed_20_ku')).';
            sigma0_dumm=double(ncread(filename_L2_ISR,'sig0_analytical_SWH_MSSfixed_20_ku')).';
            COR_dumm=double(ncread(filename_L2_ISR,'Pearson_corr_analytical_SWH_MSSfixed_20_ku')).';
        else
            %peak retracker
            SSH_dumm=double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_threshold_20_ku'))).';
            SWH_dumm=NaN(1,length(SSH_dumm));
            sigma0_dumm=double(ncread(filename_L2_ISR,'sig0_threshold_20_ku')).';
            COR_dumm=NaN(1,length(SSH_dumm));
        end
        
        
        %------------ read waveforms information ----------------------
        if ~cellfun(@isempty,input_path_L1_ISR_bs)
            ISD_i2q2_meas=double(ncread(filename_L1_ISR,'i2q2_meas_ku_l1b_echo_sar_ku').');
            ISD_waveform_scale=ncread(filename_L1_ISR,'waveform_scale_factor_l1b_echo_sar_ku');
            ISD_N_samples=length(ISD_i2q2_meas(1,:));
            ISD_i2q2_meas=squeeze(ISD_i2q2_meas).*repmat(ISD_waveform_scale,1,ISD_N_samples);
            
            % checking            
            if isempty(strfind(lower(char(name_bs(b))),'thresh'))
                max_wvfm_dumm=(max(squeeze(ISD_i2q2_meas),[],2)).';
                clear ISD_i2q2_meas ISD_waveform_scale;
                Pu_dumm=10.^(double(ncread(filename_L2_ISR,'Pu_analytical_SWH_MSSfixed_20_ku')).'./10);
                amp_fit_dumm=Pu_dumm./max_wvfm_dumm;
            else
                max_wvfm_dumm=NaN(1,length(SSH_dumm));
                Pu_dumm=NaN(1,length(SSH_dumm));
                amp_fit_dumm=NaN(1,length(SSH_dumm));
            end
            if b==1
                hr=double(ncread(filename_L1_ISR,'orb_alt_rate_l1b_echo_sar_ku')).'; %radial velocity   
                ISD_flag_type_surface=double(ncread(filename_L1_ISR,'surf_type_l1b_echo_sar_ku').');
            end
            
            if ~isempty(strfind(lower(char(name_bs(b))),'s-3 old'))
                ISD_nbeams_start_stop=double(ncread(filename_L1_ISR,'nb_stack_start_stop_l1b_echo_sar_ku')); %not exactly the same
                ISD_stack_mask_vector=double(ncread(filename_L1_ISR,'stack_mask_range_bin_l1b_echo_sar_ku')).'+1; %not exactly the same;
                for i_surf=1:length(ISD_nbeams_start_stop)
                    nb_dumm(i_surf)=length(find(ISD_stack_mask_vector(i_surf,1:ISD_nbeams_start_stop(i_surf))>1));
                end
                
                clear ISD_nbeams_start_stop ISD_stack_mask_vector;
            else
                nb_dumm=double(ncread(filename_L1_ISR,'nb_stack_l1b_echo_sar_ku')).'; %number of beams
            end
        end
    end
    ISD_num_surfaces=length(ISR_lat_surf);
    
    
    %% ----------------- Algin ESA and ISR arrays of parameters -------
    if b==1
        %---------- Filtering by mask ---------------------------------
        if ~isempty(filename_mask_KML)
            ISD_indexes_int=inpolygon(ISR_lon_surf,ISR_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));
            ISD_num_surfaces=length(find(ISD_indexes_int));
            
            ESA_indexes_int=inpolygon(ESA_lon_surf,ESA_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));
            ESA_num_surfaces=length(find(ESA_indexes_int));
            
        else
            ISD_indexes_int=logical(ones(1,ISD_num_surfaces));
            ESA_indexes_int=logical(ones(1,ESA_num_surfaces));
        end
        
        %------------- Filtering by the information in surface type flag---
        %not very high resolution but it can be trusted better
        if filter_land_surf_type
            if ~cellfun(@isempty,input_path_L1_ISR_bs)
                if ~isempty(input_path_L1_ESA)                    
                    ISD_indexes_int=ISD_indexes_int & (ISD_flag_type_surface==0);    
                    %ESA_indexes_int=ESA_indexes_int & (ESA_flag_type_surface==0);    
                end
            end
        end
        
        %consider the edges min and max latitudes (take the common parts)
        idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
        ESA_indexes_int=ESA_indexes_int & idx_not_lat_lon_zeros;
        
        
        if ~isempty(input_path_L2_ESA)
            if isempty(filename_mask_KML) && filter_land_surf_type==0
                %Forcing the number of surfaces of the ISD and ESA product be the same
                %assuming first surface not contemplated in ESA product
                idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
                indices_lat_lon_zeros=find(idx_not_lat_lon_zeros==0);
                if ISD_num_surfaces <=ESA_num_surfaces
                    ISD_indexes_int=ones(1,ISD_num_surfaces);
                    ISD_indexes_int(1)=0;
                    ESA_indexes_int=zeros(1,ESA_num_surfaces);
                    ESA_indexes_int(1:ISD_num_surfaces-1)=1;
                    for i_index=1:length(indices_lat_lon_zeros)
                        if (indices_lat_lon_zeros(i_index)+1)<=ISD_num_surfaces
                            ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
                        end
                    end
                    ISD_indexes_int=logical(ISD_indexes_int);
                    ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
                else
                    %the number of surfaces ESA limits
                    ISD_indexes_int=zeros(1,ISD_num_surfaces);
                    ISD_indexes_int(2:ESA_num_surfaces+1)=1;
                    ISD_indexes_int(1)=0;
                    ESA_indexes_int=ones(1,ESA_num_surfaces);
                    for i_index=1:length(indices_lat_lon_zeros)
                        if (indices_lat_lon_zeros(i_index)+1)<=ISD_num_surfaces
                            ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
                        end
                    end
                    ISD_indexes_int=logical(ISD_indexes_int);
                    ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
                end
            else
                if ~isempty(filename_mask_KML) && filter_land_surf_type==0
                    %Forcing same number of surfaces ESA/ISD for comparison: taking
                    %the first closest surface ESA to first ISD within the mask
                    %or index of surfaces filtered by surface land type
                    first_surface_within_mask=find(ISD_indexes_int,1,'first');
                    if first_surface_within_mask==1
                        ISD_indexes_int(1)=0; % not accounted for in the comparison
                        ISD_num_surfaces=ISD_num_surfaces-1;
                        ISD_indexes_int=logical(ISD_indexes_int);
                    end
                    
                    min_LAT_ref=min(ISR_lat_surf(ISD_indexes_int));
                    max_LAT_ref=max(ISR_lat_surf(ISD_indexes_int));
                    
                    clear ISD_indexes_int;
                    ISD_indexes_int=(ISR_lat_surf>=min_LAT_ref & ISR_lat_surf<=max_LAT_ref);
                    %search for the 1st surface closest to the first surface of
                    %isardSAT
                    [~,first_surface_ESA]=min(abs(ESA_lat_surf(ESA_indexes_int)-ISR_lat_surf(find(ISD_indexes_int,1,'first'))));
                    first_surface_ESA=first_surface_ESA+find(ESA_indexes_int,1,'first')-1;
                    ESA_indexes_int_bis=zeros(1,length(ESA_lat_surf));
                    ESA_indexes_int_bis(first_surface_ESA:end)=1;
                    ESA_indexes_int_bis=logical(ESA_indexes_int_bis);
                    ESA_indexes_int=ESA_indexes_int & ESA_indexes_int_bis;
                    
                    ISD_num_surfaces=length(find(ISD_indexes_int)==1);
                    ESA_num_surfaces=length(find(ESA_indexes_int)==1);
                    
                    %force the same number of surfaces for comparison purposes
                    if ISD_num_surfaces>ESA_num_surfaces
                        first_surface_within_mask=find(ISD_indexes_int,1,'first');
                        ISD_indexes_int=zeros(1,length(ISR_lat_surf));
                        ISD_indexes_int(first_surface_within_mask:first_surface_within_mask+ESA_num_surfaces-1)=1;
                        ISD_indexes_int=logical(ISD_indexes_int);
                    elseif ISD_num_surfaces<ESA_num_surfaces
                        first_surface_within_mask=find(ESA_indexes_int,1,'first');
                        ESA_indexes_int=zeros(1,length(ESA_lat_surf));
                        ESA_indexes_int(first_surface_within_mask:first_surface_within_mask+ISD_num_surfaces-1)=1;
                        ESA_indexes_int=logical(ESA_indexes_int);
                    end
                elseif isempty(filename_mask_KML) && filter_land_surf_type==1
                    first_surface_within_mask=find(ISD_indexes_int,1,'first');
                    if first_surface_within_mask==1
                        ISD_indexes_int(1)=0; % not accounted for in the comparison
                        ISD_num_surfaces=ISD_num_surfaces-1;
                        ISD_indexes_int=logical(ISD_indexes_int);
                    end
                    idx_int_ISD=find(ISD_indexes_int);
                    idx_int_ESA=find(ESA_indexes_int);
                    ISD_num_surfaces=length(idx_int_ISD);
                    ESA_num_surfaces=length(idx_int_ESA);
                    
                    ISD_indexes_int_bis=zeros(1,length(ISD_indexes_int));
                    ESA_indexes_int_bis=zeros(1,length(ESA_indexes_int));
                    
                    first_surface_ESA=find(ESA_indexes_int,1,'first');
                    if ISD_num_surfaces<=ESA_num_surfaces
                        for i_surf=1: ISD_num_surfaces
                            [~,dumm]=min(abs(ESA_lat_surf(ESA_indexes_int)-ISR_lat_surf(idx_int_ISD(i_surf))));
                            idx_int=first_surface_ESA+dumm(1)-1;
                            clear dumm;
                            ESA_indexes_int_bis(idx_int)=1;
                        end
                        ESA_indexes_int=logical(ESA_indexes_int_bis);
                    else                        
                        for i_surf=1: ESA_num_surfaces
                            [~,dumm]=min(abs(ESA_lat_surf(ESA_indexes_int)-ISR_lat_surf(idx_int_ISD(i_surf))));
                            idx_int=first_surface_ESA+dumm(1)-1;
                            clear dumm;
                            ESA_indexes_int_bis(idx_int)=1;
                            ISD_indexes_int_bis(idx_int_ISD(i_surf))=1;
                        end
                        ESA_indexes_int=logical(ESA_indexes_int_bis);
                        ISD_indexes_int=logical(ISD_indexes_int_bis);
                    end

                end
                
            end
            
            ISD_num_surfaces=length(find(ISD_indexes_int)==1);
            ESA_num_surfaces=length(find(ESA_indexes_int)==1);
            
%             % Checking pruposes
%             figure; geoshow('landareas.shp'); hold on;
%             geoshow(ISR_lat_surf(ISD_indexes_int),ISR_lon_surf(ISD_indexes_int),'DisplayType','Point','Marker','.','MarkerEdgeColor','red'); hold on;
%             geoshow(ESA_lat_surf(ESA_indexes_int),ESA_lon_surf(ESA_indexes_int),'DisplayType','Point','Marker','.','MarkerEdgeColor','blue'); hold on;
%             
            idx_int_ESA=find(ESA_indexes_int);
            idx_int_ISD=find(ISD_indexes_int);
            
            % ESA replacing values to meet same number surfaces as ISR
            %-------------- Geometry parameters ---------------------------------------
            ESA_lat_surf=ESA_lat_surf(ESA_indexes_int);
            ESA_lon_surf=ESA_lon_surf(ESA_indexes_int);
            if ~cellfun(@isempty,input_path_L1_ISR_bs)
                ESA_hr=ESA_hr(ESA_indexes_int);
                ESA_nb=ESA_nb(ESA_indexes_int);
            end
            
            
            %-------------- Geophysical parameters ------------------------------------
            ESA_L2_SSH_r1=ESA_L2_SSH_r1(ESA_indexes_int);
            ESA_L2_sigma0_r1=ESA_L2_sigma0_r1(ESA_indexes_int);
        end
        
        
        ISR_lat_surf=ISR_lat_surf(ISD_indexes_int);
        ISR_lon_surf=ISR_lon_surf(ISD_indexes_int);
        
        if generate_kml
            lla2kmlWaveforms_noimage(strcat(path_comparison_results,data_string,'_Geolocated_track.kml'),data_string,ISR_lat_surf,ISR_lon_surf,zeros(1,length(ISR_lon_surf)), '.');
        end
        
        num_surfaces=ISD_num_surfaces;
        
        %check STL data avialbel
        if ~isempty(input_path_L2_STL)
            STL_L2_SSH=STL_L2_SSH(ISD_indexes_int);
            STL_L2_sigma0=STL_L2_sigma0(ISD_indexes_int);
            STL_L2_SWH=STL_L2_SWH(ISD_indexes_int);
        end
        
    end
    
    %load the common surfaces of interest for ISD
    %------------- Geophysical parameters ---------------------------------
    SSH(b,:)=SSH_dumm(ISD_indexes_int);
    SWH(b,:)=SWH_dumm(ISD_indexes_int);
    sigma0(b,:)=sigma0_dumm(ISD_indexes_int);
    COR(b,:)=COR_dumm(ISD_indexes_int);
    if ~cellfun(@isempty,input_path_L1_ISR_bs)
        max_wvfm(b,:)=max_wvfm_dumm(ISD_indexes_int);
        Pu(b,:)=Pu_dumm(ISD_indexes_int);
        amp_fit(b,:)=amp_fit_dumm(ISD_indexes_int);
        
        if b==1
            hr=hr(ISD_indexes_int); %radial velocity
        end
        nb(b,:)=nb_dumm(ISD_indexes_int); %number of beams
        
        clear max_wvfm_dumm Pu_dumm amp_fit_dumm;
    end
    clear SSH_dumm SWH_dumm sigma0_dumm COR_dumm;
    
    
    
    
    %% ------------------- Outliers filtering -----------------------------
    %----------------------------------------------------------------------
    if flag_outliers_removal
        switch type_outliers_removal
            case 'percentiles'
                if ~isempty(input_path_L2_ESA) && b==1
                    %------------------- SSH ----------------------------------
                    [ESA_L2_SSH_r1,idx_outliers_SSH_ESA]=outliers_by_percentiles(ESA_L2_SSH_r1,outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_SSH_ESA=find(~idx_outliers_SSH_ESA);
                    %----------------- sigma0 ---------------------------------
                    [ESA_L2_sigma0_r1,idx_outliers_sigma0_ESA]=outliers_by_percentiles(ESA_L2_sigma0_r1,outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_sigma0_ESA=find(~idx_outliers_sigma0_ESA);
                end
                
                if ~isempty(input_path_L2_STL) && b==1
                    %------------------- SSH ----------------------------------
                    [STL_L2_SSH,idx_outliers_SSH_STL]=outliers_by_percentiles(STL_L2_SSH,outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_SSH_STL=find(~idx_outliers_SSH_STL);
                    %----------------- sigma0 ---------------------------------
                    [STL_L2_sigma0,idx_outliers_sigma0_STL]=outliers_by_percentiles(STL_L2_sigma0,outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_sigma0_STL=find(~idx_outliers_sigma0_STL);
                    %----------------- SWH ------------------------------------
                    [STL_L2_SWH,idx_outliers_SWH_STL]=outliers_by_percentiles(STL_L2_SWH,outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_SWH_STL=find(~idx_outliers_SWH_STL);
                end
                
                % compute the errors w.r.t fitting on the data using a smooth
                % function
                %------------------- SSH ----------------------------------
                [SSH(b,:),idx_outliers_SSH(b,:)]=outliers_by_percentiles(SSH(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                idx_nooutliers_SSH=find(~idx_outliers_SSH(b,:));
                %----------------- sigma0 ---------------------------------
                [sigma0(b,:),idx_outliers_sigma0(b,:)]=outliers_by_percentiles(sigma0(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                idx_nooutliers_sigma0=find(~idx_outliers_sigma0(b,:));
                %----------------- SWH ------------------------------------
                [SWH(b,:),idx_outliers_SWH(b,:)]=outliers_by_percentiles(SWH(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                idx_nooutliers_SWH=find(~idx_outliers_SWH(b,:));
                %----------------- COR ------------------------------------
                [COR(b,:),idx_outliers_COR(b,:)]=outliers_by_percentiles(COR(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                idx_nooutliers_COR=find(~idx_outliers_COR(b,:));
                
                if ~cellfun(@isempty,input_path_L1_ISR_bs)
                    %----------------- amplitude fit ----------------------
                    [amp_fit(b,:),idx_outliers_amp_fit(b,:)]=outliers_by_percentiles(amp_fit(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_amp_fit=find(~idx_outliers_amp_fit(b,:));
                    %---------------- max waveforms -----------------------
                    [max_wvfm(b,:),idx_outliers_max_wvfm(b,:)]=outliers_by_percentiles(max_wvfm(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                    idx_nooutliers_max_wvfm=find(~idx_outliers_max_wvfm(b,:));
                end
                
            case 'hampel'
                if ~isempty(input_path_L2_ESA) && b==1
                    %------------------ SSH ---------------------------
                    [~,idx_outliers_SSH_ESA] = hampel(ESA_lat_surf,ESA_L2_SSH_r1,hampel_wind,hampel_sigma);
                    ESA_L2_SSH_r1(idx_outliers_SSH_ESA)=NaN;
                    idx_nooutliers_SSH_ESA=find(~idx_outliers_SSH_ESA);
                    %----------------- sigma0 ---------------------------------
                    [~,idx_outliers_sigma0_ESA] = hampel(ESA_lat_surf,ESA_L2_sigma0_r1,hampel_wind,hampel_sigma);
                    ESA_L2_sigma0_r1(idx_outliers_sigma0_ESA)=NaN;
                    idx_nooutliers_sigma0_ESA=find(~idx_outliers_sigma0_ESA);
                end
                
                if ~isempty(input_path_L2_STL) && b==1
                    %------------------ SSH ---------------------------
                    [~,idx_outliers_SSH_STL] = hampel(ISR_lat_surf,STL_L2_SSH,hampel_wind,hampel_sigma);
                    STL_L2_SSH(idx_outliers_SSH_STL)=NaN;
                    idx_nooutliers_SSH_STL=find(~idx_outliers_SSH_STL);
                    %----------------- sigma0 ---------------------------------
                    [~,idx_outliers_sigma0_STL] = hampel(ISR_lat_surf,STL_L2_sigma0,hampel_wind,hampel_sigma);
                    STL_L2_sigma0(idx_outliers_sigma0_STL)=NaN;
                    idx_nooutliers_sigma0_STL=find(~idx_outliers_sigma0_STL);
                    %----------------- SWH ------------------------------------
                    [~,idx_outliers_SWH_STL] = hampel(ISR_lat_surf,STL_L2_SWH,hampel_wind,hampel_sigma);
                    STL_L2_SWH(idx_outliers_SWH_STL)=NaN;
                    idx_nooutliers_SWH_STL=find(~idx_outliers_SWH_STL);
                end
                
                %------------------ SSH ---------------------------
                [~,idx_outliers_SSH(b,:)] = hampel(ISR_lat_surf,SSH(b,:),hampel_wind,hampel_sigma);
                SSH(b,idx_outliers_SSH(b,:))=NaN;
                idx_nooutliers_SSH=find(~idx_outliers_SSH(b,:));
                %----------------- sigma0 ---------------------------------
                [~,idx_outliers_sigma0(b,:)] = hampel(ISR_lat_surf,sigma0(b,:),hampel_wind,hampel_sigma);
                sigma0(b,idx_outliers_sigma0(b,:))=NaN;
                idx_nooutliers_sigma0=find(~idx_outliers_sigma0(b,:));
                %----------------- SWH ------------------------------------
                [~,idx_outliers_SWH(b,:)] = hampel(ISR_lat_surf,SWH(b,:),hampel_wind,hampel_sigma);
                SWH(b,idx_outliers_SWH(b,:))=NaN;
                idx_nooutliers_SWH=find(~idx_outliers_SWH(b,:));
                %----------------- COR ------------------------------------
                [~,idx_outliers_COR(b,:)] = hampel(ISR_lat_surf,COR(b,:),hampel_wind,hampel_sigma);
                COR(b,idx_outliers_COR(b,:))=NaN;
                idx_nooutliers_COR=find(~idx_outliers_COR(b,:));
                
                if ~cellfun(@isempty,input_path_L1_ISR_bs)
                    %----------------- amplitude fit ----------------------
                    [~,idx_outliers_amp_fit(b,:)] = hampel(ISR_lat_surf,amp_fit(b,:),hampel_wind,hampel_sigma);
                    amp_fit(b,idx_outliers_amp_fit(b,:))=NaN;
                    idx_nooutliers_amp_fit=find(~idx_outliers_amp_fit(b,:));
                    
                    %----------------- max wvfm ---------------------------
                    [~,idx_outliers_max_wvfm(b,:)] = hampel(ISR_lat_surf,max_wvfm(b,:),hampel_wind,hampel_sigma);
                    max_wvfm(b,idx_outliers_max_wvfm(b,:))=NaN;
                    idx_nooutliers_max_wvfm=find(~idx_outliers_max_wvfm(b,:));
                    
                    %----------------- Pu ---------------------------
                    [~,idx_outliers_Pu(b,:)] = hampel(ISR_lat_surf,Pu(b,:),hampel_wind,hampel_sigma);
                    Pu(b,idx_outliers_Pu(b,:))=NaN;
                    idx_nooutliers_Pu=find(~idx_outliers_Pu(b,:));
                end
        end
    end
    
    %% --------------- Smoothing geophysical retrievals ---------------
    if smooth_param
        
        if ~isempty(input_path_L2_ESA) && b==1
            %----------------- SSH ------------------------------------
            ESA_L2_SSH_r1_smooth=smooth(ESA_L2_SSH_r1,win_size_detrending).';
            %----------------- SIGMA0 ---------------------------------
            ESA_L2_sigma0_r1_smooth=smooth(ESA_L2_sigma0_r1,win_size_detrending).';
        end
        
        if ~isempty(input_path_L2_STL) && b==1
            %----------------- SSH ------------------------------------
            STL_L2_SSH_smooth=smooth(STL_L2_SSH,win_size_detrending).';
            %----------------- SIGMA0 ---------------------------------
            STL_L2_sigma0_smooth=smooth(STL_L2_sigma0,win_size_detrending).';
            %----------------- SWH ---------------------------------
            STL_L2_SWH_smooth=smooth(STL_L2_SWH,win_size_detrending).';
        end
        
        %------------------ SSH ---------------------------
        SSH_smooth(b,:) = smooth(SSH(b,:),win_size_detrending).';
        %----------------- sigma0 ---------------------------------
        sigma0_smooth(b,:) = smooth(sigma0(b,:),win_size_detrending).';
        %----------------- SWH ------------------------------------
        SWH_smooth(b,:) = smooth(SWH(b,:),win_size_detrending).';
        %----------------- COR ------------------------------------
        COR_smooth(b,:) = smooth(COR(b,:),win_size_detrending).';
    end
    
    
    %% ------------------- Compute the std BLOCK-WISE ---------------------
    %----------------------------------------------------------------------
    num_boxes=floor(ISD_num_surfaces/win_size_detrending);
    for i_box=1:(num_boxes+1)
        init_sample=max([(i_box-1)*win_size_detrending+1,1]);
        last_sample=min([(i_box-1)*win_size_detrending+win_size_detrending,ISD_num_surfaces]);
        SSH_std(b,i_box)=nanstd(detrend(SSH(b,idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
        SWH_std(b,i_box)=nanstd(detrend(SWH(b,idx_nooutliers_SWH(idx_nooutliers_SWH>=init_sample & idx_nooutliers_SWH<=last_sample))));
        sigma0_std(b,i_box)=nanstd(detrend(sigma0(b,idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));
        COR_std(b,i_box)=nanstd(detrend(COR(b,idx_nooutliers_COR(idx_nooutliers_COR>=init_sample & idx_nooutliers_COR<=last_sample))));
        
        
        
        SSH_mean(b,i_box)=nanmean((SSH(b,idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
        SWH_mean(b,i_box)=nanmean((SWH(b,idx_nooutliers_SWH(idx_nooutliers_SWH>=init_sample & idx_nooutliers_SWH<=last_sample))));
        sigma0_mean(b,i_box)=nanmean((sigma0(b,idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));
        COR_mean(b,i_box)=nanmean((COR(b,idx_nooutliers_COR(idx_nooutliers_COR>=init_sample & idx_nooutliers_COR<=last_sample))));
        
        if ~cellfun(@isempty,input_path_L1_ISR_bs)
            amp_fit_std(b,i_box)=nanstd(detrend(amp_fit(b,idx_nooutliers_amp_fit(idx_nooutliers_amp_fit>=init_sample & idx_nooutliers_amp_fit<=last_sample))));
            max_wvfm_std(b,i_box)=nanstd(detrend(max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample))));
            max_wvfm_dB_std(b,i_box)=nanstd(detrend(10*log10(max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample)))));
            Pu_std(b,i_box)=nanstd(detrend(Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample))));
            Pu_dB_std(b,i_box)=nanstd(detrend(10*log10(Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample)))));
            Pu_mean(b,i_box)=nanmean((Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample))));
            Pu_dB_mean(b,i_box)=nanmean(10*log10(Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample))));
            
            
            
            amp_fit_mean(b,i_box)=nanmean((amp_fit(b,idx_nooutliers_amp_fit(idx_nooutliers_amp_fit>=init_sample & idx_nooutliers_amp_fit<=last_sample))));
            max_wvfm_mean(b,i_box)=nanmean((max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample))));
            max_wvfm_dB_mean(b,i_box)=nanmean(10*log10(max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample))));
            
            if b==1
                hr_mean(i_box)=nanmean((hr(init_sample:last_sample)));
            end
            nb_mean(b,i_box)=nanmean((nb(b,init_sample:last_sample)));
        end
        
        if b==1
            LAT_median(i_box)=nanmedian((ISR_lat_surf(init_sample:last_sample)));
            LON_median(i_box)=nanmedian((ISR_lon_surf(init_sample:last_sample)));
        end
        
        
        if ~isempty(input_path_L2_ESA) && b==1
            ESA_SSH_std(i_box)=nanstd(detrend(ESA_L2_SSH_r1(idx_nooutliers_SSH_ESA(idx_nooutliers_SSH_ESA>=init_sample & idx_nooutliers_SSH_ESA<=last_sample))));
            ESA_sigma0_std(i_box)=nanstd(detrend(ESA_L2_sigma0_r1(idx_nooutliers_sigma0_ESA(idx_nooutliers_sigma0_ESA>=init_sample & idx_nooutliers_sigma0_ESA<=last_sample))));
            ESA_SSH_mean(i_box)=nanmean((ESA_L2_SSH_r1(idx_nooutliers_SSH_ESA(idx_nooutliers_SSH_ESA>=init_sample & idx_nooutliers_SSH_ESA<=last_sample))));
            ESA_sigma0_mean(i_box)=nanmean((ESA_L2_sigma0_r1(idx_nooutliers_sigma0_ESA(idx_nooutliers_sigma0_ESA>=init_sample & idx_nooutliers_sigma0_ESA<=last_sample))));
        end
        
        
        
        if ~isempty(input_path_L2_STL) && b==1
            STL_SSH_std(i_box)=nanstd(detrend(STL_L2_SSH(idx_nooutliers_SSH_STL(idx_nooutliers_SSH_STL>=init_sample & idx_nooutliers_SSH_STL<=last_sample))));
            STL_sigma0_std(i_box)=nanstd(detrend(STL_L2_sigma0(idx_nooutliers_sigma0_STL(idx_nooutliers_sigma0_STL>=init_sample & idx_nooutliers_sigma0_STL<=last_sample))));
            STL_SWH_std(i_box)=nanstd(detrend(STL_L2_SWH(idx_nooutliers_SWH_STL(idx_nooutliers_SWH_STL>=init_sample & idx_nooutliers_SWH_STL<=last_sample))));
            
            STL_SSH_mean(i_box)=nanmean((STL_L2_SSH(idx_nooutliers_SSH_STL(idx_nooutliers_SSH_STL>=init_sample & idx_nooutliers_SSH_STL<=last_sample))));
            STL_sigma0_mean(i_box)=nanmean((STL_L2_sigma0(idx_nooutliers_sigma0_STL(idx_nooutliers_sigma0_STL>=init_sample & idx_nooutliers_sigma0_STL<=last_sample))));
            STL_SWH_mean(i_box)=nanmean((STL_L2_SWH(idx_nooutliers_SWH_STL(idx_nooutliers_SWH_STL>=init_sample & idx_nooutliers_SWH_STL<=last_sample))));
        end
        
        if ~isempty(input_path_L1_ESA) && b==1
            ESA_hr_mean(i_box)=nanmean((ESA_hr(init_sample:last_sample)));
            ESA_nb_mean(i_box)=nanmean((ESA_nb(init_sample:last_sample)));
        end
        
        
        
    end
    
    %% --------- Compute the RMSE FITTING ---------------------------------
    if smooth_param
        
        if ~isempty(input_path_L2_ESA) && b==1
            ESA_SSH_RMSE=sqrt(nanmean((ESA_L2_SSH_r1-ESA_L2_SSH_r1_smooth).^2));
            ESA_sigma0_RMSE=sqrt(nanmean((ESA_L2_sigma0_r1-ESA_L2_sigma0_r1_smooth).^2));
        end
        
        if ~isempty(input_path_L2_STL) && b==1
            STL_SSH_RMSE=sqrt(nanmean((STL_L2_SSH-STL_L2_SSH_smooth).^2));
            STL_sigma0_RMSE=sqrt(nanmean((STL_L2_sigma0-STL_L2_sigma0_smooth).^2));
            STL_SWH_RMSE=sqrt(nanmean((STL_L2_SWH-STL_L2_SWH_smooth).^2));
        end
        
        SSH_RMSE(b,1)=sqrt(nanmean((SSH(b,:)-SSH_smooth(b,:)).^2));
        SWH_RMSE(b,1)=sqrt(nanmean((SWH(b,:)-SWH_smooth(b,:)).^2));
        sigma0_RMSE(b,1)=sqrt(nanmean((sigma0(b,:)-sigma0_smooth(b,:)).^2));
        COR_RMSE(b,1)=sqrt(nanmean((COR(b,:)-COR_smooth(b,:)).^2));
    end
    
    
    %% --------- Compute the mean std over track equivalent to RMSE--------
    if ~isempty(input_path_L2_ESA) && b==1
        ESA_SSH_std_mean=nanmean(ESA_SSH_std);
        ESA_sigma0_std_mean=nanmean(ESA_sigma0_std);
    end
    
    if ~isempty(input_path_L2_STL) && b==1
        STL_SSH_std_mean=nanmean(STL_SSH_std);
        STL_sigma0_std_mean=nanmean(STL_sigma0_std);
        STL_SWH_std_mean=nanmean(STL_SWH_std);
    end
    
    SSH_std_mean(b,1)=nanmean(SSH_std(b,:));
    SWH_std_mean(b,1)=nanmean(SWH_std(b,:));
    sigma0_std_mean(b,1)=nanmean(sigma0_std(b,:));
    COR_std_mean(b,1)=nanmean(COR_std(b,:));
    
    
    
    
    %% ------------------ Compute the mean errors ESA/SLT-ISD--------------
    if ~isempty(input_path_L2_ESA)
        SSH_mean_error_ESA_ISR(b,1)=nanmean(ESA_SSH_mean-SSH_mean(b,:));
        sigma0_mean_error_ESA_ISR(b,1)=nanmean(ESA_sigma0_mean-sigma0_mean(b,:));
    end
    
    if ~isempty(input_path_L2_STL)
        SSH_mean_error_STL_ISR(b,1)=nanmean(STL_SSH_mean-SSH_mean(b,:));
        SWH_mean_error_STL_ISR(b,1)=nanmean(STL_SWH_mean-SWH_mean(b,:));
        sigma0_mean_error_STL_ISR(b,1)=nanmean(STL_sigma0_mean-sigma0_mean(b,:));
    end
    clear idx_nooutliers_SSH idx_nooutliers_SWH idx_nooutliers_sigma0 idx_nooutliers_COR;
    clear idx_nooutliers_SSH_ESA idx_nooutliers_sigma0_ESA;
    clear idx_nooutliers_SSH_STL idx_nooutliers_sigma0_STL idx_nooutliers_SWH_STL;
    clear idx_nooutliers_amp_fit idx_nooutliers_max_wvfm idx_nooutliers_Pu;
end %loop over the channels

clear ISD_i2q2_meas;
%% ---------------- PLOTING RESULTS -----------------------------------
%----------------------------------------------------------------------
%--------------------- SSH --------------------------------------------
if generate_plots==1
    if  mod(i_files_valid,plot_fits_downsampling)==0 || i_files_valid==1
        f1=figure;
        legend_text={''};
        text_in_textbox={''};
        text_errors_ESA_ISR={''};
        text_errors_STL_ISR={''};
        if ~isempty(input_path_L2_ESA)
            plot(ESA_lat_surf,ESA_L2_SSH_r1,'Marker',marker_ESA,'Color',color_ESA);
            hold on;
            legend_text=[legend_text,'ESA'];
            text_in_textbox=[text_in_textbox,...
                strcat({'ESA: RMSE [m]: '},num2str(ESA_SSH_RMSE,'%.4g'),{' , std [m]: '},num2str(ESA_SSH_std_mean,'%.4g'))];
        end
        
        if ~isempty(input_path_L2_STL)
            
            plot(ISR_lat_surf,STL_L2_SSH,'Marker',marker_STL,'Color',color_STL);
            hold on;
            legend_text=[legend_text,'STL'];
            text_in_textbox=[text_in_textbox,...
                strcat({'STL: RMSE [m]: '},num2str(STL_SSH_RMSE,'%.4g'),{' , std [m]: '},num2str(STL_SSH_std_mean,'%.4g'))];
        end
        
        for b=1:N_baselines
            plot(ISR_lat_surf,SSH(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
            hold on;
            legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
            text_in_textbox=[text_in_textbox,...
                strcat({'ISR '},name_bs(b),{': RMSE [m]: '},num2str(SSH_RMSE(b,1),'%.4g'),{' , std [m]: '},num2str(SSH_std_mean(b,1),'%.4g'))];
            if ~isempty(input_path_L2_ESA)
                text_errors_ESA_ISR=[text_errors_ESA_ISR,...
                    strcat({'Mean error ESA-ISR '},name_bs(b),{' [m]: '},num2str(SSH_mean_error_ESA_ISR(b,1),'%.4g'))];
            end
            if ~isempty(input_path_L2_STL)
                text_errors_STL_ISR=[text_errors_STL_ISR,...
                    strcat({'Mean error STL-ISR '},name_bs(b),{' [m]: '},num2str(SSH_mean_error_STL_ISR(b,1),'%.4g'))];
            end
        end
        title(strcat(upper(sh_name_nc),': Baselines comparison'));
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
        pos_leg=get(leg,'Position');
        xlabel('Latitude [deg.]'); ylabel(strcat(upper(sh_name_nc),' [m]'));
        text_in_textbox=[text_in_textbox,text_errors_ESA_ISR,text_errors_STL_ISR];
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        
        if annotation_box_active
            pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),pos_leg(4)],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
        end
        print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SSH',file_ext))
        save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SSH','.mat'),...
            'input_path_L2_ESA','ESA_lat_surf','ESA_L2_SSH_r1','marker_ESA','color_ESA',...
            'input_path_L2_STL','ISR_lat_surf','STL_L2_SSH','marker_STL','color_STL',...
            'N_baselines','SSH','marker_bs','color_bs',...
            'name_bs','sh_name_nc','legend_text','legend_fontsize','text_in_textbox','textbox_fontsize',...
            'path_comparison_results','data_string','file_ext');
        close(f1)
        %--------------------- SWH --------------------------------------------
        f1=figure;
        legend_text={''};
        text_in_textbox={''};
        text_errors_STL_ISR={''};
        if ~isempty(input_path_L2_STL)
            plot(ISR_lat_surf,STL_L2_SWH,'Marker',marker_STL,'Color',color_STL);
            hold on;
            legend_text=[legend_text,'STL'];
            text_in_textbox=[text_in_textbox,...
                strcat({'STL: RMSE [m]: '},num2str(STL_SWH_RMSE,'%.4g'),{' , std [m]: '},num2str(STL_SWH_std_mean,'%.4g'))];
        end
        for b=1:N_baselines
            plot(ISR_lat_surf,SWH(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
            hold on;
            legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
            text_in_textbox=[text_in_textbox,...
                strcat({'ISR '},name_bs(b),{': RMSE [m]: '},num2str(SWH_RMSE(b,1),'%.4g'),{' , std [m]: '},num2str(SWH_std_mean(b,1),'%.4g'))];
            if ~isempty(input_path_L2_STL)
                text_errors_STL_ISR=[text_errors_STL_ISR,...
                    strcat({'Mean error STL-ISR '},name_bs(b),{' [m]: '},num2str(SWH_mean_error_STL_ISR(b,1),'%.4g'))];
            end
        end
        title(strcat('SWH (H_s): Baselines comparison'));
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
        pos_leg=get(leg,'Position');
        xlabel('Latitude [deg.]'); ylabel('SWH [m]');
        text_in_textbox=[text_in_textbox,text_errors_STL_ISR];
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        if annotation_box_active
            pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
        end
        print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SWH',file_ext))
        save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SWH','.mat'),...
            'input_path_L2_STL','ISR_lat_surf','STL_L2_SWH','marker_STL','color_STL',...
            'N_baselines','SWH','marker_bs','color_bs',...
            'name_bs','sh_name_nc','legend_text','legend_fontsize','text_in_textbox','textbox_fontsize',...
            'path_comparison_results','data_string','file_ext');
        close(f1);
        
        %--------------------- sigma0 --------------------------------------------
        f1=figure;
        legend_text={''};
        text_in_textbox={''};
        text_errors_ESA_ISR={''};
        text_errors_STL_ISR={''};
        if ~isempty(input_path_L2_ESA)
            
            plot(ESA_lat_surf,ESA_L2_sigma0_r1,'Marker',marker_ESA,'Color',color_ESA);
            hold on;
            legend_text=[legend_text,'ESA'];
            text_in_textbox=[text_in_textbox,...
                strcat({'ESA: RMSE [dB]: '},num2str(ESA_sigma0_RMSE,'%.4g'),{' , std [dB]: '},num2str(ESA_sigma0_std_mean,'%.4g'))];
        end
        
        if ~isempty(input_path_L2_STL)
            
            plot(ISR_lat_surf,STL_L2_sigma0,'Marker',marker_STL,'Color',color_STL);
            hold on;
            legend_text=[legend_text,'STL'];
            text_in_textbox=[text_in_textbox,...
                strcat({'STL: RMSE [dB]: '},num2str(STL_sigma0_RMSE,'%.4g'),{' , std [dB]: '},num2str(STL_sigma0_std_mean,'%.4g'))];
        end
        
        for b=1:N_baselines
            plot(ISR_lat_surf,sigma0(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
            hold on;
            legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
            text_in_textbox=[text_in_textbox,...
                strcat({'ISR '},name_bs(b),{': RMSE [dB]: '},num2str(sigma0_RMSE(b,1),'%.4g'),{' , std [dB]: '},num2str(sigma0_std_mean(b,1),'%.4g'))];
            
            if ~isempty(input_path_L2_ESA)
                text_errors_ESA_ISR=[text_errors_ESA_ISR,...
                    strcat({'Mean error ESA-ISR '},name_bs(b),{' [dB]: '},num2str(sigma0_mean_error_ESA_ISR(b,1),'%.4g'))];
            end
            if ~isempty(input_path_L2_STL)
                text_errors_STL_ISR=[text_errors_STL_ISR,...
                    strcat({'Mean error STL-ISR '},name_bs(b),{' [dB]: '},num2str(sigma0_mean_error_STL_ISR(b,1),'%.4g'))];
            end
            
        end
        title(strcat('\sigma^0: Baselines comparison'));
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
        pos_leg=get(leg,'Position');
        xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]');
        text_in_textbox=[text_in_textbox,text_errors_ESA_ISR,text_errors_STL_ISR];
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        
        if annotation_box_active
            pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
        end
        print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_sigma0',file_ext))
        save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_sigma0','.mat'),...
            'input_path_L2_ESA','ESA_lat_surf','ESA_L2_sigma0_r1','marker_ESA','color_ESA',...
            'input_path_L2_STL','ISR_lat_surf','STL_L2_sigma0','marker_STL','color_STL',...
            'N_baselines','sigma0','marker_bs','color_bs',...
            'name_bs','sh_name_nc','legend_text','legend_fontsize','text_in_textbox','textbox_fontsize',...
            'path_comparison_results','data_string','file_ext');
        close(f1)
        
        
        %--------------------- nb beams --------------------------------------------
        f1=figure;
        legend_text={''};
        text_in_textbox={''};
        text_errors_ESA_ISR={''};
        text_errors_STL_ISR={''};
        if ~isempty(input_path_L2_ESA)
            plot(ESA_lat_surf,ESA_nb,'Marker',marker_ESA,'Color',color_ESA);
            hold on;
            legend_text=[legend_text,'ESA'];
        end
        
        for b=1:N_baselines
            plot(ISR_lat_surf,nb(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
            hold on;
            legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
            
            if ~isempty(input_path_L2_ESA)
                text_errors_ESA_ISR=[text_errors_ESA_ISR,...
                    strcat({'Mean error ESA-ISR '},name_bs(b),{' : '},num2str(nanmean(ESA_nb-nb(b,:)),'%.4g'))];
            end
            
        end
        title(strcat('Contributing Beams: Baselines comparison'));
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
        pos_leg=get(leg,'Position');
        xlabel('Latitude [deg.]'); ylabel('n_b');
        text_in_textbox=[text_in_textbox,text_errors_ESA_ISR];
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        
        if annotation_box_active
            pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
        end
        print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_nb',file_ext))
        save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_nb','.mat'),...
            'input_path_L2_ESA','ESA_lat_surf','ESA_nb','marker_ESA','color_ESA',...
            'N_baselines','nb','marker_bs','color_bs',...
            'name_bs','sh_name_nc','legend_text','legend_fontsize','text_in_textbox','textbox_fontsize',...
            'path_comparison_results','data_string','file_ext');
        close(f1)
        
        %--------------------- COR --------------------------------------------
        f1=figure;
        legend_text={''};
        text_in_textbox={''};
        
        for b=1:N_baselines
            plot(ISR_lat_surf,COR(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
            hold on;
            legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
            text_in_textbox=[text_in_textbox,...
                strcat({'ISR '},name_bs(b),{' --> RMSE [%]: '},num2str(COR_RMSE(b,1),'%.4g'),{' , std [%]: '},num2str(COR_std_mean(b,1),'%.4g'))];
        end
        title(strcat('\rho (Pearson correlation coeff.): Baselines comparison'));
        leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
        pos_leg=get(leg,'Position');
        xlabel('Latitude [deg.]'); ylabel('\rho [%]');
        text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
        if annotation_box_active
            pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
        end
        print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_COR',file_ext))
        save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_COR','.mat'),...
            'N_baselines','COR','marker_bs','color_bs',...
            'name_bs','sh_name_nc','legend_text','legend_fontsize','text_in_textbox','textbox_fontsize',...
            'path_comparison_results','data_string','file_ext');
        close(f1);
    end
end




end