function write_L2_product(file,out,cnf_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and geophysical correctiosn to be applied to the
% retracked range based on the info available in the L1B product and
% depending on the surface being observed
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 17/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -file    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -out = structure of the output data for the L2 product
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       
% RESTRICTIONS: 
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 

global alpha_gr_chd alpha_ga_chd;
%% ------------- NAME DEFINITION ------------------------------------------
%-------------- Define the output name for the output file ----------------
date_creation = datestr(now, '_yyyymmddTHHMMSS_');
aux=strsplit(date_creation,'_');
out.date_creation=aux(2);
clear aux;
switch cnf_p.mission    
    case 'CS2'
    %--------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                name_L2_product=strcat(strrep(file.filename_L1B_nopath,'1B','2_'),date_creation,'isd');
            case 'ISD'
                %following a la Sentinel-3 format as per SEOMs projects
                aux=file.filename_L1B_nopath(1:47);
                aux(10:15)=strcat(cnf_p.Product_type,'___');
                name_L2_product=strcat(strrep(aux(1:47),'_1_','_2_'),date_creation,'isd');
        end
    case 'S3'
    %--------------------- Sentinel-3 -----------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                % TBD
                %[data]=readL1B_S3_ESA(filename_L1B);
            case 'ISD'
                  name_L2_product=strcat(strrep(file.filename_L1B_nopath,'_1_SRA',strcat('_2_',cnf_p.Product_type)));
%                 aux=file.filename_L1B_nopath(1:47);
%                 aux(10:15)=cnf_p.Product_type+'___';
%                 name_L2_product=strcat(ststrrep(aux(1:47),'_1_','_2_'),date_creation,'isd');
        end
    case 'S6'
    %--------------------- Sentinel-6 -----------------------------------        
        switch cnf_p.L1proc
            case 'ISD'
                %needs to be defined as using the final file naming convention
                name_L2_product=strcat(strrep(file.filename_L1B_nopath,'RAW','L2_'),date_creation,'isd');
        end 
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end

if cnf_p.optional_ext_file_flag
    name_L2_product=strrep(name_L2_product,'isd',strcat(cnf_p.file_ext_string,'_isd'));
end

%% ---------------------- GENERATION OF OUTPUT PRODUCT --------------------
switch cnf_p.output_product_format
    case 'mat' % Matlab output
        % .mat file
        save([file.resultPath 'data/' name_L2_product '.mat'],'out');        
    case 'nc'  % NetCDF output product
        prepare_NetCDF_L2([file.resultPath 'data/' name_L2_product '.nc'],out,cnf_p); %can have a switch to define the variables if mission is S3 or S6
    otherwise
        error('The type of product output is not supported')        
end

%% Include the processing options product
fid_proc_opts=fopen([file.resultPath 'data/' name_L2_product,'_proc_opt.txt'],'w');
fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
fprintf(fid_proc_opts,'--------------- General Parameters & Options retrackers ----------------------------\n');
fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
fprintf(fid_proc_opts,'%s\n',char(strcat('Filtering Geographical mask',{': '},num2str(cnf_p.mask_ROI_flag))));
fprintf(fid_proc_opts,'%s\n',char(strcat('Filtering Number looks stack',{': '},num2str(cnf_p.mask_looks_flag))));
if cnf_p.mask_looks_flag
    fprintf(fid_proc_opts,'%s\n',char(strcat('Minimum Number looks stack (filtering)',{': '},num2str(cnf_p.Neff_thres))));
end
fprintf(fid_proc_opts,'%s\n',char(strcat('IFmask',{': '},num2str(cnf_p.IFmask_N))));
fprintf(fid_proc_opts,'%s\n',char(strcat('Geophysical corrections applied',{': '},num2str(cnf_p.geo_corr_application_flag))));
fprintf(fid_proc_opts,'%s\n',char(strcat('force_geocorr_surf_type',{': '},num2str(cnf_p.force_geocorr_surf_type))));
if cnf_p.force_geocorr_surf_type
    fprintf(fid_proc_opts,'%s\n',char(strcat('product_type_surface',{': '},cnf_p.product_type_surface)));
end

for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'--------------- Parameters & Options of Analytical retracker -----------------------\n');
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'%s\n',char(strcat('Use zeros in multilooking',{': '},num2str(cnf_p.use_zeros_cnf))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Zero Padding in range',{': '},num2str(cnf_p.ZP))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('window_type_a',{': '},cnf_p.window_type_a)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('window_type_r',{': '},cnf_p.window_type_r)));            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Range PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_gr_chd))))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_ga_chd))))));
            if cnf_p.fit_noise
                fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold noise',{': '},num2str(cnf_p.threshold_noise))));
            else
                fprintf(fid_proc_opts,'%s\n',char(strcat('First noise sample',{': '},num2str(cnf_p.Thn_w_first))));
                fprintf(fid_proc_opts,'%s\n',char(strcat('Width noise window',{': '},num2str(cnf_p.Thn_w_width))));
            end
            
%             disp(strcat('Sign pitch',{': '},num2str(cnf_p.sign_pitch)))
%             disp(strcat('Sign roll',{': '},num2str(cnf_p.sign_roll)))

            fprintf(fid_proc_opts,'%s\n',char(strcat('Indexation method',{': '},(cnf_p.looks_index_method))));
            switch cnf_p.looks_index_method
                case 'Look_angle'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Look angle method',{': '},cnf_p.look_ang_method)));
                case 'Doppler_freq'
                    fprintf(fid_proc_opts,'%s\n',char(strcat('Doppler freq. method',{': '},cnf_p.fd_method)));
            end
            fprintf(fid_proc_opts,'%s\n',char(strcat('Power wvfm model',{': '},cnf_p.power_wfm_model)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('LUT flag',{': '},num2str(cnf_p.lut_flag))));
            if cnf_p.pre_processing
                fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold (leading edge pre-processing)',{': '},num2str(cnf_p.percent_leading_edge))));
            end
            fprintf(fid_proc_opts,'%s\n',char(strcat('Two step fitting',{': '},num2str(cnf_p.two_step_fitting))));
            if cnf_p.two_step_fitting
                fprintf(fid_proc_opts,'%s\n',char(strcat('Two step fitting COR threshold',{': '},num2str(cnf_p.two_step_fitting_COR_threshold_rou))));
            end
            fprintf(fid_proc_opts,'%s\n',char(strcat('initial_param_fit_feedback_flag',{': '},num2str(cnf_p.initial_param_fit_feedback_flag))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('ini_Hs_rou_sliding_win_opt',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('ini_Hs_rou_sliding_win_opt_discard_std_threshold',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('ini_Hs_rou_sliding_win_size',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_size))));
            
            fprintf(fid_proc_opts,'%s\n',char(strcat('Range indexation method',{': '},cnf_p.range_index_method)));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Fitting method',{': '},cnf_p.fitting_fun_type)));
            opti_cell=strcat(fieldnames(cnf_p.fitting_options),{': '},cellfun(@num2str,struct2cell(cnf_p.fitting_options),'UniformOutput',false));
            for i_field=1:length(opti_cell)
                fprintf(fid_proc_opts,'%s\n',char(opti_cell(i_field)));
            end
            if ~isempty(cnf_p.fitting_options_lb)
                fprintf(fid_proc_opts,'%s\n',char(strcat('Lower bounds fitting: ',sprintf('%d',cnf_p.fitting_options_lb))));
            end
            if ~isempty(cnf_p.fitting_options_ub)
                fprintf(fid_proc_opts,'%s\n',char(strcat('Upper bounds fitting: ',sprintf('%d',cnf_p.fitting_options_ub))));
            end 
            
        case 'THRESHOLD'
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'--------------- Parameters & Options of Threshold retracker ------------------------\n');
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'%s\n',char(strcat('Threshold value',{': '},num2str(cnf_p.th_retracker.percentage_peak))));

        case 'OCOG'
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'--------------- Parameters & Options of OCOG retracker -----------------------------\n');
            fprintf(fid_proc_opts,'------------------------------------------------------------------------------------\n');
            fprintf(fid_proc_opts,'%s\n',char(strcat('Percentage OCOG Amplitude',{': '},num2str(cnf_p.OCOG_retracker.percentage_pow_OCOG))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('First ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n1))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Last ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n2))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('Offset',{': '},num2str(cnf_p.OCOG_retracker.offset))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('OCOG method',{': '},num2str(cnf_p.OCOG_retracker.implementation_method))));
            fprintf(fid_proc_opts,'%s\n',char(strcat('OCOG param computation method',{': '},num2str(cnf_p.OCOG_retracker.param_comp_method))));
    end
end
fclose(fid_proc_opts);


end

