function [out] = output_data_generation(file,retrackers_results,data,cnf_p)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and generates the output L2 product
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 16/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -file    =   structure with info of input path and output path,
%       name of the original L1B product processed
%       to process the data (including the L1B as well as configuration/characterization files)
%       -retrackers_results = structure of the fitting procedures for different retrackers and 
%                           different records
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       out        =   structure of output data
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% -geo_corrections_computation: computes the corresponding geophysical
% corrections to be applied as a surface dependent correction 
% (based on a record by record surface type or assuming the same type of surface for all records)
% -write_L2_product: generate the L2 product with the necessary and
% required output product information on a .mat or as a netcdf file
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
%  
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 

%% -------------- Global Variables ----------------------------------------
global c_cst fs_clock_ku_chd bw_rx_ku_chd
%delta_rb=c_cst/(2.0*fs_clock_ku_chd*cnf_p.ZP);%sampling spacing in the range dimension


%% -------------- ORGANIZE THE OUTPUT DATA --------------------------------
%--------------------------------------------------------------------------
%strcture to be used for .mat validation
out.N_records   =   data.N_records; 
out.TAI         =   data.GEO.TAI; %seconds since first janaury 2000, days and seconds in a day structure
out.lat         =   data.GEO.LAT;
out.lon         =   data.GEO.LON;
out.H_orb       =   data.GEO.H; %include the height information
out.range       =   data.MEA.win_delay*c_cst/2.0;
out.s0_sf       =   data.HRM.s0_sf; %sigma scale factor copied from L1B in dBs
%--------------------------------------------------------------------------
% ---------------------- GEOPHYSICAL CORRECTIONS --------------------------
%--------------------------------------------------------------------------
%initialize the different geo corrections
if isfield(data,'COR')
    %generate the output structure with the corresponding geophysical
    %corrections applied as well as the total
    out.COR = geo_corrections_computation(data,cnf_p);
end

%--------------------------------------------------------------------------
%------------------------- RETRACKERS OUTPUTS -----------------------------
%--------------------------------------------------------------------------
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)    
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            i_index_analytical=i_index_analytical+1;
            switch cnf_p.range_index_method
                case 'conventional'
                    delta_rb=c_cst/(2.0*bw_rx_ku_chd*cnf_p.ZP);
                case 'resampling'
                    delta_rb=c_cst/(2.0*bw_rx_ku_chd);
            end
            if cnf_p.two_step_fitting
                %include both SWH and ROUGHNESS
                %Depending on the fitting conditions one of them will be
                %NaN while the other contain the fitted value
                out.RETRACKER.ANALYTICAL(i_index_analytical).Hs              =   retrackers_results.ANALYTICAL(i_index_analytical).Hs;
                out.RETRACKER.ANALYTICAL(i_index_analytical).rou              =   retrackers_results.ANALYTICAL(i_index_analytical).rou;
            else
                switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                    case 'SWH'
                        out.RETRACKER.ANALYTICAL(i_index_analytical).Hs              =   retrackers_results.ANALYTICAL(i_index_analytical).Hs;
                    case 'MSS'
                        out.RETRACKER.ANALYTICAL(i_index_analytical).rou              =   retrackers_results.ANALYTICAL(i_index_analytical).rou;
                end

            end
            
            out.RETRACKER.ANALYTICAL(i_index_analytical).Epoch           =   retrackers_results.ANALYTICAL(i_index_analytical).Epoch;
            out.RETRACKER.ANALYTICAL(i_index_analytical).Pu              =   retrackers_results.ANALYTICAL(i_index_analytical).Pu;
            out.RETRACKER.ANALYTICAL(i_index_analytical).sigma0          =   out.RETRACKER.ANALYTICAL(i_index_analytical).Pu+data.HRM.s0_sf;
            out.RETRACKER.ANALYTICAL(i_index_analytical).corr_coeff      =   retrackers_results.ANALYTICAL(i_index_analytical).COR.*1e2; % This is the Pearson correlation coefficient                                    
            switch cnf_p.range_index_method
                case 'conventional'
                    out.RETRACKER.ANALYTICAL(i_index_analytical).retracking_cor  =   -(cnf_p.ref_sample_wd - out.RETRACKER.ANALYTICAL(i_index_analytical).Epoch)*delta_rb; 
                case 'resampling'
                    out.RETRACKER.ANALYTICAL(i_index_analytical).retracking_cor  =   -(cnf_p.ref_sample_wd/cnf_p.ZP - out.RETRACKER.ANALYTICAL(i_index_analytical).Epoch)*delta_rb; 
            end                       
            out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range   =   out.range + out.RETRACKER.ANALYTICAL(i_index_analytical).retracking_cor;
            if cnf_p.geo_corr_application_flag
                out.RETRACKER.ANALYTICAL(i_index_analytical).SSH             =   out.H_orb - (out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range+out.COR.total_geo_corr);
            else
                out.RETRACKER.ANALYTICAL(i_index_analytical).SSH             =   out.H_orb - (out.RETRACKER.ANALYTICAL(i_index_analytical).tracker_range);
            end
            
            out.RETRACKER.ANALYTICAL(i_index_analytical).flag_fitting    =   retrackers_results.ANALYTICAL(i_index_analytical).flag;
            if i_retracker==1
                out.flag_L1B        =   retrackers_results.ANALYTICAL(i_index_analytical).flag_validity_L1B; % this flag is to indicate whether L1B processing provides a valid Waveform
            end
            %disp(strcat('Analytical Retracker correction [m]',num2str(mean([out.RETRACKER.ANALYTICAL.retracking_cor]))))
        case 'THRESHOLD'
            %TBD [retracker_results.THRESHOLD]=threshold_retracker(data,cnf_p);
            delta_rb=c_cst/(2.0*bw_rx_ku_chd*cnf_p.ZP);
            out.RETRACKER.THRESHOLD.Epoch           =   retrackers_results.THRESHOLD.Epoch;
            out.RETRACKER.THRESHOLD.Pu              =   retrackers_results.THRESHOLD.Pu;
            out.RETRACKER.THRESHOLD.retracking_cor  =   -(cnf_p.ref_sample_wd - out.RETRACKER.THRESHOLD.Epoch)*delta_rb;            
            out.RETRACKER.THRESHOLD.tracker_range   =   out.range + out.RETRACKER.THRESHOLD.retracking_cor;
            if cnf_p.geo_corr_application_flag
                out.RETRACKER.THRESHOLD.SSH             =   out.H_orb - (out.RETRACKER.THRESHOLD.tracker_range+ out.COR.total_geo_corr);
            else
                out.RETRACKER.THRESHOLD.SSH             =   out.H_orb - (out.RETRACKER.THRESHOLD.tracker_range);
            end            
            out.RETRACKER.THRESHOLD.sigma0          =   out.RETRACKER.THRESHOLD.Pu+data.HRM.s0_sf;
            if i_retracker==1
                out.flag_L1B        =   retrackers_results.THRESHOLD.flag_validity_L1B; % this flag is to indicate whether L1B processing provides a valid Waveform
            end
            %disp(strcat('Threshold Retracker correction [m]',num2str(mean([out.RETRACKER.THRESHOLD.retracking_cor]))))
        case 'OCOG'
            delta_rb=c_cst/(2.0*bw_rx_ku_chd*cnf_p.ZP);
            out.RETRACKER.OCOG.Epoch           =   retrackers_results.OCOG.Epoch;
            out.RETRACKER.OCOG.Pu              =   retrackers_results.OCOG.Pu; % in dBs;
            out.RETRACKER.OCOG.retracking_cor  =   -(cnf_p.ref_sample_wd - out.RETRACKER.OCOG.Epoch)*delta_rb;            
            out.RETRACKER.OCOG.tracker_range   =   out.range + out.RETRACKER.OCOG.retracking_cor;
            if cnf_p.geo_corr_application_flag
                out.RETRACKER.OCOG.SSH             =   out.H_orb - (out.RETRACKER.OCOG.tracker_range+ out.COR.total_geo_corr);
            else
                out.RETRACKER.OCOG.SSH             =   out.H_orb - (out.RETRACKER.OCOG.tracker_range);
            end            
            out.RETRACKER.OCOG.sigma0          =   out.RETRACKER.OCOG.Pu+data.HRM.s0_sf;
            out.RETRACKER.OCOG.COG             =   retrackers_results.OCOG.COG*delta_rb; %in seconds
            out.RETRACKER.OCOG.A               =   retrackers_results.OCOG.A; % in dB
            out.RETRACKER.OCOG.W               =   retrackers_results.OCOG.W*delta_rb; % in seconds
            if i_retracker==1
                out.flag_L1B        =   retrackers_results.OCOG.flag_validity_L1B; % this flag is to indicate whether L1B processing provides a valid Waveform
            end
            %TBD [retracker_results.OCOG]=OCOG_retracker(data,cnf_p);
    end
    
end

% -------------------------------------------------------------------------
% ------------------------ SURFACE TYPE -----------------------------------
% -------------------------------------------------------------------------
out.surf_type_flag=data.surf_type_flag;

% -------------------------------------------------------------------------
% ------------------------ PROCESSING OPTIONS -----------------------------
% -------------------------------------------------------------------------
out.PROC_CNF=cnf_p;

% -------------------------------------------------------------------------
% ------------------------ GLOBAL ATTRIBUTES ------------------------------
% -------------------------------------------------------------------------
% Imported directly from L1B: Including both data file related attributes &
% orbital info
if isfield(data,'GLOBAL_ATT')
    out.GLOBAL_ATT=data.GLOBAL_ATT;
    %include the original L1B file used
    out.GLOBAL_ATT.DATA_FILE_INFO.xref_altimeter_level1B=[file.filename_L1B_nopath file.fileext_L1B];
end

%% --------------- WRITE OUTPUT PRODUCT -----------------------------------
%--------------------------------------------------------------------------
if cnf_p.write_output
    write_L2_product(file,out,cnf_p);
end

   


end

