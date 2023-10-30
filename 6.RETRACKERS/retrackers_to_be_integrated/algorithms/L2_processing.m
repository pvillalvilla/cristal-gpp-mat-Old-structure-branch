function exit_flag=L2_processing (filesBulk, i_fileL1B_input,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs the L2 processing
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1.2 21/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -filesBulk    =   structure of input files within the folder where
%       to process the data (including the L1B as well as
%       configuration/characterization files):
%        filesBulk.inputPath        --> full path to the input L1B folder
%        filesBulk.resultPath       --> full path to results folder save L2 prod
%        filesBulk.inputPath_L1BS   --> full path to input L1BS folder
%        filesBulk.cnf_chd_cst_path --> full path to cnf/chd/cst/LUTs folder
%        filesBulk.CNF_file         --> full filename cnf config file
%        filesBulk.CHD_file         --> full filename chd charac file
%        filesBulk.CST_file         --> full filename cst const file
%        filesBulk.LUT_f0_file      --> full filename LUT for f0 function
%        filesBulk.LUT_f1_file      --> full filename LUT for f1 function
%        filesBulk.nFilesL1B        --> total number of L1B to be processed
%        filesBulk.nFilesL1BS       --> total number of L1BS to be processed
%        filesBulk.L1BFiles         --> information of the L1B files to be
%                                       processed: fields--> {name,date,bytes,isdir,datenum}
%        filesBulk.L1BSFiles        --> information of the L1BS files to be
%                                       processed: fields--> {name,date,bytes,isdir,datenum}
%       -i_fileL1B_input = index of the L1B within the filesBulk.LBFiles
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       Will write the corresponding L2 product in the specified resultPath
%       folder whenever the related flag cnf_p.write_output (n the cnf_file.m) is activated 
%       exit_flag = flag indicating whether the processing has been
%       successful or not o there is an error value to -1
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%   - preparing_inputFiles_L2processing: prepare and set the L1B, L1Bs and
%   KMl files when applicable from the input filesBulk data and cnf_p
%   - read_alt_data_EM: read input data from any combination of mission and
%   processor either in netcdf mat or DBL and save it in a common data
%   structure to be used in the L2 processing
%   - analytical_retracker: performs the analytical retracking based on the
%   original model developed by Chirs Ray et al. in IEEE TGRS "SAR
%   Altimeter Backscattered Waveform Model": DOI:
%   10.1109/TGRS.2014.23330423
%   - output_data_generation: acomodate the data of interest at the output
%   of the L2 processor in a given output structure which can be optionally 
%   saved (depending on the cnf_p.write_output) in .mat or  in netCDF 
%   depending on flag cnf_p.output_product_format
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
% v1.1: Inclusion of option to run several retrackers: 
%   all the functions and processing related to the analytical or SAMOSA
%   retracker have been moved within a single function "analytical_retracker.m"
% v1.2: Accomodation and preparation of the input files into a separate
% auxiliar function "preparing_inputfiles_L2processing.m" in order to keep
% a more compact definition of L2_processing.m function
time_init=tic;
%% ------------ RUN global variables --------------------------------------
global cnf_p
run(filesBulk.CNF_file); %CNF --> generate a cnf_p structure
run(filesBulk.CST_file); %CST --> provide the global constant variables
run(filesBulk.CHD_file); %CHD --> provide the global characterization variables

%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>(3+3*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.addParamValue('targz_option_active_L1B',0);
p.addParamValue('targz_option_active_L1BS',0);
p.parse(varargin{:});
filename_mask_KML=char(p.Results.filename_mask_KML);
targz_option_active_L1B=(p.Results.targz_option_active_L1B);
targz_option_active_L1BS=(p.Results.targz_option_active_L1BS);
clear p;

%% ------------------ Ploting formating -----------------------------------
set_default_plot;
if cnf_p.visible_figures
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

%% ------------ Prepare/organize the input files definition ---------------
%--------------------------------------------------------------------------
[filename_L1B,filename_L1B_nopath,fileext_L1B,filename_L1BS]=preparing_inputFiles_L2processing(filesBulk, i_fileL1B_input, cnf_p,...                                               
                                               'targz_option_active_L1B',targz_option_active_L1B,...
                                               'targz_option_active_L1BS',targz_option_active_L1BS);


disp(strcat('Processing L1B: ',filename_L1B_nopath,fileext_L1B))

%% --------------------- LOAD DATA FROM INPUT FILE L1B --------------------
%--------------------------------------------------------------------------
%reading the data from the specific L1B & load it in a specific "data"
%structure & specific filtering (geographical and/or number of looks is further applied)
[data,flag] = read_alt_data_EM (filename_L1B, cnf_p,'filename_L1BS',filename_L1BS,'filename_mask_KML',filename_mask_KML);
if flag==-1
    exit_flag=flag;
    return;
end

%% --------------------- RETRACKER RUN ------------------------------------
%--------------------------------------------------------------------------
i_index_analytical=0;
for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            disp('Processing with Analytical retracker ....')
            %Patch to be able to run the analytical more than once (including the results in the same file) with
            %different options on the type of fitting : SWH, ROUGH or BOTH
            %Since we are using a single CNF file in this case for all of
            %them: other option is to create a different CNF for each
            %retracker
            i_index_analytical=i_index_analytical+1;
            if ~cnf_p.two_step_fitting
                switch char(cnf_p.analytical_type_of_fitting(i_index_analytical))
                    case 'SWH'
                        cnf_p.rou_flag=0;
                    case 'MSS'
                        cnf_p.rou_flag=1;
                end
            end
            [retracker_results.ANALYTICAL(i_index_analytical)]=analytical_retracker(data,cnf_p,...
                                                                'LUT_f0_file',filesBulk.LUT_f0_file,...
                                                                'LUT_f1_file',filesBulk.LUT_f1_file,...
                                                                'path_Results',filesBulk.resultPath,...
                                                                'L1B_filename',filename_L1B_nopath);
        case 'THRESHOLD'
            disp('Processing with Threshold retracker ....')
            [retracker_results.THRESHOLD]=threshold_retracker(data,cnf_p);
        case 'OCOG'
            disp('Processing with OCOG retracker ....')
            [retracker_results.OCOG]=OCOG_retracker(data,cnf_p);
    end
end

%% -------------- ORGANIZE THE OUTPUT DATA --------------------------------
%--------------------------------------------------------------------------
%generating an output data structure and writing the corresponding product
%(either using .nc or .mat according to the configuration)
file.filename_L1B_nopath=filename_L1B_nopath;
file.fileext_L1B=fileext_L1B;
file.inputPath=filesBulk.inputPath;
file.resultPath=filesBulk.resultPath;
[out]=output_data_generation(file,retracker_results,data,cnf_p);


time_end=toc(time_init);
minutes_processing = floor(time_end/60);
secs_processing = time_end - minutes_processing*60;
disp(['Processing time for ',file.filename_L1B_nopath,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
   
    


end

