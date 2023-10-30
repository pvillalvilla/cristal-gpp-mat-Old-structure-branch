function L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code runs a bulk L2 processing either in parallel or sequentially
% over the different input L1B files defined in a given folder, using a given
% set of configuration and characterization files provided in a given
% folder, storing the results in a specified folder
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 20/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -input_path_L1B    =   path containing the L1B file/s to be
%       processed either in parallel or sequentially
%       -output_path_L2 = output folder to save the different L2 products
%       (will be written whenever the corresponding flag cnf_p.write_output (in the cnf_file.m) is activated )
%       -cnf_chd_cst_path: path containing the different configuration
%       cnf_file.m
%       (configuration processing parameters for L2), characterization
%       chd_file.m
%       (different missions parameters) and the constant definition
%       cst_file.m; the Look Up Tables (LUTs) for the f0 and f1 functions
%       to be used by the single look waveform model are also included in this folder
%     OPTIONAL
%       -input_path_L1BS: path to the related L1BS products associated to
%       L1B products in input_path_L1B
%       -filename_mask_KML: fullfilename of a KML file containing the ROI
%       in order to filter out the records of interest: it is used only if
%       the related flag is activated cnf_p.mask_ROI_flag (in cnf_file.m)
%       -num_pools: indicates the number of pools or threads to be run for
%       parallel processing of different tracks: by default is set to 1 if
%       no paralelization is desired (sequential processing)
%       -targz_option_active_L1B: flag to indicate whether the available
%       L1B products within the folder to be processed are compressed and need
%       to be untar
%      -targz_option_active_L1BS: flag to indicate whether the available
%      L1BS products within the folder to be processed are compressed and need
%       to be untar
%       
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%  - L2_processing: runs the L2 processing; geophysical retrievals
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% RESTRICTIONS:
% - It is assumed all the L1B input files in the inputh_path_L1B correspond
% to the same mission and have been processed with the same configuration
% parameters
% - LUTs files are saved as .mat files, need to be updated for a more
% generic case of binary or netcdf files
% - The optional L1B-S files included in a different folder 
% - Further updates to have the possibility of having compressed and
% uncompressed files in the same folder
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 

clear global
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');

%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>(3+6*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('input_path_L1BS',{''},@(x)ischar(x));
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.addParamValue('num_pools',1);
p.addParamValue('targz_option_active_L1B',0);
p.addParamValue('targz_option_active_L1BS',0);
p.addParamValue('proc_bsln_id',{''},@(x)ischar(x));
p.parse(varargin{:});
input_path_L1BS=char(p.Results.input_path_L1BS);
filename_mask_KML=char(p.Results.filename_mask_KML);
num_pools=p.Results.num_pools;
targz_option_active_L1B=(p.Results.targz_option_active_L1B);
targz_option_active_L1BS=(p.Results.targz_option_active_L1BS);
proc_bsln_id=char(p.Results.proc_bsln_id);
clear p;

global alpha_gr_chd alpha_ga_chd;

version_matlab=version;
ttotal=tic;


%% ----- Include the different subfolders in the search path --------------
FolderInfo=dir;
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); %into a cell array where first row is 
folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
clear aux;
folders=strcat('./',folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end


%% --------- Define Paths -------------------------------------------------
filesBulk.inputPath       =   input_path_L1B;
filesBulk.resultPath      =   output_path_L2;
mkdir(filesBulk.resultPath);
mkdir([filesBulk.resultPath 'data/']);
mkdir([filesBulk.resultPath 'plots/']);
mkdir([filesBulk.resultPath 'plots/fitted_waveforms/']);



%% --------------- Configuration/characterization/LUTS --------------------
%assume all the files in the folder correspond to the same mission and
%processed with the same L1B processor: use the same
%configuration/characterization for Level-2 processing
filesBulk.cnf_chd_cst_path        =   cnf_chd_cst_path;
inputFiles      =   dir(filesBulk.cnf_chd_cst_path);
aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the
filesBulk.CNF_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,char(['cnf_file','_',proc_bsln_id])))).name];
disp(filesBulk.CNF_file);
filesBulk.CHD_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,'chd_file'))).name];
filesBulk.CST_file=[filesBulk.cnf_chd_cst_path inputFiles(~cellfun(@isempty,strfind(aux,'cst_file'))).name];
run(filesBulk.CNF_file); %CNF --> generate a cnf_p structure
% run(filesBulk.CST_file); %CST --> provide the global constant variables
run(filesBulk.CHD_file); %CHD --> provide the global characterization variables

%---------------------- LUTs files-----------------------------------------
%assuming .mat files
idx_int=~cellfun(@isempty,strfind(aux,'LUT_f0'));
if any(idx_int)
    filesBulk.LUT_f0_file=[filesBulk.cnf_chd_cst_path inputFiles(idx_int).name];
end
idx_int=~cellfun(@isempty,strfind(aux,'LUT_f1'));
if any(idx_int)
    filesBulk.LUT_f1_file=[filesBulk.cnf_chd_cst_path inputFiles(idx_int).name];
end


%% ------------------ Input L1B files filtering ---------------------------
inputFiles      =   dir(filesBulk.inputPath);
aux=struct2cell(inputFiles); aux=aux(1,:); %Keep the

if targz_option_active_L1B
    filterDATAFILES=(~cellfun(@isempty,strfind(aux,'TGZ')));
    indexFilesL1B=find(filterDATAFILES);
else
    switch cnf_p.mission
        case {'CS2','CR2'}
            % --------------------- CroySAT-2 ---------------------------------
            switch cnf_p.L1proc
                case 'ESA'
                    filter='.DBL';
                case 'ISD'
                    filter='.nc';
            end
        case {'S6','JCS','S3'}
            % --------------------- Sentinel-3/-6 -----------------------------
            filter='.nc';
        otherwise
            error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
    end
    filterDATAFILES=(~cellfun(@isempty,strfind(aux,filter)));
    indexFilesL1B=find(filterDATAFILES);
    
end

filesBulk.nFilesL1B=length(indexFilesL1B);
filesBulk.L1BFiles=inputFiles(indexFilesL1B);

disp('Total number of L1B files to be processed');
disp(num2str(filesBulk.nFilesL1B))


%% --------------- Input L1B-S files --------------------------------------
if ~isempty(input_path_L1BS)
    filesBulk.inputPath_L1BS        =   input_path_L1BS;
    inputFiles      =   dir(filesBulk.inputPath_L1BS);
    aux=struct2cell(files_L1BS.inputFiles); aux=aux(1,:); %Keep the
    
    if targz_option_active_L1BS
        filterL1BSDATAFILES=(~cellfun(@isempty,strfind(aux,'TGZ')));
        indexFilesL1BS=find(filterL1BSDATAFILES);
    else
        filterL1BSDATAFILES=(~cellfun(@isempty,strfind(aux,filter)));
        indexFilesL1BS=find(filterL1BSDATAFILES);        
    end
    
    filesBulk.nFilesL1BS=length(indexFilesL1BS);
    filesBulk.L1BSFiles=inputFiles(indexFilesL1BS);
    clear files_L1BS indexFilesL1BS filterL1BSDATAFILES;    
end

%% --------------- Displaying options of retrackers -----------------------
disp('------------------------------------------------------------------------------------')
disp('--------------- General Parameters & Options retrackers ----------------------------')
disp('------------------------------------------------------------------------------------')
disp(strcat('Filtering Geographical mask',{': '},num2str(cnf_p.mask_ROI_flag)))
disp(strcat('Filtering Number looks stack',{': '},num2str(cnf_p.mask_looks_flag)))
if cnf_p.mask_looks_flag
    disp(strcat('Minimum Number looks stack (filtering)',{': '},num2str(cnf_p.Neff_thres)))
end
disp(strcat('IFmask',{': '},num2str(cnf_p.IFmask_N)))
disp(strcat('Geophysical corrections applied',{': '},num2str(cnf_p.geo_corr_application_flag)))
disp(strcat('force_geocorr_surf_type',{': '},num2str(cnf_p.force_geocorr_surf_type)))
if cnf_p.force_geocorr_surf_type
    disp(strcat('product_type_surface',{': '},cnf_p.product_type_surface))
end

for i_retracker=1: length(cnf_p.retracker_name)
    switch char(cnf_p.retracker_name(i_retracker))
        case {'ANALYTICAL','SAMOSA'}
            disp('------------------------------------------------------------------------------------')
            disp('--------------- Parameters & Options of Analytical retracker -----------------------')
            disp('------------------------------------------------------------------------------------')
            disp(strcat('Use zeros in multilooking',{': '},num2str(cnf_p.use_zeros_cnf)))
            disp(strcat('Zero Padding in range',{': '},num2str(cnf_p.ZP)))
            disp(strcat('window_type_a',{': '},cnf_p.window_type_a));
            disp(strcat('window_type_r',{': '},cnf_p.window_type_r));
            disp(strcat('Range PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_gr_chd)))))
            disp(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_ga_chd)))))
            if cnf_p.fit_noise
                disp(strcat('Threshold noise',{': '},num2str(cnf_p.threshold_noise)))
            else
                disp(strcat('First noise sample',{': '},num2str(cnf_p.Thn_w_first)))
                disp(strcat('Width noise window',{': '},num2str(cnf_p.Thn_w_width)))
            end
            
%             disp(strcat('Sign pitch',{': '},num2str(cnf_p.sign_pitch)))
%             disp(strcat('Sign roll',{': '},num2str(cnf_p.sign_roll)))

            disp(strcat('Indexation method',{': '},(cnf_p.looks_index_method)))
            switch cnf_p.looks_index_method
                case 'Look_angle'
                    disp(strcat('Look angle method',{': '},cnf_p.look_ang_method))
                case 'Doppler_freq'
                    disp(strcat('Doppler freq. method',{': '},cnf_p.fd_method))
            end
            disp(strcat('Power wvfm model',{': '},cnf_p.power_wfm_model))
            disp(strcat('LUT flag',{': '},num2str(cnf_p.lut_flag)))
            if cnf_p.pre_processing
                disp(strcat('Threshold (leading edge pre-processing)',{': '},num2str(cnf_p.percent_leading_edge)))
            end
            disp(strcat('Two step fitting',{': '},num2str(cnf_p.two_step_fitting)))
            if cnf_p.two_step_fitting
                disp(strcat('Two step fitting COR threshold',{': '},num2str(cnf_p.two_step_fitting_COR_threshold_rou)))
            end
            disp(strcat('initial_param_fit_feedback_flag',{': '},num2str(cnf_p.initial_param_fit_feedback_flag)))
            disp(strcat('ini_Hs_rou_sliding_win_opt',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt)))
            disp(strcat('ini_Hs_rou_sliding_win_opt_discard_std_threshold',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_opt_discard_std_threshold)))
            disp(strcat('ini_Hs_rou_sliding_win_size',{': '},num2str(cnf_p.ini_Hs_rou_sliding_win_size)))
            disp(strcat('Range indexation method',{': '},cnf_p.range_index_method))
            disp(strcat('Fitting method',{': '},cnf_p.fitting_fun_type))
            cnf_p.fitting_options
            
        case 'THRESHOLD'
            disp('------------------------------------------------------------------------------------')
            disp('--------------- Parameters & Options of Threshold retracker ------------------------')
            disp('------------------------------------------------------------------------------------')
            disp(strcat('Threshold value',{': '},num2str(cnf_p.th_retracker.percentage_peak)))

        case 'OCOG'
            disp('------------------------------------------------------------------------------------')
            disp('--------------- Parameters & Options of OCOG retracker -----------------------------')
            disp('------------------------------------------------------------------------------------')
            disp(strcat('Percentage OCOG Amplitude',{': '},num2str(cnf_p.OCOG_retracker.percentage_pow_OCOG)))
            disp(strcat('First ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n1)))
            disp(strcat('Last ZP sample',{': '},num2str(cnf_p.OCOG_retracker.n2)))
            disp(strcat('Offset',{': '},num2str(cnf_p.OCOG_retracker.offset)))            
            disp(strcat('OCOG method',{': '},num2str(cnf_p.OCOG_retracker.implementation_method)))
            disp(strcat('OCOG param computation method',{': '},num2str(cnf_p.OCOG_retracker.param_comp_method)))
    end
end
%% --------------- Run parallel processing --------------------------------
if num_pools~=1
    %create pools
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_pools);
    else
        matlabpool('open',num_pools);
    end
    %% ------------- Loop per file to be processed ------------------------
    parfor i_fileL1B_input=1:filesBulk.nFilesL1B
%         try 
            L2_processing (filesBulk, i_fileL1B_input,...
                'filename_mask_KML',filename_mask_KML,...
                'targz_option_active_L1B',targz_option_active_L1B,...
                'targz_option_active_L1BS',targz_option_active_L1BS)

%         catch
%             continue
%         end
    end
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end
else
    for i_fileL1B_input=1:filesBulk.nFilesL1B
    %   try
        L2_processing (filesBulk, i_fileL1B_input,...
            'filename_mask_KML',filename_mask_KML,...
            'targz_option_active_L1B',targz_option_active_L1B,...
            'targz_option_active_L1BS',targz_option_active_L1BS)
%        catch
%            continue;
%        end
    end
end

time = toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time: ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);
%exit
end
    

