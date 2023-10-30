%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VALIDATION L2/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this algorithm is to cross-check the L2 products:
% isardSAT and ESA for CryoSat-2 data runs the validation and saves the
% results into a .mat with a given structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_L2_validation(filesBulk,i_files_valid,input_path_L2_ISD,...
                                  input_path_L2_ESA,...
                                  path_results_comparison,varargin)

%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<5 || nargin>(5+8*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('input_path_L2_STL',{''},@(x)ischar(x));
p.addParamValue('retrackers',{''},@(x)iscellstr(x)); %is a cell array with different names of different retrackers in L2 product
p.addParamValue('figures_visible',0);
p.addParamValue('flag_outliers_removal',0);
p.addParamValue('type_outliers_removal','percentiles');
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('geo_mask',[]);
p.addParamValue('annotation_box_active',1);
p.parse(varargin{:});
input_path_L2_STL=char(p.Results.input_path_L2_STL);
retrackers=p.Results.retrackers;
figures_visible=p.Results.figures_visible;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=p.Results.type_outliers_removal;
sh_name_nc=p.Results.sh_name_nc;
geo_mask=p.Results.geo_mask;
annotation_box_active=p.Results.annotation_box_active;
clear p;
                              
filename_L2_ISD=char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name);
data_string=filename_L2_ISD(17:17+30);
filename_L2_ISD=strcat(input_path_L2_ISD,filename_L2_ISD);

inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' data_string(1:15) '*C001.DBL']));
if isempty(inputL2ESAFiles)
    %add one second to initial time acquisition
    init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
    inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
    if isempty(inputL2ESAFiles)
        %add one second to initial time acquisition
        end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
        inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
        if isempty(inputL2ESAFiles)
            inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
        end
    end
    
end
if isempty(input_path_L2_STL)
    inputL2STLFiles   = 1;
else
    inputL2STLFiles   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
end

disp(char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name))
indexL2ESA   =   not([inputL2ESAFiles.isdir]);
filename_L2_ESA=strcat(input_path_L2_ESA,char(inputL2ESAFiles(indexL2ESA).name));
if ~isempty(input_path_L2_STL)
    indexL2STL   =   not([inputL2STLFiles.isdir]);
    filename_L2_STL=strcat(input_path_L2_STL,char(inputL2STLFiles(indexL2STL).name));
    [~]=L2_validation_filtering(filename_L2_ISD,filename_L2_ESA,path_results_comparison,'retrackers',retrackers,'figures_visible',figures_visible,...
        'filename_L2_STL',filename_L2_STL,'flag_outliers_removal',flag_outliers_removal,...
        'type_outliers_removal',type_outliers_removal,...
        'sh_name_nc',sh_name_nc,...
        'geo_mask',geo_mask,...
        'annotation_box_active',annotation_box_active);
else
    [~]=L2_validation_filtering(filename_L2_ISD,filename_L2_ESA,path_results_comparison,'retrackers',retrackers,'figures_visible',figures_visible,...
        'flag_outliers_removal',flag_outliers_removal,...
        'type_outliers_removal',type_outliers_removal,...
        'sh_name_nc',sh_name_nc,...
        'geo_mask',geo_mask,...
        'annotation_box_active',annotation_box_active);
    
end

end