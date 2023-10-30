%profile on;
%% -------------------- Number of Pools -----------------------------------
num_pools=1;
%% -------------------- INCLUSION OF SEARCH CODE PATH ---------------------
code_folder_full_path='C:\Users\albert\Documents\WORK\SEOM\DeDop\Matlab\L2_GPP_albert\';

cd(code_folder_full_path);
FolderInfo=dir(code_folder_full_path);
FolderInfo = FolderInfo(~cellfun('isempty', {FolderInfo.date})); 
aux=struct2cell(FolderInfo); %into a cell array where first row is 
folders=(aux(1,[FolderInfo.isdir]))'; %name is the first row
clear aux;
folders=strcat(code_folder_full_path,folders(~strcmp(folders,['.'])&~strcmp(folders,['..'])&~strcmp(folders,['.svn'])&~strcmp(folders,['inputs'])));
for i_folder=1:length(folders)
    addpath(genpath(char(folders(i_folder))));
end


%% Agulhas
%% S-3 new
proc_bsln_id='S3';
input_path_L1B='C:\Users\albert\Documents\WORK\SEOM\SHAPE\DATA\SARIn\Danube\results_v106\L1b\';
disp(input_path_L1B);
output_path_L2='C:\Users\albert\Documents\WORK\SEOM\SHAPE\DATA\SARIn\Danube\results_v106\L2\';
cnf_chd_cst_path=strcat(code_folder_full_path,'inputs/');
filename_mask_KML='';
L2_bulk_processing_paralelization(input_path_L1B,output_path_L2,cnf_chd_cst_path,'num_pools',num_pools,'filename_mask_KML',filename_mask_KML,'proc_bsln_id',proc_bsln_id)

