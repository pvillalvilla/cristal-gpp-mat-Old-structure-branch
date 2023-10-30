% Move the files of the already processed ones to a given subdirectory
resultsDir    = 'C:/Users/eduard.makhoul/isardSAT/conferences_meetings/OSTST2016/processing/data/inputs/L2_STL_S3/central_pacific/2013_new/';
originalData   = 'C:/Users/eduard.makhoul/isardSAT/publications/ARS_special_issue_CS2/processing/results/central_pacific/L2_DATA/L2_ISD_S3/data/';
originalData_processed='C:/Users/eduard.makhoul/isardSAT/publications/ARS_special_issue_CS2/processing/results/central_pacific/L2_DATA/L2_ISD_S3/data/downloaded/';
res_type_file_Ext='nc';
directory_active=0; %directories with the name of a given file
original_type_file_Ext='txt';
mkdir(originalData_processed)

if directory_active
    inputFiles=dir(resultsDir);
    inputFiles = inputFiles(~cellfun('isempty', {inputFiles.date}));
    aux=struct2cell(inputFiles); %into a cell array where first row is
    folders=(aux(1,[inputFiles.isdir]))'; %name is the first row
    clear aux;
    indexFiles      =   find(~strcmp(folders,['.'])&~strcmp(folders,['..']));        
else
    inputFiles      =   dir([resultsDir strcat('*.',res_type_file_Ext)]);
    indexFiles      =   find(not([inputFiles.isdir]));
end
nFiles          =   length(indexFiles);
iFile_moved=0;

disp(strcat('Total number of processed files ',num2str(nFiles)))
for iFile=1:nFiles

    %progressbar(iFile/nFiles);
    input_filename = [resultsDir inputFiles(indexFiles(iFile)).name];
    switch res_type_file_Ext
        case {'DBL','HDR'}
            timeacquisition_file=inputFiles(indexFiles(iFile)).name(20:20+30);
        case {'nc'}
            timeacquisition_file=inputFiles(indexFiles(iFile)).name(17:17+30);
        case {'dir'}
            %assuming the name of the directory is the time_acquisition
            timeacquisition_file=inputFiles(indexFiles(iFile)).name;
    end    
    inputFileprocessed   = dir(fullfile(originalData,['*' timeacquisition_file strcat('*.',original_type_file_Ext)]));
    if(~isempty(inputFileprocessed))
       indexFileprocessed   =   not([inputFileprocessed.isdir]);
       filename_Fileprocessed=char(inputFileprocessed(indexFileprocessed).name);
	   disp(filename_Fileprocessed);
       full_filename_Fileprocessed=strcat(originalData,filename_Fileprocessed);
       %copyfile(full_filename_Fileprocessed,originalData_processed) 
       movefile(strrep(full_filename_Fileprocessed,'DBL','HDR'),originalData_processed) 
	   iFile_moved=iFile_moved+1;
    end
end
disp(strcat('Total number of moved files',num2str(iFile_moved)))
%exit