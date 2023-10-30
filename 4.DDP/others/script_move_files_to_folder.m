% Move the files of the already processed ones to a given subdirectory
resultsDir    = '\\X8DTL\data\DATA\SCOOP\CS2_FBR_xcross_zero_vr\L1B\';
originalData   = '\\X8DTL\data\DATA\SCOOP\copy_disk_ESA_2013\central_pacific\';
originalData_processed='\\X8DTL\data\DATA\SCOOP\CS2_FBR_xcross_zero_vr\FBR\';
res_type_file_Ext='nc';
directory_active=0; %directories with the name of a given file
original_type_file_Ext='DBL';
copy_move='c'; % 'c': copy or 'm': move
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
        case {'mat'} %performance file: temporal solution
            timeacquisition_file=inputFiles(indexFiles(iFile)).name(6:6+30);
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
       switch copy_move
           case 'c'
               copyfile(full_filename_Fileprocessed,originalData_processed)
           case 'm'
               movefile(strrep(full_filename_Fileprocessed,'DBL','HDR'),originalData_processed)
       end
	   iFile_moved=iFile_moved+1;
    end
end
disp(strcat('Total number of moved files',num2str(iFile_moved)))
%exit