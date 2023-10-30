function download_CR2_ftp (ftpPath,user,pswrd,type_product_out,type_product_in,inputDir,outputDir,type_input_file)
% BULK CR2 in your ass a FULL
% Download FBR and L1B products from a given list of L2 files or whatever.
% Inputs type of data needed and type of data provided
% download_CR2(downloaded,input)
% ftp= 'SIR_SAR_FR/2015/01/';

mkdir(outputDir);
current_dir=cd;
cd(outputDir);
mw = ftp(ftpPath,user,pswrd);


inputFiles      =   dir([inputDir strcat('*.',type_input_file)]);
indexFiles      =   find(not([inputFiles.isdir]));
nFiles          =   length(indexFiles);

if ~isempty(strfind(type_product_in,'SAR'))
    L2_string_in='SIR_SAR_L2';
elseif ~isempty(strfind(type_product_in,'SIN'))
    L2_string_in='SIR_SIN_L2';
end

if ~isempty(strfind(type_product_out,'SAR'))
    L2_string_out='SIR_SAR_L2';
    FR_string_out='SIR_SAR_FR';
elseif ~isempty(strfind(type_product_out,'SIN'))
    L2_string_out='SIR_SIN_L2';
    FR_string_out='SIR_SIN_FR';
end

disp(strcat('Total number files source folder:', num2str(nFiles)))
i_valid_files=0;
for iFile=1:nFiles
    input_filename = [inputDir inputFiles(indexFiles(iFile)).name];
    switch type_input_file
        case {'DBL','HDR'}
            year    = inputFiles(indexFiles(iFile)).name(20:23);
            year_end    = inputFiles(indexFiles(iFile)).name(36:39);
            month   = inputFiles(indexFiles(iFile)).name(24:25);
            month_end   = inputFiles(indexFiles(iFile)).name(40:41);
            day     = inputFiles(indexFiles(iFile)).name(26:27);
            day_end     = inputFiles(indexFiles(iFile)).name(42:43);
            timeT1  = inputFiles(indexFiles(iFile)).name(29:34);
            timeT2  = inputFiles(indexFiles(iFile)).name(45:50);
        case 'nc'            
            year    = inputFiles(indexFiles(iFile)).name(17:20);
            year_end    = inputFiles(indexFiles(iFile)).name(33:36);
            month   = inputFiles(indexFiles(iFile)).name(21:22);
            month_end   = inputFiles(indexFiles(iFile)).name(37:38);
            day     = inputFiles(indexFiles(iFile)).name(23:24);
            day_end     = inputFiles(indexFiles(iFile)).name(39:40);
            timeT1  = inputFiles(indexFiles(iFile)).name(26:31);
            timeT2  = inputFiles(indexFiles(iFile)).name(42:47);
    end
    %progressbar(iFile/nFiles);
    cd(mw,['/' type_product_out '/' year '/' month '/']);

    inputftpFiles   = dir(mw,['/' type_product_out '/' year '/' month '/*' year month day 'T' timeT1 '*']);
    if (isempty(inputftpFiles) && strcmp(type_product_out,L2_string_out)) %&& ~strcmp(type_product_in,L2_string_in)
%         timeT1_bis=num2str(str2double(timeT1)+1);
%         inputftpFiles   = dir(mw,['/' type_product '/' year '/' month '/*' year month day 'T' timeT1_bis '*']);
%         if isempty(inputftpFiles) && strcmp(type_product,'SIR_SAR_L2')
%             timeT1_bis_bis=num2str(str2double(timeT2)+1);
%             inputftpFiles   = dir(mw,['/' type_product '/' year '/' month '/*' year month day 'T' timeT1 '*' timeT1_bis_bis '*']);
%         end
%         if isempty(inputftpFiles) && strcmp(type_product,'SIR_SAR_L2')            
%             inputftpFiles   = dir(mw,['/' type_product '/' year '/' month '/*' year month day 'T' timeT1_bis '*' timeT1_bis_bis '*']);
%         end
        init_acq_time=datestr(datenum(([year month day 'T' timeT1]),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
        end_acq_time=[year_end month_end day_end 'T' timeT2];
        inputftpFiles   = dir(mw,['/' type_product_out '/' year '/' month '/*' init_acq_time '_' end_acq_time '*']);
        if (isempty(inputftpFiles) && strcmp(type_product_out,L2_string_out))
            end_acq_time=datestr(datenum(end_acq_time,'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
            inputftpFiles   = dir(mw,['/' type_product_out '/' year '/' month '/*' [year month day 'T' timeT1] '_' end_acq_time '*']);
            if isempty(inputftpFiles) && strcmp(type_product_out,L2_string_out)
                inputftpFiles   =  dir(mw,['/' type_product_out '/' year '/' month '/*' init_acq_time '_' end_acq_time '*']);
            end
        end

    elseif (isempty(inputftpFiles) && ~strcmp(type_product_out,L2_string_out)) %&& ~strcmp(type_product_in,L2_string_in)
        init_acq_time=datestr(datenum(([year month day 'T' timeT1]),'yyyymmddTHHMMSS')-1/24/60/60,'yyyymmddTHHMMSS');
        end_acq_time=[year_end month_end day_end 'T' timeT2];
        inputftpFiles   = dir(mw,['/' type_product_out '/' year '/' month '/*' init_acq_time '_' end_acq_time '*']);
        if (isempty(inputftpFiles) && strcmp(type_product_in,L2_string_in))
            end_acq_time=datestr(datenum(end_acq_time,'yyyymmddTHHMMSS')-1/24/60/60,'yyyymmddTHHMMSS');
            inputftpFiles   = dir(mw,['/' type_product_out '/' year '/' month '/*' [year month day 'T' timeT1] '_' end_acq_time '*']);
            if (isempty(inputftpFiles) && strcmp(type_product_in,L2_string_in))
                inputftpFiles   =  dir(mw,['/' type_product_out '/' year '/' month '/*' init_acq_time '_' end_acq_time '*']);
            end
        end

    end
	
    if(~isempty(inputftpFiles))
        i_valid_files   = i_valid_files+1;
        indexftpFiles   =   find(not([inputftpFiles.isdir]));
        nftpFiles       =   length(indexftpFiles);
        
        for i_ftpfiles=1:length(indexftpFiles)
            try
                mget(mw, inputftpFiles(indexftpFiles(i_ftpfiles)).name);
                disp(inputftpFiles(indexftpFiles(i_ftpfiles)).name)
                if strcmp(type_product_out,FR_string_out)
                    untar(inputftpFiles(indexftpFiles(i_ftpfiles)).name);
                    delete(inputftpFiles(indexftpFiles(i_ftpfiles)).name);
                end
            catch
                continue;
            end
        end
    end
    
    
end
disp(strcat('Total of valid-downloaded files',num2str(i_valid_files)));
close(mw);
cd(current_dir);
end
