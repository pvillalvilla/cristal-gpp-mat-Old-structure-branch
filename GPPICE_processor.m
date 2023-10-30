%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
%
% ---------------------------------------------------------
% Objective: Read input parameters and start the chains
%
% Calling: GPPICE (input_folder, aux_folder, output_folder, options)
% INPUTs:
%
%
% OUTPUTs:
%
%
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GPPICE_processor (input_folder, aux_folder, output_folder, options)

if(nargin < 4)
    options.plotting_flag = [0 0 0]; % Stacks, L1B waveforms, track
    options.axes = [];
    options.wd_axes = [];
    options.GUI_flag = 0;
    log_flag = 1;
end

filesBulk.inputPath         = fix_paths (input_folder);
filesBulk.auxPath           = fix_paths (aux_folder);
filesBulk.outputPath        = fix_paths (output_folder);

filesBulk.inputFiles        =   dir(filesBulk.inputPath);
filesBulk.auxFiles          =   dir(filesBulk.auxPath);
filesBulk.indexaDirs        =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles        =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles            =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.auxFiles); aux=aux(1,:);

if(log_flag)
    filesBulk.fid_log = fopen([filesBulk.outputPath 'LogError.txt'],'w');
    filesBulk.LogMsg = '';
end
%% READ CONSTANTS
filesBulk.CST_file=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'cst_file'))).name];
run(filesBulk.CST_file); % cst struct output
%% READ CHD
filesBulk.CHD_file=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'.xml'))).name];
chd = read_chd(filesBulk.CHD_file, cst);
%% READ CONFIG
filesBulk.CNF_file=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'cnf_file_L1'))).name];
run(filesBulk.CNF_file); % cnf struct output
filesBulk.Success = [cnf.run_L1A cnf.run_DDP cnf.run_FFtime cnf.run_FFtime cnf.run_FFfreq cnf.run_LROS]*-1;
%% RUN L1 CHAINS
 try
    aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:);
    filesBulk.filename_L0 = [filesBulk.inputPath filesBulk.inputFiles(~cellfun(@isempty,strfind(aux,'S3NGT_SIRS_'))).name];
    [chd] = check_dimensions(filesBulk,chd,'L0');

    if (cnf.run_L1A) %% ------------------------- L1A ---------GO!!!----------------
        if(cnf.verify_L0)
            plots_L0(filesBulk, cnf, chd, cst);
        end
       
        % ---------------------------------------------------------------
        filesBulk = L1A_chain (filesBulk, chd, cnf, cst, options);
        % ---------------------------------------------------------------
        filesBulk.LogMsg = [filesBulk.LogMsg,'L1A processed SUCCESSFULLY -- '];
        filesBulk.Success(1) = 1;
    else
        filesBulk.LogMsg = [filesBulk.LogMsg,'L1A not processed -- '];
        filesBulk.Success(1) = 0;
        L0_date_find = strfind(filesBulk.filename_L0,'_201');
        L0_date_find = strfind(filesBulk.filename_L0,'_202');
        L0_mode_find = strfind(filesBulk.filename_L0, 'S3NGT_SIRS___');
        
        filesBulk.filename_L1A = strcat(filesBulk.outputPath,'S3NGT_SIRS_', ...
            filesBulk.filename_L0((L0_mode_find(1)+11):(L0_mode_find(end)+15)), '_1A_',...
            filesBulk.filename_L0((L0_date_find(1)+1):(L0_date_find(1)+30)),...
            '_isd','.nc');
    end
    
    filesBulk.outputFiles        =   dir(filesBulk.outputPath);
    aux=struct2cell(filesBulk.outputFiles); aux=aux(1,:);
    
    if (~cnf.writting_flag(1))
        L1A_name = dir([output_folder '/*_1A_20*.nc']);
        filesBulk.filename_L1A=[output_folder '/' L1A_name.name];
        %filesBulk.filename_L1A=[filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'SR_1_SRA__A_0180824T120034_20180824T120041_0001_isd.nc'))).name];
    end
    
    if(~cnf.writting_flag(3)||~cnf.run_DDP)
        L1A_date_find = strfind(filesBulk.filename_L1A,'_201');
        L1A_date_find = strfind(filesBulk.filename_L1A,'_202');
        L1A_mode_find = strfind(filesBulk.filename_L1A, 'S3NGT_SIRS_');
        
        
        filesBulk.filename_L1B = strcat(filesBulk.outputPath,'S3NGT_SIRS_',...
            filesBulk.filename_L1A((L1A_mode_find(1)+10):(L1A_mode_find(end)+14)), '_1B_',...
            filesBulk.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
            '_isd','.nc');
    end
    if(~cnf.writting_flag(4)||~cnf.run_FFtime)
        L1A_date_find = strfind(filesBulk.filename_L1A,'_201');
        L1A_date_find = strfind(filesBulk.filename_L1A,'_202');
        L1A_mode_find = strfind(filesBulk.filename_L1A, 'S3NGT_SIRS_');
    end
    
    filesBulk.filename_L1BFF_SL = strcat(filesBulk.outputPath,'S3NGT_SIRS_FF_SL_',...
        filesBulk.filename_L1A((L1A_mode_find(1)+9):(L1A_mode_find(end)+13)), '_1B_',...
        filesBulk.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
        '_isd_long','.nc');
    filesBulk.filename_L1BFF_ML = strcat(filesBulk.outputPath,'S3NGT_SIRS_FF_ML_',...
        filesBulk.filename_L1A((L1A_mode_find(1)+10):(L1A_mode_find(end)+14)), '_1B_',...
        filesBulk.filename_L1A((L1A_date_find(1)+1):(L1A_date_find(1)+30)),...
        '_isd_long','.nc');

    if(cnf.verify_L1A)
        plots_L1A_file(filesBulk, cnf, chd, cst);
        plots_compare_L0_L1A(filesBulk, cnf, chd, cst);
    end
    %% ------------------------- DDP ---------GO!!!----------------
    if (cnf.run_DDP) 
        [chd] = check_dimensions(filesBulk,chd,'L1A');
        % ---------------------------------------------------------------
        filesBulk = DDP_chain (filesBulk, chd, cnf, cst, options);
        % ---------------------------------------------------------------
        filesBulk.LogMsg = [filesBulk.LogMsg,'DDP processed SUCCESSFULLY -- '];
        filesBulk.Success(2) = 1;
        if(cnf.verify_L1B)
            plots_L1B_file(filesBulk, cnf, chd, cst);   %plot(((L1B_written.data.range_ku_l1b_echo(2:115))*L1B_written.attributes.range_ku_l1b_echo.scale_factor)+L1B_written.attributes.range_ku_l1b_echo.add_offset)
        end                                             %mesh( (double(L1B_written.data.waveform_scale_factor_l1b_echo(2:115)) * ones(1,256)).' .* double(L1B_written.data.i2q2_meas_ku_l1b_echo(:,2:115)))
    else
        filesBulk.LogMsg = [filesBulk.LogMsg,'DDP not processed -- '];
        filesBulk.Success(2) = 0;
        
    end
    %% ------------------------- FFt ---------GO!!!----------------
    if (cnf.run_FFtime)
		filesBulk = FFt_chain (filesBulk, chd, cnf, cst, options);
        filesBulk.LogMsg = [filesBulk.LogMsg,'FF Time processed SUCCESSFULLY -- '];
        filesBulk.Success(3) = 1;
    else
        filesBulk.LogMsg = [filesBulk.LogMsg,'FF Time not processed -- '];
        filesBulk.Success(3) = 0;
    end    
    %% ------------------------- FFf ---------GO!!!----------------
    if (cnf.run_FFfreq) 
        filesBulk.LogMsg = [filesBulk.LogMsg,'FF Freq processed SUCCESSFULLY -- '];
        filesBulk.Success(4) = 1;
    else
        filesBulk.LogMsg = [filesBulk.LogMsg,'FF Freq not processed -- '];
        filesBulk.Success(4) = 0;
    end    
    %% ------------------------- LROS ---------GO!!!----------------
    if (cnf.run_LROS) 
        filesBulk.LogMsg = [filesBulk.LogMsg,'LROS processed SUCCESSFULLY -- '];
        filesBulk.Success(5) = 1;
    else
        filesBulk.LogMsg = [filesBulk.LogMsg,'LROS not processed -- '];
        filesBulk.Success(5) = 0;
        
    end
    %% ------------------------- L2 ---------GO!!!----------------
    if (cnf.run_L2GR) %JPLZ: a 'GR' was missing after L2?
        %filesBulk.LUT_f0_file = 'C:\Users\albert\Documents\WORK\CRISTAL/CODE/6.RETRACKERS/auxiliar_inputs/LUT_f0.mat';
        %filesBulk.LUT_f1_file = 'C:\Users\albert\Documents\WORK\CRISTAL/CODE/6.RETRACKERS/auxiliar_inputs/LUT_f1.mat';
        filesBulk.LUT_f0_file = [aux_folder 'LUT_f0.mat'];
        filesBulk.LUT_f1_file = [aux_folder 'LUT_f1.mat'];
        
        [chd] = check_dimensions(filesBulk,chd,'L1A');
        
        
        
        filesBulk = L2GR_chain(filesBulk, cnf, chd, cst);
        filesBulk.LogMsg = [filesBulk.LogMsg,'L2 Time processed SUCCESSFULLY -- '];
        filesBulk.Success(3) = 1;
    else
        filesBulk.LogMsg = [filesBulk.LogMsg,'L2 Time not processed -- '];
        filesBulk.Success(3) = 0;
    end 
    
    
    
    
    Logtime = datestr(now, 'yyyymmddTHHMMSS');
    fprintf(filesBulk.fid_log,'%s -> EVERYTHING OK!: \nConfiguration Options>%s \n%s \n',Logtime, num2str(filesBulk.Success), filesBulk.LogMsg);
	fclose (filesBulk.fid_log);
catch err
    if(log_flag)
        
        Logtime = datestr(now, 'yyyymmddTHHMMSS');
        fprintf(filesBulk.fid_log,'%s -> ERROR!: \nConfiguration Options>%s \n%s \n',Logtime, num2str(filesBulk.Success), filesBulk.LogMsg);
        fprintf(filesBulk.fid_log, '%s\n', err.getReport('extended', 'hyperlinks','off'));
    end
	fclose (filesBulk.fid_log);
    
    if  filesBulk.Success(1) == -1        
        fclose all;
    end
    
end