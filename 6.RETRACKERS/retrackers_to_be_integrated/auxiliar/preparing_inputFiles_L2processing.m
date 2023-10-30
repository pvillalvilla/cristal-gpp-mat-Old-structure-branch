function  [filename_L1B,filename_L1B_nopath,fileext_L1B,filename_L1BS]=preparing_inputFiles_L2processing(filesBulk, i_fileL1B_input, cnf_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code sets and prepares the input files definition and untar extraction for L1B
% and L1Bs file used in the L2 processing (used to have a more clear L2_processing.m function)
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 21/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -filesBulk    =   structure of input files within the folder where
%       to process the data (including the L1B as well as configuration/characterization files)
%       -i_fileL1B_input = index of the L1B
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       -filename_L1B        =   full filename of the input L1B including
%                                the path
%       -filename_L1B_nopath =   filename of the input L1B with neither
%                                full path nor file extension
%       -fileext_L1B         =   file extension of the input L1B file
%       -filename_L1BS       =   full filename of the input L1BS including
%                                the path 
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS:
%  It is assumed that within each input folder either for the L1B or for
%  the LBS there is only a single file with the same unique_identifier 
% (inital-final measurement/validty time), i.e., it must be avoided to
% have same file with .nc or .DBL and the same compressed file .TGZ
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 
%% ---------------- Handling input variables ------------------------------
if(nargin<3 || nargin>(3+3*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('targz_option_active_L1B',0);
p.addParamValue('targz_option_active_L1BS',0);
p.parse(varargin{:});
targz_option_active_L1B=(p.Results.targz_option_active_L1B);
targz_option_active_L1BS=(p.Results.targz_option_active_L1BS);
clear p;

%% --------------------- Prepare the input files definition ---------------
%--------------------------------------------------------------------------
% ------------- L1B file --------------------------------------------------
% -------------------------------------------------------------------------

filename_L1B=char([filesBulk.inputPath filesBulk.L1BFiles(i_fileL1B_input).name]);
if targz_option_active_L1B
    %untar the file    
    extracted_files=untar(filename_L1B,filesBulk.inputPath);
    if any(~cellfun(@isempty,strfind(extracted_files,'.DBL')))
        %DBL like files
        filename_L1B=strrep(filename_L1B,'TGZ','DBL');
    elseif any(~cellfun(@isempty,strfind(extracted_files,'.nc')))
        filename_L1B=strrep(filename_L1B,'TGZ','.nc');
    elseif any(~cellfun(@isempty,strfind(extracted_files,'.mat')))
        filename_L1B=strrep(filename_L1B,'TGZ','.mat');
    else
        error('Extracted untar files do not contain any valid L1B input file (.DBL or .nc)');
    end    
end

[~,filename_L1B_nopath,fileext_L1B]=fileparts(filename_L1B);
%-------------------- File unqiue identifier extraction -------------------
%identifier of a given file is the initial and final acq times
switch cnf_p.mission
    case 'CS2'
        %--------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                switch fileext_L1B
                    case '.DBL'
                        file_unique_id=filename_L1B_nopath(20:50);
                    case '.nc'
                        %TBD
                end
            case 'ISD'
                %following a la Sentinel-3 format as per SEOMs projects
                file_unique_id=filename_L1B_nopath(17:47);
        end
    case 'S3'
        %--------------------- Sentinel-3 -----------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                % TBD
            case 'ISD'
                file_unique_id=filename_L1B_nopath(17:end);                
        end
    case 'S6'
        %--------------------- Sentinel-6 -----------------------------------
        switch cnf_p.L1proc
            case 'ISD'
                switch fileext_L1B
                    case '.nc'
                        % according to name convention in IODD "JC-ID-ESA-GP-0175" issue 1.3
                        file_unique_id=filename_L1B_nopath;%filename_L1B_nopath(17:17+35);
                        %for the .mat file we assume the L1BS info is already
                        %saved in the .mat workspace
                    case '.mat'
                        file_unique_id = ''; %a way to overpass the extraction of 
                                             %the L1BS files as for
                                             %Sentinel-6 .mat we wil have
                                             %stack in the .mat workspace
                end
        end
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end

%--------------------------------------------------------------------------
%------------------------- L1BS file --------------------------------------
%--------------------------------------------------------------------------
if isfield(filesBulk,'L1BSFiles') && (~isempty(file_unique_id))
    %if it is not empty means we have provided some specific folder
    % we should look for the specific file linked to the L1B to be
    % processed
    %------------ Check whether file exist ----------------------
    %assume there is a single file definition:
    idx_int=find(~cellfun(@isempty,strfind({filesBulk.L1BSFiles(:).name},file_unique_id)), 1);
    
    if ~isempty(idx_int)
        filename_L1BS=char([filesBulk.inputPath filesBulk.L1BSFiles(idx_int).name]);
        if targz_option_active_L1BS
            %untar the file
            extracted_files=untar(filename_L1BS,filesBulk.inputPath);
            if any(~cellfun(@isempty,strfind(extracted_files,'.DBL')))
                %DBL like files
                filename_L1BS=strrep(filename_L1BS,'TGZ','DBL');
            elseif any(~cellfun(@isempty,strfind(extracted_files,'.nc')))
                filename_L1BS=strrep(filename_L1BS,'TGZ','.nc');
            elseif any(~cellfun(@isempty,strfind(extracted_files,'.mat')))
                filename_L1BS=strrep(filename_L1BS,'TGZ','.mat');
            else
                error('Extracted untar files do not contain any valid L1B input file (.DBL or .nc)');
            end
        end
        
    else
        filename_L1BS=''; %define the file as empty since there is no such L1BS related to the L1B
    end
else
    filename_L1BS='';
end


end

