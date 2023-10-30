function [data,exit_flag] = read_alt_data_EM (filename_L1B, cnf_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows for reading altimetry data from any L1 processor we work
% with and store it in a common structure for the L2 processor
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 13/06/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       filename_L1B    =   L1B filename with the fullpath information
%       cnf_p       =   configuration parameters
% OUTPUT:
%       data        =   structure of data as defined by our L2 processor
%       exit_flag   =   flag indicating whether processing succesful 1 or
%       there is an error -1
%  
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - readL1B_CS2_ESA: read Cryosat-2 L1B products produced by ESA 
%                    (currently only developed for .DBL, takes into account netCDF .nc in future)
% - readL1B_CS2_ISD: read Cryosat-2 L1B products produced by ISD 
%                    (currently developed for netCDF .nc files using a la Sentinel-3 format)
% - readL1B_S3_ISD: read Sentinel-3 L1B products produced by ISD 
%                    (exactly the same as the one for CS-2 a la Sentinel-3 format)
% - readL1B_S6_ISD: read Sentinel-6 L1B products produced by ISD 
%                    (reading either the .mat files or the .nc with the final format for Sentinel-6 using PSD issue 1.3)
% - filter_L1B_data: filter the data by regions of interest and_/or number
%                    of looks within the stack 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
%
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: Based on read_alt_data.m defined different reading functions for
% the different missions and different potential processors 
% (based on the formating of the output products, .nc, DBL. and .mat)

%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(2+2*2))
    error('Wrong number of input parameters');   
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_L1BS',{''},@(x)ischar(x));
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.parse(varargin{:});
filename_L1BS=char(p.Results.filename_L1BS);
filename_mask_KML=char(p.Results.filename_mask_KML);
clear p;
exit_flag=1;

%% ------------------------------------------------------------------------- 
% Loading data L1B
% ------------------------------------------------------------------------- 
switch cnf_p.mission
    
    case {'CS2','CR2'}
    %% --------------------- CroySAT-2 ------------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                [data]=readL1B_CS2_ESA(filename_L1B);
            case 'ISD'
                [data]=readL1B_CS2_ISD(filename_L1B,cnf_p,'filename_L1BS',filename_L1BS);
        end
    case 'S3'
    %% --------------------- Sentinel-3 -----------------------------------
        switch cnf_p.L1proc
            case 'ESA'
                % TBD
                %[data]=readL1B_S3_ESA(filename_L1B);
            case 'ISD'
                [data]=readL1B_S3_ISD(filename_L1B,cnf_p,'filename_L1BS',filename_L1BS);
        end
    case {'S6','JCS'}
    %% --------------------- Sentinel-6 -----------------------------------        
        switch cnf_p.L1proc
            case 'ISD'
                [data]=readL1B_S6_ISD(filename_L1B,cnf_p,'filename_L1BS',filename_L1BS);
        end 
    otherwise
        error(strcat('Mission ',cnf_p.mission,' is not currently contemplated or not valid'));
end


%% -------------------- Inclusion of the seed information -----------------


%% -------------------- Filter data ---------------------------------------
%Bye geographic location using an external .kml or/and depending on the
%size of the associated stack
[data,flag]=filter_L1B_data (data,cnf_p,'filename_mask_KML',filename_mask_KML);
if flag==-1
    %track is not within the limits of the 
    exit_flag=flag;
    return
end

end

