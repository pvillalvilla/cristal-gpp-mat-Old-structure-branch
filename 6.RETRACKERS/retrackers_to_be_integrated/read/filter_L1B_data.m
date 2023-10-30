function [data,exit_flag] = filter_L1B_data (data,cnf_p,varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L.
% -------------------------------------------------------------------------
% This code allows for filtering out those L1B records specified either as
% by a KML file over a specific ROI and/or because the number of looks in stack
% is below a given threshold
%
% -------------------------------------------------------------------------
%
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         M�nica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 15/06/2016
% This software is built with internal funding
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       MANDATORY:
%           -data            =   data structure for L2 processing
%           -cnf_p           =   Structure with configuration parameters of the
%                                L2 processor
%       OPTIONAL:
%           -filename_mask_KML   =   Geograhpical mask kml file with the fullpath name
%
%
% OUTPUT:
%       data        =   filtered structure of data as defined by our L2 processor
%       exit_flag   =   indicating whether the processing was successful (1) or
%       not (-1)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - kml2lla: reads a given geographical mask in a kml file and provides the longitude, latitude and altitude information
%            to be used for filtering purposes
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


%% ---------------- Handling input variables ------------------------------
if(nargin<2 || nargin>(2+1*2))
    error('Wrong number of input parameters');
end
%option to include the L1B_S product to read the exact info from Look to
%Look within the stack
p = inputParser;
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.parse(varargin{:});
filename_mask_KML=char(p.Results.filename_mask_KML);
clear p;
exit_flag=1;

%% ----------------- GEOGRAPHICAL MASKING ---------------------------------
original_num_records=length(data.GEO.LAT);
if cnf_p.mask_ROI_flag
    product_mask  = kml2lla(filename_mask_KML);    
    ISD_lon_surf_bis=data.GEO.LON;
%     idx_lt_0= ISD_lon_surf_bis<0;
%     %longitudes +- values (-180,180)
%     if any(idx_lt_0)
%         ISD_lon_surf_bis(idx_lt_0)=ISD_lon_surf_bis(idx_lt_0)+360.0;
%     end    
    idx_int=inpolygon(ISD_lon_surf_bis,data.GEO.LAT,product_mask.coord(:,1),product_mask.coord(:,2));
    clear ISD_lon_surf_bis;
    [~,name_file,ext_file]=fileparts(filename_mask_KML);
    data.GLOBAL_ATT.DATA_FILE_INFO.geographical_mask_kml=[name_file ext_file];
    if ~any(idx_int)
        disp('Track outside the limits of the geographical mask')
        exit_flag=-1;
        return
    end
else
    idx_int=ones(1,original_num_records);
end
%% ----------------- Look number masking ----------------------------------
% mask out those records with number of looks below a given threshold
if cnf_p.mask_looks_flag
    idx_int=idx_int & (data.HRM.Neff >= cnf_p.Neff_thres);
end
idx_int=find(idx_int);
data.N_records=length(idx_int);


%% --------------------- FILTER DATA ---------------------------------------
% -----------------------------------------------------------------
% GEO: geographical information
% -----------------------------------------------------------------
data.GEO.TAI.total             =   data.GEO.TAI.total(idx_int);
data.GEO.TAI.days              =   data.GEO.TAI.days(idx_int);
data.GEO.TAI.secs              =   data.GEO.TAI.secs(idx_int);
data.GEO.TAI.microsecs         =   data.GEO.TAI.microsecs(idx_int);
data.GEO.LAT                   =   data.GEO.LAT(idx_int);
data.GEO.LON                   =   data.GEO.LON(idx_int); %defined between [0,360]
data.GEO.H_rate                =   data.GEO.H_rate(idx_int); % m/s
data.GEO.V                     =   data.GEO.V(idx_int);
data.GEO.H                     =   data.GEO.H(idx_int);
data.GEO.pitch                 =   data.GEO.pitch(idx_int);
data.GEO.roll                  =   data.GEO.roll(idx_int);
data.GEO.yaw                   =   data.GEO.yaw(idx_int);

% -----------------------------------------------------------------
% MEA: measurements
% -----------------------------------------------------------------
data.MEA.win_delay = data.MEA.win_delay(idx_int);

% ---------------------------------------------------------------------
% COR: Geophysical corrections
% ---------------------------------------------------------------------
if isfield(data,'COR')
    %-------- load individual corrections ---------------------------------
    % will be replicated for each data block of 20 surfaces
    data.COR.dry_trop                   =   data.COR.dry_trop(idx_int);
    data.COR.wet_trop                   =   data.COR.wet_trop(idx_int);
    data.COR.inv_bar                    =   data.COR.inv_bar(idx_int);
    data.COR.dac                        =   data.COR.dac(idx_int);
    data.COR.gim_ion                    =   data.COR.gim_ion(idx_int);
    data.COR.model_ion                  =   data.COR.model_ion(idx_int);
    data.COR.ocean_equilibrium_tide     =   data.COR.ocean_equilibrium_tide(idx_int);
    data.COR.ocean_longperiod_tide      =   data.COR.ocean_longperiod_tide(idx_int);
    data.COR.ocean_loading_tide         =   data.COR.ocean_loading_tide(idx_int);
    data.COR.solidearth_tide            =   data.COR.solidearth_tide(idx_int);
    data.COR.geocentric_polar_tide      =   data.COR.geocentric_polar_tide(idx_int);
    %---------- Combined corrections --------------------------------------
    data.COR.prop_GIM_ion   =   data.COR.prop_GIM_ion(idx_int); % not clear??
    data.COR.prop_Mod_ion   =   data.COR.prop_Mod_ion(idx_int); % not clear ??
    data.COR.surf_dac       =   data.COR.surf_dac(idx_int);
    data.COR.surf_invb      =   data.COR.surf_invb(idx_int);
    data.COR.geop           =   data.COR.geop(idx_int);
end
if isfield(data,'surf_type_flag')
    data.surf_type_flag=data.surf_type_flag(idx_int);
end



% -----------------------------------------------------------------
% HRM: High-resolution mode: Waveforms
% -----------------------------------------------------------------
% ------------------ Waveforms ------------------------------------
data.HRM.power_wav=data.HRM.power_wav(:,idx_int);
data.HRM.Neff       =   data.HRM.Neff(idx_int); %effective number of beams that form the stack including possuible looks that are set entirely to zero
data.HRM.FLAG.mlQ   =   data.HRM.FLAG.mlQ(idx_int); % if 0 no error if 1 a error ocurred in the stack or multilook
data.HRM.FLAG.pQ    =   data.HRM.FLAG.pQ(idx_int); % if 1 then ok, 0 error in the power
if isfield(data.HRM,'ThN')
    data.HRM.ThN = data.HRM.ThN(idx_int);
end
if isfield(data.HRM,'wfm_count')
    data.HRM.wfm_count  = data.HRM.wfm_count(idx_int);
end


% ----sigma0 scaling factor for conversion from Power to sigma0
%units
if isfield(data.HRM,'s0_sf')
    data.HRM.s0_sf=data.HRM.s0_sf(idx_int);
end

%--------------------------------------------------------------
%------------- Stack characterization parameters --------------
%--------------------------------------------------------------
%----- Dopppler mask ------------------------------------------
if isfield(data.HRM,'Doppler_mask')
    data.HRM.Doppler_mask   = data.HRM.Doppler_mask(:,idx_int);
end
%------------- Geometry-related parameters ------------------------
%exact
if all(isfield(data.HRM,{'pri_stack','V_stack'}))
    data.HRM.pri_stack = data.HRM.pri_stack(:,idx_int);
    data.HRM.V_stack = data.HRM.V_stack(:,idx_int);
end
if isfield(data.HRM,'beam_ang_stack')
    data.HRM.beam_ang_stack = data.HRM.beam_ang_stack(:,idx_int);
elseif isfield(data.HRM,'look_ang_stack')
    data.HRM.look_ang_stack = data.HRM.look_ang_stack(:,idx_int);
end
%approximate
if isfield(data.HRM,'pri_surf')
    data.HRM.pri_surf = data.HRM.pri_surf(idx_int); %in seconds
end
if all(isfield(data.HRM,{'pri_surf','pointing_ang_start_surf','pointing_ang_stop_surf','doppler_angle_start_surf','doppler_angle_stop_surf'}))
    data.HRM.pointing_ang_start_surf = data.HRM.pointing_ang_start_surf(idx_int); % in radians
    data.HRM.pointing_ang_stop_surf = data.HRM.pointing_ang_stop_surf(idx_int); % in radians
    data.HRM.doppler_angle_start_surf = data.HRM.doppler_angle_start_surf(idx_int); % in radians
    data.HRM.doppler_angle_stop_surf = data.HRM.doppler_angle_stop_surf(idx_int); % in radians
elseif all(isfield(data.HRM,{'look_ang_start_surf','look_ang_stop_surf'}))
    data.HRM.look_ang_start_surf = data.HRM.look_ang_start_surf(idx_int); % in radians
    data.HRM.look_ang_stop_surf = data.HRM.look_ang_stop_surf(idx_int); % in radians
end


 
end

