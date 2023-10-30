function nf_p = gen_nonfit_params_EM (data, cst, chd, cnf_L2) 
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This function initializates the non fitting parameters' structure
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Monica Roca / isardSAT
%
% Last revision:    Albert Garcia / isardSAT v1.1 10/06/2019
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       data    =   input data structure for the L2 processor
%       cnf_L2           =   configuration parameters
% OUTPUT:
%       ini_p           =   Fitting parameters seed

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
% - 
% Versions control:
% v1.0:
% v1.1: 2019/06/10 Updated struct and variable name
%% -------------- Global variables ----------------------------------------
% ---------------- Configuration parametres -------------------------------


chd.A_s2Ga=1.0196;
chd.alpha_ga=1.0/(2*(0.36012).^2);

%% ------------------- Height statistical distribution --------------------
%if ~cnf_L2.rou_flag
    nf_p.waveskew   =   0;
    nf_p.EMbias     =   0;
    nf_p.rou    =   cnf_L2.rou; % from Yuguan Liu
%end

%% ------------------ Array values ----------------------------------------
%% ------------------ Geometry/Projections --------------------------------
% Height
nf_p.h       =   data.alt;
% Spherical earth term
R_Earth =sqrt((cst.semi_major_axis^4*cos(data.lat.*pi/180).^2+cst.semi_minor_axis^4*sin(data.lat*pi/180).^2)./((cst.semi_major_axis*cos(data.lat.*pi/180)).^2+(cst.semi_minor_axis*sin(data.lat.*pi/180)).^2));
alpha        =   1 + (nf_p.h)./(R_Earth);

% Pitch and roll angles projection as distances
% assuming no change of sign is performed in the lecture
%take into account that the reference system for Model (x,y,z) is having z
%normal to the surface and x in the plane formed by satellite velocity and
%z : positive y component on the left_side; x positive forward movement
%satellite
switch cnf_L2.mission
    case {'CS2','CR2'}
        % Satellite-frame coordinates: x-along-track, z-opposite direction
        % to nadir and positive y on left-side of satellite
        %pitch is defined as positive/negative nose down/up respectively        
        nf_p.xp      =   data.alt.*(-1.0*data.pitch);
        %nf_p.xp      =   data.alt.*(cnf_L2.sign_pitch*1.0*data.pitch);
        %roll-angle is defined as positive/negative antenna right down/up        
        nf_p.yp      =   data.alt.*(data.roll);
        %nf_p.yp      =   data.alt.*(cnf_L2.sign_roll*data.roll);
    case 'S3'
        %TBD
        %pitch is defined as positive/negative nose up/down respectively (clockwise w.r.t Y-axis looking from origin )
        nf_p.xp      =   data.alt.*(data.pitch);
        %roll-angle is defined as positive/negative satellite-right-side
        %down/up (clockwise angle w.r.t x-axis looking from origin)
        nf_p.yp      =   data.alt.*(-1.0*data.roll);
    case 'S6'
        % Satellite-frame coordinates: x-along-track, z direction
        % to nadir and positive y on right-side of satellite
        %pitch is defined as positive/negative nose up/down respectively (clockwise w.r.t Y-axis looking from origin )
        nf_p.xp      =   data.alt.*(1.0*data.pitch);
        %roll-angle is defined as positive/negative satellite-right-side
        %down/up (clockwise angle w.r.t x-axis looking from origin)
        nf_p.yp      =   data.alt.*(data.roll);        
end


%Equivalent antenna beamwidth projection on ground
nf_p.alphax  =   (8*log(2))./(data.alt.*chd.antenna_beamwidth_along_track).^2;
nf_p.alphay  =   (8*log(2))./(data.alt.*chd.antenna_beamwidth_across_track).^2; %8 comes from fact that using the two way antenna pattern

%% --------------- Sampling spacing ---------------------------------------
% along-track
nf_p.Lx      =   cst.c*data.alt./(2*data.vel.*chd.freq*chd.N_ku_pulses_burst.*data.pri);
% across-track
%nf_p.Ly      =   sqrt(cst.c*data.alt./(alpha*chd.bw*cnf_L2.ZP));
% nf_p.Ly      =   sqrt(cst.c*data.alt./(alpha*chd.bw));
nf_p.Ly      =   sqrt(cst.c*data.alt./(alpha*chd.alt_freq*cnf_L2.ZP));
% vertical
%nf_p.Lz      =   0.5*cst.c/(chd.bw*cnf_L2.ZP);
% nf_p.Lz      =   0.5*cst.c/(chd.bw);
 nf_p.Lz      =   0.5*cst.c/(chd.alt_freq*cnf_L2.ZP);

nf_p.oversampling=chd.alt_freq/chd.bw; %due to differences in BW and fs (like in Sentinel-6)

% not clear from where this term is coming
nf_p.Lgamma  =   alpha./(2*data.alt.*nf_p.alphay);

%% --------------- Number of looks that form the stack --------------------
% it would be probably more convenient to define a unique structure value of Neff
% no need to use data.nb_stack_l1b_echo (will be kept as originally proposed just in case there are some future elaborations on the SARin data)
switch  cnf_L2.mode
    case 'SARin'
        nf_p.Neff    =   data.nb_stack_l1b_echo;
    case {'SAR','HR','RMC'}        
        nf_p.Neff    =   data.nb_stack_l1b_echo;
    otherwise
        error('No valid operational mode')
end

%% ------------------- Constant Values ------------------------------------
nf_p.Npulses    =   chd.N_ku_pulses_burst;
nf_p.alphag_a   =   chd.alpha_ga;
nf_p.alphag_r   =   chd.alpha_gr;
nf_p.A_s2Ga     =   chd.A_s2Ga;
nf_p.A_s2Gr     =   chd.A_s2Gr;
nf_p.Nbcycle    =   chd.N_bursts_cycle_sar;
nf_p.bw_Rx      =   chd.bw;
nf_p.Nsamples   =   chd.N_samples_sar*cnf_L2.ZP;


end