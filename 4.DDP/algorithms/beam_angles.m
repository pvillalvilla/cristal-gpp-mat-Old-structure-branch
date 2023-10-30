%% HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% Polar ICE topography mission
% aligned with isardSAT_GPPICE_ATBD_v0a
% ---------------------------------------------------------
% Objective: The purpose of the beam angles is to compute, for every burst, 
% the angles between the nadir direction and the direction defined by 
% the satellite location and each surface location under the satellite's boresight. 
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT
%            Eduard Makhoul / isardSAT
%            JP López-Zaragoza / isardSAT
% ----------------------------------------------------------
% Version  record
% 1.0 2018/07/18 First version imported from Dedop rev 125
% 1.1 2019/01/22 introduced the floor to force to 33 when having 65 pulses for beam_ang_index
% (JP López-Zaragoza), 2.0 2023/03/15, Changed the beam steering method. Now we don't calculate the angle 
%                                      just for the beams between qmin and qmax. We calculate the angle between 
%                                      all the surfaces created and the current burst satellite position and then take the ones in which we are interested.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [L1A]       = beam_angles (L1A,L1BS, N_total_surf_loc,i_surf,i_burst_focussed,N_bursts, cnf, chd, cst)

L1A.beam_ang_nadir_index = 0;
i_surf_stacked=i_surf;
beam_ang                        = 1./zeros(1,chd.N_pulses_burst*cnf.zp_fact_azimut);
L1A.beam_ang                    = 1./zeros(1,chd.N_pulses_burst*cnf.zp_fact_azimut);
surf_loc_index                  = zeros(1,chd.N_pulses_burst*cnf.zp_fact_azimut);
L1A.surf_loc_index              = zeros(1,chd.N_pulses_burst*cnf.zp_fact_azimut);

vel_sat = [L1A.x_vel_sat_sar,L1A.y_vel_sat_sar,L1A.z_vel_sat_sar];
bursts_indexes=1:chd.N_pulses_burst;

% For all the ground surfaces created so far, we compute the angle to the satellite position in the current burst
for i_surf=i_surf_stacked:N_total_surf_loc
    
    v(i_surf,:) = [L1BS(i_surf).x_surf, L1BS(i_surf).y_surf, L1BS(i_surf).z_surf] - [L1A.x_sar_sat, L1A.y_sar_sat, L1A.z_sar_sat];
    w = vel_sat;
    
    b(i_surf) = acos ((v(i_surf,:)*w') / (norm(v(i_surf,:)) * norm(w)));
end

% We compute the angular azimuth resolution and look for the surface which
% points closer to the nadir satellite position in the current burst (surf_doppler0)
ang_az_beam_res = asin (cst.c / (2 * chd.freq * norm(vel_sat) * (chd.burst_duration_sar(1))))/cnf.zp_fact_azimut;
[doppler_diff,surf_doppler0] = min(abs((b(:)-cst.pi/2)));
if(doppler_diff>ang_az_beam_res)
    doppler_shift_beams=round(doppler_diff/ang_az_beam_res);
    surf_doppler0=surf_doppler0-+doppler_shift_beams;
end
shift_backward = 0;
shift_forward = 0;

% We calculate which are the initial and final surfaces the beams can see
init_surf = max([surf_doppler0-32,1]);
end_surf = min([surf_doppler0+31,N_total_surf_loc]);

% % when at burst 72 the 1st one is a zero
% if b(init_surf)==0
%     firstnonzero_init_surf=find(b);
%     init_surf=firstnonzero_init_surf(1);
% end

% We look for where the pi/2 surface is, and if outside of the angles the
% beams see, calculate how many positions/beams it would correspond, to do a shift
if (surf_doppler0-32<i_surf_stacked)
    %     init_surf = 1;
    shift_backward=surf_doppler0-32-i_surf_stacked;
end
if(surf_doppler0+31>N_total_surf_loc)
    %     end_surf = N_total_surf_loc;
    shift_forward=surf_doppler0+31-N_total_surf_loc;
end

% We save the angles and surfaces the beams can see in an array with the equivalent beam position within the burst
try
    beam_ang(1-shift_backward:64-shift_forward)= b(init_surf:end_surf);
    surf_loc_index(1-shift_backward:64-shift_forward)= init_surf:end_surf;
end

% Having problems now when in last bursts when 1-shift_backward=-1.
% Small trap here
if sum(isinf(beam_ang))==64
    beam_ang(min(1-shift_backward,1):64-shift_forward)= b(init_surf:end_surf);
    surf_loc_index(min(1-shift_backward,1):64-shift_forward)= init_surf:end_surf;
end
% % Small trap to shift the beams of the first bursts (1:N_total_surf_loc), since they all point to the same surfaces even when the satellite is advancing
% if i_burst_focussed<N_total_surf_loc+10
%     ang_az_beam_res = asin (cst.c / (2 * chd.freq * norm(vel_sat) * (chd.burst_duration_sar(1))))/cnf.zp_fact_azimut;
%     shift_val=floor( abs(beam_ang(33)-cst.pi/2) / ang_az_beam_res );
%     clear beam_ang_aux
%     beam_ang_aux=circshift(beam_ang,shift_val);
%     beam_ang_aux(1:32)=Inf;
%     beam_ang=beam_ang_aux;
%     clear surf_loc_index_aux
%     surf_loc_index_aux=circshift(surf_loc_index,shift_val);
%     surf_loc_index_aux(1:32)=0;
%     surf_loc_index=surf_loc_index_aux;
% end

L1A.beam_ang=beam_ang;
L1A.surf_loc_index=surf_loc_index;

noninf_indexes=~isinf(beam_ang);
noninf_indexes=find(noninf_indexes);

L1A.beam_ang_index=33;
%L1A.start_beam      = 1-shift_backward;
L1A.start_beam      = noninf_indexes(1);
%L1A.end_beam        = 64-shift_forward;
L1A.end_beam        = noninf_indexes(end);
L1A.beam_index      = L1A.start_beam:L1A.end_beam;
L1A.N_beams         = length(L1A.beam_index);

%If beam 33 is empty, we assign it the value it would have, even if it does not see any surface. We need it to
%have a value if we want in the azimuth processing to focus the burst correctly using the approximate method
if(~isfinite(beam_ang(L1A.beam_ang_index)))
    [doppler_diff,surf_doppler0] = min(abs((b(:)-cst.pi/2)));
    L1A.beam_ang(L1A.beam_ang_index)=b(surf_doppler0)+doppler_shift_beams*ang_az_beam_res;
    
end

end %end function



