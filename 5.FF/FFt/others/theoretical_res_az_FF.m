 

c_cst = 299792458;
PRI=1/16.258e3;% 1/18e3 CB, 1/16.258e3 OB
H_orb=[688.652e3 698.950e3 717.227e3];
v_sat=[7547.718512 7534.8250441 7514.015129];
for kk=1:3
disp(kk)

% % ------------ KU --------------------------------
 wv_length_ku = 0.022084158968692; %3e8/13.575e9 %
%  v_sat=7.593524114339879e+03;%7.530281573680895e+03%
%  H_orb=7.200000000174579e+05;%717e3%
 delta_beam3db=1.05*pi/180;%1.095*pi/180%
 R_earth=6367.444657e3;
 
 
alpha_r=1+H_orb(kk)/R_earth;
%ground velocity of the beam
vg=v_sat(kk)/alpha_r;
v_r = sqrt(vg*v_sat(kk));
%time of exposure
Texp=delta_beam3db*H_orb(kk)/(vg);
%synthetic angle (different from aperture angle due to earth curvature)
Delta_syn=v_sat(kk)/vg*delta_beam3db;
%theoretical resolution
res_az=0.886*wv_length_ku/(2*Delta_syn);
BW_A=2*vg/wv_length_ku*Delta_syn; 
disp(['res_az FF SL Ku ' num2str(res_az) ' meters']);
% disp(['Tint Ku ' num2str(Texp) ' seconds']);
res_az_UF=0.886*wv_length_ku*H_orb(kk)/(2*v_sat(kk)*64*PRI);
disp(['res_az UF Ku ' num2str(res_az_UF) ' meters']);
% --------------- Ka ----------------------------
wv_length_ku = 0.008385803020979; %3e8/13.575e9 %
% v_sat=7.593524114339879e+03;%7.530281573680895e+03%
% H_orb=7.200000000174579e+05;%717e3%
delta_beam3db=0.42*pi/180;%1.095*pi/180%
R_earth=6367.444657e3;


alpha_r=1+H_orb(kk)/R_earth;
%ground velocity of the beam
vg=v_sat(kk)/alpha_r;
v_r = sqrt(vg*v_sat(kk));
%time of exposure
Texp=delta_beam3db*H_orb(kk)/(vg);
%synthetic angle (different from aperture angle due to earth curvature)
Delta_syn=v_sat(kk)/vg*delta_beam3db;
%theoretical resolution
res_az=0.886*wv_length_ku/(2*Delta_syn);
BW_A=2*vg/wv_length_ku*Delta_syn;
disp(['res_az FF SL Ka ' num2str(res_az) ' meters']);
% disp(['Tint Ka ' num2str(Texp) ' seconds']);


res_az_UF=0.886*wv_length_ku*H_orb(kk)/(2*v_sat(kk)*64*PRI);
disp(['res_az UF Ka ' num2str(res_az_UF) ' meters']);


tau_c=1/600e6;  
res_az_lr=2*sqrt(R_earth/(R_earth+H_orb(kk))*c_cst*tau_c*H_orb(kk));
disp(['res_az LR Ku and Ka ' num2str(res_az_lr) ' meters']);
end

% if res_az is given

%  Texp = 0.886 * wv_length_ku/(2*res_az)*H_orb/v_sat(kk);
 
 
