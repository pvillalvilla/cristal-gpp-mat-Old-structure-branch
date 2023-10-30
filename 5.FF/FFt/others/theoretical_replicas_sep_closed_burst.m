PRF=18e3; %18.181e3%
BRI = 0.0125;
N_pulses=64;
% %KU
% wv_length_ku = 0.022084158968692; %3e8/13.575e9 %
% v_sat=7.593524114339879e+03;%7.530281573680895e+03%
% H_orb=7.200000000174579e+05;%717e3%
% delta_beam3db=1.06*pi/180;%1.095*pi/180%
% R_earth=6367.444657e3;
%Ka
wv_length_ku = 0.008385803020979; %3e8/13.575e9 %
v_sat=7.593524114339879e+03;%7.530281573680895e+03%
H_orb=7.200000000174579e+05;%717e3%
delta_beam3db=1.06*pi/180;%1.095*pi/180%
R_earth=6367.444657e3;



alpha_r=1+H_orb/R_earth;
%ground velocity of the beam
vg=v_sat/alpha_r;
v_r = sqrt(vg*v_sat);

%Doppler bandwidth giving resolution
Delta_syn=v_sat/vg*delta_beam3db;
%theoretical resolution
res_az=0.886*wv_length_ku/(2*Delta_syn);
BW_A=2*vg/wv_length_ku*Delta_syn;

%Azimuth rate
FM = 2*v_r^2/(wv_length_ku*H_orb);
%spectral burst cycle period
W_p = FM*BRI;

% Period of ghost replicas
Period_replicas = 1/W_p*vg;

% Burst bandwidth
BW_B = FM*(N_pulses/PRF);



