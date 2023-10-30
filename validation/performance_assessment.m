length_dataset=100000;

range_uncertainty_flag=1;
orbit_uncertainty_flag=1;
theta_uncertainty_flag=1;
c_cst = 299792458;

x_cog_ant= 1.67;
z_cog_ant = 0.55;
y_cog_ant = 1.67/2;
%time tag
pitch_uncert=0.000555556/180*pi; %2 arcseconds
roll_uncert=0.000555556/180*pi;%2 arcseconds
z_cog_ant_corr_pitch    = - x_cog_ant * sin(pitch_uncert) - z_cog_ant * cos(pitch_uncert);
        delta_pitch             = 2 * z_cog_ant_corr_pitch / c_cst;
        z_cog_ant_corr_roll     = y_cog_ant * sin(roll_uncert);
        delta_roll              = 2 * z_cog_ant_corr_roll / c_cst;
%         L1A.win_delay = L1A.win_delay + delta_pitch + delta_roll;



Mission='CRISTAL'; %[CS2,CRISTAL]

switch Mission
    case 'CS2'
        snow_uncertainty = 6e-2;
        orbit_uncertainty = 2e-2;
        range_retrack_uncertainty = 4.9e-2;
        theta_uncertainty = 0.0112/180*pi; %40 arcsec
        
    case 'CRISTAL'
        snow_uncertainty = 1.2e-2;
        range_retrack_uncertainty = 0.8e-2;
        orbit_uncertainty = 1.9e-2;
        theta_uncertainty = 0.00638889/180*pi; %23 arcsec
        
end


Re                  = 6367.444657e3;%[m]
orbit_uncertainty   = orbit_uncertainty_flag*orbit_uncertainty; %[m]
mean_orbit          = 7.200000000174579e+05;
h                   = mean_orbit+orbit_uncertainty*randn(1,length_dataset); %[m]
geocorr_uncertainty = range_uncertainty_flag*sqrt(1e-2^2+3e-2^2+1.1e-2^2+snow_uncertainty); %[cm] dry, wet, iono, snow
geocorr_mean        = -2.3;
delta_r_li          = geocorr_mean+geocorr_uncertainty*randn(1,length_dataset);

ret_uncertainty     = range_uncertainty_flag*range_retrack_uncertainty; %[m] the uncertainty from the retracker (it should be 0 with Swath)
theta_uncertainty   = theta_uncertainty_flag*theta_uncertainty;%0.00638889/180*pi; %0.0112º for 40 arcsec, for 0.00638889º for  23 arcsec
range_variation     = [10 100 200 500 800]; %[m]
theta_min = -0.9/180*pi;
theta_max = 0.9/180*pi;

% PropError((h-(Rret+delta_r_li).*cos(theta)+Re.*(1-cos(h/Re.*theta))).*cos(h/Re.*theta),[h Rret theta],[mean_orbit mean_orbit-100 0],[orbit_uncertainty ret_uncertainty theta_uncertainty])

% loop over different slopes by modifying the max_elevation, keeping the
% theta and the Xtrack distance
% figure;
range_bias=linspace(000,1000,10);
for i_case=1:length(range_bias)
    
    min_range = mean_orbit -500 - range_bias(i_case);
    max_range = min_range -500 + 256 - range_bias(i_case);
   
    Rret = linspace(max_range,min_range,length_dataset)+ ret_uncertainty*randn(1,length_dataset); %[m]
    theta = (linspace(theta_min,theta_max,length_dataset) + theta_uncertainty*randn(1,length_dataset));  % [rad]
    
    LIH =(h-(Rret+delta_r_li).*cos(theta)+Re.*(1-cos(h/Re.*theta))).*cos(h/Re.*theta); %[m]
    
    
    %ideal case
    Rret = linspace(max_range,min_range,length_dataset); %[m]
    theta = (linspace(theta_min,theta_max,length_dataset)); % [rad]
    x_track_distance=  sin(theta) .* Rret;
    
    LIH_ideal=(h-(Rret+geocorr_mean).*cos(theta)+Re.*(1-cos(h/Re.*theta))).*cos(h/Re.*theta); %[m]
    
    slope(i_case)=atan((max(LIH_ideal)- min(LIH_ideal))/(max(x_track_distance)- min(x_track_distance)))*180/pi;
    
    
%     subplot(2,1,1); 
    figure(24); plot(x_track_distance/1000,LIH); hold all;% plot(x_track_distance(1:1000:end)/1000,LIH_ideal(1:1000:end),'k.');
    figlabels('Across track distance [km]','Elevation [m]','','Combined Uncertainity ',20)
%     subplot(2,1,2); 
%     plot(x_track_distance/1000,LIH-LIH_ideal,'.'); hold all;
%     figlabels('Across track distance [km]','Elevation uncertainity [m]','',' ',20)
%     
    uncertainity(i_case)=std(LIH-LIH_ideal);
    
    for kk=1:100
        
        uncert(kk)=std(LIH((kk-1)*1000+1:kk*1000)-LIH_ideal((kk-1)*1000+1:kk*1000));
        
    end
    
    
    figure(25); hold all;plot(x_track_distance(500:1000:end),uncert);
    p(i_case,:)=polyfit(x_track_distance(500:1000:100000/2)/1000,uncert(1:length(x_track_distance(500:1000:100000/2))),1);
end

legend([num2str(slope(1),'%0.02f') 'º'], ...
[num2str(slope(2),'%0.02f') 'º' ], ...
[num2str(slope(3),'%0.02f') 'º' ], ...
[num2str(slope(4),'%0.02f') 'º' ], ...
[num2str(slope(5),'%0.02f') 'º' ], ...
[num2str(slope(6),'%0.02f') 'º' ], ...
[num2str(slope(7),'%0.02f') 'º' ], ...
[num2str(slope(8),'%0.02f') 'º' ], ...
[num2str(slope(9),'%0.02f') 'º' ], ...
[num2str(slope(10),'%0.02f') 'º' ] )

disp(['slopes: ' num2str(slope)])
disp(['uncertainity: ' num2str(uncertainity)])
disp(['uncertainity: ' num2str(mean(uncertainity))])

%%%%

cs=c_cst/1.281 ; % corresponds to a snow layer with a density of 320 kg m-3, https://tc.copernicus.org/preprints/tc-2018-164/tc-2018-164.pdf
ros= 320; %kg/m3
uncert_ros=0.01;
uncert_hs=0.012;
hs=0.2;
hs=linspace(0,1,1000);
uncert_delta_hs=sqrt(((c_cst/cs-1)^2*uncert_hs^2)+(-hs*c_cst/cs^2).^2*(-c_cst*(14*ros+17)/(20*(7/(10*ros^2)+17/(10*ros)+1)^(3/2)))^2*uncert_ros^2)*100;
freeboard_uncert=sqrt(3.3/4^2+2.88/4^2+1.13^2*720e5*((23e-6)^2*(0*pi/180).^2+(20e-6)^2*(0*pi/180).^2)+uncert_delta_hs/4.^2);
figure;plot(hs,freeboard_uncert);
figlabels('Snow depth [m]','Sea Ice Freeboard Uncertainty [cm]','','Freeboard Uncertainty ',20)
mean(freeboard_uncert)
