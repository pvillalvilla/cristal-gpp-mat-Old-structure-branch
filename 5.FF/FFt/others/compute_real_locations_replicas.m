

%% Computation of the ghost peak replicas
%load the data
load('C:\Users\eduard.makhoul\isardSAT\projects\PICE\data\DATA - PICE\unit_test\results\KA_pm_300_around_TRP\PICE_1B_VHRt_KA_20180824T120034_20180824T12003_isd_interp.mat')
zp_fact_az = 16;
[arclen,az]=distance(lat_surf,lon_surf,ones(1,length(lat_surf))*chd.lat_trp,ones(1,length(lat_surf))*chd.lon_trp,[cst.semi_major_axis sqrt(1-(cst.semi_minor_axis/cst.semi_major_axis).^2)]);
[~,min_dist]=min(arclen);
distances_over_arclengths=arclen;
distances_over_arclengths(1:min_dist-1)=-1.0*distances_over_arclengths(1:min_dist-1);
distances_over_arclengths_interp=interp(distances_over_arclengths,zp_fact_az);


% call theoretical script to get the value of period of replicas
theoretical_replicas_sep_closed_burst;

N_replicas = floor((max(distances_over_arclengths_interp))./(Period_replicas));

cut_along_track = abs(wfm_AC_interp(1:length(distances_over_arclengths_interp),pos_max_across_focused)).^2;
[~,pos_max]=max(cut_along_track);

distance_replicas_left = [];
pos_replicas_left = [];
distance_replicas_right = [];
pos_replicas_right = [];

for i_replica=1:N_replicas
    %right
    idx_int=distances_over_arclengths_interp>(Period_replicas*(i_replica-1)+Period_replicas/2);
    [~,pos_replica]=max(cut_along_track(idx_int));
    if ~isempty(pos_replica)
        pos_replica =pos_replica+find(idx_int,1,'first')-1;
        pos_replicas_right = [pos_replicas_right,...
            pos_replica];
        distance_replicas_right=[distance_replicas_right,...
            abs(distances_over_arclengths_interp(pos_max)- ...
            distances_over_arclengths_interp(pos_replica))];
    end
    %left
    idx_int=distances_over_arclengths_interp<(-1.0*(Period_replicas*(i_replica-1)+Period_replicas/2));
    [~,pos_replica]=max(cut_along_track(idx_int));
    if ~isempty(pos_replica)
        pos_replica =pos_replica+find(idx_int,1,'first')-1;
        pos_replicas_left = [pos_replicas_left,...
            pos_replica];
        distance_replicas_left=[distance_replicas_left,...
            abs(distances_over_arclengths_interp(pos_max)- ...
            distances_over_arclengths_interp(pos_replica))];
    end
end

figure;
IRF=10*log10(cut_along_track./max(cut_along_track));
plot(distances_over_arclengths_interp,IRF,'-b')
ylim([min_image,max_image]);
title('PTR Along-track cut')
xlabel('Along-track distance w.r.t TRP location [m]'); ylabel('[dB]');
for i_replica=1:N_replicas
    text(distances_over_arclengths_interp(pos_replicas_left(i_replica)),...
        IRF(pos_replicas_left(i_replica)),strcat({'\Delta_X = '},num2str(distance_replicas_left(i_replica)),{' [m]'}));
    
    text(distances_over_arclengths_interp(pos_replicas_right(i_replica)),...
        IRF(pos_replicas_right(i_replica)),strcat({'\Delta_X = '},num2str(distance_replicas_right(i_replica)),{' [m]'}));
end
print('-dpng','-r150',strcat(plotsPath,'CUTs_PTR_interp_',filename_L1B_FFt,'.png'));



%IRF of the burst theoretical
t=(distances_over_arclengths_interp-distances_over_arclengths_interp(pos_max))/vg;%linspace(-T_exp/2,T_exp/2,N_points_representation);
% Burst PTR response (assuming no antenna weighting)
H_B = (sinc(BW_B*t));

H_A = zeros(1,length(t));
burst_idx=(0:1:N_replicas*2)-N_replicas*2/2;
for i_burst=1:N_replicas*2+1
H_A =H_A+abs(H_B).*abs(sinc(BW_A*(t-burst_idx(i_burst)/W_p)));
end
H_B = abs(H_B).^2;
H_A = abs(H_A).^2;

figure;
IRF=10*log10(cut_along_track./max(cut_along_track));
plot(distances_over_arclengths_interp,10*log10(H_B),'-r');
hold on;
plot(distances_over_arclengths_interp,10*log10(H_A),'-g');
plot(distances_over_arclengths_interp,IRF,'-k','LineWidth',2)
[res_az]=cc_ldb(20*log10(abs(wfm_AC_interp(:,pos_max_across_focused))),distances_over_arclengths_interp);
xlim(min_max_along); ylim([min_image,max_image]);
title(strcat('\delta_{al}=',num2str(res_az),'[m]'));
xlabel('Along-track distance w.r.t TRP location [m]'); ylabel('[dB]');
for i_replica=1:N_replicas
    text(distances_over_arclengths_interp(pos_replicas_left(i_replica)),...
        IRF(pos_replicas_left(i_replica)),strcat({'\Delta_X = '},num2str(distance_replicas_left(i_replica)),{' [m]'}));
    
    text(distances_over_arclengths_interp(pos_replicas_right(i_replica)),...
        IRF(pos_replicas_right(i_replica)),strcat({'\Delta_X = '},num2str(distance_replicas_right(i_replica)),{' [m]'}));
end
% 
% 
% figure;
% plot(distances_over_arclengths_interp,(abs(wfm_AC_interp(1:length(distances_over_arclengths_interp),pos_max_across_focused))./(10^(max_data/20))))
% [res_az]=cc_ldb(20*log10(abs(wfm_AC_interp(:,pos_max_across_focused))),distances_over_arclengths_interp);
% xlim(min_max_along); ylim([0,1]);
% title(strcat('\delta_{al}=',num2str(res_az),'[m]'));
% xlabel('Along-track distance w.r.t TRP location [m]'); ylabel('[dB]');