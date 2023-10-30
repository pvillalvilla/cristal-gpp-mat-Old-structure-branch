%% PLOTTING
set(0,'defaultFigureVisible','on');
zp_fact_trp=1;
h=figure;
subplot(2,2,1);
k=surf(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS_buffer(i_surf_stacked).N_beams_stack,abs((fft(L1BS_buffer(i_surf_stacked).beams_surf(:,:).', chd.N_samples_sar *zp_fact_trp))).'.^2);
set(k, 'edgecolor','none');view(10,50);
set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',12);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',12);
figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' before alignment'],12)

subplot(2,2,3);
plot(L1BS_buffer(i_surf_stacked).shift-1,1:L1BS_buffer(i_surf_stacked).N_beams_stack);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',12);
set(gca,'XLim',[- chd.N_samples_sar   chd.N_samples_sar ],'FontSize',12);
figlabels('Shift [samples]','Beams','','Range corrections',12)
subplot(2,2,2);
k=surf(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS_buffer(i_surf_stacked).N_beams_stack,abs((fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).', chd.N_samples_sar *zp_fact_trp))).'.^2);
set(k, 'edgecolor','none');
%view(10,50);
set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',12);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',12);
figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' aligned'],12)

subplot(2,2,4);
zp_fact_trp=512;
[max_val,max_pos]=max(abs((fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).', chd.N_samples_sar *zp_fact_trp))).^2);
max_pos = max_pos/zp_fact_trp-1/zp_fact_trp; %reference from [0 to N] IMPORTANT!!



% subplot(2,1,1);
% k=surf(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS_buffer(i_surf_stacked).N_beams_stack,abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).', chd.N_samples_sar *zp_fact_trp),1)).'.^2);
% set(k, 'edgecolor','none');
% %view(10,50);
% set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',12);
% set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',12);
% figlabels('Samples','Beams','Power',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' aligned'],12)
% subplot(2,1,2);
% zp_fact_trp=512;
% % [max_val,max_pos]=max(abs(fftshift(ifft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).', chd.N_samples_sar *zp_fact_trp),1)).^2);
% [max_val,max_pos]=max(abs(fftshift(fft(L1BS_buffer(i_surf_stacked).beam_geo_corr(:,:).', chd.N_samples_sar *zp_fact_trp),1)).^2);
% max_pos = max_pos/zp_fact_trp-1/zp_fact_trp; %reference from [0 to N] IMPORTANT!!



%% GEOPHYSICAL CORRECTIONS

alt_rate_sat=L1BS_buffer(i_surf_stacked).alt_rate_sat;


%% IN CODE %%
geophysical_correction=wet+dry+iono+solid_earth+geocentric_tide+ocean_loading; %09/04/2016


%%
[central, valid_beams] = TRP_beam_analysis(L1BS_buffer(i_surf_stacked));

% init_beam=25; end_beam=217; % %07/11/2015 CR2
% init_beam=29; end_beam=222; % %30/04/2012 CR2
% init_beam=21; end_beam=216; % %29/05/2012 CR2
% init_beam=37; end_beam=217; % %10/05/2013 CR2
% valid_beams=init_beam:end_beam;
% % valid_beams = valid_beams(find(valid_beams~=150));
% valid_beams = valid_beams(find(valid_beams~=115));

plot(max_pos,1:L1BS_buffer(i_surf_stacked).N_beams_stack);
set(gca,'YLim',[1 L1BS_buffer(i_surf_stacked).N_beams_stack],'FontSize',12);
set(gca,'XLim',[max_pos(floor(L1BS_buffer(i_surf_stacked).N_beams_stack/2))-5 max_pos(floor(L1BS_buffer(i_surf_stacked).N_beams_stack/2))+5],'FontSize',12);
% set(gca,'XLim',[1  chd.N_samples_sar ],'FontSize',12);

[Slope_coef_A1,SS,MUMU] = polyfit(valid_beams,max_pos(valid_beams),1);
Slope_A1 = polyval(Slope_coef_A1, valid_beams,SS,MUMU); %Slope_coef_A1(1) related with datation error
Stack_noise_A1=std(max_pos(valid_beams)-Slope_A1)/chd.bw*cst.c/2; % units [m]
Stack_alignment_A1=(Slope_A1(end)-Slope_A1(1))/length(valid_beams)/chd.bw*cst.c/2; % units [m/beam]
hold all;
plot(Slope_A1,(valid_beams))
figlabels('Samples','Beams','',['Stack #' num2str(L1BS_buffer(i_surf_stacked).surf_counter) ' alignment: ' num2str(Stack_alignment_A1) ' [m/beam] and noise: ' num2str(Stack_noise_A1) ' [m]; zero padding: ' num2str(zp_fact_trp)],12)
hold off
saveas (h,[filesBulk.resultPath 'Stack_alignment_' num2str(L1BS_buffer(i_surf_stacked).surf_counter) '.png']);
saveas (h,[filesBulk.resultPath 'Stack_alignment_' num2str(L1BS_buffer(i_surf_stacked).surf_counter) '.fig']);


%% COMPUTING RESULTS
% [L1B]=readanyNETCDF_V1('./S3A_SR_1_SRA____20160822T193828_20160822T194828_20160822T204103_0599_008_013______SVL_O_NR_001.SEN3/measurement.nc');
% [L2]=readanyNETCDF_V1('./L2/standard_measurement_NT.nc');
% time_1Hz = double(L2.data.time_01).';
% time_20Hz = double(L2.data.time_20_ku).';
% [~,pos] = min(abs(time_1Hz-L1BS_buffer(i_surf_stacked).time_surf));
% [~,pos20] = min(abs(time_20Hz-L1BS_buffer(i_surf_stacked).time_surf));
% [geophysical_correction] = compute_L2_corrections(L2,pos,'linear');

if chd.lat_trp>78
    TRP_int_delay = 9.88; % Svalbard
elseif chd.lat_trp<35.4
    TRP_int_delay = 4.954; % Crete
else
    TRP_int_delay = 6.533; % University of Crete, Old location
end

%  TRP_int_delay   = 6.533; % University of Crete, Old location
% TRP_int_delay   = 5.160; % Provided by Stelios 31/05/2016
%TRP_int_delay   = 4.954; % Provided by Stelios 22/02/2016
 %TRP_int_delay   = 9.88;  %Svalbard
           %TRP insitu  Wet,  Dry ,   Iono 
% 
% geophysical_correction=-2.555; %10/05/2013 CR2 FBR
geophysical_correction=-2.374681; %07/11/2015 CR2 FBR
% geophysical_correction=-2.385; %30/04/2012 CR2 FBR
% geophysical_correction=-2.355; %29/05/2012 CR2 FBR
% geophysical_correction=(-52.1-2057.4-55.24)*1e-3; %07/11/2015 CR2 FBR


p = lla2ecef([chd.lat_trp,chd.lon_trp,chd.alt_trp],cst.flat_coeff,cst.semi_major_axis);
    x_TRP = p(:,1).';
    y_TRP = p(:,2).';
    z_TRP = p(:,3).';
TRP_coord = [x_TRP,y_TRP,z_TRP];
% TRP_coord = [L1BS_buffer(i_surf_stacked).x_surf,L1BS_buffer(i_surf_stacked).y_surf,L1BS_buffer(i_surf_stacked).z_surf];
window_corr = L1BS_buffer(i_surf_stacked).win_delay_surf;

if cnf.window_delay_source == 0 || strcmp(mission,'CR2') %L0
    time_stack = window_corr +((max_pos-(chd.N_samples_sar/2))/chd.bw);
elseif cnf.window_delay_source == 1 %L1A
    time_stack = window_corr +((max_pos-(44))/chd.bw);
end
range_stack = time_stack *cst.c/2  + geophysical_correction-TRP_int_delay/2;
range_without_slant =range_stack-(L1BS_buffer(i_surf_stacked).slant_range_corr(:).')/chd.bw*cst.c/2;

for i_beam=1:length(range_without_slant)
    time_beam(i_beam)= L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(i_beam)).time_sar_ku;
    SAT_coord = [L1BS_buffer(i_surf_stacked).x_sar_sat_beam(i_beam),L1BS_buffer(i_surf_stacked).y_sar_sat_beam(i_beam),L1BS_buffer(i_surf_stacked).z_sar_sat_beam(i_beam)];
    dif       = (SAT_coord(1)-TRP_coord(1))^2 + (SAT_coord(2)-TRP_coord(2))^2 + (SAT_coord(3)-TRP_coord(3))^2;
    theo_range(i_beam) =  sqrt(dif);
end
theo_range_aligned=theo_range+(L1BS_buffer(i_surf_stacked).slant_range_corr(:).')/chd.bw*cst.c/2;

% time interpolation:
interp_fact=0.001;
delta_time=mean(diff(time_beam));
delta_time_interp= delta_time*interp_fact;
times_STACK = min(time_beam):delta_time_interp:max(time_beam);

[range_theo_coef,SS,MUMU] = polyfit(time_beam,theo_range,2);
range_theo_interp = polyval(range_theo_coef, times_STACK,SS,MUMU); 
range_theo_interp_2 = interp1(time_beam, theo_range,times_STACK,'splines');
[range_meas_coef,SS,MUMU] = polyfit(time_beam(valid_beams),range_without_slant(valid_beams),2);
range_meas_interp = polyval(range_meas_coef, times_STACK,SS,MUMU); 
range_meas_interp_2 = interp1(time_beam(valid_beams), range_without_slant(valid_beams),times_STACK);

deltaR = fitfunctions_yaxis (range_meas_interp, range_theo_interp);
deltaT = -1*fitfunctions_xaxis ((range_meas_interp-deltaR), range_theo_interp); %--20110711-- reduce Y diff
if deltaT==0;
  deltaT = fitfunctions_xaxis (range_theo_interp, range_meas_interp-deltaR);
end   
range_error_A1_v2 = deltaR; %[m] 
datation_error_A1_v2 = deltaT*delta_time_interp *1e6;%[micros]
[range_theo,pos_theo]=min(range_theo_interp);
[range_meas,pos_meas]=min(range_meas_interp);
datation_error_A1_v1 = (times_STACK(pos_meas)-times_STACK(pos_theo))*1e6; %[micros]
range_error_A1_v1=range_meas-range_theo;
disp(['Range error fitting: ' num2str(range_error_A1_v2*1e3) ' mm & minimum value: ' num2str(range_error_A1_v1*1e3) ' mm & aligned ranges: ' num2str(mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3) ' mm']);
disp(['Datation error fitting: ' num2str(datation_error_A1_v2) ' microseconds & minimum value: ' num2str(datation_error_A1_v1) ' microseconds']);
disp(['alignment: ' num2str(Stack_alignment_A1*1000) ' [mm/beam] noise: ' num2str(Stack_noise_A1*1000) ' [mm]; zero padding: ' num2str(zp_fact_trp)]);

%% View side parameter%%
switch mission
    case 'S3'
        if ( sum([L1BS_buffer(1:5).lat_surf]) > sum([L1BS_buffer(6:10).lat_surf]) ) %Descending
            if (L1BS_buffer(i_surf_stacked).lon_surf > chd.lon_trp)
                sat_view = 'Right';
            else
                sat_view = 'Left';
            end
        else %Ascending
            if (L1BS_buffer(i_surf_stacked).lon_surf > chd.lon_trp)
                sat_view = 'Left';
            else
                sat_view = 'Right';
            end
        end
    case 'CR2'
        CR2_mode = filename_L1A(13:15);
        %{
        if strcmp(CR2_mode,'SIN')
            if ( sum([L1A_buffer(1:5).lat_surf]) > sum([L1A_buffer(6:10).lat_surf]) ) %Descending
                if (L1A_buffer(pos).lon_surf > chd.lon_trp)
                    sat_view = 'Right';
                else
                    sat_view = 'Left';
                end
            else %Ascending
                if (L1A_buffer(pos).lon_surf > chd.lon_trp)
                    sat_view = 'Left';
                else
                    sat_view = 'Right';
                end
            end
        elseif strcmp(CR2_mode,'SAR')
            if ( sum([L1A_buffer(1:5).lat_sar_surf]) > sum([L1A_buffer(6:10).lat_sar_surf]) ) %Descending
                if (L1A_buffer(pos).lon_sar_surf > chd.lon_trp)
                    sat_view = 'Right';
                else
                    sat_view = 'Left';
                end
            else %Ascending
                if (L1A_buffer(pos).lon_sar_surf > chd.lon_trp)
                    sat_view = 'Left';
                else
                    sat_view = 'Right';
                end
            end    
        end
        %}
end



%% Writing 
fidResults=fopen('../results_Dedop.csv','a+');
cycle='';
%{
if chd.lat_trp>78
    results_file=dir('./../*Svalbard.csv');%Svalbard
else
    results_file=dir('./../*Crete.csv');%Crete
end

fidResults=fopen(['../'  results_file.name],'a+');
%}
fprintf(fidResults,'Cycle;Data; Range error fitting [mm]; Range error minimum value [mm]; Range error aligned ranges [mm]; Datation error fitting [microseconds]; Datation error minimum value[microseconds]; alignment [mm/beam]; noise [mm]; wet[mm];dry[mm];iono[mm]; solid_earth[mm]; geocentric_tide[mm]; ocean_loading[mm];geophysical_correction [mm]; zero padding; n_beams; dist;\n');
%fprintf(fidResults,'%d; %s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; %s; %d; %d; %d;  \n',...
switch mission
    case 'S3_'
        fprintf(fidResults,'%d; %s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; %s; %d; %d; %d;\n',...
            cycle,currenfolder,...
            (range_error_A1_v2*1e3),...
            (range_error_A1_v1*1e3),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos), 
            (wet), (dry), (iono),...
            (solid_earth),(geocentric_tide),(ocean_loading),(geophysical_correction),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),L1BS_buffer(i_surf_stacked).dist);%,...
            %(mode_L1),(mode_L2),(TRP_int_delay*1000));%,L1BS_buffer(i_surf_stacked).dist); %,instrument_range_correction_tx_rx,USO_correction)
        %lla2kml(['./'  currenfolder '_L1BS'], [L1BS_buffer.lat_surf], [L1BS_buffer.lon_surf]);
        %lla2kml_tour(['./'  currenfolder '_tour'],L1BS_buffer);
        %lla2kml(['./' currenfolder '_wet_min'], L1A_buffer(wet_min_pos).lat_sar_sat, L1A_buffer(wet_min_pos).lon_sar_sat, L1A_buffer(wet_min_pos).alt_sar_surf,'.','yellow');
        %lla2kml(['./' currenfolder '_dry_min'], L1A_buffer(dry_min_pos).lat_sar_sat, L1A_buffer(dry_min_pos).lon_sar_sat, L1A_buffer(dry_min_pos).alt_sar_surf,'.','green');
        %lla2kml(['./' currenfolder '_pos'], L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(pos)).lat_sar_sat, L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(pos)).lon_sar_sat, L1A_buffer(L1BS_buffer(i_surf_stacked).burst_index(pos)).alt_sar_surf,'.','red');
    case 'CR2'
        fprintf(fidResults,'%s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; \n',...
            currenfolder,...
            (range_error_A1_v2),...
            (range_error_A1_v1),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos), (wet), (wet_min_pos), (dry), (dry_min_pos), (iono),...
            (pos), (wet), (dry), (iono),...
            (solid_earth),(geocentric_tide),(ocean_loading),(geophysical_correction),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),...
            (TRP_int_delay*1000),L1BS_buffer(i_surf_stacked).dist);%,sat_view); %,instrument_range_correction_tx_rx,USO_correction)
        
    case 'LRM'
        fprintf(fidResults,'%s; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %d; %s; \n',...
            currenfolder,...
            (range_error_A1_v2),...
            (range_error_A1_v1),...
            (mean(range_stack(valid_beams)-theo_range_aligned(valid_beams))*1e3),...
            (datation_error_A1_v2),...
            (datation_error_A1_v1),...
            (Stack_alignment_A1*1000),...
            (Stack_noise_A1*1000),...%(pos), (wet), (wet_min_pos), (dry), (dry_min_pos), (iono),...
            (pos), (wet), (dry), (iono),...
            (solid_earth),(geocentric_tide),(ocean_loading),(geophysical_correction),(alt_rate_sat),...
            (zp_fact_trp),(valid_beams(1)),(valid_beams(end)),(length(valid_beams)),...
            (TRP_int_delay*1000),L1BS_buffer(i_surf_stacked).dist);

end
fclose(fidResults);


% h=figure;
% subplot(1,2,1);
% plot(range_without_slant'); hold all; plot(theo_range');
% figlabels('Beams','Range [m]','','Range' ,14);
% legend('Range meas', 'Range theo');
% subplot(1,2,2);plot(range_without_slant-theo_range);
% figlabels('Beams','Range [m]','',['Error = Range meas - Range theo, zp:' num2str(zp_fact_trp)] ,14);
%                                                %L1BS_buffer(i_surf_stacked).win_delay_surf-TRP_int_delay/cst.c*2;
% elevation = L1BS_buffer(i_surf_stacked).alt_sat-L1BS_buffer(i_surf_stacked).win_delay_surf*cst.c/2-((max_pos-(chd.N_samples_sar/2))/chd.bw*cst.c/2) + TRP_int_delay/2- geophysical_correction;
%    
% % disp(elevation(floor(length(elevation)/2))-chd.alt_trp);
% saveas (h,[filesBulk.resultPath 'Stack_range_bias_' num2str(L1BS_buffer(i_surf_stacked).surf_counter) '.png']);
% close(h);
% 
% 
