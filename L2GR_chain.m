
% v1.0 First version
% v2.0 added multiple retrackers and plotting at the end
% 

function filesBulk = L2GR_chain(filesBulk, cnf, chd, cst)



font_size = 12;

set_default_plot;
    
%% READ CONFIG
aux=struct2cell(filesBulk.auxFiles); aux=aux(1,:);
filesBulk.CNF_file_L2=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'cnf_file_L2'))).name];
run(filesBulk.CNF_file_L2); % cnf struct output


plotting_retrackers = 1;

[~,outputPath_size] = size(filesBulk.outputPath);

% filesBulk.filename_L1B = 'G:/My Drive/S3NGT/codes/S3NGT_TOO4/results/S3NGT_SIRS_HR_1B_0280101T234415_20280101T23442__isd_long.nc';
plots_folder = [filesBulk.outputPath 'plots/'];

%reading L1B parameters
cnf.RMC_RMC_on=0;
if(cnf.RMC_RMC_on)
    [L1B]= reading_L1B_SIN(filesBulk, chd, cnf, cst);
else
    [L1B]= reading_L1B (filesBulk, chd, cnf, cst);
%     [L1B_FF]= reading_L1B_FF (filesBulk, chd, cnf, cst);
end
L1B_written = readanyNETCDF_V1(filesBulk.filename_L1B);

%JPLZ: Temporary load by hand the FF-ML variables
%TSI2.ku_ML_1
%load('C:\Users\juanpedro\isardSAT\Projects\CRISTAL\Test_cases_IRIS_SVALP\script_callers\TSI2_SARin_OB_Ku\results\FF_048mresSL_10mresML_300tint\PICE_1B_FF-ML_variables_KU_20280101T193201_20280101T19320_isd.mat');
%TSI.2ka_ML
%load('C:\Users\juanpedro\isardSAT\Projects\CRISTAL\Test_cases_IRIS_SVALP\script_callers\TSI2_SAR_OB_Ka\results\FF_048mresSL_10mresML_250tint\PICE_1B_FF-ML_variables_KA_20280101T193201_20280101T19320_isd.mat');
%TSI2v2.ku_ML_1
%load('C:\Users\juanpedro\isardSAT\Projects\CRISTAL\Test_cases_IRIS_SVALP\script_callers\TSI2v2_SARin_OB_Ku\results\FF_045mresSL_3mresML_200tint\PICE_1B_FF-ML_variables_KU_20280101T193201_20280101T19320_isd.mat');
%TSIv2.2ka_ML
%load('C:\Users\juanpedro\isardSAT\Projects\CRISTAL\Test_cases_IRIS_SVALP\script_callers\TSI2v2_SAR_OB_Ka\results\FF_045mresSL_3mresML_100tint\PICE_1B_FF-ML_variables_KA_20280101T193201_20280101T19320_isd.mat');
%TSI8.ku_ML_1
%load('C:\Users\juanpedro\isardSAT\Projects\CRISTAL\Test_cases_IRIS_SVALP\script_callers\TSI8_SARin_OB_Ku\results\FF_048mresSL_10mresML_300tint\PICE_1B_FF-ML_variables_KU_20280101T193201_20280101T19320_isd.mat');
%TSI8.ka_ML_1
%load('C:\Users\juanpedro\isardSAT\Projects\CRISTAL\Test_cases_IRIS_SVALP\script_callers\TSI8_SAR_OB_Ka\results\FF_048mresSL_10mresML_250tint\PICE_1B_FF-ML_variables_KA_20280101T193201_20280101T19320_isd.mat');
%L1B.scaled_waveforms = (((L1B.waveform_scale_factor) * ones(1,chd.N_samples_sar*cnf.zp_fact_range)).' .* double(L1B.i2q2_meas_ku_l1b_echo).').';

% L1B.scaled_waveforms_ML=pow_ML_FF;

%waveform_scale_factor_ML = single(max(pow_ML_FF.')/ (2^16-2));
%i2q2_meas_ku_l1b_echo_ML = double(round(pow_ML_FF.' ./ waveform_scale_factor_ML));

%scaled_waveforms_ML = (((waveform_scale_factor_ML.') * ones(1,chd.N_samples_sar*cnf.zp_fact_range)) .* double(i2q2_meas_ku_l1b_echo_ML).');

% if(strcmp(cnf.processing_mode,'SIN'))
% L1B.scaled_waveforms_ML_2=pow_ML_FF_2;
% end

% % load ML win delay for computing range
% L1B.alt_sat_surf_ML=alt_sat_surf_ML.';
% L1B.win_delay_surf_ML=win_delay_surf_ML.';
% %L1B.range_ku_l1b_echo_ML=double((L1B.win_delay_surf_ML*cst.c/2-700000)*1e4);
% L1B.alt_sat_surf_ML=alt_sat_surf_ML.';
% L1B.range_ku_l1b_echo_ML=double((L1B.win_delay_surf_ML*cst.c/2));
% 
% range_scale_factor = L1B_written.attributes.range_ku_l1b_echo.scale_factor;
% range_offset = L1B_written.attributes.range_ku_l1b_echo.add_offset;

%%for TSI2v2Ku, cut out of scenario zones
% L1B.scaled_waveforms=L1B.scaled_waveforms(25044:39999,:);
% L1B.alt=L1B.alt(25044:39999);
% L1B.range_ku_l1b_echo=L1B.range_ku_l1b_echo(25044:39999);
% 
% L1B.scaled_waveforms_ML=L1B.scaled_waveforms_ML(3130:5002,:);
% L1B.alt_sat_surf_ML=L1B.alt_sat_surf_ML(3130:5002);
% L1B.range_ku_l1b_echo_ML=L1B.range_ku_l1b_echo_ML(3130:5002);

%%

delta_h = 0;
                                                      
%call analytical retracker with the common strcuture from the retracker repository
set_default_plot;
[data, cst_p, chd_p]=read_adapt_cristal2retracker(L1B, cnf_L2, cst, chd);
[ocean_sar] =analytical_retracker_commonstruct(data, cnf_L2, chd_p, cst_p, 'LUT_f0_file',filesBulk.LUT_f0_file,...
                                                         'LUT_f1_file',filesBulk.LUT_f1_file,...
                                                         'path_Results',filesBulk.outputPath, ...
                                                         'L1B_filename',filesBulk.filename_L1B);   
                                                      
height_ssh        = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - ocean_sar.Epoch.')*chd.T0_nom * cst.c / 2))+ delta_h.');

[ocog] = OCOG_retracker(L1B.scaled_waveforms,L1B.alt,cnf_L2);
height_ocog        = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - ocog.COG.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');


figure; plot(L1B.lat(1:44),height_ssh(1:44),'bo')
hold on
plot(L1B.lat(56:58),height_ssh(56:58),'ro')
grid
xlabel('Latitude')
ylabel('Height [m]')
legend('Sea-ice','Ocean lead')
ylim([0 1])

figure; plot(L1B.lat,height_ocog,'bo')
grid
xlabel('Latitude')
ylabel('Elevation [m]')
% legend('')
ylim([240 250])



% Extract the coordinates of the first position
first_position = ecef(1, :);

% Calculate the Euclidean distances
distances = pdist2(first_position, ecef);
figure; plot(L1B.lat(10:120),height_ocog(10:120),'bo')
grid
xlabel('Distance')
ylabel('Elevation [m]')

range_random_error = std(height_ssh(1:29*3)) / sqrt(29)
% hold on
% m = mean(height_ssh(1:80));
% plot(1:142, ones(1,142)*m,'k--')
% title('TOO4 range random error')

    %% Correct for window delay
    % DDP-HR waveforms
    geoid_correction = geoidheight(L1B.lat,L1B.lon);
    elev_HR = (L1B.alt-L1B.range_ku_l1b_echo)-geoid_correction';
    % We will now take the first record as reference to build the elevation axis. This can be reviewed in the future if needed.
    top_elev_HR  = elev_HR(1)+size(L1B.scaled_waveforms,2)/2/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2;
    bottom_elev_HR    = elev_HR(1)-(size(L1B.scaled_waveforms,2)/2)/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2;
    elev_axis_HR = (top_elev_HR-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:bottom_elev_HR);
    
    clear wvf_raw_shifted_HR
    wvf_raw_shifted_HR=cat(1,0.*L1B.scaled_waveforms.',0.*L1B.scaled_waveforms.',0.*L1B.scaled_waveforms.');
    wvf_raw_shifted_HR=wvf_raw_shifted_HR.';
    if(strcmp(cnf.processing_mode,'SIN'))
        ph_diff_shifted_HR=cat(1,0.*L1B.phase_diff_meas_ku_l1b_echo.',0.*L1B.phase_diff_meas_ku_l1b_echo.',0.*L1B.phase_diff_meas_ku_l1b_echo.');
        ph_diff_shifted_HR=ph_diff_shifted_HR.';
    end

    clear shift_mat
    for i_surf=1:size(L1B.scaled_waveforms,1)
        %shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/c_cst*2*2*Bw); %This extra 2 is for zeropadding
        shift_mat(i_surf)=round((elev_HR(1)-elev_HR(i_surf))/cst.c*2*(1/chd.T0_nom)*cnf.zp_fact_range);
        %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
        wvf_raw_shifted_HR(i_surf, size(L1B.scaled_waveforms,2) + 1 +shift_mat(i_surf):(size(L1B.scaled_waveforms,2)*2+shift_mat(i_surf)))= L1B.scaled_waveforms(i_surf,:) ;
        if(strcmp(cnf.processing_mode,'SIN'))
            ph_diff_shifted_HR(i_surf, size(L1B.phase_diff_meas_ku_l1b_echo,2) + 1 +shift_mat(i_surf):(size(L1B.phase_diff_meas_ku_l1b_echo,2)*2+shift_mat(i_surf)))= L1B.phase_diff_meas_ku_l1b_echo(i_surf,:) ;
        end
    end
    L1B.scaled_waveforms_aligned=wvf_raw_shifted_HR(:,257:512);
    if(strcmp(cnf.processing_mode,'SIN'))
        L1B.phase_diff_aligned=ph_diff_shifted_HR(:,257:512);
    end
    
    % FF-SL waveforms 
    geoid_correction_SL = geoidheight(L1B_FF.lat,L1B_FF.lon);  
    elev_SL = (L1B_FF.alt-L1B_FF.range_ku_l1b_echo)-geoid_correction_SL';
    top_elev_SL  = elev_SL(1)+size(L1B_FF.scaled_waveforms,2)/2/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2;
    bottom_elev_SL    = elev_SL(1)-(size(L1B_FF.scaled_waveforms,2)/2)/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2;
    elev_axis_SL = (top_elev_SL-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:bottom_elev_SL);
    
    clear wvf_raw_shifted_SL
    wvf_raw_shifted_SL=cat(1,0.*L1B_FF.scaled_waveforms.',0.*L1B_FF.scaled_waveforms.',0.*L1B_FF.scaled_waveforms.');
    wvf_raw_shifted_SL=wvf_raw_shifted_SL.';
    clear shift_mat
    for i_surf=1:size(L1B_FF.scaled_waveforms,1)
        %shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/c_cst*2*2*Bw); %This extra 2 is for zeropadding
        shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/cst.c*2*(1/chd.T0_nom)*cnf.zp_fact_range);
        %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
        wvf_raw_shifted_SL(i_surf, size(L1B_FF.scaled_waveforms,2) + 1 +shift_mat(i_surf):(size(L1B_FF.scaled_waveforms,2)*2+shift_mat(i_surf)))= L1B_FF.scaled_waveforms(i_surf,:) ;
    end
    L1B_FF.scaled_waveforms_aligned_SL=wvf_raw_shifted_SL(:,257:512);

   
    if(strcmp(cnf.processing_mode,'SIN'))
     %2nd Ku wavefom   
    end
    
    % FF-ML waveforms
    geoid_correction_ML = geoidheight(L1B_FF.lat_surf_ML,L1B_FF.lon_surf_ML);
    elev_ML = (L1B_FF.alt_surf_ML-L1B_FF.range_ML)-geoid_correction_ML';
    top_elev_ML  = elev_ML(1)+size(L1B_FF.power_scaled_waveforms_ML,2)/2/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2;
    bottom_elev_ML    = elev_ML(1)-(size(L1B_FF.power_scaled_waveforms_ML,2)/2)/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2;
    elev_axis_ML = (top_elev_ML-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:bottom_elev_ML);
    
    clear wvf_raw_shifted_ML
    wvf_raw_shifted_ML=cat(1,0.*L1B_FF.power_scaled_waveforms_ML.',0.*L1B_FF.power_scaled_waveforms_ML.',0.*L1B_FF.power_scaled_waveforms_ML.');
    wvf_raw_shifted_ML=wvf_raw_shifted_ML.';
    if(strcmp(cnf.processing_mode,'SIN'))
    ph_diff_shifted_ML=cat(1,0.*L1B_FF.phase_difference_ML.',0.*L1B_FF.phase_difference_ML.',0.*L1B_FF.phase_difference_ML.');
    ph_diff_shifted_ML=ph_diff_shifted_ML.';
    end

    clear shift_mat
    for i_surf=1:size(L1B_FF.power_scaled_waveforms_ML,1)
        %shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/c_cst*2*2*Bw); %This extra 2 is for zeropadding
        shift_mat(i_surf)=round((elev_ML(1)-elev_ML(i_surf))/cst.c*2*(1/chd.T0_nom)*cnf.zp_fact_range);
        %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
        wvf_raw_shifted_ML(i_surf, size(L1B_FF.power_scaled_waveforms_ML,2) + 1 +shift_mat(i_surf):(size(L1B_FF.power_scaled_waveforms_ML,2)*2+shift_mat(i_surf)))= L1B_FF.power_scaled_waveforms_ML(i_surf,:) ;
        if(strcmp(cnf.processing_mode,'SIN'))
            ph_diff_shifted_ML(i_surf, size(L1B_FF.phase_difference_ML,2) + 1 +shift_mat(i_surf):(size(L1B_FF.phase_difference_ML,2)*2+shift_mat(i_surf)))= L1B_FF.phase_difference_ML(i_surf,:) ;
        end
    end
    L1B_FF.power_scaled_waveforms_ML_aligned=wvf_raw_shifted_ML(:,257:512);
    if(strcmp(cnf.processing_mode,'SIN'))
        L1B_FF.phase_difference_ML_aligned=ph_diff_shifted_ML(:,257:512);
    end
    
    % ***TO DO***JPLZ: implement method to cut the aligned waveforms to the
    % N_sar_samples length optimally, for now we do it manually
    % TSI2ku (:,253:508)
    %L1B.scaled_waveforms_aligned_ML=wvf_raw_shifted_ML(:,253:508);
    
    L1B_FF.power_scaled_waveforms_ML_aligned=wvf_raw_shifted_ML(:,257:512);

    if(strcmp(cnf.processing_mode,'SIN'))
        %2nd Ku waveform
    end
    
    %% method updated FF-ML
%     geoid_correction_ML = geoidheight(L1B_FF.lat_surf_ML,L1B_FF.lon_surf_ML);
%     elev_ML = (L1B_FF.alt_surf_ML-L1B_FF.range_ML)-geoid_correction_ML';
%     max_elev_ML  = max(elev_ML);
%     min_elev_ML    = min(elev_ML);
%     %     elev_axis_ML_2 = (top_elev_ML-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:-1/(1/chd.T0_nom)/cnf.zp_fact_range*cst.c/2:bottom_elev_ML);
%     [~,reference_elevation_surf_ML] = max(elev_ML);
%     %select the elevation to use as a reference, we will shift all the others and keep this one fixed
%     N_windows_ML = floor((max_elev_ML-min_elev_ML)/(1/chd.bw/cnf.zp_fact_range*cst.c/2)/size(L1B_FF.power_scaled_waveforms_ML,2))  +1;
%     
%     clear wvf_raw_shifted_ML 
%     
%     wvf_raw_shifted_ML = zeros(size(L1B_FF.power_scaled_waveforms_ML,1),(N_windows+1)*size(L1B_FF.power_scaled_waveforms_ML,2));
%      
%     top_elev_ML = elev_ML(reference_elevation_surf_ML)+size(L1B_FF.power_scaled_waveforms_ML,2)/2/chd.bw/cnf.zp_fact_range*cst.c/2;
%     bottom_elev_ML = elev_ML(reference_elevation_surf_ML)-(size(L1B_FF.power_scaled_waveforms_ML,2))/chd.bw/cnf.zp_fact_range*cst.c/2/2-N_windows*size(L1B_FF.power_scaled_waveforms_ML,2)/chd.bw/cnf.zp_fact_range*cst.c/2;
%     elev_axis_ML = (top_elev_ML:-1/chd.bw/cnf.zp_fact_range*cst.c/2:bottom_elev_ML+1/chd.bw/cnf.zp_fact_range*cst.c/2);
%     
%     clear shift_mat
% 
% for i_surf=1:length(L1B_FF.range_ML)
%     shift_mat(i_surf)=round((elev_ML(reference_elevation_surf_ML)-elev_ML(i_surf))/cst.c*cnf.zp_fact_range*1*chd.bw); %This extra 2 is for zeropadding JPL: I put a 1 since we use zp_fact_rng=1 now
%     %shift_mat(i_surf)=round((elev_DDP(1)-elev_DDP(i_surf))/c_cst*2*Bw);
%     %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
%     wvf_raw_shifted_ML(i_surf, 1+shift_mat(i_surf):size(L1B_FF.power_scaled_waveforms_ML,2)+shift_mat(i_surf))= L1B_FF.power_scaled_waveforms_ML(i_surf,:) ;
% %     cs2_ph_diff_shifted_ML(1+shift_mat(i_surf):size(cs2_ph_diff,1)+shift_mat(i_surf),i_surf)= cs2_ph_diff(:,i_surf) ;
% %     cs2_coh_shifted_ML(1+shift_mat(i_surf):size(cs2_coh,1)+shift_mat(i_surf),i_surf)= cs2_coh(:,i_surf) ;
% end
%     
%    %Plot waveforms for all samples, re-aligned with the window shift
% figure; 
%     %imagesc(cs2_lat,elev_axis_ML,20*log10(wvf_raw_shifted_ML));
%     imagesc(L1B_FF.lat_surf_ML,elev_axis_ML,20*log10(wvf_raw_shifted_ML/max(wvf_raw_shifted_ML(:))).'); %normalised
%     %caxis([0 100]) %changes colour scale limits (originally set to 65 100)
%     %figlabels('Latitude','Elevation','','L1B waveforms(dB) ',23);set(gca,'Fontsize',23);
%     figlabels('Latitude','Elevation [m]','',['L1B waveforms(dB): ' ],23);set(gca,'Fontsize',23);
%     set(gca,'YDir', 'normal'); %corrects Y axis if upside down
%     ylim([-50 20]) 
% %     xlim([58.89 59.12])
%     caxis([-100 0]) %if normalised
%     colormap(colormap_blues);
%     hcb=colorbar;
%     hcb.Label.String = 'Rx Power, Normalised';
    
    %% Blue radargram plot
    %HR
    figure;
    imagesc(L1B.lat, elev_axis_HR, 10.*log10(L1B.scaled_waveforms_aligned/max(L1B.scaled_waveforms_aligned(:))).');
    title('HR')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Rx Power Normalised [dB]';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')
    load('colormap_blues.mat');
    colormap(colormap_blues);
    
    %FF-SL
    figure;
    imagesc(L1B_FF.lat, elev_axis_SL, 10.*log10(abs(L1B_FF.scaled_waveforms_aligned_SL)/max(abs(L1B_FF.scaled_waveforms_aligned_SL(:)))).');
    title('FF-SL')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Rx Power Normalised [dB]';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')
    load('colormap_blues.mat');
    colormap(colormap_blues);
    
    %FF-ML
    figure;
    imagesc(L1B_FF.lat_surf_ML, elev_axis_ML, 10.*log10(L1B_FF.power_scaled_waveforms_ML_aligned/max(L1B_FF.power_scaled_waveforms_ML_aligned(:))).');
    title('FF-ML')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Rx Power Normalised [dB]';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')
    load('colormap_blues.mat');
    colormap(colormap_blues);
    
    for i=1:size(L1B.phase_diff_meas_ku_l1b_echo,1)
        smoothed_phase(i,:)=smooth(L1B.phase_diff_meas_ku_l1b_echo(i,:));
        smoothed_phase_aligned(i,:)=smooth(L1B.phase_diff_aligned(i,:));
        unwraped_phase(i,:)=unwrap(L1B.phase_diff_meas_ku_l1b_echo(i,:));
        unwraped_phase_aligned(i,:)=unwrap(L1B.phase_diff_aligned(i,:),1.5*pi);     
    end
   
    %Phase diff radargram plot -HR
    figure;
    %imagesc(L1B.lat, elev_axis_HR, L1B.phase_diff_meas_ku_l1b_echo.');
    imagesc(L1B.lat, elev_axis_HR,  L1B.phase_diff_aligned.');
    %imagesc(L1B.lat, elev_axis_HR,  smoothed_phase_aligned.');
    %imagesc(L1B.lat, elev_axis_HR,  unwraped_phase_aligned.');
    title('HR')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Phase difference [rad]';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')
%     load('colormap_blues.mat');
%     colormap(colormap_blues);

    %Coherence radargram plot -HR
    figure;
    imagesc(L1B.lat, elev_axis_HR+10, L1B.phase_diff_aligned.');
    title('HR')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Phase difference';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')
%     load('colormap_blues.mat');
%     colormap(colormap_blues);
    
%Phase diff radargram plot -FF-ML
 figure;
    %imagesc(L1B.lat, elev_axis_HR, L1B.phase_diff_meas_ku_l1b_echo.');
    %imagesc(L1B.lat, elev_axis_HR,  L1B.phase_diff_aligned.');
    %imagesc(L1B.lat, elev_axis_HR,  smoothed_phase_aligned.');
    imagesc(L1B_FF.lat_surf_ML, elev_axis_ML,  L1B_FF.phase_difference_ML_aligned.');
    title('FF-ML')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Phase difference [rad]';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')

    
%% Iceberg detection method
if cnf_L2.iceberg_detection
    disp('------------Entering iceberg detection algorithm------------')
    [icebergs_properties_L1B]=iceberg_detection(filesBulk,L1B,cnf,cnf_L2,chd,cst);
else
    disp('------------Skipping iceberg detection algorithm------------')
end
%% SWATH
if(cnf.run_swath==1)
    
%     [SWATH]= swath_processing (filesBulk,L1B,cnf,chd,cst);    
    [SWATH,filesBulk]= swath_processing (filesBulk,L1B,cnf,chd,cst);
    
    %write swath product TBD
end
%% ----------------------- RETRACKERS -----------------------------

cnf.range_index1 = 50;% ranges of noise floor in the aligned ML waveforms
cnf.range_index2 = 54;
cnf.range_index1 = 10;
cnf.range_index2 = 20;

PP = peakiness(L1B.scaled_waveforms_aligned.',cnf, chd);
PP_ML = peakiness(L1B.scaled_waveforms_aligned_ML.',cnf, chd);
delta_h = compute_snow_depth_correction (PP,chd);
delta_h_ML = compute_snow_depth_correction (PP_ML,chd);

%ML delta_h filter and plot
windowSize = 20;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
PP_ML_filtered=filter(b,a,PP_ML);
delta_h_ML_filtered=filter(b,a,delta_h_ML);

x1=1:length(delta_h_ML);

figure;
plot(PP);
hold;
plot(x1,PP_ML,x1,PP_ML_filtered);
hold off;
figure;
plot(delta_h);
hold;
plot(x1,delta_h_ML,x1,delta_h_ML_filtered);
hold off;

delta_h =0;
delta_h_ML =0;

%% Classification Sea Ice / Open Ocean / Leads Done using the peakiness although it should be done using aux data.
% Leads, the most specular ones
% Sea Ice, medium peaky, depending on the snow depth and the size
% Open Ocean, least peaky.


% index_sea_ice = 82:240;
% index_ocean = [1:80 242:length(PP)];
index_sea_ice = 34:90;
index_ocean = [1:31 93:length(PP)];


%% 0 PEAK 
[max_value,max_pos] = peak_retracker(L1B.scaled_waveforms.');
max_pos=max_pos/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
height_peak     = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - max_pos.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');

clear max_value max_pos
% 0 PEAK ML
[max_value,max_pos] = peak_retracker(L1B_FF.power_scaled_waveforms_ML.');
max_pos=max_pos/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
height_peak_ML     = ((L1B_FF.alt_surf_ML - (L1B_FF.range_ML-(chd.N_samples_sar/2 - max_pos.'+1)*chd.T0_nom * cst.c / 2))+ delta_h_ML.');

%% 0 GAUSSIAN
if(chd.simulation  == 'PT')
    % DDP
    parab     = fitgaus(L1B.scaled_waveforms);
    parab_max = parab(:,1)/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
    maxpos_gaus = parab_max';
    height_gaussian = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - maxpos_gaus.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');
    
    clear parab parab_max maxpos_gaus
    % ML
    parab     = fitgaus(L1B_FF.power_scaled_waveforms_ML);
    parab_max = parab(:,1)/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
    maxpos_gaus = parab_max';
    height_gaussian_ML = ((L1B_FF.alt_surf_ML - (L1B_FF.range_ML-(chd.N_samples_sar/2 - maxpos_gaus.'+1)*chd.T0_nom * cst.c / 2))+ delta_h_ML.');
end
%% 1 TFMRA 
thresholds  = [0.95 0.87 0.6 0.4];
for i_th_index=1:length(thresholds)
    i_th=thresholds(i_th_index);
    [tfmr_val(i_th_index,:),tfmr_pos(i_th_index,:)] = threshold_first_maximum_restracker(L1B.scaled_waveforms.',i_th,cnf);
    tfmr_pos(i_th_index,:)=tfmr_pos(i_th_index,:)/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
    height_tfmr(i_th_index,:)        = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - tfmr_pos(i_th_index,:).'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');
end

clear tfmr_val tfmr_pos
% 1 TFMRA ML
for i_th_index=1:length(thresholds)
    i_th=thresholds(i_th_index);
    [tfmr_val(i_th_index,:),tfmr_pos(i_th_index,:)] = threshold_first_maximum_restracker(L1B_FF.power_scaled_waveforms_ML.',i_th,cnf);
    tfmr_pos(i_th_index,:)=tfmr_pos(i_th_index,:)/cnf.zp_fact_range  -1/cnf.zp_fact_range  +1;
    height_tfmr_ML(i_th_index,:)        = ((L1B_FF.alt_surf_ML - (L1B_FF.range_ML-(chd.N_samples_sar/2 - tfmr_pos(i_th_index,:).'+1)*chd.T0_nom * cst.c / 2))+ delta_h_ML.');
    %height_tfmr_ML(i_th_index,:)        = ((L1B.alt_sat_surf_ML - (L1B.range_ku_l1b_echo_ML-(chd.N_samples_sar/2 - tfmr_pos(i_th_index,:).'+1)*chd.T0_nom * cst.c / 2)));
end

% %ML filter and plot
% windowSize = 50;
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% height_tfmr_ML_filtered=filter(b,a,height_tfmr_ML(1,:));
% x1=1:length(height_tfmr_ML);

figure;
scatter(1:length(height_tfmr_ML(1,:)), height_tfmr_ML(1,:));


%% snow and Ku-Ka band comparison plots if both loaded


%plot snow depth if both files loaded
snow_depth_095=2.4/3*(-height_tfmr_ML_ku_PPfilter(1,:)+height_tfmr_ML_ka_PPfilter(1,:));
coincident_indxs=find(~isnan(snow_depth_095));

figure;
scatter(1:length(height_tfmr_ML_ku_PPfilter(1,:)), height_tfmr_ML_ku_PPfilter(1,:), ...
        'MarkerFaceColor', [0 0.45 0.74], 'MarkerEdgeColor', 'black');
xlabel('L2 records') 
ylabel('Height [m]')
legend('Ku band freeboard')
xlim([0 length(height_tfmr_ML_ku_PPfilter(1,:))]) 
ylim([8.5 11])
grid on
box on

figure;
scatter(1:length(height_tfmr_ML_ka_PPfilter(1,:)), height_tfmr_ML_ka_PPfilter(1,:), ...
        'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor', 'black');
xlabel('L2 records') 
ylabel('Height [m]')
legend('Ka band freeboard')
xlim([0 length(height_tfmr_ML_ka_PPfilter(1,:))]) 
ylim([9 11])
grid on
box on

figure;
scatter(1:length(snow_depth_095),snow_depth_095, ...
        'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', 'black');
xlabel('L2 records') 
ylabel('Height [m]')
legend('Snow depth')
xlim([0 length(snow_depth_095)]) 
ylim([-0.5 1.5])
grid on
box on



%ML delta_h filter and plot
windowSize = 20;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
height_tfmr_ML_ku_filtered=filter(b,a,height_tfmr_ML_ku(2,:));
height_tfmr_ML_ka_filtered=filter(b,a,height_tfmr_ML_ka(2,:));
snow_depth_087_filtered=2.4/3.*(height_tfmr_ML_ku_filtered-height_tfmr_ML_ka_filtered);
x2=1:length(height_tfmr_ML_ku);

figure;
plot(x1,delta_h_ML,x1,delta_h_ML_filtered);

%PP Ku-Ka
figure;
plot(x1,PP_ML_ku,x1,PP_ML_ka);

%delta_h Ku-Ka
figure;
plot(x1,delta_h_ML_ku,x1,delta_h_ML_ka);

% Heights-snow depth Ku-Ka
x1=1:length(height_tfmr_ML_ku(2,:));
figure;
plot(x1,height_tfmr_ML_ku_nodeltah(1,:),x1,height_tfmr_ML_ka_nodeltah(1,:),x1,snow_depth_087);
figure; %filters
plot(x2,height_tfmr_ML_ku_filtered,x2,height_tfmr_ML_ka_filtered,x2,snow_depth_087_filtered);


%% 2 OCOG 
[ocog] = OCOG_retracker(L1B.scaled_waveforms,L1B.alt,cnf_L2);
height_ocog        = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - ocog.COG.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');

% [ocog] = OCOG_retracker(L1B,cnf_L2);

% OCOG ML

[ocog_ML] = OCOG_retracker(L1B_FF.power_scaled_waveforms_ML,L1B_FF.alt_surf_ML,cnf_L2);
height_ocog_ML        = ((L1B_FF.alt_surf_ML - (L1B_FF.range_ML-(chd.N_samples_sar/2 - ocog_ML.COG.'+1)*chd.T0_nom * cst.c / 2))+ delta_h_ML.');

%[ocog_ML] = OCOG_retracker(L1B,cnf_L2);


%% 3 TCOG

[tcog] = TCOG_retracker(L1B.scaled_waveforms,L1B.alt,cnf_L2, chd);
height_tcog        = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - tcog.Epoch.'+1)*chd.T0_nom * cst.c / 2))+ delta_h.');

[tcog_ML] = TCOG_retracker(L1B_FF.power_scaled_waveforms_ML,L1B_FF.alt_surf_ML,cnf_L2, chd);
height_tcog_ML        = ((L1B_FF.alt_surf_ML - (L1B_FF.range_ML-(chd.N_samples_sar/2 - tcog_ML.Epoch.'+1)*chd.T0_nom * cst.c / 2))+ delta_h_ML.');

%% 4 LED * particular case for TFMRA when 0.5


%% 5 Physical OCEAN SAR

[ocean_sar] =analytical_retracker(L1B, cst, chd, cnf_L2, 'LUT_f0_file',filesBulk.LUT_f0_file,...
                                                          'LUT_f1_file',filesBulk.LUT_f1_file,...
                                                          'path_Results',filesBulk.outputPath);
                                                      
%call analytical retracker with the common strcuture from the retracker repository
set_default_plot;
[data, cst_p, chd_p]=read_adapt_cristal2retracker(L1B, cnf_L2, cst, chd);
[ocean_sar] =analytical_retracker_commonstruct(data, cnf_L2, chd_p, cst_p, 'LUT_f0_file',filesBulk.LUT_f0_file,...
                                                         'LUT_f1_file',filesBulk.LUT_f1_file,...
                                                         'path_Results',filesBulk.outputPath, ...
                                                         'L1B_filename',filesBulk.filename_L1B);
    
                                                      
height_ssh        = ((L1B.alt - (L1B.range_ku_l1b_echo-(chd.N_samples_sar/2 - ocean_sar.Epoch.')*chd.T0_nom * cst.c / 2))+ delta_h.');


% i_th_ocean= 0.87; % Ocean
% % 87% of peak, first value
% [max_value,max_pos]=max(L1B.scaled_waveforms.');
% 
% [~,records] = size(max_value);
% pc87_value=zeros(1,records);
% pc87_pos=zeros(1,records);
% for i=1:records
%     [pc87_value(i),pc87_pos(i)] = min( abs(L1B.scaled_waveforms(i,1:max_pos)-(i_th_ocean*max_value(i))) );
% end


%% 6 Lead Detection
if cnf.RMC_RMC_on    
    for b = 1:length(cnf.RMC.window_sizes);

    window_size  = cnf.RMC.window_sizes(b);

    peaks_all = [];
    lead_pks1 = [];
    lead_pks2 = [];


    cells_idx = find(~cellfun(@isempty,{filesBulk.SURF.N_surfs}));


    for index = cells_idx(1):cells_idx(end)
        pks = [];

        %% 1.EXTRACT L1B Parameters
        if isempty(filesBulk.SURF(index).surf_ini_sample)
            k = k + 1;
        else
             k= 0;
        end

        % waveform scaling
        scale_factor    = L1B.echo_scale_factor(index-k);
        scale_power     = L1B.echo_scale_power(index-k);

        surf_length = filesBulk.SURF(index-k).surf_end_sample - filesBulk.SURF(index-k).surf_ini_sample+1;
        [~,max_s_pos] = max(surf_length);
        start_idx = filesBulk.SURF(index-k).surf_ini_sample(max_s_pos);
        end_idx = filesBulk.SURF(index-k).surf_end_sample(max_s_pos);

        % scaling (we already have it done?)
        plot_Wvf        = 10*log10(L1B.averaged_power_echo_waveform(:,index-k).*(scale_factor).*2.^scale_power);
        % find peaks
        [~,~,~,pks(:,4)] = findpeaks(plot_Wvf);

        plot_Wvf        = L1B.averaged_power_echo_waveform(:,index-k).*(scale_factor).*2.^scale_power;

        [pks(:,1),pks(:,2),pks(:,3),~] = findpeaks(plot_Wvf);



        pks(:,5) = L1B.coherence(pks(:,2),index-k);
        pks(:,6) = length(pks(pks(:,1)>max(pks(:,1))/20));
        pks(:,7) = index*ones(length(pks(:,1)),1);

        pks = sortrows(pks);

        pks = flip(pks);


        pwr_thresh = max(pks(:,1))/cnf.RMC.pwr_thresh;

        pks_idx = find(((pks(:,1)>pwr_thresh) & (pks(:,3)>cnf.RMC.width_thresh) & (pks(:,4)>cnf.RMC.prominence_thresh) & (pks(:,5)>cnf.RMC.coherence_thresh) & pks(:,1)>=pks(1,1)/5) & ((pks(:,2)>end_idx) | (pks(:,2)<start_idx)));
        pks_idx2 = find(((pks(:,1)>pwr_thresh) & (pks(:,3)>cnf.RMC.width_thresh) & (pks(:,4)>cnf.RMC.prominence_thresh) & (pks(:,5)>cnf.RMC.coherence_thresh) & pks(:,1)>=pks(1,1)/5) & ((pks(:,2)>end_idx) | (pks(:,2)<start_idx) | (pks(:,2)>-10)) );


        peaks = pks(pks_idx,:);
        if(~isempty(peaks))
            [maxpkht,maxpkloc] = max(peaks(:,1));
            max_peaks(index) = peaks(maxpkloc,2);
        else
            max_peaks(index) = max_peaks(index-1);
        end

        if(~isempty(peaks))
            [~,mainpkloc] = min(peaks(peaks(:,1)>maxpkht/4,2));
            main_peaks(index) = peaks(mainpkloc,2);
        else
            main_peaks(index) = main_peaks(index-1);
        end

        pks_idx3 = find(((pks(:,1)>pwr_thresh) & (pks(:,3)>cnf.RMC.width_thresh) & (pks(:,4)>cnf.RMC.prominence_thresh) & (pks(:,5)>cnf.RMC.coherence_thresh) & pks(:,1)>=pks(1,1)/5) & ((pks(:,2)>end_idx) | (pks(:,2)<start_idx)) & (pks(:,2)>main_peaks(index)));% & (abs(SSHA_c)<=SSHA_thresh)));% | ((pks(:,3)>5) & (pks(:,1)>pwr_thresh/2))& pks(:,2)>=pks(1,2));


        peaks_all = [peaks_all ; peaks];


        if ~isempty(peaks) 
            lead1 = pks(pks_idx2(pks(pks_idx2,1)>cnf.RMC.sea_ice_power_thresh),:);
            lead2 = pks(pks_idx3,:);
        end

        lead_pks1 = [lead_pks1 ; lead1];
        lead_pks2 = [lead_pks2 ; lead2];

    end

        lead_pks = [lead_pks1 ; lead_pks2];


    if(RMC_on)
        peak_points = [];
        peak_leads = [];

        lower = ceil(cnf.RMC.tracking_gate/N_samples_sar_chd*(window_size*4))-1;
        upper = ceil((1 - cnf.RMC.tracking_gate/N_samples_sar_chd)*(window_size*4));


        for a = cells_idx(1):cells_idx(end)

            pk_pts = find((main_peaks(a)-lower)<peaks_all(peaks_all(:,7) == a,2)& peaks_all(peaks_all(:,7) == a,2)<(main_peaks(a)+upper));
            peak_points = [peak_points ; pk_pts];
            pk_leads = find((main_peaks(a)-lower)<lead_pks(lead_pks(:,7) == a,2)& lead_pks(lead_pks(:,7) == a,2)<(main_peaks(a)+upper));
            peak_leads = [peak_leads ; pk_leads];

            kept_peak_points = length(peak_points)/length(peaks_all)*100;
            removed_peak_points = length(peaks_all) - length(peak_points);
            kept_peak_leads = length(peak_leads)/length(lead_pks)*100;
            removed_peak_leads = length(lead_pks) - length(peak_leads);

        end

        cum_len = 0;
        cum_max_len = 0;
        start_idx = 0;  
        end_idx = 0;

        for a = cells_idx(1):cells_idx(end)
            if filesBulk.SURF(a).N_surfs_OK == 0
                continue
            elseif isempty(filesBulk.SURF(a).POCA)
                continue
            end
            surf_length = filesBulk.SURF(a).surf_end_sample - filesBulk.SURF(a).surf_ini_sample+1;
            [max_s_len,max_s_pos] = max(surf_length);
            start_idx(a) = filesBulk.SURF(a).surf_ini_sample(max_s_pos);
            end_idx(a) = filesBulk.SURF(a).surf_end_sample(max_s_pos);
            surf_idx = filesBulk.SURF(a).surf_ini_sample(max_s_pos):filesBulk.SURF(a).surf_end_sample(max_s_pos);
            cum_len = cum_len + length(surf_idx(surf_idx<main_peaks(a) + upper & surf_idx>main_peaks(a) - lower));
            cum_max_len = cum_max_len + max_s_len;
        end


        kept_swath_points = cum_len/cum_max_len * 100;
        removed_swath_points = cum_max_len - cum_len;


        RMC_out{b}(ia).Kept_Swath_Perc = kept_swath_points;
        RMC_out{b}(ia).Kept_Peak_Points_Perc = kept_peak_points;
        RMC_out{b}(ia).Kept_Peak_Leads_Perc = kept_peak_leads;
        RMC_out{b}(ia).No_of_Samples = length(SURF);
        RMC_out{b}(ia).Removed_Swath_Points = removed_swath_points;
        RMC_out{b}(ia).Removed_Peak_Points = removed_peak_points;
        RMC_out{b}(ia).Removed_Peak_Leads = removed_peak_leads;
    end
end
end
%% 7 SSHA
records_L1b=length(L1B.alt);

for index_L1b=1: records_L1b
    
    far_Range       = L1B.range_ku_l1b_echo(index_L1b)+chd.N_samples_sar/2*chd.T0_nom/cnf.zp_fact_range*cst.c/2;
    short_Range     = L1B.range_ku_l1b_echo(index_L1b)-(chd.N_samples_sar/2)*chd.T0_nom/cnf.zp_fact_range*cst.c/2;
    Range = (short_Range:chd.T0_nom/cnf.zp_fact_range*cst.c/2:far_Range-chd.T0_nom/cnf.zp_fact_range*cst.c/2);
    AoA    = (smooth(chd.wv_length * unwrap(L1B.phase_diff_meas_ku_l1b_echo(index_L1b,:))./(2*pi*B)-roll,'moving')/pi*180).';   %[deg]
    jump_AoA     = (chd.wv_length * 2*pi./(2*pi*B))/pi*180;   %[deg] 360deg to AoA
    AoA = AoA - jump_AoA* round(mean(AoA(270:274))/jump_AoA);
    Xdist  =  sin(AoA /180*pi) .* Range;

    if(index_L1b>1)
        [~,az_v] = distance('gc',L1B.lat(index_L1b-1),L1B.lon(index_L1b-1),L1B.lat(index_L1b),L1B.lon(index_L1b),[cst.semi_major_axis Eccentricity]);
    else
        [~,az_v] = distance('gc',L1B.lat(index_L1b),L1B.lon(index_L1b),L1B.lat(index_L1b+1),L1B.lon(index_L1b+1),[cst.semi_major_axis Eccentricity]);
    end

    az_d=az_v-90;

    [geoid_lat,geoid_lon]   = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b), -10000:100:10000 ,az_d,[cst.semi_major_axis Eccentricity]);
    geoid_simple = geoidheight(geoid_lat,geoid_lon);           
    geoid_x = -10000:100:10000;
    geoid_p = polyfit(geoid_x,geoid_simple,1);

    beta = geoid_p(1);
    alpha = beta/eta;
    window_retracking_shift = (512-272)*cst.c/2*chd.T0_nom/2;
    dh = (eta*(L1B.range_ku_l1b_echo(index_L1b)*cst.c/2-window_retracking_shift))/2.*(AoA.*AoA - 2*alpha*AoA);

    Range = Range - dh;

    SSH = cos(AoA .* Range/cst.semi_major_axis).*(L1B.alt(index_L1b) - Range.*cos(AoA) + cst.semi_major_axis.*(1 - cos(Range/cst.semi_major_axis .* AoA)));

    [lat_across(index_L1b),lon_across(index_L1b)] = reckon(L1B.lat(index_L1b),L1B.lon(index_L1b),Xdist,az_d,[cst.semi_major_axis Eccentricity]);

    geoid_across = geoidheight(lat_across(index_L1b),lon_across(index_L1b),'EGM2008'); % Must download EGM2008 using the aeroDataPackage function. Select Download from Internet and Aerospace Geoid Data. The file required is geoidegm2008grid.mat. Must add wherever it is downloaded to the Path. https://www.mathworks.com/help/aerotbx/ug/aerodatapackage.html for more info.

    SSHA = SSH - geoid_across';

    filesBulk.SSHA(index_L1b) = SSHA;
end         


if (plotting_retrackers)
    
    
    indexes = 40:length(height_peak);
%     indexes=59:86;
    simulation_tag = '';
    for i_surf=1:2 % ocean and ice results splitted
        if(strcmp(chd.simulation,'SI'))
            indexes= index_sea_ice;
            simulation_tag = 'seaice';
            if(i_surf==2)
                indexes=index_ocean;
                simulation_tag = 'ocean';
            end
        end
        figure; 
        plot(indexes,height_peak(indexes),   'DisplayName',['Peak > Mean; ' num2str(mean(height_peak(indexes))) '; m , Std: ' num2str(std(height_peak(indexes))) '; m, Range Noise 1Hz performance ' num2str(std(height_peak(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        hold all;
        disp(['Peak > Mean; ' num2str(mean(height_peak(indexes))) '; m , Std; ' num2str(std(height_peak(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_peak(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))])        
        plot((indexes),height_ssh(indexes),    'DisplayName',['Physical > Mean; ' num2str(mean(height_ssh(indexes))) '; m , Std: ' num2str(std(height_ssh(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_ssh(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        disp(['Physical > Mean; ' num2str(mean(height_ssh(indexes))) '; m , Std; ' num2str(std(height_ssh(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_ssh(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))])        
        plot((indexes),height_ocog(indexes),   'DisplayName',['OCOG > Mean; ' num2str(mean(height_ocog(indexes))) '; m , Std: ' num2str(std(height_ocog(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_ocog(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        disp(['OCOG > Mean; ' num2str(mean(height_ocog(indexes))) '; m , Std; ' num2str(std(height_ocog(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_ocog(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        plot((indexes),height_tcog(indexes),   'DisplayName',['TCOG > Mean; ' num2str(mean(height_tcog(indexes))) '; m , Std: ' num2str(std(height_tcog(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_tcog(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        disp(['TCOG > Mean; ' num2str(mean(height_tcog(indexes))) '; m , Std; ' num2str(std(height_tcog(indexes))) '; m, Range Noise 1Hz performance; ' num2str(std(height_tcog(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        for i_th_index=1:length(thresholds)
            plot((indexes),height_tfmr(i_th_index,(indexes)), 'DisplayName',['TFMR ' num2str(thresholds(i_th_index)) ' > Mean: ' num2str(mean(height_tfmr(i_th_index,(indexes)))) '; m, Std: ' num2str(std(height_tfmr(i_th_index,:))) '; m,' 'Range Noise 1Hz performance ' num2str(std(height_tfmr(i_th_index,(indexes)))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
            disp(['TFMR ' num2str(thresholds(i_th_index)) ' > Mean: ' num2str(mean(height_tfmr(i_th_index,(indexes)))) '; m, Std: ' num2str(std(height_tfmr(i_th_index,:))) '; m,' 'Range Noise 1Hz performance ' num2str(std(height_tfmr(i_th_index,(indexes)))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))]);
        end
        legend('Location','southoutside')
        if(chd.simulation  == 'PT')
            plot((indexes),height_gaussian(indexes),'DisplayName',['Gaussian > Mean: ' num2str(mean(height_gaussian(indexes))) '; m , Std: ' num2str(std(height_gaussian(indexes))) '; m']);

        end
        set(gca,'YLim',[-1+(0) 1+(0)],'FontSize',font_size);
        if(strcmp(chd.simulation,'SI'))
            set(gca,'YLim',[-3+(mean(height_peak)) 3+(mean(height_peak))],'FontSize',font_size);
        end
        set(gca,'XLim',[1 length(height_peak)],'FontSize',font_size);
        figlabels('L2 Records','Elevation [m]','',simulation_tag,font_size);
        figName = [plots_folder '7_Retrackers_' simulation_tag];
        print(gcf, figName,'-dpng');
        saveas (gcf,[figName,'.png']);
        
        
        if(strcmp(chd.simulation,'OC')||strcmp(simulation_tag,'ocean'))
            figure;
            subplot(2,2,1); plot((indexes),height_ssh(indexes),     'DisplayName',['SSH , Mean: ' num2str(mean(height_ssh(indexes))) '; m , Std: ' num2str(std(height_ssh(indexes)),3) '; m']);
            set(gca,'XLim',[1 length(height_peak)],'FontSize',font_size); figlabels('L2 Records','SSH [m]','',['SSH , Mean: ' num2str(mean(height_ssh(indexes)),3) '; m , Std: ' num2str(std(height_ssh(indexes)),3) '; m, Noise 1Hz performance ' num2str(std(height_ssh(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))],font_size);
            subplot(2,2,2); plot((indexes),ocean_sar.Hs(indexes),   'DisplayName',['SWH , Mean: ' num2str(mean(ocean_sar.Hs(indexes)),3) '; m , Std: ' num2str(std(ocean_sar.Hs(indexes)),3) '; m']);
            set(gca,'XLim',[1 length(height_peak)],'FontSize',font_size);figlabels('L2 Records','SWH [m]','',['SWH , Mean: ' num2str(mean(ocean_sar.Hs(indexes)),3) '; m , Std: ' num2str(std(ocean_sar.Hs(indexes)),3) '; m, Noise 1Hz performance ' num2str(std(ocean_sar.Hs(indexes))/sqrt(length(L1B.time)/(L1B.time(end)-L1B.time(1))))],font_size);
            subplot(2,2,3); plot((indexes),ocean_sar.Pu(indexes),   'DisplayName',['Pu , Mean: ' num2str(mean(ocean_sar.Pu(indexes)),3) ' dB , Std: ' num2str(std(ocean_sar.Pu(indexes)),3) ' dB']);
            set(gca,'XLim',[1 length(height_peak)],'FontSize',font_size); figlabels('L2 Records','Pu [dB]','',['Pu , Mean: ' num2str(mean(ocean_sar.Pu(indexes)),3) ' dB , Std: ' num2str(std(ocean_sar.Pu(indexes)),3) ' dB'],font_size);
            subplot(2,2,4); plot((indexes),ocean_sar.COR(indexes),  'DisplayName',['r , Mean: ' num2str(mean(ocean_sar.COR(indexes)),3) ' % , Std: ' num2str(std(ocean_sar.COR(indexes)),3) ' %']);
            set(gca,'XLim',[1 length(height_peak)],'FontSize',font_size); figlabels('L2 Records','Correlation coefficient [%]','',['r , Mean: ' num2str(mean(ocean_sar.COR(indexes)),3) ' % , Std: ' num2str(std(ocean_sar.COR(indexes)),3) ' %'],font_size);
           
            figName = [plots_folder 'Ocean_Retracker_results' ];
            print(gcf, figName,'-dpng');
            saveas (gcf,[figName,'.png']);    
        end
        if~(strcmp(chd.simulation,'SI'))
            break;
        end

    end
    
    
    close all
    %% -----------------------------
    figure; plot(max_value);
    figlabels('Peak Power','Power','','',font_size);
    figName = [plots_folder '8_Peak_power_L1B'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
   
    %% -----------------------------
    close all;
end