%% Matlab script to retreive results from transponder

%% Zero padding of the Tranponders Surfaces
global cnf.zp_fact_range chd.N_samples_sar_sar_chd cst.flat_coeff cst.semi_major_axis
global N_transponders_chd transponder_latitude_list_chd 
global transponder_longitude_list_chd transponder_altitude_list_chd
global cst.c chd.bw T0_nom


if cnf.force_exact_method==1
    file_specific_ext='_exact'
else
    file_specific_ext='_approximate'
end


SAVE = 1;
VISIBLE = 1;
set_default_plot
pt_index = transponder_surf_idx;

figure_format='jpg';
res_fig='-r150';

switch lower(figure_format)
    case 'pdf'
        file_ext='.pdf';
        print_file='-dpdf';
    case 'eps'
        file_ext='.eps';
        print_file='-depsc';
    case 'png'
        file_ext='.png';
        print_file='-dpng';
    case 'jpg'
        file_ext='.jpg';
        print_file='-djpeg';        
end

set(0,'defaultFigureVisible','off');
file_ext=strcat(file_specific_ext,file_ext);

switch(Process_ID)
    case 59
        mode_string=strcat('RAW',file_specific_ext);
    case 60
        mode_string=strcat('RMC',file_specific_ext);
end

input_filename=strrep(filename_ISP,'.DBL','');
result_Path ='C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\AR\results\L1B\TRP_performance\';

%N_transponders_chd=1;

for i_pt = 1:N_transponders_chd
    clear max_val max_pos beams_rng_cmpr_PT beams_rng_PT center_stack error_elevation range_theo range_meas max_val
    
    % Range compression of the stacks after geometry corrections
    [beams_rng_cmpr_PT]  = range_compression(1,N_beams_stack(pt_index(i_pt)),...
        wfm_geo_corr(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt)),:),w_rg);
    beams_rng_PT(:,:) = beams_rng_cmpr_PT;
    % Compute the positions of the across-track maximums within stack
    % [retracker]
    [max_val,max_pos] = max(beams_rng_PT.');    
    % Locate the central beam (as beam with maximum power)
    [~,center_stack] = max(max_val);
    center_stack_pos(i_pt) = center_stack;
    % Compute the elevation error within stack
    error_elevation = alt_sat(pt_index(i_pt))-(win_delay_surf(pt_index(i_pt))*cst.c/2+...
        ((max_pos)-(((chd.N_samples_sar_sar_chd/2+1)-1)*cnf.zp_fact_range+1))./cnf.zp_fact_range.*cst.c/2.*T0_sar_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))))-...
        transponder_altitude_list_chd(i_pt);
%     index_pt=max(1,center_stack_pos(i_pt)-125):min(center_stack_pos(i_pt)+124,N_beams_stack(pt_index(i_pt)));
%     index_pt=max(1,center_stack_pos(i_pt)-75):min(min(446,center_stack_pos(i_pt)+74),N_beams_stack(pt_index(i_pt)));
    % Select the central beams
    index_pt=max(1,center_stack_pos(i_pt)-125):min(min(N_beams_stack(pt_index(i_pt)),center_stack_pos(i_pt)+124),N_beams_stack(pt_index(i_pt)));
%       index_pt=1:N_beams_stack(pt_index(i_pt));
    % Compute the measured range for each beam (including the slant range variation)
    range_meas = win_delay_surf(pt_index(i_pt))*cst.c/2+...
        ((max_pos)-(((chd.N_samples_sar_sar_chd/2+1)-1)*cnf.zp_fact_range+1))./cnf.zp_fact_range.*cst.c/2.*T0_sar_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt)))-...
        slant_range_corr(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*cst.c/2.*T0_sar_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt)));
    xyz_pt = lla2ecef([transponder_latitude_list_chd(i_pt) transponder_longitude_list_chd(i_pt) transponder_altitude_list_chd(i_pt)], cst.flat_coeff, cst.semi_major_axis);
    burst_sel = (burst_index(pt_index(i_pt),find(burst_index(pt_index(i_pt),:))));
    xyz_sat = lla2ecef([lat_sar_sat(burst_sel) lon_sar_sat(burst_sel) alt_sar_sat(burst_sel)], cst.flat_coeff, cst.semi_major_axis);
    
    % Compute the theoretical range as distance between satellite & PTR 
    for i_beam=1:N_beams_stack(pt_index(i_pt))
        range_theo(i_beam) = sqrt((xyz_sat(i_beam,1)-xyz_pt(1))^2+(xyz_sat(i_beam,2)-xyz_pt(2))^2+(xyz_sat(i_beam,3)-xyz_pt(3))^2);
    end
    range_theo_coef = polyfit(burst_sel(index_pt)',range_theo(index_pt)',2);
    range_meas_coef = polyfit(burst_sel(index_pt)',range_meas(index_pt)',2);
    times_STACK = min((burst_sel)):1e-3:max((burst_sel));
    range_theo_STACK = range_theo_coef(3)+ range_theo_coef(2)*times_STACK+range_theo_coef(1)*(times_STACK).^2;
    range_meas_STACK = range_meas_coef(3) + range_meas_coef(2) *times_STACK+range_meas_coef(1) *(times_STACK).^2;
    [a_value,a_time] = min(range_theo_STACK);  
    [b_value,b_time] = min(range_meas_STACK);
    range_error_v1 = b_value-a_value; %[mm]
    datation_error_v1 = (times_STACK(b_time)-times_STACK(a_time))*mean(diff(time)); %[seconds]

    figure;
    plot((burst_sel-burst_sel(1)).*mean(diff(time)),range_meas/1e3,'-r');
    hold all;
    plot((burst_sel-burst_sel(1)).*mean(diff(time)),range_theo/1e3,'-b');
    plot((times_STACK-times_STACK(1)).*mean(diff(time)),range_meas_STACK/1e3,'-g');
    plot((times_STACK-times_STACK(1)).*mean(diff(time)),range_theo_STACK/1e3,'-m');
    ylab = get(gca,'YLabel');set(ylab,'String','Slant range [Km]','FontSize',26,'FontName','Arial');
    xlab = get(gca,'XLabel');set(xlab,'String','Time stack [s]','FontSize',26,'FontName','Arial');  
    set(gca,'XLim', [min((burst_sel-burst_sel(1)).*mean(diff(time)))...
        max((burst_sel-burst_sel(1)).*mean(diff(time)))],'FontSize',26) ;
    leg=legend ('Measured ','Expected','Polyfit measured','Polyfit expected','Location','NorthEast');
    title(strcat('Stack #',num2str(pt_index(i_pt)),' LAT. ',{' '},num2str(lat_surf(pt_index(i_pt))),'[deg] (',mode_string,')'),'Interpreter','None')
    pos_leg=get(leg,'Position');
    annotation('textbox',[pos_leg(1)-0.05,pos_leg(2)-pos_leg(4)-0.05,pos_leg(3),pos_leg(4)],...    
        'String',{strcat('Range bias [cm]=',num2str(range_error_v1*100)),...
        strcat('Datation error [\mus]=',num2str(datation_error_v1*1e6))},...
        'FitBoxToText','on','FontSize',26);
    print(print_file,res_fig,strcat(result_Path,input_filename,'_Slant_range_variation_measured_expected_',mode_string,'_PTR_',num2str(i_pt),file_ext)); 
    
%     plot_beam(:,:) = beams_surf(pt_index(i_pt),:,:);
%     figure; imagesc(abs(fftshift(fft(plot_beam.'),1)).');
%     plot_beam(:,:) = wfm_geo_corr(pt_index(i_pt),:,:);
%     figure; imagesc(abs(fftshift(fft(plot_beam.'),1)).');
%     figure; h=surf(1:N_beams_stack(pt_index(i_pt)),1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,beams_rng_PT.');
%     ylab = get(gca,'YLabel');set(ylab,'String','Samples','FontSize',40,'FontName','Arial');xlab = get(gca,'XLabel');set(xlab,'String','Beams','FontSize',40,'FontName','Arial');
%     set(h, 'edgecolor','none');colormap(colormapATDD);view(45,30)
    %Gaussian fit
    f = fit(index_pt',max_val(index_pt)'./max(max_val).','gauss1');
    gauss =  f.a1*exp(-((index_pt-f.b1)/f.c1).^2);
    skew(i_pt) = skewness(gauss);
    kurt(i_pt) = kurtosis(gauss);
    width(i_pt) = f.c1*sqrt(2*log(2));
    center_stack_pos_gauss(i_pt)=f.b1;
    
    %Gaussian fitting & statistics over the averaged-stack across-track +
    %(stats as per S-3 DPM: kurtosis, skewness and mean value are statistical parameters)
    %computed using the look angle as reference
    beam_power=sum(beams_rng_PT,2).';
    beam_power=beam_power/max(beam_power);
    f = fit(look_ang_surf(pt_index(i_pt),index_pt).',beam_power(index_pt).','gauss1');
    gauss =  f.a1*exp(-((look_ang_surf(pt_index(i_pt),index_pt)-f.b1)/f.c1).^2);
    [mu_stack(i_pt),std_stack(i_pt),skew(i_pt),kurt(i_pt)]=...
        stack_statistical_param(beam_power(index_pt),look_ang_surf(pt_index(i_pt),index_pt));
    
    h=figure; plot(look_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*180/pi,10*log10(beam_power),'Linewidth',3);
    hold all; plot(look_ang_surf(pt_index(i_pt),index_pt).*180/pi,10*log10(beam_power(index_pt)),'Linewidth',3);
    plot(look_ang_surf(pt_index(i_pt),index_pt).*180/pi,10*log10(gauss),'.-','Linewidth',2);  
    ylab = get(gca,'YLabel');set(ylab,'String','Norm. power [dB]','FontSize',26,'FontName','Arial');
    xlab = get(gca,'XLabel');set(xlab,'String','Look angle [deg.]','FontSize',26,'FontName','Arial');  
    set(gca,'XLim', [min(look_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*180/pi) max(look_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*180/pi)],'FontSize',26) 
    legend ('Norm. power ',[num2str(length(index_pt)) ' central Beams'],....
         ['Gauss fitting: Skewness ' num2str(skew(i_pt),3) ' Kurtosis ' num2str(kurt(i_pt),3)],'Location','SouthEast');
    title(strcat('Stack #',num2str(pt_index(i_pt)),' LAT. ',{' '},num2str(lat_surf(pt_index(i_pt))),'[deg] (',mode_string,')'),'Interpreter','None')
    print(print_file,res_fig,strcat(result_Path,input_filename,'Stack_fitting_',mode_string,'_PTR_',num2str(i_pt),file_ext));
%     h=figure; plot(max_val/max((max_val)),'Linewidth',3);
%     hold all; plot(index_pt,max_val(index_pt)/max(max_val),'Linewidth',3);
%     plot(index_pt,gauss,'.-','Linewidth',2);  ylab = get(gca,'YLabel');set(ylab,'String','Beam max power','FontSize',26,'FontName','Arial');xlab = get(gca,'XLabel');set(xlab,'String','Beam index','FontSize',26,'FontName','Arial');  set(gca,'XLim', [1 length(error_elevation)],'FontSize',26) 
%      legend ('Max power normalised',[num2str(length(index_pt)) ' central Beams'],['Gauss fitting: Skewness ' num2str(skew(i_pt),3) ' Kurtosis ' num2str(kurt(i_pt),3)],'Location','SouthEast');
%     title(strcat('Stack #',num2str(pt_index(i_pt)),' LAT. ',{' '},num2str(lat_surf(pt_index(i_pt))),'[deg] (',mode_string,')'))
% %     saveas (h, ['./ATTD/PT_',num2str(i_pt),'_stack_power_',TEST,'_',num2str(length(index_pt)),'.jpg']);close (h);

    
    
    % linear regression on stack alignment
    p = polyfit(index_pt,error_elevation(index_pt),1);
               
     h=figure; plot(look_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*180/pi,error_elevation,'Linewidth',3,'DisplayName','Stack alignment'); ylab = get(gca,'YLabel');
     set(ylab,'String','Elevation error [m]','FontSize',26,'FontName','Arial');
     xlab = get(gca,'XLabel');set(xlab,'String','Look angle [deg.]','FontSize',26,'FontName','Arial');  %set(gca,'YLim',[-0.01 0.01],'XLim', [1 length(error_elevation)],'FontSize',26) 
    hold all; plot(look_ang_surf(pt_index(i_pt),index_pt).*180/pi,error_elevation(index_pt),'Linewidth',3,'DisplayName',[num2str(length(index_pt)) ' central Beams']);
    %     figure; plot((max_pos /cnf.zp_fact_range) * T0_nom * cst.c/2); ylab = get(gca,'YLabel');set(ylab,'String','metres','FontSize',26,'FontName','Arial'); % set(gca,'YLim',[41 43],'FontSize',26) 
    set(gca,'YLim',[mean(error_elevation(index_pt))-0.01 mean(error_elevation(index_pt))+0.01])
    set(gca,'XLim', [min(look_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*180/pi) max(look_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt))).*180/pi)],'FontSize',26); 
    yfit =  p(1) * index_pt + p(2);
    plot(look_ang_surf(pt_index(i_pt),index_pt).*180/pi,yfit,'.-','Linewidth',2,'DisplayName',['linear regression: slope ' num2str(p(1)*1000,3) ' mm/beam' '--> range variation ' num2str(p(1)*1e3*(index_pt(end)-index_pt(1)),3) ' mm']);
    leg=legend('show');
    pos_leg=get(leg,'Position');
    title(strcat('Stack #',num2str(pt_index(i_pt)),' LAT. ',{' '},num2str(lat_surf(pt_index(i_pt))),'[deg] (',mode_string,')'),'Interpreter','None')
    annotation('textbox',[pos_leg(1)-0.05,pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),pos_leg(4)],...
        'String',{strcat('Mean elevation error [mm]=',num2str(mean((error_elevation(index_pt)))*1e3)),...
        strcat('Std elevation error [mm]=',num2str(std(detrend(error_elevation(index_pt)))*1e3))},...
        'FitBoxToText','on','Fontsize',26);
    print(print_file,res_fig,strcat(result_Path,input_filename,'Elevation_error_alignment_',mode_string,'_PTR_',num2str(i_pt),file_ext));

    %compute the selection of T0 for the surface (brust right above surface)
    [~,beam_index_nadir]=min(abs(-pi/2+beam_ang_surf(pt_index(i_pt),1:N_beams_stack(pt_index(i_pt)))));
    T0_surf = T0_sar_surf(pt_index(i_pt),beam_index_nadir(1));

    for i_sample = 1:chd.N_samples_sar_sar_chd*cnf.zp_fact_range 
        wfm_cor_i2q2_PT(i_pt,i_sample) = sum(beams_rng_cmpr_PT(1,:,i_sample))/N_beams_stack(pt_index(i_pt));
    end
%     figure; plot(wfm_cor_i2q2_PT(i_pt,:));                                                                                               
    [max_val_L1,max_pos_L1] = max(wfm_cor_i2q2_PT(i_pt,:));
    [~,left_half]= min(abs(wfm_cor_i2q2_PT(i_pt,1:max_pos_L1)-max_val_L1/2));
    [~,right_half]= min(abs(wfm_cor_i2q2_PT(i_pt,max_pos_L1:end)-max_val_L1/2));
    error_elevation_L1(i_pt) = alt_sat(pt_index(i_pt))-(win_delay_surf(pt_index(i_pt))*cst.c/2+((max_pos_L1)-(((chd.N_samples_sar_sar_chd/2+1)-1)*cnf.zp_fact_range+1))./cnf.zp_fact_range*cst.c/2*T0_surf)-transponder_altitude_list_chd(i_pt);
%     samples=62600:67400;
%     f_L1 = fit((samples)', wfm_cor_i2q2_PT(i_pt,samples)','gauss1');
%     gauss_L1 =  f_L1.a1*exp(-(((samples)-f_L1.b1)/f_L1.c1).^2);
    
    %save(strcat(resultPath,'point_target_analysis.mat'));
    
    
    
    Req.std_error_stack(i_pt) = std(detrend(error_elevation(index_pt)))*1e3; % [mm]
    Req.range_variation_stack(i_pt) = p(1)*(index_pt(end)-index_pt(1))*1e3; % [mm]
    Req.mean_error_stack(i_pt) = mean((error_elevation(index_pt)))*1e3; % [mm]
    Req.error_L1B(i_pt) =  error_elevation_L1(i_pt)*1e3; % [mm]
    %Req.5(i_pt) =  ((right_half+max_pos_L1-left_half)/cnf.zp_fact_range-0.866)*cst.c/2*T0_surf*1000; % [mm]
    Req.datation_error(i_pt) = datation_error_v1*1e6; %[seconds]
%     fprintf(fidResults,'%d;%d;%d;%d;%d;%d;%d\n',Req1(i_pt),Req2(i_pt),Req3(i_pt),Req4(i_pt),Req5(i_pt),Req6(i_pt),length(index_pt));
end

save(strcat(result_Path,input_filename,'_point_target_analysis_',file_specific_ext,'.mat'),'Req');


%fclose (fidResults);
%       
% figure; plot(max_val_72);
% figure; plot(max_val_130);
% figure; plot(max_val_189);
% figure; plot(max_val_247);
% 
% figure; plot(max_pos_72 /cnf.zp_fact_range); ylab = get(gca,'YLabel');set(ylab,'String','Retracked Samples','FontSize',26,'FontName','Arial'); set(gca,'YLim',[41 43],'FontSize',26) 
% figure; plot(max_pos_130/cnf.zp_fact_range); ylab = get(gca,'YLabel');set(ylab,'String','Retracked Samples','FontSize',26,'FontName','Arial'); set(gca,'YLim',[41 43],'FontSize',26)
% figure; plot(max_pos_189/cnf.zp_fact_range); ylab = get(gca,'YLabel');set(ylab,'String','Retracked Samples','FontSize',26,'FontName','Arial'); set(gca,'YLim',[41 43],'FontSize',26)
% figure; plot(max_pos_247/cnf.zp_fact_range); ylab = get(gca,'YLabel');set(ylab,'String','Retracked Samples','FontSize',26,'FontName','Arial'); set(gca,'YLim',[41 43],'FontSize',26)

% max_image=max([max(10*log10(beams_rng_PT_CPP(:,:))),max(10*log10(beams_rng_PT(:,:)))]);
% min_image=max_image-40;
% figure;
% subplot(1,2,1)
% imagesc(10*log10(beams_rng_PT)); colormap('jet'); c=colorbar;
% caxis([min_image,max_image]); ylabel(c,'[dBW]');
% title('TRP #1 MAT')
% xlabel('Samples'); ylabel('Beams');
% 
% subplot(1,2,2)
% imagesc(10*log10(beams_rng_PT_CPP)); colormap('jet'); c=colorbar;
% caxis([min_image,max_image]); ylabel(c,'[dBW]');
% title('TRP #1 C++')
% xlabel('Samples'); ylabel('Beams');

% wfm_0_rng_Dopp = fftshift(fft(wfm_0,[],1),1); 
% 
% wfm_2_rng_Dopp = fftshift(fft(ifft(fftshift(wfm_2,2),[],2),[],1),1);
% 
% wfm_3_rng_Dopp = fftshift(fft(ifft(fftshift(wfm_3,2),[],2),[],1),1);
% 
% wfm_4_rng_Dopp = fftshift(fft(ifft(fftshift(wfm_4,2),[],2),[],1),1);
% 
% 
% max_image = 10*log10(max([max(abs(wfm_0_rng_Dopp(:)).^2),...
%                           max(abs(wfm_2_rng_Dopp(:)).^2),...
%                           max(abs(wfm_3_rng_Dopp(:)).^2),...
%                           max(abs(wfm_4_rng_Dopp(:)).^2),...
%                           max(abs(wfm_7).^2)]));
% min_image =max_image-40.0;
% 
% figure; 
% imagesc(10*log10(abs(wfm_0_rng_Dopp).^2));
% colormap('jet'); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% xlabel('Samples'); ylabel('Dopp. bins');
% title(strcat('Input waveform - range & Dopp domains #',num2str(i_burst)))
% 
% figure; 
% imagesc(10*log10(abs(wfm_2_rng_Dopp).^2));
% colormap('jet'); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% xlabel('Samples'); ylabel('Dopp. bins');
% title(strcat('After az. weights - range & Dopp domains #',num2str(i_burst)))
% 
% figure; 
% imagesc(10*log10(abs(wfm_3_rng_Dopp).^2));
% colormap('jet'); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% xlabel('Samples'); ylabel('Dopp. bins');
% title(strcat('After CAI/FAI canc. & radial velocity - range & Dopp domains #',num2str(i_burst)))
% 
% figure; 
% imagesc(10*log10(abs(wfm_4_rng_Dopp).^2));
% colormap('jet'); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% xlabel('Samples'); ylabel('Dopp. bins');
% title(strcat('After Dopp. corr. - range & Dopp domains #',num2str(i_burst)))
% 
% figure; 
% imagesc(10*log10(abs(wfm_7).^2));
% colormap('jet'); c=colorbar; ylabel(c,'[dB]');
% caxis([min_image,max_image]);
% xlabel('Samples'); ylabel('Dopp. bins');
% title(strcat('After RMC correction. - range & Dopp domains #',num2str(i_burst)))
% 
% 
% 
% 
% 
% RMC=load('C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\results\S6_PT12_0001_20170918_lat46_fixedH0\RMC\no_hr_applied_reversion\wks_L1A_PTR_RMC.mat');
% RAW_RMC=load('C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\results\S6_PT12_0001_20170918_lat46_fixedH0\RMC\from_RAW_RMC\dig_gain_corr_begin_quan_no_hr_rever\wks_L1A_RAW_RMC_var_gain_begin_quantization_delta_tau_RMC.mat');
% N_total_bursts_isp=RMC.N_total_bursts_isp;
% chd.N_pulses_burst=RMC.chd.N_pulses_burst;
% chd.N_samples_sar_sar_chd=RMC.chd.N_samples_sar_sar_chd;
% pri_sar_pre_dat=RMC.pri_sar_pre_dat;
% prf_mean=mean(1./(pri_sar_pre_dat));
% resultPath=RAW_RMC.resultPath;
% Doppler_freq=-prf_mean/2+(0:1:chd.N_pulses_burst-1).*prf_mean/chd.N_pulses_burst;
% 
% figure_format='jpg';
% res_fig='-r300';
% 
% switch lower(figure_format)
%     case 'pdf'
%         file_ext='.pdf';
%         print_file='-dpdf';
%     case 'eps'
%         file_ext='.eps';
%         print_file='-depsc';
%     case 'png'
%         file_ext='.png';
%         print_file='-dpng';
%     case 'jpg'
%         file_ext='.jpg';
%         print_file='-djpeg';        
% end
% 
% %% compute the waveforms in range-Pulse and range-Doppler domains
% i_echo=0;
% wvfms_rng_pulse_rearranged_RMC=zeros(N_total_bursts_isp*chd.N_pulses_burst,chd.N_samples_sar_sar_chd);
% wvfms_rng_Doppler_RMC=zeros(N_total_bursts_isp,chd.N_pulses_burst,chd.N_samples_sar_sar_chd);
% wvfms_rng_pulse_rearranged_RAW_RMC=zeros(N_total_bursts_isp*chd.N_pulses_burst,chd.N_samples_sar_sar_chd);
% wvfms_rng_Doppler_RAW_RMC=zeros(N_total_bursts_isp,chd.N_pulses_burst,chd.N_samples_sar_sar_chd);
% for i_burst=1:N_total_bursts_isp
%     dumm_RMC=fft(fftshift(squeeze(RMC.wfm_cal_gain_corrected(i_burst,:,:)),2),[],2);
%     wvfms_rng_Doppler_RMC(i_burst,:,:)=fftshift(fft(dumm_RMC,[],1),1)./sqrt(chd.N_pulses_burst);
%     dumm_RAW_RMC=fft(fftshift(squeeze(RAW_RMC.wfm_cal_gain_corrected(i_burst,:,:)),2),[],2);
%     wvfms_rng_Doppler_RAW_RMC(i_burst,:,:)=fftshift(fft(dumm_RAW_RMC,[],1),1)./sqrt(chd.N_pulses_burst);
%     for i_pulse=1:chd.N_pulses_burst
%         i_echo=i_echo+1;
%         wvfms_rng_pulse_rearranged_RMC(i_echo,:)=dumm_RMC(i_pulse,:);
%         wvfms_rng_pulse_rearranged_RAW_RMC(i_echo,:)=dumm_RAW_RMC(i_pulse,:);
%     end
% end
% %waveforms in range-Pulse domain
% max_image=10*log10(max([max(abs(wvfms_rng_pulse_rearranged_RMC(:)).^2),...
% max(abs(wvfms_rng_pulse_rearranged_RAW_RMC(:)).^2)]));
% min_image=max_image-40.0;
% figure;
% subplot(2,2,1);
% imagesc(10*log10(abs(wvfms_rng_pulse_rearranged_RMC).^2).');
% colormap('jet'); c=colorbar; caxis([min_image,max_image]);
% ylabel(c,'[dBW]'); xlabel('Pulses'); ylabel('Range samples');
% title('RMC');
% subplot(2,2,2);
% imagesc(10*log10(abs(wvfms_rng_pulse_rearranged_RAW_RMC).^2).');
% colormap('jet'); c=colorbar; caxis([min_image,max_image]);
% ylabel(c,'[dBW]'); xlabel('Pulses'); ylabel('Range samples');
% title('RAW --> RMC');
% subplot(2,2,3);
% imagesc(abs(abs(wvfms_rng_pulse_rearranged_RAW_RMC).^2-abs(wvfms_rng_pulse_rearranged_RMC).^2).');
% colormap('jet'); c=colorbar;
% ylabel(c,'[W]'); xlabel('Pulses'); ylabel('Range samples');
% title('Absolute error');
% subplot(2,2,4);
% imagesc(abs((abs(wvfms_rng_pulse_rearranged_RAW_RMC).^2-abs(wvfms_rng_pulse_rearranged_RMC).^2)./(abs(wvfms_rng_pulse_rearranged_RAW_RMC).^2).*100).');
% colormap('jet'); c=colorbar; caxis([0,100])
% ylabel(c,'[%]'); xlabel('Pulses'); ylabel('Range samples');
% title('Relative error');
% [axT,hT]=suplabel('Waveforms in range-Pulse domain','t');

% % Waveforms in the range-Doppler domain
% max_image=10*log10(max([max(abs(wvfms_rng_Doppler_RMC(:)).^2),...
% max(abs(wvfms_rng_Doppler_RAW_RMC(:)).^2)]));
% min_image=max_image-40.0;
% for i_burst=1:N_total_bursts_isp
%     
%     f1=figure;
%     subplot(2,2,1);
%     imagesc(Doppler_freq/1e3,1:chd.N_samples_sar_sar_chd,10*log10(abs(squeeze(wvfms_rng_Doppler_RMC(i_burst,:,:))).^2).');
%     colormap('jet'); c=colorbar; caxis([min_image,max_image]);
%     ylabel(c,'[dBW]'); xlabel('Doppler [KHz]'); ylabel('Range samples');
%     title('RMC');
%     subplot(2,2,2);
%     imagesc(Doppler_freq/1e3,1:chd.N_samples_sar_sar_chd,10*log10(abs(squeeze(wvfms_rng_Doppler_RAW_RMC(i_burst,:,:))).^2).');
%     colormap('jet'); c=colorbar; caxis([min_image,max_image]);
%     ylabel(c,'[dBW]'); xlabel('Doppler [KHz]'); ylabel('Range samples');
%     title('RAW --> RMC');
%     subplot(2,2,3);
%     imagesc(Doppler_freq/1e3,1:chd.N_samples_sar_sar_chd,abs((abs(squeeze(wvfms_rng_Doppler_RMC(i_burst,:,:))).^2).' - ...
%         (abs(squeeze(wvfms_rng_Doppler_RAW_RMC(i_burst,:,:))).^2).'));
%     colormap('jet'); c=colorbar;
%     ylabel(c,'[W]'); xlabel('Doppler [KHz]'); ylabel('Range samples');
%     title('Absolute error');
%     subplot(2,2,4);
%     imagesc(Doppler_freq/1e3,1:chd.N_samples_sar_sar_chd,abs((abs(squeeze(wvfms_rng_Doppler_RMC(i_burst,:,:))).^2).' - ...
%         (abs(squeeze(wvfms_rng_Doppler_RAW_RMC(i_burst,:,:))).^2).')./((abs(squeeze(wvfms_rng_Doppler_RMC(i_burst,:,:))).^2).')*100);
%     colormap('jet'); c=colorbar; caxis([0,100])
%     ylabel(c,'[%]'); xlabel('Doppler [KHz]'); ylabel('Range samples');
%     title('Relative error');
%     [axT,hT]=suplabel('Waveforms in range-Doppler domain','t');
%     print(print_file,res_fig,strcat(resultPath,'wvfms_Doppler_domain_RMC_RAW_',num2str(i_burst),file_ext));
%     close(f1);
%     
% end
