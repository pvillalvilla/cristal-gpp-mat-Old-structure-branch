SAVE      = 1;
VISIBLE   = 1;
set_default_plot;
plots_folder = [filesBulk.outputPath filesBulk.filename_L0(end-56:end-3)];
if(~exist (plots_folder, 'dir'))
    mkdir(plots_folder);
end

font_size=20;
%% PLOT Time 
figure; plot([L1A_buffer(:).time],'.');
set(gca,'XLim',[1 N_bursts],'FontSize',font_size) % --> eixos
xlab = get(gca,'XLabel'); set(xlab,'String','Record indexes','FontSize',font_size)
ylab = get(gca,'YLabel'); set(ylab,'String','seconds','FontSize',font_size);
if SAVE == 1
    figName = [plots_folder '/11_Orbit_time_for_burst'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
end

figure; plot(diff([L1A_buffer(:).time]),'.')
set(gca,'XLim',[1 N_bursts-1],'FontSize',font_size) % --> eixos
ylab = get(gca,'YLabel'); set(ylab,'String','seconds','FontSize',font_size)
if SAVE == 1
    figName = [plots_folder '/B12_Spacing_between_L1_orbit_times'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end


% 1. Average waveforms within a burst, containing all the data set, with AGC
% attenuations applied.
zp = 1;
wfm_isp_av              = zeros(N_bursts,chd.N_samples_sar);
wfm_cal_gain_av         = zeros(N_bursts,chd.N_samples_sar*zp);
wd_shift                = zeros(1,N_bursts);
wfm_cal_gain_wdcorr     = zeros(N_bursts,chd.N_samples_sar*zp);

clear wfmAUX wfmAUX2 max_isp max_cal_gain
for i_burst_plot = 1:N_bursts
    % FFT & power
    %         wfmAUX = abs(fftshift(fft(squeeze(wfm_iq_sar_ku_isp(i_burst,:,:)).',chd.N_samples_sar*zp),1).').^2/(chd.N_samples_sar);
    wfmAUX2 = abs(fftshift(fft(squeeze(L1A_buffer(i_burst_plot).wfm_cal_gain_corrected(:,:)).',chd.N_samples_sar*zp),1).').^2/(chd.N_samples_sar);
    
    % Average
    for i_sample = 1:chd.N_samples_sar*zp
        %             wfm_isp_av(i_burst,i_sample) = sum(wfmAUX(:,i_sample))/chd.N_pulses_burst;
        wfm_cal_gain_av(i_burst_plot,i_sample) = sum(wfmAUX2(:,i_sample))/chd.N_pulses_burst;
    end
end
%     max_isp = max(wfm_isp_av.');
max_burst = max(wfm_cal_gain_av.');
max_range = max(wfm_cal_gain_av);
%******************************************************
% ISP waveforms
figure; imagesc(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_cal_gain_av.')
set(gca,'XLim',[1 N_bursts])
set(gca,'YLim',[1 chd.N_samples_sar*zp])
xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
if SAVE == 1
    figName = [plots_folder '/A12_isp_waveforms'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.png'])
end


%******************************************************
% ISP waveforms with H0 and L1_WinDelay applied
wfm_cal_gain_av_aux = wfm_cal_gain_av;
for i_burst_plot = 1:N_bursts
    wd_shift(i_burst_plot) = L1A_buffer(i_burst_plot).win_delay/ chd.T0_nom + 20;
    wfm_cal_gain_wdcorr(i_burst_plot,:) = circshift(wfm_cal_gain_av_aux(i_burst_plot,:),[0,round(wd_shift(i_burst_plot)*zp)]);
end
figure; mesh(1:N_bursts, 1:chd.N_samples_sar*zp, wfm_cal_gain_wdcorr.');
set(gca,'XLim',[1 N_bursts]) % --> eixos
set(gca,'YLim',[1 chd.N_samples_sar*zp]) % --> eixos
xlab = get(gca,'XLabel'); set(xlab,'String','Bursts')
ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
if SAVE == 1
    figName = [TEST,'_A13_isp_waveforms_with_ATT_CAL_H0_WD_applied'];
    print(gcf, figName,'-dpng')
    saveas (gcf,[figName,'.fig'])
end



%% PLOT locations

figure; plot3([L1A_buffer(:).lat_sar_surf].',[L1A_buffer(:).lon_sar_surf].',[L1A_buffer(:).alt_sar_surf].');

