path=pwd;
file=dir('*.nc');
% lat=ncread(file(1).name,'lat_85_ku');
index=1:1969;

% tot_gain_ch1_85_ku=ncread(file(1).name,'tot_gain_ch1_85_ku');
% noise_power_85_ku=ncread(file(1).name,'noise_power_85_ku');
% instr_cor_gain_tx_rx_85_ku=ncread(file(1).name,'instr_cor_gain_tx_rx_85_ku');
N_samples_sar_chd=256;
burst=zeros(N_samples_sar_chd,64,length(index));
burst_phase=zeros(N_samples_sar_chd,64,length(index));

zp=1;

%     sig0_cal_ku_l1a_echo_sar_ku =ncread(file(1).name,'sig0_cal_ku_l1a_echo_sar_ku',index,[1]);
%     iq_scale  = ncread(file(1).name,'scale_factor_ku_l1a_echo_sar_ku',[index(1) index(end)],[1]);


for i_burst=1:10:length(index)
    i_samples=ncread(file(1).name,'rx1_complex_waveforms_i_samples',[1 1 index(i_burst)],[N_samples_sar_chd 64 1]);
    q_samples=ncread(file(1).name,'rx1_complex_waveforms_q_samples',[1 1 index(i_burst)],[N_samples_sar_chd 64 1]);
    
    burst(:,:,i_burst)   = ((i_samples+1i*q_samples));
%     burst_phase(:,:,i_burst)   = angle(fft(squeeze(burst(:,:,i_burst)),N_samples_sar_chd*zp));
    
end

% for i_burst = 1:length(index)
%     % FFT & power
%     %         wfmAUX = abs(fftshift(fft(squeeze(wfm_iq_sar_ku_isp(i_burst,:,:)).',N_samples_sar_chd*zp),1).').^2/(N_samples_sar_chd);
%     wfmAUX2 = abs(fftshift(fft(squeeze(burst(:,:,i_burst)),N_samples_sar_chd*zp),1).').^2/(N_samples_sar_chd);
%     %         wfmAUX2         = abs(squeeze(wfm_iq_sar_ku_isp_RAW2RMC(i_burst,:,:)));
%     %         wfmAUX2_power_scaled  = 10*log10(abs(squeeze(wfmAUX2.*2)));
%     wfm_cal_gain(i_burst,:)   = sum(wfmAUX2,1)/64;
%     % Average
%     
% end

for i_burst = 1:10:length(index)
% 2D FFT
h=figure(1);imagesc(1/64:1/64:64-1/64,1/64:1/64:N_samples_sar_chd-1/64,20*log10(abs((fft(fftshift(fft(squeeze(burst(:,:,i_burst)).',64*64).',2),N_samples_sar_chd*64))).^2/(N_samples_sar_chd)));
xlab = get(gca,'XLabel'); set(xlab,'String','Doppler Beams')
ylab = get(gca,'YLabel'); set(ylab,'String','Range Bins')
title(['Burst # ' num2str(i_burst)]);% ' Index ' num2str((index(i_burst))) ]);
hcb=colorbar; colormap(colormap_blues)
hcb.Label.String = 'Unscaled Power [dB]';
caxis([-400 -200]);
    pause(0.1);
    saveas (h,['Burst_' num2str(i_burst, '%03d') '.png']);
end 
    %     max_isp = max(wfm_isp_av.');
%     max_cal_gain = max(wfm_cal_gain.');
%     disp(' Generating 1 2D averaged bursts for RAW2RMC');
%     figure;
%     subplot(2,2,1);
%     
%     imagesc(1:64,1:N_samples_rmc_chd*zp, 10*log10(abs(squeeze(wfm_iq_sar_ku_isp_RMC(burst_reference-init_burst+1,:,:).*2.^rx_asic_shift_s_ku_isp_RMC(burst_reference-init_burst+1)))).')
%     xlab = get(gca,'XLabel'); set(xlab,'String','Pulses')
%     ylab = get(gca,'YLabel'); set(ylab,'String','Samples')
%     title(['Burst # ' num2str(burst_reference)  ]);
%     hcb=colorbar; colormap(jet)
%     hcb.Label.String = 'Amplitude [counts]';