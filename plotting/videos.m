%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code CREATE PLOTS FOR VIDEOS
% ---------------------------------------------------------
% Objective: Plo. 
% 
% INPUTs : - L1B
% OUTPUTs: - L1A, L1BS, L1B, and the corresponding breakpoints
%
% ----------------------------------------------------------
% Author:    Albert Garcia / isardSAT
% Reviewer:  M�nica Roca   / isardSAT
% Last rev.: M�nica Roca   / isardSAT (11/09/2013)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mida = get(0,'ScreenSize');


%% backup RMC variables
% Compare
beams_surf_RMC                  = beams_surf;
beams_rng_cmpr_RMC              = beams_rng_cmpr;
% beams_focused_shifted_surf_RMC  = beams_focused_shifted_surf;
wfm_cor_i2q2_RMC         = wfm_cor_i2q2;
wfm_iq_isp_RMC           = wfm_iq_isp;
wfm_sar_reversed_RMC            = wfm_sar_reversed;
wfm_geo_corr_RMC                = wfm_geo_corr;

save('./results/RMC_4Bwaveforms.mat', 'beams_surf_RMC', 'beams_rng_cmpr_RMC', 'wfm_cor_i2q2_RMC', 'wfm_iq_isp_RMC', 'wfm_sar_reversed_RMC','wfm_geo_corr_RMC');


%% RMC REVERSED VS SAR RAW in the RANGE WINDOW

h = figure('Position',mida);
cnf.zp_fact_range = 1;
x_max = 60;
alt_sar_surf(150)= alt_trp(1);
alt_sar_surf(610)= alt_trp(2);
alt_sar_surf(1068)= alt_trp(3);
alt_sar_surf(1526)= alt_trp(4);

for i_burst=1: N_total_bursts_isp-800

   clear plot_RMC plot_RAW fft_RMC fft_RAW
   plot_RMC=zeros(64,chd.N_samples_sar_sar_chd); 
   plot_RMC(:,1:128) = abs(wfm_iq_isp_RMC(i_burst,:,:));
   plot_RAW(:,:) = wfm_iq_isp(i_burst,:,:);
   plot_RAW2(:,:) = wfm_iq2_isp(i_burst,:,:);
   if(i_burst> 830)
       x_max = 128;
   end
%    fft_RMC(:,:)  = abs(fftshift(fft(plot_RMC.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
%    fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
%    fft_RAW2(:,:)  = abs(fftshift(fft(plot_RAW2.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
%    
   subplot(3,1,1);plot(alt_sar_surf(1:i_burst),'k', 'LineWidth',3); figlabels('Burst Index','Window Position [m]','',['Tracking Window Position for Burst ',num2str(i_burst)],18);set(gca,'XLim',[1 N_total_bursts_isp-800],'FontSize',12);set(gca,'YLim',[0 max(max(alt_sar_surf),max(alt_trp))+3],'FontSize',12) ;
   hold on; plot(i_burst,alt_sar_surf(i_burst),'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
   plot(150,alt_trp(1),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');
   if(i_burst>300)  plot(610,alt_trp(2),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
   if(i_burst>800) plot(1068,alt_trp(3),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
   if(i_burst>1200)plot(1526,alt_trp(4),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end 
   hold off  
    
   subplot(3,1,2);imagesc(1:chd.N_samples_sar_sar_chd,1:64,plot_RMC);figlabels('','Pulses','','RMC reversed',12);colormap(hot);freezeColors;set(gca,'XLim',[1 x_max],'FontSize',12) ;
%    subplot(3,2,4);imagesc(1:chd.N_samples_sar_sar_chd,1:64,fft_RAW);figlabels('','Pulses','','RAW ',12);colormap(hot);freezeColors;set(gca,'XLim',[1 256],'FontSize',12) ;
   subplot(3,1,3);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(plot_RMC),'LineWidth',3);figlabels('Samples','','','',12);colormap(hot);freezeColors;set(gca,'XLim',[0 x_max],'FontSize',12) ; set(gca,'YLim',[0 8000],'FontSize',12) ;
%    subplot(3,2,6);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW),'LineWidth',3);figlabels('Samples','','',' ',12);colormap(hot);freezeColors;set(gca,'XLim',[1 256],'FontSize',12);
%    hold all;subplot(3,2,6);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW2),'LineWidth',3);figlabels('Samples','','',' ',12);colormap(hot);freezeColors;set(gca,'XLim',[1 256],'FontSize',12);
   
   pause(0.01)
%    saveas (h, ['./results/reversion/',num2str(i_burst, '%04d'),'_RMCreverted_vs_RAW.png']);
    
   
%    subplot(3,2,5);plot(abs(sum(plot_RMC,1))); figlabels('Samples','','','RMC ISP Multilooked ',12);set(gca,'XLim',[1 256],'FontSize',12); 
    
    
end

i_burst=1540;
zp_factor=8;
for i_beam=1:64
    clear plot_exact plot_aprox fft_exact fft_aprox fft_exact_zp fft_aprox_zp
    plot_exact = squeeze(beams_focused_shifted_exact(i_burst,i_beam,:)); 
    plot_aprox = squeeze(beams_focused_shifted_aprox(i_burst,i_beam,:));
%     fft_exact= plot_exact.';
%     fft_aprox= plot_aprox.';
     fft_exact(:,:)  = abs(fftshift(fft(plot_exact,256),1));
     fft_aprox(:,:)  = abs(fftshift(fft(plot_aprox,256),1));
%     [~,max_pos]= max(fft_exact);
%     fft_exact_zp(:,:)  = abs(fftshift(fft(plot_exact,256*zp_factor),1));
%     fft_aprox_zp(:,:)  = abs(fftshift(fft(plot_aprox,256*zp_factor),1));
    rms(i_beam)=sqrt(sum((fft_aprox(isfinite(fft_aprox))-fft_exact(isfinite(fft_aprox))).^2));
%     rms_zp(i_beam)=sqrt(sum((fft_aprox_zp-fft_exact_zp).^2));
    h=figure; plot(fft_exact,'linewidth',4); hold all; plot(fft_aprox,'g','linewidth',2); 
%     plot(1:1/zp_factor:(257-1/zp_factor),fft_exact_zp,'b.','linewidth',1); plot(1:1/zp_factor:(257-1/zp_factor),fft_aprox_zp,'g.','linewidth',1); 
    legend ('exact method', 'approx method');%, 'exact method ZP 32', 'approx method ZP 32');
    title(['Beam ' num2str(i_beam, '%02d'),'  RMS= ', num2str(rms(i_beam))]);
%     set(gca,'XLim',[max_pos-5 max_pos+5]);
    set(gca,'XLim',[1 256]);
    saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_burst_',num2str(i_beam, '%02d'),'_beam_EXACT_vs_APPOXIMATE.png']);
    close(h)
end
h=figure; plot(rms);
figlabels('L1B index','RMS','','',30);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_burst_slant_doppler_WD_RMS_EXACT_vs_APPOXIMATE.fig']);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_burst_slant_doppler_WD_RMS_EXACT_vs_APPOXIMATE.png']);
close(h);
% h=figure; imagesc(abs(fftshift(fft(squeeze(beams_wdcorr_exact(i_burst,:,:)).',256*zp_factor),1)).');colormap(colormapATDD);colorbar;
h=figure; imagesc((((squeeze(beams_wdcorr_exact(i_burst,:,:)).'))).');colormap(colormapATDD);colorbar;
figlabels('Samples','Beam','','',30);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_burst_slant_doppler_WD_zp1_EXACT.png']);
% h=figure;imagesc(abs(fftshift(fft(squeeze(beams_wdcorr_aprox(i_burst,:,:)).',256*zp_factor),1)).');colormap(colormapATDD);colorbar;
h=figure; imagesc((((squeeze(beams_wdcorr_aprox(i_burst,:,:)).'))).');colormap(colormapATDD);colorbar;
figlabels('Samples','Beam','','',30);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_burst_slant_doppler_WD_zp1_APPROXIMATE.png']);
close(h);


for i_beam=1:459
    clear plot_exact plot_aprox fft_exact fft_aprox
    plot_exact = squeeze(wfm_geo_corr_exact(i_burst,i_beam,:)); 
    plot_aprox = squeeze(wfm_geo_corr_aprox(i_burst,i_beam,:)); 
    fft_exact(:,:)  = abs(fftshift(fft(plot_exact,256*8),1));
    fft_aprox(:,:)  = abs(fftshift(fft(plot_aprox,256*8),1));
    rms(i_beam)=sqrt(sum((fft_aprox-fft_exact).^2));
    h=figure; plot(1:1/zp_factor:(256-1/zp_factor),fft_exact,'linewidth',4); hold all; plot(1:1/zp_factor:(256-1/zp_factor),fft_aprox,'g','linewidth',2); legend ('exact method', 'approx method');
    title(['Beam ' num2str(i_beam, '%02d'),'  RMS= ', num2str(rms(i_beam))]);
    set(gca,'XLim',[1 256]);
    saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_stack_aligned_',num2str(i_beam, '%02d'),'_beam_aligned_EXACT_vs_APPOXIMATE.png']);
    close(h)
end
h=figure; plot(rms);
figlabels('Beam index','RMS','','',30);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_stack_RMS_EXACT_vs_APPOXIMATE.fig']);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_stack_RMS_EXACT_vs_APPOXIMATE.png']);
close(h);
h=figure; imagesc(abs(fftshift(fft(squeeze(wfm_geo_corr_exact(i_burst,:,:)).',256*zp_factor),1)).');colormap(colormapATDD);colorbar;
figlabels('Samples','Beam','','',30);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_stack_zp1_aligned_EXACT.png']);
h=figure; imagesc(abs(fftshift(fft(squeeze(wfm_geo_corr_aprox(i_burst,:,:)).',256*zp_factor),1)).');colormap(colormapATDD);colorbar;
figlabels('Samples','Beam','','',30);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_stack_zp1_aligned_APPROXIMATE.png']);
close(h);


figure; plot(wfm_cor_i2q2_exact(i_burst,:));hold all; plot(wfm_cor_i2q2_aprox(i_burst,:),'g','linewidth',2); legend ('exact method', 'approx method');
set(gca,'XLim',[1 512]);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_L1B_EXACT_vs_APPROXIMATE.png']);
set(gca,'XLim',[50 100]);
saveas (h, ['./results/Exact vs Approximate/',num2str(i_burst, '%04d'),'_L1B_EXACT_vs_APPROXIMATE_zoom.png']);


for i_burst=1: N_total_bursts_isp-800

   clear plot_RMC plot_RMC_az fft_RMC fft_RMC_az
   
   plot_RMC=zeros(64,chd.N_samples_sar_sar_chd); 
   plot_RMC(:,:) = wfm_sar_reversed_RMC(i_burst,:,:);
   
   plot_RAW(:,:) = wfm_iq_isp(i_burst,:,:);
   
   if(i_burst> 830)
       x_max = 128;
   end
   fft_RMC(:,:)  = abs(fftshift(fft(plot_RMC.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
   
   fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
  
   
   subplot(3,1,1);plot(alt_sar_surf(1:i_burst),'k', 'LineWidth',3); figlabels('Burst Index','Window Position [m]','',['Tracking Window Position for Burst ',num2str(i_burst)],18);set(gca,'XLim',[141 168],'FontSize',12);set(gca,'YLim',[0 max(max(alt_sar_surf),max(alt_trp))+3],'FontSize',12) ;
   hold on; plot(i_burst,alt_sar_surf(i_burst),'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
   plot(150,alt_trp(1),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');
   if(i_burst>300)  plot(610,alt_trp(2),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
   if(i_burst>800) plot(1068,alt_trp(3),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
   if(i_burst>1200)plot(1526,alt_trp(4),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end 
   hold off  
    
   subplot(3,2,3);imagesc(1:chd.N_samples_sar_sar_chd,1:64,fft_RMC);figlabels('','Pulses','','RAW',12);colormap(hot);freezeColors;set(gca,'XLim',[1 x_max],'FontSize',12) ;
   subplot(3,2,4);imagesc(1:chd.N_samples_sar_sar_chd,1:64,fft_RAW);figlabels('','Pulses','','RAW I2Q2 ',12);colormap(hot);freezeColors;set(gca,'XLim',[1 x_max],'FontSize',12) ;
   subplot(3,2,5);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC),'LineWidth',3);figlabels('Samples','','',' ',12);colormap(hot);freezeColors;set(gca,'XLim',[1 x_max],'FontSize',12);
   subplot(3,2,6);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW),'LineWidth',3);figlabels('Samples','','',' ',12);colormap(hot);freezeColors;set(gca,'XLim',[1 x_max],'FontSize',12);
   
   pause(0.01)
%    saveas (h, ['./results/reversion/',num2str(i_burst, '%04d'),'_RAW_vs_RAWI2Q2.png']);
    
   
%    subplot(3,2,5);plot(abs(sum(plot_RMC,1))); figlabels('Samples','','','RMC ISP Multilooked ',12);set(gca,'XLim',[1 256],'FontSize',12); 
    
    
end


%% RMC ISP VS RMC reversed Azimuth


cnf.zp_fact_range = 4;
x_min=1;
x_max = 256;
alt_sar_surf(150)= alt_trp(1);
alt_sar_surf(610)= alt_trp(2);
alt_sar_surf(1068)= alt_trp(3);
alt_sar_surf(1526)= alt_trp(4);
h = figure('Position',mida);

for i_burst=1: N_total_bursts_isp

   clear plot_RMC plot_RAW fft_RMC fft_RAW
   plot_RMC=zeros(64,chd.N_samples_sar_sar_chd); 
   plot_RMC(:,:) = wfm_sar_reversed_RMC(i_burst,:,:);
   plot_RAW(:,:) = wfm_iq_isp(i_burst,:,:);
   
%    plot_RMC_az(:,:) = beams_focused_shifted_surf_RMC(i_burst,:,:);
%    if(i_burst> 830)
%        x_max = chd.N_samples_sar_sar_chd;
%    end
   fft_RMC(:,:)  = abs(fftshift(fft(plot_RMC.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
   fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
  
%    fft_RMC_az(:,:)  = abs(fftshift(fft(plot_RMC_az.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
   
   subplot(3,1,1);plot(41.743253 + alt_sar_surf(1:i_burst),'k', 'LineWidth',3); figlabels('Burst Index','Window Position [m]','',['Tracking Window Position for Burst ',num2str(i_burst)],18);set(gca,'XLim',[1 N_total_bursts_isp],'FontSize',12); set(gca,'YLim',[min(alt_sar_surf)+38.743253, max(alt_sar_surf)+44.743253],'FontSize',12) ;
%    hold on; plot(i_burst,41.743253 + alt_sar_surf(i_burst),'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
%    plot(150,alt_trp(1),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');
%    if(i_burst>300)  plot(610,alt_trp(2),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
%    if(i_burst>800) plot(1068,alt_trp(3),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
%    if(i_burst>1200)plot(1526,alt_trp(4),'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end 
%    hold off  
    
   subplot(3,2,3);imagesc(1:chd.N_samples_sar_sar_chd,1:64,fft_RMC);figlabels('','Pulses','','RMC reversed',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
   subplot(3,2,4);imagesc(1:chd.N_samples_sar_sar_chd,1:64,fft_RAW);figlabels('','Pulses','','RAW ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
   subplot(3,1,3);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC),'LineWidth',3);figlabels('Samples','','','',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ; set(gca,'YLim',[0 4000],'FontSize',12) ;
   hold all; subplot(3,1,3);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW)/sqrt(256)*2,'LineWidth',3);figlabels('Samples','','',' ',12);colormap(hot);%set(gca,'XLim',[40 6],'FontSize',12); set(gca,'YLim',[0 8000],'FontSize',12) ;
   hold off
   pause(0.01)
   saveas (h, ['./results/reversion/',num2str(i_burst, '%04d'),'_RMCreverted_vs_RAW_S3.png']);
    
   
%    subplot(3,2,5);plot(abs(sum(plot_RMC,1))); figlabels('Samples','','','RMC ISP Multilooked ',12);set(gca,'XLim',[1 256],'FontSize',12); 
    
    
end
   
    
%% RMC VS SAR RAW before and after GEOMETRY CORRECTIONS

%%SC2
h = figure('Position',mida);
max_beams=max(N_beams_stack);
for i_surf=1: N_total_surf_loc
    plot_RAW=zeros(max_beams,chd.N_samples_sar_sar_chd); 
    plot_RAW(:,:) = wfm_geo_corr(i_surf,:,:);
    fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       kk=surf(1:chd.N_samples_sar_sar_chd*cnf.zp_fact_range,1:max_beams,fft_RAW);figlabels('Samples','Doppler Beams','','',30);colormap(colormapATDD);colorbar;
       set(kk, 'edgecolor','none');colormap(colormapATDD);view(45,30)
%            imagesc(1:chd.N_samples_sar_sar_chd*cnf.zp_fact_range,1:max_beams,fft_RAW);figlabels('Samples','Doppler Beams','','',30);colormap(colormapATDD);colorbar;
     saveas (h, ['./results/SC2/RMC/3-2/',num2str(i_surf, '%04d'),'_RAW_STACKS_3D_geo_S1.png']);   
end

figure; plot(wfm_cor_i2q2_wdcorr(50,:));figlabels('Samples','','','',30);set(gca,'XLim',[1 chd.N_samples_sar_sar_chd*cnf.zp_fact_range],'FontSize',30);
figure; plot(wfm_cor_i2q2_wdcorr(80,:));figlabels('Samples','','','',30);set(gca,'XLim',[1 chd.N_samples_sar_sar_chd*cnf.zp_fact_range],'FontSize',30);
figure; plot(wfm_cor_i2q2_wdcorr(70,:));figlabels('Samples','','','',30);set(gca,'XLim',[1 chd.N_samples_sar_sar_chd*cnf.zp_fact_range],'FontSize',30);

figure; imagesc(wfm_cor_i2q2_wdcorr);figlabels('Samples','L1 record','','',30);
figure; imagesc(wfm_cor_i2q2);figlabels('Samples','L1 record','','',30);
colormap(colormapATDD);colorbar;

if(trp_flag)
    h = figure('Position',mida);
    max_beams=max(N_beams_stack);
    x_min =1;
    x_max=256;
    y_max = 40;
    for i_surf=1: N_total_surf_loc

       clear plot_RMC plot_RAW fft_RMC fft_RAW plot_RMC_ plot_RAW_ fft_RMC_ fft_RAW_
       plot_RMC=zeros(max_beams,chd.N_samples_sar_sar_chd); plot_RMC_=zeros(max_beams,chd.N_samples_sar_sar_chd);
       plot_RAW=zeros(max_beams,chd.N_samples_sar_sar_chd); plot_RAW_=zeros(max_beams,chd.N_samples_sar_sar_chd);

       plot_RMC(:,:) = beams_surf_RMC(i_surf,:,:); plot_RMC_(:,:) = wfm_geo_corr_RMC(i_surf,:,:);
       plot_RAW(:,:) = beams_surf(i_surf,:,:);     plot_RAW_(:,:) = wfm_geo_corr(i_surf,:,:);
       fft_RMC(:,:)  = abs(fftshift(fft(plot_RMC.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       fft_RMC_(:,:)  = abs(fftshift(fft(plot_RMC_.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       fft_RAW_(:,:)  = abs(fftshift(fft(plot_RAW_.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');


       subplot(1,3,3);plot(alt_surf(1:i_surf),1:i_surf,'k', 'LineWidth',3); figlabels('Window Position [m]','Surface Index','',['Tracking Window Position ',num2str(i_surf)],18);set(gca,'YLim',[1 N_total_surf_loc],'FontSize',12);set(gca,'XLim',[0 max(max(alt_sar_surf),max(alt_trp))+3],'FontSize',12) ;set(gca,'YDir','reverse');
       hold on; plot(alt_surf(i_surf),i_surf,'ro','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','k'); 
       plot(alt_trp(1),30,'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');
%        if(i_surf>60)  plot(alt_trp(2),90,'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
%        if(i_surf>120) plot(alt_trp(3),150,'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end
%        if(i_surf>180) plot(alt_trp(4),210,'ro','MarkerSize',4,'MarkerFaceColor','g','MarkerEdgeColor','k');end 
       hold off

       subplot(3,3,1);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RMC);figlabels('','Doppler Beams','','Stack RMC ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,2);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RAW);figlabels('','Doppler Beams','','Stack RAW ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,4);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RMC_);figlabels('','Doppler Beams','','Stack aligned RMC ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,5);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RAW_);figlabels('','Doppler Beams','','Stack aligned RAW ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,7);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC_)            ,'LineWidth',3);figlabels('Samples','','','Stack Sum RMC',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);
       subplot(3,3,8);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW_)/sqrt(256)*2,'LineWidth',3);figlabels('Samples','','','Stack Sum RAW',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);

       pause(0.01)
       saveas (h, ['./results/stacks/',num2str(i_surf, '%04d'),'_RMC_vs_RAW_STACKS_S1.png']);


    end
else
    h = figure('Position',mida);
    max_beams=max(N_beams_stack);
    x_min =1;
    x_max=256;
    y_max= 12e-5;
    for i_surf=1: N_total_surf_loc

       clear plot_RMC plot_RAW fft_RMC fft_RAW plot_RMC_ plot_RAW_ fft_RMC_ fft_RAW_
       plot_RMC=zeros(max_beams,chd.N_samples_sar_sar_chd); plot_RMC_=zeros(max_beams,chd.N_samples_sar_sar_chd);
       plot_RAW=zeros(max_beams,chd.N_samples_sar_sar_chd); plot_RAW_=zeros(max_beams,chd.N_samples_sar_sar_chd);

       plot_RMC(:,:) = beams_surf_RMC(i_surf,:,:); plot_RMC_(:,:) = wfm_geo_corr_RMC(i_surf,:,:);
       plot_RAW(:,:) = beams_surf(i_surf,:,:);     plot_RAW_(:,:) = wfm_geo_corr(i_surf,:,:);
       fft_RMC(:,:)  = abs(fftshift(fft(plot_RMC.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       fft_RMC_(:,:)  = abs(fftshift(fft(plot_RMC_.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
       fft_RAW_(:,:)  = abs(fftshift(fft(plot_RAW_.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');


       subplot(1,3,3);plot(41.743253 + alt_surf(1:i_surf),1:i_surf,'k', 'LineWidth',3); figlabels('Window Position [m]','Surface Index','',['Tracking Window Position ',num2str(i_surf)],18);set(gca,'YLim',[1 N_total_surf_loc],'FontSize',12);set(gca,'XLim',[min(alt_sar_surf)+37 max(alt_sar_surf)+45],'FontSize',12) ;set(gca,'YDir','reverse');
       
       subplot(3,3,1);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RMC);figlabels('','Doppler Beams','','Stack RMC ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,2);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RAW);figlabels('','Doppler Beams','','Stack RAW ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,4);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RMC_);figlabels('','Doppler Beams','','Stack aligned RMC ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,5);imagesc(1:chd.N_samples_sar_sar_chd,1:max_beams,fft_RAW_);figlabels('','Doppler Beams','','Stack aligned RAW ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ;
       subplot(3,3,7);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC_)/max_beams           ,'LineWidth',3);figlabels('Samples','','','Stack Multilooked RMC',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);
       hold on; subplot(3,3,7);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC_(100:(max(N_beams_stack)-100),:)/(max_beams-200))            ,'g','LineWidth',2);figlabels('Samples','','','Stack Sum RMC',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);hold off
       legend('Ambiguities','No Ambiguities')
       subplot(3,3,8);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW_)/max_beams,'LineWidth',3);figlabels('Samples','','','Stack Multilooked RAW',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);
       hold on; subplot(3,3,8);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW_(100:(max(N_beams_stack)-100),:)/(max_beams-200)),'g','LineWidth',2);figlabels('Samples','','','Stack Sum RMC',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);hold off
       legend('Ambiguities ','No Ambiguities')
       pause(0.01)
       saveas (h, ['./results/stacks/',num2str(i_surf, '%04d'),'_RMC_vs_RAW_STACKS_S3.png']);
    end

     i_surf=100;
     plot_RMC(:,:) = beams_surf_RMC(i_surf,:,:); plot_RMC_(:,:) = wfm_geo_corr_RMC(i_surf,:,:);
     plot_RAW(:,:) = beams_surf(i_surf,:,:);     plot_RAW_(:,:) = wfm_geo_corr(i_surf,:,:);
     fft_RMC(:,:)  = abs(fftshift(fft(plot_RMC.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
     fft_RAW(:,:)  = abs(fftshift(fft(plot_RAW.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
     fft_RMC_(:,:)  = abs(fftshift(fft(plot_RMC_.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');
     fft_RAW_(:,:)  = abs(fftshift(fft(plot_RAW_.',chd.N_samples_sar_sar_chd*cnf.zp_fact_range),1).');

     for i_beam_amb=0: max_beams/2

       clear plot_RMC plot_RAW  plot_RMC_ plot_RAW_ 
       plot_RMC=zeros(max_beams,chd.N_samples_sar_sar_chd); plot_RMC_=zeros(max_beams,chd.N_samples_sar_sar_chd);
       plot_RAW=zeros(max_beams,chd.N_samples_sar_sar_chd); plot_RAW_=zeros(max_beams,chd.N_samples_sar_sar_chd);

       
       subplot(2,3,1);imagesc(1:chd.N_samples_sar_sar_chd,1+i_beam_amb:max_beams-i_beam_amb,fft_RMC_(1+i_beam_amb:max_beams-i_beam_amb,:));figlabels('','Doppler Beams','','Stack RMC ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ; set(gca,'YLim',[1 max_beams],'FontSize',12);
       subplot(2,3,4);imagesc(1:chd.N_samples_sar_sar_chd,1+i_beam_amb:max_beams-i_beam_amb,fft_RAW_(1+i_beam_amb:max_beams-i_beam_amb,:));figlabels('','Doppler Beams','','Stack RAW ',12);colormap(hot);set(gca,'XLim',[x_min x_max],'FontSize',12) ; set(gca,'YLim',[1 max_beams],'FontSize',12);
       subplot(2,3,2);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC_)/max_beams           ,'LineWidth',3);figlabels('Samples','','','Stack Multilooked RMC',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);
       hold on; subplot(2,3,2);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RMC_(1+i_beam_amb:max_beams-i_beam_amb,:)/(max_beams-2*i_beam_amb))            ,'g','LineWidth',2);figlabels('Samples','','','Stack Sum RMC',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);hold off
       legend('Ambiguities','No Amb.')
       subplot(2,3,5);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW_)/max_beams/sqrt(256)*2,'LineWidth',3);figlabels('Samples','','','Stack Multilooked RAW',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);
       hold on; subplot(2,3,5);plot(1/cnf.zp_fact_range:1/cnf.zp_fact_range:chd.N_samples_sar_sar_chd,sum(fft_RAW_(1+i_beam_amb:max_beams-i_beam_amb,:)/(max_beams-2*i_beam_amb)/sqrt(256)*2),'g','LineWidth',2);figlabels('Samples','','','Stack Sum RAW',12);set(gca,'XLim',[1 x_max],'FontSize',12); set(gca,'YLim',[0 y_max],'FontSize',12);hold off
       legend('Ambiguities ','No Amb.')
       error_RMC(i_beam_amb+1) = norm(sum(fft_RMC_)/max_beams-sum(fft_RMC_(1+i_beam_amb:max_beams-i_beam_amb,:)/(max_beams-2*i_beam_amb)));
       subplot(2,3,3);plot(0:i_beam_amb,error_RMC(1:i_beam_amb+1),'LineWidth',3);figlabels('Beams filtered','','','Error RMC',12);set(gca,'XLim',[1 max_beams/2],'FontSize',12); set(gca,'YLim',[0 0.1],'FontSize',12);
       error_RAW(i_beam_amb+1) = norm(sum(fft_RMC_)/max_beams-sum(fft_RMC_(1+i_beam_amb:max_beams-i_beam_amb,:)/(max_beams-2*i_beam_amb)));
       subplot(2,3,6);plot(0:i_beam_amb,error_RAW(1:i_beam_amb+1),'LineWidth',3);figlabels('Beams filtered','','','Error RAW',12);set(gca,'XLim',[1 max_beams/2],'FontSize',12); set(gca,'YLim',[0 0.1],'FontSize',12);
       
       pause(0.01)
       saveas (h, ['./results/stacks/',num2str(i_beam_amb, '%04d'),'_RMC_vs_RAW_STACKS_S3_100.png']);
    end
    
    
end
