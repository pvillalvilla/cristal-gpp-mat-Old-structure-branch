%% FREQ ABS
figure;
subplot(3,2,1);
mesh(abs(L1BS.beams_surf(1:L1BS.N_beams_stack,:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude before alignment (freq domain)',20);
view(62,30);
subplot(3,2,2);
plot(abs(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude before alignment (freq domain)',20);
subplot(3,2,3);
mesh(abs(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples)));
view(62,30);
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Amplitude shift (radians)',20);
subplot(3,2,4);
plot(abs(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Amplitude shift Central beam (radians)',20);
subplot(3,2,5);
mesh(abs(((L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:)))));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude after alignment (time domain)',20);
view(62,30);
colormap('jet');
subplot(3,2,6);
plot(abs(((L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:))))); hold all
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Amplitude Shift of the central beam',20);
legend('Post', 'Pre');

%% FREQ PHASE
figure;
subplot(3,2,1);
mesh(angle(L1BS.beams_surf(1:L1BS.N_beams_stack,:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack phase before alignment (freq domain)',20);
view(62,30);
subplot(3,2,2);
plot(angle(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase before alignment (freq domain)',20);
subplot(3,2,3);
mesh(angle(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples)));
view(62,30);
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Phase shift (radians)',20);
subplot(3,2,4);
plot(angle(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Phase shift Central beam (radians)',20);
subplot(3,2,5);
mesh(angle(fft(fftshift(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack phase after alignment (freq domain)',20);
view(62,30);
colormap('jet');
subplot(3,2,6);
plot(angle(((L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:))))); hold all
plot(angle(((L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)))));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Phase of  the central beam',20);
legend('Post', 'Pre');


%% TIME
figure;
subplot(3,2,1);
mesh(abs(fft(fftshift(L1BS.beams_surf(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude before alignment (time domain)',20);
view(62,30);
subplot(3,2,2);
plot(abs(fft(fftshift(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude before alignment (time domain)',20);
subplot(3,2,3);
mesh(abs(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples),2),[],2)));
view(62,30);
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Amplitude shift (radians)',20);
subplot(3,2,4);
plot(abs(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Amplitude shift Central beam (radians)',20);
subplot(3,2,5);
mesh(abs(fft(fftshift(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude after alignment (time domain)',20);
view(62,30);
colormap('jet');
subplot(3,2,6);
plot(abs(fft(fftshift(L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:),2),[],2))); hold all
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude after alignment (time domain)',20);

figure;
subplot(3,2,1);
mesh(angle(fft(fftshift(L1BS.beams_surf(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack phase before alignment (time domain)',20);
view(62,30);
subplot(3,2,2);
plot(angle(fft(fftshift(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase before alignment (time domain)',20);
subplot(3,2,3);
mesh(angle(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples),2),[],2)));
view(62,30);
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Phase shift (radians)',20);
subplot(3,2,4);
plot(angle(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Phase shift Central beam (radians)',20);
subplot(3,2,5);
mesh(angle(fft(fftshift(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack phase after alignment (time domain)',20);
view(62,30);
colormap('jet');
subplot(3,2,6);
plot(angle(fft(fftshift(L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:),2),[],2))); hold all
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase after alignment (time domain)',20);

%% FFT no Shift

figure;
subplot(3,2,1);
mesh(abs(fft((L1BS.beams_surf(1:L1BS.N_beams_stack,:)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude before alignment (time domain no fftshift)',20);
view(62,30);
subplot(3,2,2);
plot(abs(fft((L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude before alignment (time domain no fftshift)',20);
subplot(3,2,3);
mesh(abs(fft((exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples)),[],2)));
view(62,30);
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Amplitude shift (time domain no fftshift)',20);
subplot(3,2,4);
plot(abs(fft((exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Amplitude shift central beam (time domain no fftshift)',20);
subplot(3,2,5);
mesh(abs(fft((L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude after alignment (time domain no fftshift)',20);
view(62,30);
colormap('jet');
subplot(3,2,6);hold off;
plot(abs(fft((L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:)),[],2))); hold all
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude after alignment (time domain no fftshift)',20);



figure;
subplot(3,2,1);
mesh(angle(fft((L1BS.beams_surf(1:L1BS.N_beams_stack,:)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack phase before alignment (time domain no fftshift)',20);
view(62,30);
subplot(3,2,2);
plot(angle(fft((L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase before alignment (time domain no fftshift)',20);
subplot(3,2,3);
mesh(angle(fft((exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples)),[],2)));
view(62,30);
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Phase shift (time domain no fftshift)',20);
subplot(3,2,4);
plot(angle(fft((exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Phase shift Central beam (time domain no fftshift)',20);
subplot(3,2,5);
mesh(angle(fft((L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:)),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack phase after alignment (time domain no fftshift)',20);
view(62,30);
colormap('jet');
subplot(3,2,6);hold off;
plot(angle(fft((L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:)),[],2))); hold all
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase after alignment (time domain no fftshift)',20);



%% Input signal FREQ
figure;
subplot(2,2,1);
imagesc(abs(L1BS.beams_surf(1:L1BS.N_beams_stack,:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude before alignment (freq domain)',20); colormap('jet'); 

subplot(2,2,2);
imagesc(angle(L1BS.beams_surf(1:L1BS.N_beams_stack,:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','Stack phase before alignment (freq domain)',20); 

subplot(2,2,3);
plot(abs(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude before alignment (freq domain)',20);

subplot(2,2,4);
plot(angle(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase before alignment (freq domain)',20);

%% Input signal Time
figure;
subplot(2,2,1);
mesh(abs(fft(fftshift(L1BS.beams_surf(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);colormap('jet'); 
figlabels('Samples','Beam index','','Stack amplitude before alignment (time domain)',20);

subplot(2,2,2);
mesh(angle(fft(fftshift(L1BS.beams_surf(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','Stack phase before alignment (time domain)',20); 

subplot(2,2,3);
plot(abs(fft(fftshift(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude before alignment (time domain)',20);

subplot(2,2,4);
plot(angle(fft(fftshift(L1BS.beams_surf(floor(L1BS.N_beams_stack/2),:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase before alignment (time domain)',20);


%% Shift FREQ
figure;
subplot(2,2,1);
imagesc(abs(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','EXP amplitude (freq domain)',20); colormap('jet'); 

subplot(2,2,2);
imagesc(angle(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','EXP phase (freq domain)',20); 

subplot(2,2,3);
plot(abs(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','EXP Central amplitude (freq domain)',20);

subplot(2,2,4);
plot(angle(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','EXP Central phase (freq domain)',20);

%% Shift time
figure;
subplot(2,2,1);
mesh(abs(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);colormap('jet'); 
figlabels('Samples','Beam index','','EXP amplitude  (time domain)',20);

subplot(2,2,2);
mesh(angle(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(1:L1BS.N_beams_stack)'*i_samples),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','EXP phase  (time domain)',20); 

subplot(2,2,3);
plot(abs(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','EXP Central amplitude (time domain)',20);

subplot(2,2,4);
plot(angle(fft(fftshift(exp(2i.*cst.pi./chd.N_samples_sar.*L1BS.shift(floor(L1BS.N_beams_stack/2))'*i_samples),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','EXP Central phase  (time domain)',20);

%% Output signal FREQ
figure;
subplot(2,2,1);
imagesc(abs(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);
figlabels('Samples','Beam index','','Stack amplitude after alignment (freq domain)',20); colormap('jet'); 

subplot(2,2,2);
imagesc(angle(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','Stack phase after alignment (freq domain)',20); 

subplot(2,2,3);
plot(abs(L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude after alignment (freq domain)',20);

subplot(2,2,4);
plot(angle(L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase after alignment (freq domain)',20);


%% Output signal Time
zp_fact_trp=1;
figure;
subplot(2,2,1);
mesh(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS.N_beams_stack,abs(fft(fftshift(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);colormap('jet'); 
figlabels('Samples','Beam index','','Stack amplitude after alignment (time domain)',20);

subplot(2,2,2);
mesh(angle(fft(fftshift(L1BS.beam_geo_corr(1:L1BS.N_beams_stack,:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','Stack phase after alignment (time domain)',20); 

subplot(2,2,3);
plot(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,abs(fft(fftshift(L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam amplitude after alignment (time domain)',20);

subplot(2,2,4);
plot(angle(fft(fftshift(L1BS.beam_geo_corr(floor(L1BS.N_beams_stack/2),:),2),[],2)));
set(gca,'XLim',[1 chd.N_samples_sar ],'FontSize',20);
figlabels('Samples','','','Central beam phase after alignment (time domain)',20);



zp_fact_trp=16;

non_centered_spectra=fftshift(squeeze(L1BS.beam_geo_corr(:,:))...
                    ,2);
    beams_zp_fft = fft([non_centered_spectra(:,1:chd.N_samples_sar/2),...
                       zeros(max(L1BS.N_beams_stack),(zp_fact_trp-1)*chd.N_samples_sar),...
                       non_centered_spectra(:,chd.N_samples_sar/2+1:chd.N_samples_sar)],...
                       chd.N_samples_sar * zp_fact_trp,2);

%% Output signal Time
figure;
subplot(2,2,1);
mesh(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,1:L1BS.N_beams_stack,abs(beams_zp_fft));
set(gca,'XLim',[1 chd.N_samples_sar],'FontSize',20);
set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20);colormap('jet'); 
figlabels('Samples','Beam index','','Stack amplitude after alignment (time domain)',20);

subplot(2,2,2);
mesh(angle(beams_zp_fft));
set(gca,'XLim',[1 chd.N_samples_sar * 16 ],'FontSize',20); set(gca,'YLim',[1 L1BS.N_beams_stack],'FontSize',20); 
figlabels('Samples','Beam index','','Stack phase after alignment (time domain)',20); 

subplot(2,2,3);
plot(0:1/zp_fact_trp: chd.N_samples_sar -1/zp_fact_trp,abs(((beams_zp_fft(floor(L1BS.N_beams_stack/2),:)))));
set(gca,'XLim',[1 chd.N_samples_sar * 16],'FontSize',20);
figlabels('Samples','','','Central beam amplitude after alignment (time domain)',20);

subplot(2,2,4);
plot(angle(((beams_zp_fft(floor(L1BS.N_beams_stack/2),:)))));
set(gca,'XLim',[1 chd.N_samples_sar * 16 ],'FontSize',20);
figlabels('Samples','','','Central beam phase after alignment (time domain)',20);