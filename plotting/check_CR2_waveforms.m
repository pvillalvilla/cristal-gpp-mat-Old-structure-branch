%check CS2 waveforms
input_files=dir('*.nc');
for i_file=1:length(input_files)
    filename = input_files(i_file).name;
    Bw=320 * 1e6;
    c_cst = 299792458;
    zp_factor_DDP=2;
    load('colormap_blues.mat');
    
    
    
    
    pwr_waveform_20_ku=ncread(filename,'pwr_waveform_20_ku');
    window_del_20_ku=double(ncread(filename,'window_del_20_ku'))*1e-12;
    ph_diff_waveform_20_ku=ncread(filename,'ph_diff_waveform_20_ku');
    coherence_waveform_20_ku=ncread(filename,'coherence_waveform_20_ku');
    lon_20_ku=ncread(filename,'lon_20_ku');
    lat_20_ku=ncread(filename,'lat_20_ku');
    alt_20_ku=ncread(filename,'alt_20_ku');
    scale_factor = double(ncread(filename,'echo_scale_factor_20_ku'));
    scale_power = double(ncread(filename,'echo_scale_pwr_20_ku'));
    
    lon_20_ku(lon_20_ku<0)=lon_20_ku(lon_20_ku<0)+360;
    
    indexes=find (lat_20_ku<-66.5 & lat_20_ku>-68.5);
    lla2kml(['.\plots\' filename(20:32)],lat_20_ku(indexes),lon_20_ku(indexes),100+zeros(length(lon_20_ku(indexes)),1),'.');
    
%     geoid_correction_DDP=geoidheight(lat_20_ku(indexes).',lon_20_ku(indexes).');
    elev_DDP = (alt_20_ku(indexes)-window_del_20_ku(indexes)*c_cst/2);
    
    % compute the number of windows to add
    max_elev =max(elev_DDP);
    min_elev = min(elev_DDP);
    
    N_windows = floor((max_elev-min_elev)/(1/Bw/zp_factor_DDP*c_cst/2)/size(pwr_waveform_20_ku(:,indexes),1))+1;
    wvf_raw_shifted_DDP = zeros(N_windows*size(pwr_waveform_20_ku(:,indexes),1),size(pwr_waveform_20_ku(:,indexes),2));
    chr_raw_shifted_DDP = zeros(N_windows*size(pwr_waveform_20_ku(:,indexes),1),size(pwr_waveform_20_ku(:,indexes),2));
    phd_raw_shifted_DDP = zeros(N_windows*size(pwr_waveform_20_ku(:,indexes),1),size(pwr_waveform_20_ku(:,indexes),2));
    [~,reference_elevation_surf]=max(elev_DDP);
    top_elev_DDP  = elev_DDP(reference_elevation_surf)+size(pwr_waveform_20_ku(:,indexes),1)/2/Bw/zp_factor_DDP*c_cst/2;
    bottom_elev_DDP    = elev_DDP(reference_elevation_surf)-(size(wvf_raw_shifted_DDP,1))/Bw/zp_factor_DDP*c_cst/2+size(pwr_waveform_20_ku(:,indexes),1)/2/Bw/zp_factor_DDP*c_cst/2;
    elev_axis_DDP = (top_elev_DDP:-1/Bw/zp_factor_DDP*c_cst/2:bottom_elev_DDP+1/Bw/zp_factor_DDP*c_cst/2);
    clear wvf_raw_shifted_DDP chr_raw_shifted_DDP
    %     wvf_raw_shifted_DDP = cat(1,0.*fp_0928_DDP_p,0.*fp_0928_DDP_p,0.*fp_0928_DDP_p);
    clear shift_mat
    
    for i_surf=1:length(window_del_20_ku(indexes))
        shift_mat(i_surf)=round((elev_DDP(reference_elevation_surf)-elev_DDP(i_surf))/c_cst*zp_factor_DDP*2*Bw); %This extra 2 is for zeropadding
        %shift_mat(i_surf)=round((elev_DDP(1)-elev_DDP(i_surf))/c_cst*2*Bw);
        %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
        wvf_raw_shifted_DDP(size(pwr_waveform_20_ku(:,indexes),1) + 1 +shift_mat(i_surf):2*size(pwr_waveform_20_ku(:,indexes),1)+shift_mat(i_surf),i_surf)= pwr_waveform_20_ku(:,indexes(i_surf)).*scale_factor(i_surf).*2.^scale_power(i_surf) ;
        chr_raw_shifted_DDP(size(pwr_waveform_20_ku(:,indexes),1) + 1 +shift_mat(i_surf):2*size(coherence_waveform_20_ku(:,indexes),1)+shift_mat(i_surf),i_surf)= coherence_waveform_20_ku(:,indexes(i_surf));
        phd_raw_shifted_DDP(size(pwr_waveform_20_ku(:,indexes),1) + 1 +shift_mat(i_surf):2*size(ph_diff_waveform_20_ku(:,indexes),1)+shift_mat(i_surf),i_surf)= ph_diff_waveform_20_ku(:,indexes(i_surf));
        
    end
    
    [maxval,maxpos]=max(wvf_raw_shifted_DDP);
    figure; imagesc(-mean(maxpos):size(wvf_raw_shifted_DDP,1)-mean(maxpos),lat_20_ku(indexes), 10.*log10(wvf_raw_shifted_DDP).');
    hcb=colorbar;
    hcb.Label.String = 'Rx Power [dB]';
    % xlim([-59.5033  -59.4054])
    ylabel('Latitude [degrees]')
    xlabel('Range bin')
    colormap(colormap_blues);
    set(gca,'YDir', 'normal');
    set(gca,'XLim',[-300 -300+1024])
    figName = ['.\plots\' filename(20:32) '_WVF'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
    
    
    figure; imagesc(-mean(maxpos):size(wvf_raw_shifted_DDP,1)-mean(maxpos),lat_20_ku(indexes), (chr_raw_shifted_DDP).');
    hcb=colorbar;
    hcb.Label.String = 'Coherence ';
    % xlim([-59.5033  -59.4054])
    ylabel('Latitude [degrees]')
    xlabel('Range bin')
    %     colormap(colormap_blues);
    set(gca,'YDir', 'normal');
    caxis([0 1])
    set(gca,'XLim',[-300 -300+1024])
    figName = ['.\plots\' filename(20:32) '_CHR'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
    
    
%     
    
    
    %unwrapping
    figure; imagesc(-mean(maxpos):size(wvf_raw_shifted_DDP,1)-mean(maxpos),lat_20_ku(indexes),180/pi*phd_raw_shifted_DDP.' );
    hcb=colorbar;
    hcb.Label.String = 'Phase difference [degres]';
    % xlim([-59.5033  -59.4054])
    ylabel('Latitude [degrees]')
    xlabel('Range bin')
    %     colormap(colormap_blues);
    set(gca,'YDir', 'normal');
    caxis([-180 180])
    set(gca,'XLim',[-300 -300+1024])
    figName = ['.\plots\' filename(20:32) '_PHD'];
    print(gcf, figName,'-dpng');
    saveas (gcf,[figName,'.png']);
close all
end