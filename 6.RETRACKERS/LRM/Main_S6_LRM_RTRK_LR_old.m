% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Created by isardSAT S.L.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% This is the main file to run isardSAT's LRM retracker for J-CS GPP LRM
% L1b
%
% -------------------------------------------------------------------------
% 
% Author:               Cristina Martin-Puig / isardSAT
%
% Reviewer:             Cristina Martin-Puig / isardSAT
%
% Last revision:        Cristina Martin-Puig / isardSAT 01/08/2014
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% This software is built under the Jason-CS Ph3 contract 
% -------------------------------------------------------------------------

% close all
% clear
% clc

% -------------------------------------------------------------------------
% CONFIGURATION
% -------------------------------------------------------------------------

% - START EDITING ---------------------------------------------------------

% Limit of fitting routine 
% -------------------------------------------------------------------------
% We've seen that Laiba's model is fitting nicely the S6/Jason-CS waveforms
% except for the tail. Therefore, we have introduced in the configuration 
% parameters that allow for defining up to which range bin we do want to
% do the fitting

clear
% clc

%% Editing

% Thermal noise or noise floor estimation window configuration
cnf_p.TNini         =   10; % window positioning initial range bin for thermal noise estimation
cnf_p.TNlast        =   20; % window positioning last range bin for thermal noise estimation


cnf_p.ini_Epoch         =   50;   % Epoch initial values --> this can be estimated from the max peak detection or leave it like this
cnf_p.ini_Hs            =   1;  % SWH initial value
cnf_p.ini_Pu            =   1;    % Initial amplitude
% note we assume xi = misspointing = to zero , otherwise we would have to
% input it as additional model parameter
cnf_p.plot_fit          =   1; % 1 - yes, 0 -no
cnf_p.plot_out          =   1; % 1 - yes, 0 -no


%% Main path
cnf_p.MainPath          =   'I:\Work\Sentinel-6\Data\Simulations\Pre-QR2\outputs\LR\';
cnf_p.Filename          =   'l1b_lr_product_lat40.nc';



%% Characterisation

nf_p.semi_major_axis    =   6378137;            % Equatorial radius [m] ['float32']
nf_p.semi_minor_axis    =   6356752.3142;
nf_p.c                  =   299792458;
nf_p.beamwidth          =   (1.35*pi/180);
nf_p.BWClock            =   395e6;

% cnf_p.ZP                =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'range_oversampling_factor_ku'));  % specify if Zero Padding


nf_p.rt                 =   1/(nf_p.BWClock);
nf_p.Ns                 =   256;
cnf_p.ZP_sample         =   nf_p.Ns/2 + 1; % Mid window location to compensate Epoch and range
nf_p.pi                 =   3.141592653589790;  % pi_cst from the CST file


%% DEFINITION OF STRUCTURE
% -----------------------------------------------------------------
fit_params_ini          =   [cnf_p.ini_Epoch cnf_p.ini_Hs cnf_p.ini_Pu];



%% DATA STRUCTURE
% -----------------------------------------------------------------
data.LAT        =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'latitude_20_ku'));
data.LON        =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'longitude_20_ku'));
data.H      	=   double(ncread([cnf_p.MainPath,cnf_p.Filename],'com_altitude_20_ku'));
data.WD         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'tracker_range_calibrated_20_ku'));

data.Re         =   sqrt(nf_p.semi_major_axis^2*cos(data.LAT).^2+nf_p.semi_minor_axis^2*sin(data.LAT).^2);
data.pitch      =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'off_nadir_pitch_angle_pf_20_ku'));
data.xi         =   data.pitch*pi/180; %mispointing
data.alpha_c    =   1 + (data.H)./data.Re;
data.h          =   data.H.*data.alpha_c;


wfm_lrm         =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'lr_power_waveform_20_ku'));
wfm_sc_ft       =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'waveform_scale_factor_20_ku'));
sigma0_sc_ft    =   double(ncread([cnf_p.MainPath,cnf_p.Filename],'sigma0_scaling_factor_20_ku'));


%t           =   1:cnf_p.Nfit_lim;
options     =   optimset('Algorithm', 'levenberg-marquardt','Display','off');

out.Epoch       =   zeros(1,length(data.LAT));
out.Hs          =   zeros(1,length(data.LAT));
out.Pu          =   zeros(1,length(data.LAT));
out.xi          =   zeros(1,length(data.LAT));
out.flag        =   zeros(1,length(data.LAT));
out.cor         =   zeros(1,length(data.LAT));
out.SSH         =   zeros(1,length(data.LAT));
out.sigma0      =   zeros(1,length(data.LAT));
wfm_lrm_norm    =   zeros(nf_p.Ns,length(data.LAT));


for m = 1:length(data.LAT)
    
    wfm_lrm_norm(:,m) = wfm_lrm(:,m).'/max(wfm_lrm(:,m));
    
%     disp(['Fitting wfm # = ',num2str(m)]);
    nf_p.xi         =   data.xi(m);
    nf_p.h          =   data.h(m);
%     index           =   find(wfm_lrm(:,m) >= 1.3*mean(wfm_lrm(1:10,m)));
    
    noise_th(m)     =   1.4*mean(wfm_lrm_norm(cnf_p.TNini:cnf_p.TNlast,m));
    le_index        =   find(wfm_lrm_norm(:,m) >= noise_th(m));
    
    if (~isempty(find(diff(le_index)>1,1)))
        le_index = le_index(find(diff(le_index)>1)+1:length(le_index)); %to make sure there are not some higher samples at the beginning
                                                            %of the waveform and we have an index vector like [1,2,104,105,...]
    end
    
%     plot(wfm_lrm(:,m),'LineWidth',2)
%     hold on
%     plot(index,wfm_lrm(index,m),'r','LineWidth',2)
%     if strcmp(MODE,'RMC') == 0
        le_index           =   le_index(1);
%     end
%     plot(1:index+cnf_p.Nfit_width*cnf_p.ZP,wfm_lrm(1:index+cnf_p.Nfit_width*cnf_p.ZP,m),'g','LineWidth',2)
%     hold off

    
    
%%  Waveform fitting
    
%     index = le_index - 20;
    index = le_index;
    wfm_init = 1;
    wfm_last = nf_p.Ns;
    
    cnf_p.Nfit_width   =   wfm_last - wfm_init + 1;
    
    fit_wfm             =   wfm_lrm_norm(wfm_init:wfm_last,m).';
    
    
    nf_p.TN     =   noise_th(m)/1.4;
    
    %     nf_p.TN         =   mean(fit_wfm(cnf_p.TNini*cnf_p.ZP:cnf_p.TNlast*cnf_p.ZP));
    
    
    mpfun           =   @(fit_params,t)sl_lrm_wave_gen_LA(t, fit_params, nf_p);
    
    [fit_params,~,res,flagfit]     =   lsqcurvefit(mpfun,fit_params_ini,(1:length(fit_wfm)),fit_wfm,[],[],options);          

    fit_params_ini  = fit_params;
    out.Epoch(m)    = fit_params(1);
%     out.Epoch(m)    =   fit_params(1)*395/320;
    out.Hs(m)       =   fit_params(2);
%     out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(1:index+cnf_p.Nfit_width*cnf_p.ZP,m))*wfm_sc_ft(m)); %when starting at 1
    out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(wfm_init:wfm_last,m))*wfm_sc_ft(m)); %when starting at 1, ending at cnf_p.Nfit_widht
%     out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(index:cnf_p.Nfit_width*cnf_p.ZP,m))*wfm_sc_ft(m)); %when starting at index
%     out.Pu(m)       =   10*log10(fit_params(3)*max(wfm_lrm(1:index+cnf_p.Nfit_width*cnf_p.ZP,m))); %wfm_lrm, in this case, has the scaling factor already applied.
    out.sigma0(m)   =   out.Pu(m)+sigma0_sc_ft(m);
    
    
    
    % xi(m)       =   fit_params(4);
    out.flag(m)     =   flagfit;
    sl_wfm          =   sl_lrm_wave_gen_LA((1:length(fit_wfm)),fit_params,nf_p);
    correlation_fit =   corrcoef(fit_wfm,sl_wfm);
    out.COR(m)      =   correlation_fit(1,2);
    
    if cnf_p.plot_fit
        plot(fit_wfm,'k.-','MarkerSize',8)
        hold on
%         plot((1:length(fit_wfm))*320/395,sl_wfm, 'r.-', 'MarkerSize',8)
        plot(sl_wfm, 'g.-','LineWidth',2,'MarkerSize',10)
        grid on
        xlabel('Samples','FontName','Helvetica','FontSize',16);
        ylabel('Normalised power','FontName','Helvetica','FontSize',16);
        title(['wfm #', num2str(m)],'FontName','Helvetica','FontSize',11,'FontWeight','bold');
        annotation('textbox', [0.7 0.15 0.2 0.2],...
                            'String',{['Epoch = ', num2str(out.Epoch(m),4), ' [r.b]'],['Hs = ' ,num2str(abs(out.Hs(m)),4), ' [m]'],...
                            ['Pu = ', num2str(out.Pu(m),4), ' [ - ]'], ['r = ', num2str(out.COR(m)*100,5), '[%]'], ...
                            ['LSFlag = ', num2str(out.flag(m))] },...
                            'FontSize',20,...
                            'FontName','Helvetica',...
                            'BackgroundColor',[1 1 1]);
        axis([1 length(fit_wfm) 0 1])
%         print('-deps',['./results/',TEST,'/LR_from_HR/retracker_fits/',MODE,'/fit_S6GPP_LRM_UPTON_',num2str(cnf_p.Nfit_width),'_wfmCOUNT_',num2str(m),'.eps']);
%         saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_fits/',MODE,'/fit_S6GPP_LRM_UPTON_',num2str(cnf_p.Nfit_lim),'_wfmCOUNT_',num2str(m),'.jpg']);
        hold off
    end

end


%% PLOTS SSH,Hs,Pu,Sigma0
out.retracking_cor  = (out.Epoch'-cnf_p.ZP_sample)*nf_p.c*0.5/nf_p.BWClock;


% if strcmp(MODE,'LRM') == 1
%     out.SSH     = 1E-4 * (double(L1B_LRM.data.com_altitude_ku)-double(L1B_LRM.data.altimeter_range_calibrated_ku)) - (out.Epoch'-cnf_p.ZP_sample)*nf_p.c*0.5/nf_p.BWClock;
%     %No need to add the 1.3e6 to the variables because it's been removed from both
%     
% %     out.range   = double(L1B_LRM.data.altimeter_range_calibrated_ku)*1e-4         - (out.Epoch.'-cnf_p.ZP_sample)*nf_p.c*0.5/nf_p.BWClock ;
%     out.range   = double(L1B_LRM.data.altimeter_range_calibrated_ku)*1e-4 + 1.3e6 - (out.Epoch.'-cnf_p.ZP_sample)*nf_p.c*0.5/nf_p.BWClock ;
% else
    out.SSH = data.H - data.WD - out.retracking_cor;
%     out.SSH = (com_altitude_ku-altimeter_range_calibrated_ku - (out.Epoch-cnf_p.ZP_sample+index-1)/cnf_p.ZP*nf_p.c*0.5/nf_p.BWClock).';
% end


disp('SSH (m)'); disp(mean(out.SSH))
disp('SWH (m)'); disp(mean(abs(out.Hs)))
disp('Sigma0 (dB)'); disp(mean(out.sigma0))
disp('STD SSH (cm)'); disp(std(out.SSH)*100)
disp('STD SWH (cm)'); disp(std(abs(out.Hs))*100)
disp('STD Sigma0 (dB)'); disp(std(out.sigma0))
disp('Corr coef (%)'); disp(mean(out.COR(m)*100))

if cnf_p.plot_out
    figure;
    plot(data.LAT,out.SSH,'.-','Markersize',8);
    grid on
    title(['SSH [m] with mean value = ',num2str(mean(out.SSH))])
    ylabel('[m]')
    xlabel('Latitude [deg]')
%     saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,'/OUTPUT_L2A_S6GPP_MODE_LRM_SSH_',num2str(cnf_p.Nfit_width),'.fig']);
%     saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,'/OUTPUT_L2A_S6GPP_MODE_LRM_SSH_',num2str(cnf_p.Nfit_width),'.jpg']);
    
    figure;
    plot(data.LAT,abs(out.Hs),'.-','Markersize',8);
    grid on
    title(['Hs [m] with mean value = ',num2str(mean(abs(out.Hs)))])
    ylabel('[m]')
    xlabel('Latitude [deg]')
%     saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,'/OUTPUT_L2A_S6GPP_MODE_LRM_Hs_',num2str(cnf_p.Nfit_width),'.fig']);
%     saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,'/OUTPUT_L2A_S6GPP_MODE_LRM_Hs_',num2str(cnf_p.Nfit_width),'.jpg']);
    
    figure
    subplot(3,1,1)
    plot(out.Pu,'.-','Markersize',8);
    grid on
    title(['Pu [dB] with mean value = ',num2str(mean(out.Pu))])
    ylabel('[dB]')
    xlabel('Latitude [deg]')
    
    subplot(3,1,2)
    plot(sigma0_sc_ft,'.-','Markersize',8);
    grid on
    title(['sigma0-scaling-factor [dB] with mean value = ',num2str(mean(sigma0_sc_ft))])
    ylabel('[dB]')
    xlabel('Latitude [deg]')
    
    subplot(3,1,3)
    plot(out.sigma0,'.-','Markersize',8);
    grid on
    title(['sigma0 [dB] with mean value = ',num2str(mean(out.sigma0))])
    ylabel('[dB]')
    xlabel('Latitude [deg]')
%     saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,'/OUTPUT_L2A_S6GPP_MODE_LRM_Pu_Sigma0_',num2str(cnf_p.Nfit_width),'.fig']);
%     saveas(gcf,['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,'/OUTPUT_L2A_S6GPP_MODE_LRM_Pu_Sigma0_',num2str(cnf_p.Nfit_width),'.jpg']);
    
end
% save(['./results/',TEST,'/LR_from_HR/retracker_outputs/',MODE,...
%     '/OUTPUT_L2A_S6GPP_MODE_LRM_MOVINGWINDOW_',num2str(cnf_p.Nfit_width),...
%     '_ZP_',num2str(cnf_p.ZP),'.mat'],'out')


