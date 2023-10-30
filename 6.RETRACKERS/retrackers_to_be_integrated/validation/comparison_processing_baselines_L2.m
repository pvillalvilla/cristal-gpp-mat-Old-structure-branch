function res=comparison_processing_baselines_L2(input_path_L2_validation_b1,b1_id,input_path_L2_validation_b2,b2_id,path_comparison_results,varargin)
time_init=tic;
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
version_matlab=version;
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<5 || nargin>(5+4*2))
    error('Wrong number of input parameters');
end
p = inputParser;
p.addParamValue('retrackers',{''},@(x)iscellstr(x)); %is a cell array with different names of different retrackers in L2 product
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('figures_visible',0);
p.addParamValue('STL_results_active',1);
p.parse(varargin{:});
retrackers=p.Results.retrackers;
sh_name_nc=p.Results.sh_name_nc;
figures_visible=p.Results.figures_visible;
STL_results_active=p.Results.STL_results_active;
clear p;

%----------------- Define linstyles for bulk comparison -------------------
linestyle_ESA='or';
linestyle_analytical_SWH_MSSfixed_b1='*b';
linestyle_analytical_SWH_MSSfixed_b2='^c';
linestyle_STL_b1='sg';
linestyle_STL_b2='+m';

linestyle_analytical_MSS_SWHfixed_b1='^g';
linestyle_analytical_MSS_SWHfixed_b2='^g';
linestyle_threshold_b1='+k';
linestyle_threshold_b2='+k';
linestyle_OCOG_b1='dc';
linestyle_OCOG_b2='dc';
fontsize_xlabel_tracks=8;
size_marker=7;
title_name_SWH_MSSfixed='analytical';
title_name_MSS_SWHfixed='analytical-fit-MSS-and-SWH-fixed';

min_rmse_SSH=0.0; %meters 
max_rmse_SSH=1.0;
min_rmse_fit_SSH=0.0;
max_rmse_fit_SSH=1.5;
min_mean_SSH=0.0; %meters 
max_mean_SSH=1.0;
step_RMSE_fit_SSH_hist2=0.025; %meters

min_rmse_SIG0=0.0; %dB 
max_rmse_SIG0=1.0;
min_rmse_fit_SIG0=0.0;
max_rmse_fit_SIG0=1.0;
min_mean_SIG0=-1.5; %dB 
max_mean_SIG0=1.5;
step_RMSE_fit_SIG0_hist2=0.005; %dB

min_rmse_SWH=0.0; %M 
max_rmse_SWH=1.0;
min_rmse_fit_SWH=0.0;
max_rmse_fit_SWH=1.0;
min_mean_SWH=-1.0; %M 
max_mean_SWH=1.0;
step_RMSE_fit_SWH_hist2=0.005; %m

min_rmse_COR=0.0;  
max_rmse_COR=0.6;
min_rmse_fit_COR=0.0;
max_rmse_fit_COR=1.0;
min_mean_COR=95.0;  
max_mean_COR=100.0;
step_RMSE_fit_COR_hist2=0.0075; 

nanclr=[1 1 1];


%------------------- Check the available retrackers to be processed ------
idx_int_analytical_SWH_MSSfixed=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL_SWH')), 1);
idx_int_analytical_MSS_SWHfixed=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL_MSS')), 1);
idx_int_thres=find(~cellfun(@isempty,strfind(retrackers,'THRESHOLD')), 1);
idx_int_ocog=find(~cellfun(@isempty,strfind(retrackers,'OCOG')), 1);


mkdir(path_comparison_results);

%% -------------COMPARISON OF TRACKS --------------------------------------
%--------------------------------------------------------------------------
%loading and reordering the data into a single array of structures
%filesBulk.inputFilesEvaluation      =   dir(fullfile(path_comparison_results,'*_L2_Evaluation.mat'));

% ---- load results validation Baseline -1 --------------------------------
load(strcat(input_path_L2_validation_b1,'L2_Bulk_validation_information.mat'));
SIGMA0_b1=SIGMA0;
SSH_b1=SSH;
SWH_b1=SWH;
COR_b1=COR;
clear SIGMA0 SSH SWH COR;
% ---- load results validation Baseline -2 --------------------------------
load(strcat(input_path_L2_validation_b2,'L2_Bulk_validation_information.mat'));
SIGMA0_b2=SIGMA0;
SSH_b2=SSH;
SWH_b2=SWH;
COR_b2=COR;
clear SIGMA0 SSH SWH COR;
min_track=1;
max_track=length(SIGMA0_b1);

%% ----------  Ploting ----------------------------------------------------
if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end
set(0,'defaultLineMarkerSize',size_marker);  % set the default line marker size
%% ---------------------------  SSH ---------------------------------------
%--------------------------------------------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison ISR w.r.t ESA & STL $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SSH_b1(:).THRESHOLD];
    res.b1.SSH.THRESHOLD.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SSH.THRESHOLD.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_threshold_b1)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).THRESHOLD];
    res.b2.SSH.THRESHOLD.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SSH.THRESHOLD.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_threshold_b2)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b2_id,']'])];
end
if ~isempty(idx_int_ocog)
    results=[SSH_b1(:).OCOG];
    res.b1.SSH.OCOG.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SSH.OCOG.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_OCOG_b1)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).OCOG];
    res.b2.SSH.OCOG.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SSH.OCOG.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_OCOG_b2)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b2_id,']'])];
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SSH.ANALYTICAL_SWH_MSSfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    %comaprison with Starlab
    if STL_results_active
        results=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b1.SSH.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b2.SSH.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH_b1(:).ANALYTICAL_MSS_SWHfixed];
    res.b1.SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nantsd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_MSS_SWHfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).ANALYTICAL_MSS_SWHfixed];
    res.b2.SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SSH.ANALYTICAL_MSS_SWHfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_MSS_SWHfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    %comparison results Starlab
    if STL_results_active
        results=[SSH_b1(:).ANALYTICAL_MSS_SWHfixed];
        res.b1.SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b1.SSH.ANALYTICAL_MSS_SWHfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SSH_b2(:).ANALYTICAL_MSS_SWHfixed];
        res.b2.SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_STL_ISR=nanmean([results.RMSE_error_L2_STL]);
        res.b2.SSH.ANALYTICAL_MSS_SWHfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    end
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{',upper(sh_name_nc),{'} '},'[m]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_rmse_SSH max_rmse_SSH])
if STL_results_active    
    title(strcat('RMSE error on',{' '},upper(sh_name_nc),': isardSAT (ISR) vs ESA & Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_SSH_ESA_STL_ISR.png'))
else
    title(strcat('RMSE error on',{' '},upper(sh_name_nc),': isardSAT (ISR) vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_SSH_ESA_ISR.png'))
end

% %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SSH_b1(:).THRESHOLD];
    res.b1.SSH.THRESHOLD.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SSH.THRESHOLD.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_threshold_b1)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).THRESHOLD];
    res.b2.SSH.THRESHOLD.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SSH.THRESHOLD.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_threshold_b2)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b2_id,']'])];
end
if ~isempty(idx_int_ocog)
    results=[SSH_b1(:).OCOG];
    res.b1.SSH.OCOG.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SSH.OCOG.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_OCOG_b1)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).OCOG];
    res.b2.SSH.OCOG.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SSH.OCOG.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_OCOG_b2)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b2_id,']'])];
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SSH.ANALYTICAL_SWH_MSSfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SSH.ANALYTICAL_SWH_MSSfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    %comaprison with Starlab
    if STL_results_active
        results=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b1.SSH.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b2.SSH.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH_b1(:).ANALYTICAL_MSS_SWHfixed];
    res.b1.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SSH.ANALYTICAL_MSS_SWHfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_MSS_SWHfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).ANALYTICAL_MSS_SWHfixed];
    res.b2.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SSH.ANALYTICAL_MSS_SWHfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_MSS_SWHfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    %comparison results Starlab
    if STL_results_active
        results=[SSH_b1(:).ANALYTICAL_MSS_SWHfixed];
        res.b1.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b1.SSH.ANALYTICAL_MSS_SWHfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SSH_b2(:).ANALYTICAL_MSS_SWHfixed];
        res.b2.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b2.SSH.ANALYTICAL_MSS_SWHfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    end
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('\epsilon_{',upper(sh_name_nc),'} [m]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_mean_SSH max_mean_SSH])
if STL_results_active    
    title(strcat('Mean error on',{' '},upper(sh_name_nc),': isardSAT (ISR) vs ESA & Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_SSH_ESA_STL_ISR.png'))
else
    title(strcat('Mean error on',{' '},upper(sh_name_nc),': isardSAT (ISR) vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_SSH_ESA_ISR.png'))
end


%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
results=[SSH_b1(:).ESA_L2];
res.ESA.SSH.mean_rmse_fit=nanmean(([results.rmse_fitting]));
res.ESA.SSH.std_rmse_fit=nanstd(([results.rmse_fitting]));
plot([results.rmse_fitting],linestyle_ESA)
legend_text={'ESA'};
hold on;
grid on;
if ~isempty(idx_int_thres)
    results=[SSH_b1(:).THRESHOLD];
    res.b1.SSH.THRESHOLD.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SSH.THRESHOLD.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_threshold_b1)
    legend_text=[legend_text,strcat({'ISR threshold '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).THRESHOLD];
    res.b2.SSH.THRESHOLD.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SSH.THRESHOLD.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_threshold_b2)
    legend_text=[legend_text,strcat({'ISR threshold '},['[',b2_id,']'])];
end
if ~isempty(idx_int_ocog)
    results=[SSH_b1(:).OCOG];
    res.b1.SSH.OCOG.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SSH.OCOG.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_OCOG_b1)
    legend_text=[legend_text,strcat({'ISR OCOG '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).OCOG];
    res.b2.SSH.OCOG.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SSH.OCOG.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_OCOG_b2)
    legend_text=[legend_text,strcat({'ISR OCOG '},['[',b2_id,']'])];
end
%comaprison with Starlab
if STL_results_active
    results=[SSH_b1(:).ANALYTICAL_STL];
    res.b1.SSH.ANALYTICAL_STL.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SSH.ANALYTICAL_STL.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_STL_b1)
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'},[' [',b1_id,']'])];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SSH.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SSH.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SSH.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SSH.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH_b1(:).ANALYTICAL_MSS_SWHfixed];
    res.b1.SSH.ANALYTICAL_MSS_SWHfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SSH.ANALYTICAL_MSS_SWHfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SSH_b2(:).ANALYTICAL_MSS_SWHfixed];
    res.b2.SSH.ANALYTICAL_MSS_SWHfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SSH.ANALYTICAL_MSS_SWHfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{',upper(sh_name_nc),'} [m]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_rmse_fit_SSH max_rmse_fit_SSH])
if STL_results_active    
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR) vs ESA & Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SSH_ESA_STL_ISR.png'))
else
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR) vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SSH_ESA_ISR.png'))
end

%%---------------------------- SCATTER PLOTS ------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    
    results_ESA=[SSH_b1(:).ESA_L2];
    results_ISR_b1=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([min([results_ISR_b1.rmse_fitting]),min([results_ISR_b2.rmse_fitting])]) max([max([results_ISR_b1.rmse_fitting]),max([results_ISR_b2.rmse_fitting])])];
    line_ref_y=[min([results_ESA.rmse_fitting]) max([results_ESA.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SSH_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SSH_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_ESA.rmse_fitting].',[results_ISR_b1.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    hist2_b2=hist3([[results_ESA.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);        
    %hist2_b2(hist2_b2==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max([max(hist2_b1),max(hist2_b2)]);
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},b1_id,{'} '},'[m]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ESA} [m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_rmse_FIT_SSH_ESA_ISR.png'))
    
    %---------------------- baseline 2 ------------------------------------
    figure;
    %surf(flipud(X),flipud(Y),flipud(hist2_b2),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]);
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b2)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},b2_id,{'} '},'[m]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ESA} [m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)',{' ['},b2_id,{'] '},'vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B2_rmse_FIT_SSH_ESA_ISR.png'))
    clear ctrs X Y edges_x edges_y hist2_b1 hist2_b2;
end
if STL_results_active
    
    results_STL=[SSH_b1(:).ANALYTICAL_STL];
    results_ISR_b1=[SSH_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SSH_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([min([results_ISR_b1.rmse_fitting]),min([results_ISR_b2.rmse_fitting])]) max([max([results_ISR_b1.rmse_fitting]),max([results_ISR_b2.rmse_fitting])])];
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SSH_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SSH_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_STL.rmse_fitting].',[results_ISR_b1.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    hist2_b2=hist3([[results_STL.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);        
    %hist2_b2(hist2_b2==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max([max(hist2_b1),max(hist2_b2)]);
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},b1_id,{'} '},'[m]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{STL} [m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_rmse_FIT_SSH_STL_ISR.png'))
    
    %---------------------- baseline 2 ------------------------------------
    figure;
    %surf(flipud(X),flipud(Y),flipud(hist2_b2),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]);
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b2)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},b2_id,{'} '},'[m]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{STL} [m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)',{' ['},b2_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B2_rmse_FIT_SSH_STL_ISR.png'))
    clear ctrs X Y edges_x edges_y hist2_b1 hist2_b2;
end
%% ---------------------------  SIGMA0 ------------------------------------
%--------------------------------------------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison ISR w.r.t ESA & STL $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SIGMA0_b1(:).THRESHOLD];
    res.b1.SIGMA0.THRESHOLD.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SIGMA0.THRESHOLD.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_threshold_b1)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).THRESHOLD];
    res.b2.SIGMA0.THRESHOLD.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SIGMA0.THRESHOLD.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_threshold_b2)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b2_id,']'])];
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0_b1(:).OCOG];
    res.b1.SIGMA0.OCOG.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SIGMA0.OCOG.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_OCOG_b1)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).OCOG];
    res.b2.SIGMA0.OCOG.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SIGMA0.OCOG.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_OCOG_b2)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b2_id,']'])];
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    %comaprison with Starlab
    if STL_results_active
        results=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0_b1(:).ANALYTICAL_MSS_SWHfixed];
    res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_MSS_SWHfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).ANALYTICAL_MSS_SWHfixed];
    res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
    res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
    plot([results.RMSE_error_L2],linestyle_analytical_MSS_SWHfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    %comparison results Starlab
    if STL_results_active
        results=[SIGMA0_b1(:).ANALYTICAL_MSS_SWHfixed];
        res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SIGMA0_b2(:).ANALYTICAL_MSS_SWHfixed];
        res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    end
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{','\sigma^0',{'} '},'[dB]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_rmse_SIG0 max_rmse_SIG0])
if STL_results_active    
    title(strcat('RMSE error on',{' '},'\sigma^0',': isardSAT (ISR) vs ESA & Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_SIG0_ESA_STL_ISR.png'))
else
    title(strcat('RMSE error on',{' '},'\sigma^0',': isardSAT (ISR) vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_SIG0_ESA_ISR.png'))
end

% %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SIGMA0_b1(:).THRESHOLD];
    res.b1.SIGMA0.THRESHOLD.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SIGMA0.THRESHOLD.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_threshold_b1)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).THRESHOLD];
    res.b2.SIGMA0.THRESHOLD.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SIGMA0.THRESHOLD.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_threshold_b2)
    legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',b2_id,']'])];
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0_b1(:).OCOG];
    res.b1.SIGMA0.OCOG.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SIGMA0.OCOG.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_OCOG_b1)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).OCOG];
    res.b2.SIGMA0.OCOG.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SIGMA0.OCOG.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_OCOG_b2)
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',b2_id,']'])];
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    %comaprison with Starlab
    if STL_results_active
        results=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0_b1(:).ANALYTICAL_MSS_SWHfixed];
    res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_MSS_SWHfixed_b1)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).ANALYTICAL_MSS_SWHfixed];
    res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],linestyle_analytical_MSS_SWHfixed_b2)
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    %comparison results Starlab
    if STL_results_active
        results=[SIGMA0_b2(:).ANALYTICAL_MSS_SWHfixed];
        res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SIGMA0_b2(:).ANALYTICAL_MSS_SWHfixed];
        res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
    end
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('\epsilon_{','\sigma^0','} [dB]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_mean_SIG0 max_mean_SIG0])
if STL_results_active    
    title(strcat('Mean error on',{' '},'\sigma^0',': isardSAT (ISR) vs ESA & Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_SIG0_ESA_STL_ISR.png'))
else
    title(strcat('Mean error on',{' '},'\sigma^0',': isardSAT (ISR) vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_SIG0_ESA_ISR.png'))
end


%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
% results=[SIGMA0_b1(:).ESA_L2];
% res.ESA.SIGMA0.mean_rmse_fit=nanmean(([results.rmse_fitting]));
% res.ESA.SIGMA0.std_rmse_fit=nanstd(([results.rmse_fitting]));
% plot([results.rmse_fitting],linestyle_ESA)
% legend_text={'ESA'};
% hold on;
% grid on;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SIGMA0_b1(:).THRESHOLD];
    res.b1.SIGMA0.THRESHOLD.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SIGMA0.THRESHOLD.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_threshold_b1)
    legend_text=[legend_text,strcat({'ISR threshold '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).THRESHOLD];
    res.b2.SIGMA0.THRESHOLD.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SIGMA0.THRESHOLD.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_threshold_b2)
    legend_text=[legend_text,strcat({'ISR threshold '},['[',b2_id,']'])];
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0_b1(:).OCOG];
    res.b1.SIGMA0.OCOG.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SIGMA0.OCOG.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_OCOG_b1)
    legend_text=[legend_text,strcat({'ISR OCOG '},['[',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).OCOG];
    res.b2.SIGMA0.OCOG.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SIGMA0.OCOG.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_OCOG_b2)
    legend_text=[legend_text,strcat({'ISR OCOG '},['[',b2_id,']'])];
end
%comaprison with Starlab
if STL_results_active
    results=[SIGMA0_b1(:).ANALYTICAL_STL];
    res.b1.SIGMA0.ANALYTICAL_STL.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SIGMA0.ANALYTICAL_STL.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_STL_b1)
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'},[' [',b1_id,']'])];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SIGMA0.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0_b1(:).ANALYTICAL_MSS_SWHfixed];
    res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_MSS_SWHfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SIGMA0_b2(:).ANALYTICAL_MSS_SWHfixed];
    res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SIGMA0.ANALYTICAL_MSS_SWHfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_MSS_SWHfixed,[' [',b2_id,']'])];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{','\sigma^0','} [dB]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_rmse_fit_SIG0 max_rmse_fit_SIG0])
if STL_results_active    
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR) vs ESA & Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SIG0_ESA_STL_ISR.png'))
else
%     title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR) vs ESA'))
%     print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SIG0_ESA_ISR.png'))
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SIG0_ISR.png'))
end

%%---------------------------- SCATTER PLOTS ------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    
    results_ESA=[SIGMA0_b1(:).ESA_L2];
    results_ISR_b1=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([min([results_ISR_b1.rmse_fitting]),min([results_ISR_b2.rmse_fitting])]) max([max([results_ISR_b1.rmse_fitting]),max([results_ISR_b2.rmse_fitting])])];
    line_ref_y=[min([results_ESA.rmse_fitting]) max([results_ESA.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_ESA.rmse_fitting].',[results_ISR_b1.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    hist2_b2=hist3([[results_ESA.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);        
    %hist2_b2(hist2_b2==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max([max(hist2_b1),max(hist2_b2)]);
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},b1_id,{'} '},'[dB]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','\sigma^0','}^{ESA} [dB]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_rmse_FIT_SIG0_ESA_ISR.png'))
    
    %---------------------- baseline 2 ------------------------------------
    figure;
    %surf(flipud(X),flipud(Y),flipud(hist2_b2),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]);
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b2)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},b2_id,{'} '},'[dB]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','\sigma^0','}^{ESA} [dB]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},b2_id,{'] '},'vs ESA'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B2_rmse_FIT_SIG0_ESA_ISR.png'))
    clear ctrs X Y edges_x edges_y hist2_b1 hist2_b2;
end
%-------------- ISR BASELINE 1 VS ISR BASELINE 2 --------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    % Baseline 1 against baseline 2

    results_ISR_b1=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([results_ISR_b2.rmse_fitting]) max([results_ISR_b2.rmse_fitting])];
    line_ref_y=[min([results_ISR_b1.rmse_fitting]) max([results_ISR_b1.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_ISR_b1.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;

    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(max(hist2_b1));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    ylabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},b1_id,{'} '},'[dB]'),'Interpreter','Tex');
    xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},b2_id,{'} '},'[dB]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs',{' ['},b2_id,{'] '}))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_vs_B2_rmse_FIT_SIG0_ISR.png'))
end

if STL_results_active
    
    results_STL=[SIGMA0_b1(:).ANALYTICAL_STL];
    results_ISR_b1=[SIGMA0_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SIGMA0_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([min([results_ISR_b1.rmse_fitting]),min([results_ISR_b2.rmse_fitting])]) max([max([results_ISR_b1.rmse_fitting]),max([results_ISR_b2.rmse_fitting])])];
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_STL.rmse_fitting].',[results_ISR_b1.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    hist2_b2=hist3([[results_STL.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);        
    %hist2_b2(hist2_b2==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max([max(hist2_b1),max(hist2_b2)]);
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},b1_id,{'} '},'[dB]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','\sigma^0','}^{STL} [dB]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_rmse_FIT_SIG0_STL_ISR.png'))
    
    %---------------------- baseline 2 ------------------------------------
    figure;
    %surf(flipud(X),flipud(Y),flipud(hist2_b2),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]);
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b2)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},b2_id,{'} '},'[dB]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','\sigma^0','}^{STL} [dB]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},b2_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B2_rmse_FIT_SIG0_STL_ISR.png'))
    clear ctrs X Y edges_x edges_y hist2_b1 hist2_b2;
end
%% ---------------------------  SWH ---------------------------------------
%--------------------------------------------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison ISR w.r.t ESA & STL $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %comaprison with Starlab
    if STL_results_active
        results=[SWH_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.SWH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b1.SWH.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SWH_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.SWH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b2.SWH.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
        
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel(strcat('RMSE_{','SWH',{'} '},'[m]'),'Interpreter','Tex');
        xlabel(xlabel_str)
        axis([min_track max_track min_rmse_SWH max_rmse_SWH])
%         if STL_results_active
            title(strcat('RMSE error on',{' '},'SWH',': isardSAT (ISR) vs Starlab (STL)'))
            print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_SWH_STL_ISR.png'))
%         else
%             title(strcat('RMSE error on',{' '},'SWH',': isardSAT (ISR) vs ESA'))
%             print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_SWH_ESA_ISR.png'))
%         end
    end
end



% %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %comaprison with Starlab
    if STL_results_active
        results=[SWH_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b1.SWH.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[SWH_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b2.SWH.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel(strcat('\epsilon_{','SWH','} [m]'),'Interpreter','Tex');
        xlabel(xlabel_str)
        axis([min_track max_track min_mean_SWH max_mean_SWH])
%         if STL_results_active
            title(strcat('Mean error on',{' '},'SWH',': isardSAT (ISR) vs Starlab (STL)'))
             print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_SWH_STL_ISR.png'))
%         else
%             title(strcat('Mean error on',{' '},'SWH',': isardSAT (ISR) vs ESA'))
%             print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_SWH_ESA_ISR.png'))
%         end
    end
end




%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%comaprison with Starlab
if STL_results_active
    results=[SWH_b1(:).ANALYTICAL_STL];
    res.b1.SWH.ANALYTICAL_STL.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SWH.ANALYTICAL_STL.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_STL_b1)
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'},[' [',b1_id,']'])];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SWH_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.SWH.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.SWH.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[SWH_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.SWH.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.SWH.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{','SWH','} [m]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_rmse_fit_SWH max_rmse_fit_SWH])
if STL_results_active    
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR) vs Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SWH_STL_ISR.png'))
else
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_SWH_ISR.png'))
end

%%---------------------------- SCATTER PLOTS ------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    % Baseline 1 against baseline 2

    results_ISR_b1=[SWH_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SWH_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([results_ISR_b2.rmse_fitting]) max([results_ISR_b2.rmse_fitting])];
    line_ref_y=[min([results_ISR_b1.rmse_fitting]) max([results_ISR_b1.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SWH_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SWH_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_ISR_b1.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;

    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(max(hist2_b1));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    ylabel(strcat('RMSE_{','SWH','}^{ISR',{','},b1_id,{'} '},'[m]'),'Interpreter','Tex');
    xlabel(strcat('RMSE_{','SWH','}^{ISR',{','},b2_id,{'} '},'[m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs',{' ['},b2_id,{'] '}))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_vs_B2_rmse_FIT_SWH_ISR.png'))
end
if STL_results_active
    
    results_STL=[SWH_b1(:).ANALYTICAL_STL];
    results_ISR_b1=[SWH_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SWH_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([min([results_ISR_b1.rmse_fitting]),min([results_ISR_b2.rmse_fitting])]) max([max([results_ISR_b1.rmse_fitting]),max([results_ISR_b2.rmse_fitting])])];
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_STL.rmse_fitting].',[results_ISR_b1.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    hist2_b2=hist3([[results_STL.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);        
    %hist2_b2(hist2_b2==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max([max(hist2_b1),max(hist2_b2)]);
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','SWH','}^{ISR',{','},b1_id,{'} '},'[m]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','SWH','}^{STL} [m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_rmse_FIT_SWH_STL_ISR.png'))
    
    %---------------------- baseline 2 ------------------------------------
    figure;
    %surf(flipud(X),flipud(Y),flipud(hist2_b2),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]);
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b2)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    caxis([min_image-dmap,max_image]); 
    %hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','SWH','}^{ISR',{','},b2_id,{'} '},'[m]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','SWH','}^{STL} [m]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},b2_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B2_rmse_FIT_SWH_STL_ISR.png'))
    clear ctrs X Y edges_x edges_y hist2_b1 hist2_b2;
end

%% ----------------------  PEARSON CORR COEF. -----------------------------
%--------------------------------------------------------------------------

%$$$$$$$$$$$$$$$$$ Comparison ISR w.r.t ESA & STL $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %comaprison with Starlab
    if STL_results_active
        results=[COR_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.COR.ANALYTICAL_COR_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b1.COR.ANALYTICAL_COR_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_COR_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[COR_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.COR.ANALYTICAL_COR_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
        res.b2.COR.ANALYTICAL_COR_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
        plot([results.RMSE_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
        
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel(strcat('RMSE_{','gof',{'} '},'[%]'),'Interpreter','Tex');
        xlabel(xlabel_str)
        axis([min_track max_track min_rmse_COR max_rmse_COR])
%         if STL_results_active
            title(strcat('RMSE error on',{' '},'gof (goodness of fitting)',': isardSAT (ISR) vs ESA & Starlab (STL)'))
            print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_COR_ESA_STL_ISR.png'))
%         else
%             title(strcat('RMSE error on',{' '},'gof (goodness of fitting)',': isardSAT (ISR) vs ESA'))
%             print('-dpng',strcat(path_comparison_results,'Baseline_comparison_RMSE_COR_ESA_ISR.png'))
%         end
    end
end



% %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %comaprison with Starlab
    if STL_results_active
        results=[COR_b1(:).ANALYTICAL_SWH_MSSfixed];
        res.b1.COR.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b1.COR.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b1)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
        hold on;
        grid on;
        results=[COR_b2(:).ANALYTICAL_SWH_MSSfixed];
        res.b2.COR.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.b2.COR.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],linestyle_STL_b2)
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel(strcat('\epsilon_{','gof','} [%]'),'Interpreter','Tex');
        xlabel(xlabel_str)
        axis([min_track max_track min_mean_COR max_mean_COR])
%         if STL_results_active
            title(strcat('Mean error on',{' '},'gof (goodness of fitting)',': isardSAT (ISR) vs ESA & Starlab (STL)'))
            print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_COR_ESA_STL_ISR.png'))
%         else
%             title(strcat('Mean error on',{' '},'gof (goodness of fitting)',': isardSAT (ISR) vs ESA'))
%             print('-dpng',strcat(path_comparison_results,'Baseline_comparison_mean_err_COR_ESA_ISR.png'))
%         end
    end
end




%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%comaprison with Starlab
if STL_results_active
    results=[COR_b1(:).ANALYTICAL_STL];
    res.b1.COR.ANALYTICAL_STL.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.COR.ANALYTICAL_STL.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_STL_b1)
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'},[' [',b1_id,']'])];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[COR_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.COR.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b1.COR.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[COR_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.COR.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.b2.COR.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{','gof','} [%]'),'Interpreter','Tex');
xlabel(xlabel_str)
axis([min_track max_track min_rmse_fit_COR max_rmse_fit_COR])
if STL_results_active    
    title(strcat('RMSE error on fitted',{' '},'gof (goodness of fit)',': isardSAT (ISR) vs Starlab (STL)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_COR_STL_ISR.png'))
else
    title(strcat('RMSE error on fitted',{' '},'gof (goodness of fit)',': isardSAT (ISR)'))
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_rmse_FIT_COR_ISR.png'))
end

%%---------------------------- SCATTER PLOTS ------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %baseline 1 vs baseline 2 isardSAT
    results_ISR_b1=[COR_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[COR_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([results_ISR_b2.rmse_fitting]) max([results_ISR_b2.rmse_fitting])];
    line_ref_y=[min([results_ISR_b1.rmse_fitting]) max([results_ISR_b1.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_COR_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_COR_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_ISR_b1.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(max(hist2_b1));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    ylabel(strcat('RMSE_{','gof','}^{ISR',{','},b1_id,{'} '},'[%]'),'Interpreter','Tex');
    xlabel(strcat('RMSE_{','gof','}^{ISR',{','},b2_id,{'} '},'[%]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'gof',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs',{' ['},b2_id,{'] '}))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_vs_B2_rmse_FIT_COR_ISR.png'))
    
end
if STL_results_active
    
    results_STL=[SWH_b1(:).ANALYTICAL_STL];
    results_ISR_b1=[SWH_b1(:).ANALYTICAL_SWH_MSSfixed];    
    results_ISR_b2=[SWH_b2(:).ANALYTICAL_SWH_MSSfixed];
    line_ref_x=[min([min([results_ISR_b1.rmse_fitting]),min([results_ISR_b2.rmse_fitting])]) max([max([results_ISR_b1.rmse_fitting]),max([results_ISR_b2.rmse_fitting])])];
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    hist2_b1=hist3([[results_STL.rmse_fitting].',[results_ISR_b1.rmse_fitting].'],'Edges',ctrs);
    %hist2_b1(hist2_b1==0)=NaN;
    hist2_b2=hist3([[results_STL.rmse_fitting].',[results_ISR_b2.rmse_fitting].'],'Edges',ctrs);        
    %hist2_b2(hist2_b2==0)=NaN;
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max([max(hist2_b1),max(hist2_b2)]);
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    figure;    
    %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]); 
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b1)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    %caxis([min_image-dmap,max_image]); 
    hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;    
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','gof','}^{ISR',{','},b1_id,{'} '},'[%]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','gof','}^{STL} [%]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},b1_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B1_rmse_FIT_SWH_STL_ISR.png'))
    
    %---------------------- baseline 2 ------------------------------------
    figure;
    %surf(flipud(X),flipud(Y),flipud(hist2_b2),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
    %caxis([min_image, max_image]);
    imagesc(edges_x,fliplr(edges_y),flipud(hist2_b2)); colormap([nanclr; jet]); 
    set(gca,'YDir','normal')
    caxis([min_image-dmap,max_image]); hcb=colorbar; ylim(hcb,[min_image max_image]);
    hold on;
    plot(line_ref_x,line_ref_y,'--k');
    xlabel(strcat('RMSE_{','gof','}^{ISR',{','},b2_id,{'} '},'[%]'),'Interpreter','Tex');
    ylabel(strcat('RMSE_{','gof','}^{STL} [%]'),'Interpreter','Tex');
    title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},b2_id,{'] '},'vs STL SAMOSA-3 (S-3)'))
    print('-dpng',strcat(path_comparison_results,'Scatter_plot_B2_rmse_FIT_SWH_STL_ISR.png'))
    clear ctrs X Y edges_x edges_y hist2_b1 hist2_b2;
end


%$$$$$$$$$$$$$$$$$$$$$$$$$ MEAN VALUE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure;
legend_text={''};
if STL_results_active
    results=[COR_b1(:).ANALYTICAL_STL];
    res.b1.COR.ANALYTICAL_STL.mean_mean_value=nanmean(([results.mean]));
    res.b1.COR.ANALYTICAL_STL.std_mean_value=nanstd(([results.mean]));
    plot([results.rmse_fitting],linestyle_STL_b1)
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'},[' [',b1_id,']'])];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[COR_b1(:).ANALYTICAL_SWH_MSSfixed];
    res.b1.COR.ANALYTICAL_SWH_MSSfixed.mean_mean_value=nanmean(([results.mean]));
    res.b1.COR.ANALYTICAL_SWH_MSSfixed.std_mean_value=nanstd(([results.mean]));
    plot([results.mean],linestyle_analytical_SWH_MSSfixed_b1)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b1_id,']'])];
    hold on;
    grid on;
    results=[COR_b2(:).ANALYTICAL_SWH_MSSfixed];
    res.b2.COR.ANALYTICAL_SWH_MSSfixed.mean_mean_value=nanmean(([results.mean]));
    res.b2.COR.ANALYTICAL_SWH_MSSfixed.std_mean_value=nanstd(([results.mean]));
    plot([results.mean],linestyle_analytical_SWH_MSSfixed_b2)
    legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,[' [',b2_id,']'])];
end
axis([min_track max_track min_mean_COR max_mean_COR])
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','southeast');
ylabel('\rho_{pearson} [%]','Interpreter','Tex');
xlabel(xlabel_str);
if STL_results_active
    title('Mean value pearson corr. coeff. \rho_{pearson}: isardSAT (ISR) & Starlab (STL)','Interpreter','Tex')    
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_Mean_ISR_COR.png'))
else
    title('Mean value pearson corr. coeff. \rho_{pearson}: isardSAT (ISR)','Interpreter','Tex')    
    print('-dpng',strcat(path_comparison_results,'Baseline_comparison_Mean_ISR_STL_COR.png'))
end



%close all

time_end=toc(time_init);
minutes_processing = floor(time_end/60);
secs_processing = time_end - minutes_processing*60;
disp(['Validation/processing time: ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);    
end

