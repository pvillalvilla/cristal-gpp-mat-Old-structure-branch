function res=comparison_baselines_L2(input_path_L2_validation_baselines,baselines_id,path_comparison_results,varargin)
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
linestyle_analytical_SWH_MSSfixed_bs={'*b','^c','+m'};
linestyle_STL_bs={'sg','+m'};
linestyle_analytical_MSS_SWHfixed_bs={'^g','^g'};
linestyle_threshold_bs={'+k','+k'};
linestyle_OCOG_bs={'dc','dc'};
fontsize_xlabel_tracks=8;
%size_marker=7;
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
n_baselines=length(input_path_L2_validation_baselines);
for i_base=1:n_baselines
    % ---- load results validation Baseline -1 --------------------------------
    load(strcat(char(input_path_L2_validation_baselines(i_base)),'L2_Bulk_validation_information.mat'));
    bs(i_base).SIGMA0=SIGMA0;
    bs(i_base).SSH=SSH;
    bs(i_base).SWH=SWH;
    bs(i_base).COR=COR;
    clear SIGMA0 SSH SWH COR;
end
min_track=1;
max_track=length([bs(i_base).SIGMA0]);

%% ----------  Ploting ----------------------------------------------------
if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end
%set(0,'defaultLineMarkerSize',size_marker);  % set the default line marker size
%% ---------------------------  SSH ---------------------------------------
%--------------------------------------------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison ISR w.r.t ESA & STL $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.THRESHOLD];
        res.bs(i_base).SSH.THRESHOLD.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SSH.THRESHOLD.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_threshold_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_ocog)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH(:).OCOG];
        res.bs(i_base).SSH.OCOG.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SSH.OCOG.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_OCOG_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
        %comaprison with Starlab
        if STL_results_active
            results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
            res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
            plot([results.RMSE_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
    end
    
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed];
        res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.std_RMSE_ESA_ISR=nantsd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_analytical_MSS_SWHfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
        %comparison results Starlab
        if STL_results_active
            results=[bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed];
            res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
            res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
            plot([results.RMSE_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.THRESHOLD];
        res.bs(i_base).SSH.THRESHOLD.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
        res.bs(i_base).SSH.THRESHOLD.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
        plot([results.mean_error_L2],CHAR(linestyle_threshold_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_ocog)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.OCOG];
        res.bs(i_base).SSH.OCOG.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
        res.bs(i_base).SSH.OCOG.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
        plot([results.mean_error_L2],char(linestyle_OCOG_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
        res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
        plot([results.mean_error_L2],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
        %comaprison with Starlab
        if STL_results_active
            results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
            res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
            plot([results.mean_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed];
        res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
        res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
        plot([results.mean_error_L2],char(linestyle_analytical_MSS_SWHfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
        %comparison results Starlab
        if STL_results_active
            results=[bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed];
            res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
            res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
            plot([results.mean_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
results=[bs(1).SSH.ESA_L2];
res.ESA.SSH.mean_rmse_fit=nanmean(([results.rmse_fitting]));
res.ESA.SSH.std_rmse_fit=nanstd(([results.rmse_fitting]));
plot([results.rmse_fitting],linestyle_ESA)
legend_text={'ESA'};
hold on;
grid on;
if ~isempty(idx_int_thres)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.THRESHOLD];
        res.bs(i_base).SSH.THRESHOLD.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SSH.THRESHOLD.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_threshold_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR threshold '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_ocog)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.OCOG];
        res.bs(i_base).SSH.OCOG.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SSH.OCOG.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_OCOG_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR OCOG '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
%comaprison with Starlab
if STL_results_active
    results=[bs(1).ANALYTICAL_STL];
    res.STL.SSH.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.STL.SSH.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],CHAR(linestyle_STL_bs(1)))
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'})];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed];
        res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SSH.ANALYTICAL_MSS_SWHfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_analytical_MSS_SWHfixed_bs(i_base)));
        legend_text=[legend_text,strcat({'ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
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
    
    results_ESA=[bs(1).SSH.ESA_L2];
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end
    line_ref_y=[min([results_ESA.rmse_fitting]) max([results_ESA.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SSH_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SSH_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines
        hist2_bs(:,:,i_base)=hist3([[results_ESA.rmse_fitting].',[ ISR_bs(i_base).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
    end
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    for i_base=1:n_baselines
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},baselines_id(i_base),{'} '},'[m]'),'Interpreter','Tex');
        ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ESA} [m]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)',{' ['},baselines_id(i_base),{'] '},'vs ESA'))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b',num2str(i_base),'_rmse_FIT_SSH_ESA_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs ISR_bs;
end
%-------------- ISR BASELINE 1 VS ISR BASELINE OTHERS --------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    % Baseline 1 against baseline 2
    
    %All against baseline-1
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];  
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end       
    
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines-1
        hist2_bs(:,:,i_base)=hist3([[ISR_bs(1).results.rmse_fitting].',[ISR_bs(i_base+1).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
    end
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    for i_base=1:n_baselines-1
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},baselines_id(i_base+1),{'} '},'[m]'),'Interpreter','Tex');
        xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},baselines_id(1),{'} '},'[dB]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},baselines_id(1),{'] '},'vs',{' ['},baselines_id(i_base+1),{'] '}))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b1_vs_b',num2str(i_base+1),'_rmse_FIT_SSH_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs ISR_bs;
end
if STL_results_active    
    results_STL=[bs(i_base).SSH.ANALYTICAL_STL];
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SSH.ANALYTICAL_SWH_MSSfixed];
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SSH_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SSH_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines
         hist2_bs(:,:,i_base)=hist3([[results_STL.rmse_fitting].',[ISR_bs(i_base).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
    end
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    for i_base=1:n_baselines
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        xlabel(strcat('RMSE_{',upper(sh_name_nc),'}^{ISR',{','},baselines_id(i_base),{'} '},'[m]'),'Interpreter','Tex');
        ylabel(strcat('RMSE_{',upper(sh_name_nc),'}^{STL} [m]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)',{' ['},baselines_id(i_base),{'] '},'vs STL SAMOSA-3 (S-3)'))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b',num2str(i_base),'_rmse_FIT_SSH_STL_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
end
%% ---------------------------  SIGMA0 ------------------------------------
%--------------------------------------------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison ISR w.r.t ESA & STL $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.THRESHOLD];
        res.bs(i_base).SIGMA0.THRESHOLD.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SIGMA0.THRESHOLD.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_threshold_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_ocog)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.OCOG];
        res.bs(i_base).SIGMA0.OCOG.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SIGMA0.OCOG.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_OCOG_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
        %comaprison with Starlab
        if STL_results_active
            results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
            res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
            plot([results.RMSE_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed];
        res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_RMSE_ESA_ISR=nanmean(([results.RMSE_error_L2]));
        res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.std_RMSE_ESA_ISR=nanstd(([results.RMSE_error_L2]));
        plot([results.RMSE_error_L2],char(linestyle_analytical_MSS_SWHfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
        %comparison results Starlab
        if STL_results_active
            results=[bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed];
            res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
            res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
            plot([results.RMSE_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.THRESHOLD];
        res.bs(i_base).SIGMA0.THRESHOLD.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
        res.bs(i_base).SIGMA0.THRESHOLD.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
        plot([results.mean_error_L2],char(linestyle_threshold_bs(i_base)))
        legend_text=[legend_text,strcat({'ESA - ISR threshold '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_ocog)
    for i_base=1:n_baselines
    results=[bs(i_base).SIGMA0.OCOG];
    res.bs(i_base).SIGMA0.OCOG.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.bs(i_base).SIGMA0.OCOG.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],char(linestyle_OCOG_bs(i_base)))
    legend_text=[legend_text,strcat({'ESA - ISR OCOG '},['[',baselines_id(i_base),']'])];
    hold on;
    grid on;
    end
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
    results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];
    res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
    hold on;
    grid on;
    %comaprison with Starlab
    if STL_results_active
        results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],char(linestyle_STL_bs(i_base)))
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    for i_base=1:n_baselines
    results=[bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed];
    res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_ESA_ISR=nanmean(([results.mean_error_L2]));
    res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.std_error_ESA_ISR=nanstd(([results.mean_error_L2]));
    plot([results.mean_error_L2],char(linestyle_analytical_MSS_SWHfixed_bs(i_base)))
    legend_text=[legend_text,strcat({'ESA - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
    hold on;
    grid on;
    %comparison results Starlab
    if STL_results_active
        results=[bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed];
        res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
        res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
        plot([results.mean_error_L2_STL],char(linestyle_STL_bs(i_base)))
        legend_text=[legend_text,strcat({'STL - ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
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
% results=[bs(i_base).SIGMA0.ESA_L2];
% res.ESA.SIGMA0.mean_rmse_fit=nanmean(([results.rmse_fitting]));
% res.ESA.SIGMA0.std_rmse_fit=nanstd(([results.rmse_fitting]));
% plot([results.rmse_fitting],linestyle_ESA)
% legend_text={'ESA'};
% hold on;
% grid on;
legend_text={''};
if ~isempty(idx_int_thres)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.THRESHOLD];
        res.bs(i_base).SIGMA0.THRESHOLD.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SIGMA0.THRESHOLD.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_threshold_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR threshold '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_ocog)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.OCOG];
        res.bs(i_base).SIGMA0.OCOG.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SIGMA0.OCOG.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_OCOG_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR OCOG '},['[',baselines_id(i_base),']'])];
        hold on;
        grid on;
    end
end
%comaprison with Starlab
if STL_results_active
    results=[bs(1).SIGMA0.ANALYTICAL_STL];
    res.STL.SIGMA0.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.STL.SIGMA0.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],char(linestyle_STL_bs(1)))
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'})];
    hold on;
    grid on;    
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    for i_base=1:n_baselines
    results=[bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed];
    res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.bs(i_base).SIGMA0.ANALYTICAL_MSS_SWHfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],char(linestyle_analytical_MSS_SWHfixed_bs(i_base)))
    legend_text=[legend_text,strcat({'ISR '},title_name_MSS_SWHfixed,' [',baselines_id(i_base),']')];
    hold on;
    grid on;
    end
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
    
    results_ESA=[bs(i_base).SIGMA0.ESA_L2];
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];  
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end        
    line_ref_y=[min([results_ESA.rmse_fitting]) max([results_ESA.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines
        hist2_bs(:,:,i_base)=hist3([[results_ESA.rmse_fitting].',[ISR_bs(i_base).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
        %[X,Y]=meshgrid(edges_x,edges_y);
    end
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    
    for i_base=1:n_baselines
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},baselines_id(i_base),{'} '},'[dB]'),'Interpreter','Tex');
        ylabel(strcat('RMSE_{','\sigma^0','}^{ESA} [dB]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},baselines_id(i_base),{'] '},'vs ESA'))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b',num2str(i_base),'_rmse_FIT_SIG0_ESA_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
end
%-------------- ISR BASELINE 1 VS ISR BASELINE OTHERS --------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    % Baseline 1 against baseline 2
    
    %All against baseline-1
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];  
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end       
    
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines-1
        hist2_bs(:,:,i_base)=hist3([[ISR_bs(1).results.rmse_fitting].',[ISR_bs(i_base+1).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
    end
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    for i_base=1:n_baselines-1
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        ylabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},baselines_id(i_base+1),{'} '},'[dB]'),'Interpreter','Tex');
        xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},baselines_id(1),{'} '},'[dB]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},baselines_id(1),{'] '},'vs',{' ['},baselines_id(i_base+1),{'] '}))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b1_vs_b',num2str(i_base+1),'_rmse_FIT_SIG0_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
end

if STL_results_active    
    results_STL=[bs(1).SIGMA0.ANALYTICAL_STL];
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SIGMA0.ANALYTICAL_SWH_MSSfixed];  
        if i_base==1
            line_ref_x=[min(ISR_bs(i_base).results.rmse_fitting) max(ISR_bs(i_base).results.rmse_fitting)];
        else
            line_ref_x=[min([line_ref_x(1),min(ISR_bs(i_base).results.rmse_fitting)]) max([line_ref_x(2),max(ISR_bs(i_base).results.rmse_fitting)])];
        end
    end     
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines
        hist2_bs(:,:,i_base)=hist3([[results_STL.rmse_fitting].',[ISR_bs(i_base).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
    end
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    for i_base=1:n_baselines
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        xlabel(strcat('RMSE_{','\sigma^0','}^{ISR',{','},baselines_id(i_base),{'} '},'[dB]'),'Interpreter','Tex');
        ylabel(strcat('RMSE_{','\sigma^0','}^{STL} [dB]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'\sigma^0',': isardSAT (ISR)',{' ['},baselines_id(i_base),{'] '},'vs STL SAMOSA-3 (S-3)'))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b',num2str(i_base),'_rmse_FIT_SIG0_STL_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
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
        for i_base=1:n_baselines
            results=[bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
            res.bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
            plot([results.RMSE_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
        for i_base=1:n_baselines
            results=[bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
            res.bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
            plot([results.mean_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
    results=[bs(1).SWH.ANALYTICAL_STL];
    res.STL.SWH.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.STL.SWH.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],char(linestyle_STL_bs(1)))
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'})];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
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
    %All against baseline-1
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed];  
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end  

    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SWH_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SWH_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines-1
        hist2_bs(:,:,i_base)=hist3([[ISR_bs(1).results.rmse_fitting].',[ISR_bs(i_base+1).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
    end
    %[X,Y]=meshgrid(edges_x,edges_y);
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    for i_base=1:n_baselines-1
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        ylabel(strcat('RMSE_{','SWH','}^{ISR',{','},baselines_id(1),{'} '},'[m]'),'Interpreter','Tex');
        xlabel(strcat('RMSE_{','SWH','}^{ISR',{','},baselines_id(i_base+1),{'} '},'[m]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},baselines_id(1),{'] '},'vs',{' ['},baselines_id(i_base+1),{'] '}))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b1_vs_b',num2str(i_base+1),'_rmse_FIT_SWH_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
end
if STL_results_active
    
    results_STL=[bs(1).SWH.ANALYTICAL_STL];
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).SWH.ANALYTICAL_SWH_MSSfixed];
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    
    for i_base=1:n_baselines
        hist2_bs(:,:,i_base)=hist3([[results_STL.rmse_fitting].',[ISR_bs(i_base).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
        %[X,Y]=meshgrid(edges_x,edges_y);
    end
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    

    for i_base=1:n_baselines
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        xlabel(strcat('RMSE_{','SWH','}^{ISR',{','},baselines_id(i_base),{'} '},'[m]'),'Interpreter','Tex');
        ylabel(strcat('RMSE_{','SWH','}^{STL} [m]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'SWH',': isardSAT (ISR)',{' ['},baselines_id(i_base),{'] '},'vs STL SAMOSA-3 (S-3)'))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b',num2str(i_base),'_rmse_FIT_SWH_STL_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
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
        for i_base=1:n_baselines
            results=[bs(i_base).COR.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).COR.ANALYTICAL_COR_MSSfixed.mean_RMSE_STL_ISR=nanmean(([results.RMSE_error_L2_STL]));
            res.bs(i_base).COR.ANALYTICAL_COR_MSSfixed.std_RMSE_STL_ISR=nanstd(([results.RMSE_error_L2_STL]));
            plot([results.RMSE_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_COR_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
        for i_base=1:n_baselines
            results=[bs(i_base).COR.ANALYTICAL_SWH_MSSfixed];
            res.bs(i_base).COR.ANALYTICAL_SWH_MSSfixed.mean_error_STL_ISR=nanmean(([results.mean_error_L2_STL]));
            res.bs(i_base).COR.ANALYTICAL_SWH_MSSfixed.std_error_STL_ISR=nanstd(([results.mean_error_L2_STL]));
            plot([results.mean_error_L2_STL],char(linestyle_STL_bs(i_base)))
            legend_text=[legend_text,strcat({'STL - ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
            hold on;
            grid on;
        end
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
    results=[bs(1).COR.ANALYTICAL_STL];
    res.STL.COR.mean_rmse_fit=nanmean(([results.rmse_fitting]));
    res.STL.COR.std_rmse_fit=nanstd(([results.rmse_fitting]));
    plot([results.rmse_fitting],char(linestyle_STL_bs(1)))
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'})];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).COR.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).COR.ANALYTICAL_SWH_MSSfixed.mean_rmse_fit=nanmean(([results.rmse_fitting]));
        res.bs(i_base).COR.ANALYTICAL_SWH_MSSfixed.std_rmse_fit=nanstd(([results.rmse_fitting]));
        plot([results.rmse_fitting],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
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
    %All against baseline-1
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).COR.ANALYTICAL_SWH_MSSfixed];
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end
    
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_COR_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_COR_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins 
    for i_base=1:n_baselines-1
        hist2_bs(:,:,i_base)=hist3([[ISR_bs(1).results.rmse_fitting].',[ISR_bs(i_base+1).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
        %[X,Y]=meshgrid(edges_x,edges_y);
    end
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    for i_base=1:n_baselines-1
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        ylabel(strcat('RMSE_{','gof','}^{ISR',{','},baselines_id(1),{'} '},'[%]'),'Interpreter','Tex');
        xlabel(strcat('RMSE_{','gof','}^{ISR',{','},baselines_id(i_base+1),{'} '},'[%]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'gof',': isardSAT (ISR)',{' ['},baselines_id(1),{'] '},'vs',{' ['},baselines_id(i_base+1),{'] '}))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b1_vs_b',num2str(i_base+1),'_rmse_FIT_COR_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
end
if STL_results_active
    
    results_STL=[bs(i_base).COR.ANALYTICAL_STL];
    for i_base=1:n_baselines
        ISR_bs(i_base).results=[bs(i_base).COR.ANALYTICAL_SWH_MSSfixed];
        if i_base==1
        line_ref_x=[min([ISR_bs(i_base).results.rmse_fitting]) max([ISR_bs(i_base).results.rmse_fitting])];
        else
            line_ref_x=[min([line_ref_x(1),min([ISR_bs(i_base).results.rmse_fitting])]) max([line_ref_x(2),max([ISR_bs(i_base).results.rmse_fitting])])];
        end
    end
    
    line_ref_y=[min([results_STL.rmse_fitting]) max([results_STL.rmse_fitting])];
    
    line_ref_x=[min([line_ref_x,line_ref_y]) max([line_ref_x,line_ref_y])];
    line_ref_y=line_ref_x;
    
    edges_x=line_ref_x(1):step_RMSE_fit_SIG0_hist2:line_ref_x(2);
    edges_y=line_ref_y(1):step_RMSE_fit_SIG0_hist2:line_ref_y(2);

    ctrs{1}=edges_x;%centers of bins 
    ctrs{2}=edges_y;%centers of bins
    for i_base=1:n_baselines
        hist2_bs(:,:,i_base)=hist3([[results_STL.rmse_fitting].',[ISR_bs(i_base).results.rmse_fitting].'],'Edges',ctrs);
        %hist2_b1(hist2_b1==0)=NaN;
        %[X,Y]=meshgrid(edges_x,edges_y);
    end
    
    min_image=1;
    max_image=max(hist2_bs(:));
    
    %# size of colormap
    n = size(jet,1);
    %# color step
    dmap=(max_image-min_image)/n;
    
    %---------------------- baseline 1 ------------------------------------
    for i_base=1:n_baselines
        figure;
        %surf(flipud(X),flipud(Y),flipud(hist2_b1),'LineStyle','none','FaceColor','flat'); colormap('jet'); c=colorbar; ylabel(c,'[# tracks]'); view(2); grid off;
        %caxis([min_image, max_image]);
        imagesc(edges_x,fliplr(edges_y),flipud(hist2_bs(:,:,i_base))); colormap([nanclr; jet]);
        set(gca,'YDir','normal')
        %caxis([min_image-dmap,max_image]);
        hcb=colorbar; ylim(hcb,[min_image max_image]);
        hold on;
        plot(line_ref_x,line_ref_y,'--k');
        xlabel(strcat('RMSE_{','gof','}^{ISR',{','},baselines_id(i_base),{'} '},'[%]'),'Interpreter','Tex');
        ylabel(strcat('RMSE_{','gof','}^{STL} [%]'),'Interpreter','Tex');
        title(strcat('RMSE error on fitted',{' '},'gof',': isardSAT (ISR)',{' ['},baselines_id(i_base),{'] '},'vs STL SAMOSA-3 (S-3)'))
        print('-dpng',strcat(path_comparison_results,'Scatter_plot_b',num2str(i_base),'_rmse_FIT_COR_STL_ISR.png'))
    end
    clear ctrs X Y edges_x edges_y hist2_bs;
end


%$$$$$$$$$$$$$$$$$$$$$$$$$ MEAN VALUE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure;
legend_text={''};
if STL_results_active
    results=[bs(1).COR.ANALYTICAL_STL];
    res.STL.COR.mean_mean_value=nanmean(([results.mean]));
    res.STL.COR.std_mean_value=nanstd(([results.mean]));
    plot([results.rmse_fitting],char(linestyle_STL_bs(i_base)))
    legend_text=[legend_text,strcat({'STL SAMOSA-3 (S-3)'})];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    for i_base=1:n_baselines
        results=[bs(i_base).COR.ANALYTICAL_SWH_MSSfixed];
        res.bs(i_base).COR.ANALYTICAL_SWH_MSSfixed.mean_mean_value=nanmean(([results.mean]));
        res.bs(i_base).COR.ANALYTICAL_SWH_MSSfixed.std_mean_value=nanstd(([results.mean]));
        plot([results.mean],char(linestyle_analytical_SWH_MSSfixed_bs(i_base)))
        legend_text=[legend_text,strcat({'ISR '},title_name_SWH_MSSfixed,' [',baselines_id(i_base),']')];
        hold on;
        grid on;
    end
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

