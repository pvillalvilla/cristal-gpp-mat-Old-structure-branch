function [statistics]=comparison_GR_L2_files(input_files, resultPath, ID_baselines, ID_GR_file_string,varargin)
%% optional inputs managing
if(nargin<4 || nargin>(4+12*2))
    error('Wrong number of input parameters');
end
p = inputParser;
p.addParamValue('ref_SSH',NaN); %reference SSH to compute error for S6 sim
p.addParamValue('ref_SWH',NaN); %reference SWH to compute error for S6 sim
p.addParamValue('ref_sigma0',NaN); %reference sigma0 to compute error for S6
p.addParamValue('ref_SSH_whole',[]); %reference SSH whole track (used for slope simulations)
p.addParamValue('lat_SSH_whole',[]); %reference lat whole track
p.addParamValue('figure_format','jpg'); 
p.addParamValue('res_fig','-r300'); 
p.addParamValue('win_size_stat',20); 
p.addParamValue('marker_bs',{'none','none','none','none','none','none','none'}); 
p.addParamValue('color_bs',[rgb('red'); rgb('blue'); rgb('DarkGreen'); rgb('Cyan'); rgb('Green'); rgb('Magenta'); rgb('Black');]); 
p.addParamValue('line_bs',{'-','--',':','-.'}); 
p.addParamValue('plot_figure',1); 

p.parse(varargin{:});
ref_SSH       = p.Results.ref_SSH;
ref_SWH       = p.Results.ref_SWH;
ref_sigma0    = p.Results.ref_sigma0;
ref_SSH_whole = p.Results.ref_SSH_whole;
lat_SSH_whole = p.Results.lat_SSH_whole;
figure_format = p.Results.figure_format;
res_fig       = p.Results.res_fig;
win_size_stat = p.Results.win_size_stat;
marker_bs     = p.Results.marker_bs;
color_bs      = p.Results.color_bs;
line_bs       = p.Results.line_bs;
plot_figure   = p.Results.plot_figure;
clear p;

%% Formating of ploting options
set_default_plot;
%set(0,'defaultFigureVisible','off');
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


%% Local variables 
N_files = length(input_files);
cst_p.sec_in_day_cst=86400;
cst_p.c_cst=299792458.0; 

SSH_mean_std       = NaN(1,N_files);
SSH_mean_bias      = NaN(1,N_files);
SWH_mean_std       = NaN(1,N_files);
SWH_mean_bias      = NaN(1,N_files);
sigma0_mean_std    = NaN(1,N_files);
sigma0_mean_bias   = NaN(1,N_files);
COR_mean_std       = NaN(1,N_files);
COR_mean_mean      = NaN(1,N_files);


SSH_mean_std_str       = cell(1,N_files);
SSH_mean_bias_str      = cell(1,N_files);
SWH_mean_std_str       = cell(1,N_files);
SWH_mean_bias_str      = cell(1,N_files);
sigma0_mean_std_str    = cell(1,N_files);
sigma0_mean_bias_str   = cell(1,N_files);
COR_mean_std_str       = cell(1,N_files);
COR_mean_mean_str      = cell(1,N_files);



%% Comparison of different files
for i_file = 1:N_files
    % Read file    
    data = readL2_ISD (char(input_files(i_file)),{'SWH'},cst_p);
    
%     if strcmp(char(ID_baselines(i_file)),'RAW')
%         in_lat = min(data.GEO.LAT);
%         end_lat = max(data.GEO.LAT);        
%     end
%     switch char(ID_baselines(i_file))
%         case {'LR-RMC','LR-RMC RAW'}
%             idx_int = find(data.GEO.LAT>=in_lat & data.GEO.LAT<=end_lat);
%             data.GEO.LAT = data.GEO.LAT(idx_int);
%             data.GEO.LON = data.GEO.LON(idx_int);
%             data.GR.analytical_SWH.SSH=data.GR.analytical_SWH.SSH(idx_int);
%             data.GR.analytical_SWH.SWH=data.GR.analytical_SWH.SWH(idx_int);
%             data.GR.analytical_SWH.sig0=data.GR.analytical_SWH.sig0(idx_int);
%             data.GR.analytical_SWH.COR=data.GR.analytical_SWH.COR(idx_int);
%     end
    
    
    % Compute statistics
    if ~isempty(ref_SSH_whole)
        
        [~,~,SSH_mean_std(i_file),~]=block_statistics([data.GR.analytical_SWH.SSH],win_size_stat);        
        SSH_mean_std_str(i_file) = {num2str(SSH_mean_std(i_file)*100)};
        
        %look for the SSH ref closest to each surface 
        idx_ref_SSH = 1:1:length([data.GR.analytical_SWH.SSH]);
        for i_surface = 1:length([data.GR.analytical_SWH.SSH])
            [~,idx_ref_SSH(i_surface)] = min(abs(data.GEO.LAT(i_surface)-ref_SSH_lat));
        end        
        [~,~,~,SSH_mean_bias(i_file)]=block_statistics([data.GR.analytical_SWH.SSH].'-ref_SSH_whole(idx_ref_SSH),win_size_stat);
        SSH_mean_bias_str(i_file) = {num2str(SSH_mean_bias(i_file)*100)};
    else
        [~,~,SSH_mean_std(i_file),SSH_mean_bias(i_file)]=block_statistics([data.GR.analytical_SWH.SSH],win_size_stat);
        SSH_mean_bias(i_file) = SSH_mean_bias(i_file)-ref_SSH;
        SSH_mean_std_str(i_file) = {num2str(SSH_mean_std(i_file)*100)};
        SSH_mean_bias_str(i_file) = {num2str(SSH_mean_bias(i_file)*100)};
    end

    
    [~,~,SWH_mean_std(i_file),SWH_mean_bias(i_file)]=block_statistics([data.GR.analytical_SWH.SWH],win_size_stat);
    SWH_mean_bias(i_file) = SWH_mean_bias(i_file) - ref_SWH;
    SWH_mean_std_str(i_file) = {num2str(SWH_mean_std(i_file)*100)};
    SWH_mean_bias_str(i_file) = {num2str(SWH_mean_bias(i_file)*100)};
    
%     if ~isempty(strfind(char(ID_baselines(i_file)),'LR-RMC RAW'))        
%         onboard_proc_sar_raw_chd = 88.27;
%         onboard_proc_sar_rmc_chd = 70.21;
%         data.GR.analytical_SWH.sig0=data.GR.analytical_SWH.sig0-(onboard_proc_sar_raw_chd-onboard_proc_sar_rmc_chd);               
%     end
    [~,~,sigma0_mean_std(i_file),sigma0_mean_bias(i_file)]=block_statistics([data.GR.analytical_SWH.sig0],win_size_stat);
    sigma0_mean_bias(i_file) = sigma0_mean_bias(i_file) - ref_sigma0;
    sigma0_mean_std_str(i_file) = {num2str(sigma0_mean_std(i_file))};
    sigma0_mean_bias_str(i_file) = {num2str(sigma0_mean_bias(i_file))};
    
    
    [~,~,COR_mean_std(i_file),COR_mean_mean(i_file)]=block_statistics([data.GR.analytical_SWH.COR],win_size_stat);
    COR_mean_std_str(i_file) = {num2str(COR_mean_std(i_file))};
    COR_mean_mean_str(i_file) = {num2str(COR_mean_mean(i_file))};
    
    
    if plot_figure == 1
        % Plot
        if i_file == 1
            f1=figure;
        end
        
        % SSH
        ax_SSH=subplot(4,1,1);
        plot(ax_SSH,data.GEO.LAT,data.GR.analytical_SWH.SSH,...
            'Marker',char(marker_bs(i_file)),'Color',color_bs(i_file,:),'Linestyle',char(line_bs(i_file)));
        hold on;
        if i_file==N_files
            if ~isempty(ref_SSH_whole)
                plot(ax_SSH,lat_SSH_whole,ref_SSH_whole,...
                    'Marker','None','Color',rgb('magenta'),'Linestyle','-');
                legend(ax_SSH,[ID_baselines,{'ref. SSH'}],'Orientation','Horizontal','Location','Best');
            else
                legend(ax_SSH,ID_baselines,'Orientation','Horizontal','Location','Best');
            end
            xlabel(ax_SSH,'Latitude [deg.]'); ylabel(ax_SSH,'SSH [m]');
            title(ax_SSH,strcat(' std  [cm]: [',strjoin(strcat(ID_baselines,{': '},SSH_mean_std_str),'  '),']',{' || '},...
                ' bias [cm]: [',strjoin(strcat(ID_baselines,{': '},SSH_mean_bias_str),'  '),']'),'Interpreter','None')
        end
        
        
        % SWH
        ax_SWH=subplot(4,1,2);
        plot(ax_SWH,data.GEO.LAT,data.GR.analytical_SWH.SWH,...
            'Marker',char(marker_bs(i_file)),'Color',color_bs(i_file,:),'Linestyle',char(line_bs(i_file)));
        hold on;
        if i_file==N_files
            xlabel(ax_SWH,'Latitude [deg.]'); ylabel(ax_SWH,'H_s [m]');
            legend(ax_SWH,ID_baselines,'Orientation','Horizontal','Location','Best');
            title(ax_SWH,strcat(' std  [cm]: [',strjoin(strcat(ID_baselines,{': '},SWH_mean_std_str),'  '),']',{' || '},...
                ' bias [cm]: [',strjoin(strcat(ID_baselines,{': '},SWH_mean_bias_str),'  '),']'),'Interpreter','None')
        end
        
        % Sigma0
        ax_sigma0=subplot(4,1,3);
        plot(ax_sigma0,data.GEO.LAT,data.GR.analytical_SWH.sig0,...
            'Marker',char(marker_bs(i_file)),'Color',color_bs(i_file,:),'Linestyle',char(line_bs(i_file)));
        hold on;
        if i_file==N_files
            xlabel(ax_sigma0,'Latitude [deg.]'); ylabel(ax_sigma0,'\sigma^0 [dB]');
            legend(ax_sigma0,ID_baselines,'Orientation','Horizontal','Location','Best');
            title(ax_sigma0,strcat(' std  [dB]: [',strjoin(strcat(ID_baselines,{': '},sigma0_mean_std_str),'  '),']',{' || '},...
                ' bias [dB]: [',strjoin(strcat(ID_baselines,{': '},sigma0_mean_bias_str),'  '),']'),'Interpreter','None')
        end
        
        % COR
        ax_COR=subplot(4,1,4);
        plot(ax_COR,data.GEO.LAT,data.GR.analytical_SWH.COR,...
            'Marker',char(marker_bs(i_file)),'Color',color_bs(i_file,:),'Linestyle',char(line_bs(i_file)));
        hold on;
        if i_file==N_files
            xlabel(ax_COR,'Latitude [deg.]'); ylabel(ax_COR,'\rho [%]');
            legend(ax_COR,ID_baselines,'Orientation','Horizontal','Location','Best');
            title(ax_COR,strcat(' std  [%]: [',strjoin(strcat(ID_baselines,{': '},COR_mean_std_str),'  '),']',{' || '},...
                ' mean [%]: [',strjoin(strcat(ID_baselines,{': '},COR_mean_mean_str),'  '),']'),'Interpreter','None')
        end
        
        if i_file==N_files
            [axT,hT]=suplabel(ID_GR_file_string,'t');
            print(print_file,res_fig,strcat(resultPath,char(strcat('GR_comparison_',...
                strjoin(strrep(ID_baselines,'-','_'),'__'),'_',ID_GR_file_string,file_ext))))
            close(f1)
        end
        
    end
    clear data;
    
end

statistics.SSH_mean_std=SSH_mean_std;
statistics.SSH_mean_bias=SSH_mean_bias;
statistics.SWH_mean_std=SWH_mean_std;
statistics.SWH_mean_bias=SWH_mean_bias;
statistics.sigma0_mean_std=sigma0_mean_std;
statistics.sigma0_mean_bias=sigma0_mean_bias;
statistics.COR_mean_std=COR_mean_std;
statistics.COR_mean_mean=COR_mean_mean;

end
