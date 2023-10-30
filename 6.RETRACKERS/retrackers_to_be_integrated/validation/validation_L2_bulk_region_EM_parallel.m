function validation_L2_bulk_region_EM_parallel(input_path_L2_ISD,input_path_L2_ESA,path_comparison_results,varargin)
time_init=tic;
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
version_matlab=version;
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<3 || nargin>(3+11*2))
    error('Wrong number of input parameters');
end
p = inputParser;
p.addParamValue('input_path_L2_STL',{''},@(x)ischar(x));
p.addParamValue('retrackers',{''},@(x)iscellstr(x)); %is a cell array with different names of different retrackers in L2 product
p.addParamValue('figures_visible',0);
p.addParamValue('num_pools',1);
p.addParamValue('flag_outliers_removal',0);
p.addParamValue('type_outliers_removal','percentiles');
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('active_comparison_tracks',1);
p.addParamValue('active_validation_tracks',1);
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.addParamValue('annotation_box_active',1);
p.parse(varargin{:});
input_path_L2_STL=char(p.Results.input_path_L2_STL);
retrackers=p.Results.retrackers;
figures_visible=p.Results.figures_visible;
num_pools=p.Results.num_pools;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=p.Results.type_outliers_removal;
sh_name_nc=p.Results.sh_name_nc;
active_comparison_tracks=p.Results.active_comparison_tracks;
active_validation_tracks=p.Results.active_validation_tracks;
filename_mask_KML=p.Results.filename_mask_KML;
annotation_box_active=p.Results.annotation_box_active;
clear p;

%----------------- Define linstyles for bulk comparison -------------------
linestyle_ESA='or';
linestyle_analytical_SWH_MSSfixed='*b';
linestyle_analytical_MSS_SWHfixed='^g';
linestyle_threshold='+k';
linestyle_OCOG='dc';
linestyle_STL='sg';
fontsize_xlabel_tracks=8;
size_marker=7;
title_name_SWH_MSSfixed='analytical';
title_name_MSS_SWHfixed='analytical-fit-MSS-and-SWH-fixed';

%-------------- Minimum Number of Surfaces within geo mask ----------------
min_num_surf_validation=50;


%------------------- Check the available retrackers to be processed ------
idx_int_analytical_SWH_MSSfixed=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL_SWH')), 1);
idx_int_analytical_MSS_SWHfixed=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL_MSS')), 1);
idx_int_thres=find(~cellfun(@isempty,strfind(retrackers,'THRESHOLD')), 1);
idx_int_ocog=find(~cellfun(@isempty,strfind(retrackers,'OCOG')), 1);


filesBulk.inputPath       =   input_path_L2_ISD;
mkdir(path_comparison_results);

filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,'.nc')));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

i_files_valid=0;
%% -------------- check available files -----------------------------------
if isempty(filename_mask_KML)
    within_geo_mask=1;
    geo_mask  = [];
else
    geo_mask  = kml2lla(filename_mask_KML);
end
%--------------------------------------------------------------------------
for i_file=1:filesBulk.nFilesNC
    filename_L2_ISD=char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);
    data_string=filename_L2_ISD(17:17+30);
    file_char_name=filename_L2_ISD(17+30+1+17:end-3);
    
    % ---- Checking whether the track within geomask if available ---------
    if ~isempty(filename_mask_KML)
        LAT_ISD_L2=double(ncread([input_path_L2_ISD filename_L2_ISD],'lat_20_ku')).';        
        LON_ISD_L2=double(ncread([input_path_L2_ISD filename_L2_ISD],'lon_20_ku')).';
        
%         idx_lt_0= LON_ISD_L2<0;
%         %longitudes +- values (-180,180)
%         if any(idx_lt_0)
%             LON_ISD_L2(idx_lt_0)=LON_ISD_L2(idx_lt_0)+360.0;
%         end
% 
%         clear idx_lt_0;
                
        idx_int=inpolygon(LON_ISD_L2,LAT_ISD_L2,geo_mask.coord(:,1),geo_mask.coord(:,2));
        if ~any(idx_int)
            disp(strcat('Track,',{' '},filename_L2_ISD,{' '},'outside the limits of the geographical mask'))
            within_geo_mask=0;
        else
            if (length(find(idx_int)) < min_num_surf_validation)
                disp(strcat('Track,',{' '},filename_L2_ISD,{' '},'has a # surfaces within geographical mask below threshold'))
                within_geo_mask=0;
            else
                within_geo_mask=1;
            end            
        end
        
        clear LAT_ISD_L2 LON_ISD_L2 idx_int idx_lt_180 idx_lt_90;
        
    end
    
    inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' data_string(1:15) '*_C001.DBL']));
    if isempty(inputL2ESAFiles)
        %add one second to initial time acquisition
        init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
        inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
        if isempty(inputL2ESAFiles)
            %add one second to initial time acquisition
            end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
            inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
            if isempty(inputL2ESAFiles)
                inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
            end
        end
        
    end
    %--------------- Checking whether Starlab L2 file is available --------
    if isempty(input_path_L2_STL)
        inputL2STLFiles   = 1;
    else
        inputL2STLFiles   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
    end
    
    if ~isempty(inputL2ESAFiles) && ~isempty(inputL2STLFiles) && within_geo_mask
        i_files_valid=i_files_valid+1;
        filesBulk.indexFilesNC_valid(i_files_valid)=filesBulk.indexFilesNC(i_file);
        date_file_id(i_files_valid,1:(length(data_string)+length(file_char_name)+1))=strcat(strrep(data_string(1:31),'_','-'),'_',file_char_name);
        
    end

end


%% ------------------ RUN the validation for each available file ----------
%--------------------------------------------------------------------------
filesBulk.nFilesNC_valid=length(filesBulk.indexFilesNC_valid);
disp(strcat('Total number of input L2 ISR for evaluation: ',num2str(filesBulk.nFilesNC)));
disp(strcat('Total number of valid files for evaluation: ',num2str(filesBulk.nFilesNC_valid)));
if active_validation_tracks
    if num_pools~=1
        %create pools
        if str2double(version_matlab(end-5:end-2))>2013
            parpool(num_pools);
        else
            matlabpool('open',num_pools);
        end
        %% ------------- Loop per file to be processed ------------------------
        parfor i_files_valid=1:filesBulk.nFilesNC_valid
            try
                run_L2_validation(filesBulk,i_files_valid,input_path_L2_ISD,...
                    input_path_L2_ESA,...
                    path_comparison_results,...
                    'input_path_L2_STL',input_path_L2_STL,....
                    'figures_visible',figures_visible,...
                    'retrackers',retrackers,'flag_outliers_removal',flag_outliers_removal,...
                    'type_outliers_removal',type_outliers_removal,...
                    'sh_name_nc',sh_name_nc,...
                    'geo_mask',geo_mask,...
                    'annotation_box_active',annotation_box_active)
            catch
                disp(strcat('Some error in validation for file',{''},char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name)));
                continue;
            end
        end
        %close pools
        if str2double(version_matlab(end-5:end-2))>2013
            poolobj = gcp('nocreate');
            delete(poolobj);
        else
            matlabpool('close');
        end
    else
        for i_files_valid=1:filesBulk.nFilesNC_valid
            try
                run_L2_validation(filesBulk,i_files_valid,input_path_L2_ISD,...
                    input_path_L2_ESA,...
                    path_comparison_results,...
                    'input_path_L2_STL',input_path_L2_STL,....
                    'figures_visible',figures_visible,...
                    'retrackers',retrackers,'flag_outliers_removal',flag_outliers_removal,...
                    'type_outliers_removal',type_outliers_removal,...
                    'sh_name_nc',sh_name_nc,...
                    'geo_mask',geo_mask,...
                    'annotation_box_active',annotation_box_active)
            catch
                disp(strcat('Some error in validation for file',{''},char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name)));
                continue;
            end
        end
    end
end
%% -------------COMPARISON OF TRACKS --------------------------------------
%--------------------------------------------------------------------------
if active_comparison_tracks
%loading and reordering the data into a single array of structures
%filesBulk.inputFilesEvaluation      =   dir(fullfile(path_comparison_results,'*_L2_Evaluation.mat'));
for i_files_valid=1:filesBulk.nFilesNC_valid    
    %load(strcat(path_comparison_results,char(filesBulk.inputFilesEvaluation(i_files_valid).name)))
    name_file_L2_ISD=char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name);
    name_file_L2_ISD=name_file_L2_ISD(17:17+30);
%    aux=name_file_L2_ISD;
%    aux2=strcat(strrep(char(name_file_L2_ISD),'.nc','_'),'L2_Evaluation.mat');
    disp(name_file_L2_ISD)
    load(strcat(path_comparison_results,name_file_L2_ISD,'/',strcat(name_file_L2_ISD,'_L2_Evaluation.mat')));
    SIGMA0(i_files_valid)=res.SIGMA0;
    SSH(i_files_valid)=res.SSH;
    SWH(i_files_valid)=res.SWH;
    COR(i_files_valid)=res.COR;
end
save(strcat(path_comparison_results,'L2_Bulk_validation_information.mat'),'SIGMA0','SSH','SWH','COR');


%% ----------  Ploting ----------------------------------------------------
xlabels=date_file_id;
if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end
set(0,'defaultLineMarkerSize',size_marker);  % set the default line marker size
%% ---------------------------  SSH ---------------------------------------
%--------------------------------------------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison retracker w.r.t ESA $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SSH(:).THRESHOLD];
    plot([results.RMSE_error_L2],linestyle_threshold)
    legend_text=[legend_text,{'L2 ESA-ISR (threshold-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(idx_int_ocog)
    results=[SSH(:).OCOG];
    plot([results.RMSE_error_L2],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ESA-ISR (OCOG-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(input_path_L2_STL)
    results=[SSH(:).ANALYTICAL_STL];
    plot([results.RMSE_error_L2],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    hold on;
    grid on;
    title(strcat('RMSE error on',{' '},upper(sh_name_nc),': ESA - isardSAT (ISR) & ESA - Starlab (STL)'))
else
    title(strcat('RMSE error on',{' '},upper(sh_name_nc),': ESA - isardSAT'))
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.RMSE_error_L2],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.RMSE_error_L2],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('RMSE_{',upper(sh_name_nc),{'} '},'[m]'),'Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'RMSE_SSH_ESA_ISD.png'))
%&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SSH(:).THRESHOLD];
    plot([results.mean_error_L2],linestyle_threshold)
    legend_text=[legend_text,{'L2 ESA-ISR (threshold-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(idx_int_ocog)
    results=[SSH(:).OCOG];
    plot([results.mean_error_L2],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ESA-ISR (OCOG-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(input_path_L2_STL)
    results=[SSH(:).ANALYTICAL_STL];
    plot([results.mean_error_L2],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    hold on;
    grid on;
    title(strcat('Mean error on',{' '},upper(sh_name_nc),': ESA - isardSAT (ISR) & ESA - Starlab (STL)'))
else
    title(strcat('Mean error on',{' '},upper(sh_name_nc),': ESA - isardSAT'))
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean_error_L2],linestyle_analytical_SWH_MSSfixed)
   legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.mean_error_L2],linestyle_analytical_MSS_SWHfixed)
   legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel(strcat('\epsilon_{',upper(sh_name_nc),'} [m]'),'Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'Mean_error_SSH_ESA_ISD.png'))

%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
results=[SSH(:).ESA_L2];
plot([results.rmse_fitting],linestyle_ESA)
legend_text={'L2 ESA'};
hold on;
grid on;
if ~isempty(idx_int_thres)
    results=[SSH(:).THRESHOLD];
    plot([results.rmse_fitting],linestyle_threshold)
    legend_text=[legend_text,{'L2 ISR (threshold-retracker)'}];
end
if ~isempty(idx_int_ocog)
    results=[SSH(:).OCOG];
    plot([results.rmse_fitting],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ISR (OCOG-retracker)'}];
end
title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': ESA & isardSAT (ISR)'))
if ~isempty(input_path_L2_STL)
    results=[SSH(:).ANALYTICAL_STL];
    plot([results.rmse_fitting],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': ESA & isardSAT (ISR) & Starlab (STL)'))
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel('RMSE_{SSH} [m]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SSH.png'))
%&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
results=[SSH(:).ESA_L2];
plot([results.mean_error_fitting],linestyle_ESA)
legend_text={'L2 ESA'};
hold on;
grid on;
if ~isempty(idx_int_thres)
    results=[SSH(:).THRESHOLD];
    plot([results.mean_error_fitting],linestyle_threshold)
    legend_text=[legend_text,{'L2 ISR (threshold-retracker)'}];
end
if ~isempty(idx_int_ocog)
    results=[SSH(:).OCOG];
    plot([results.mean_error_fitting],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ISR (OCOG-retracker)'}];
end
title(strcat('Mean error on fitted',{' '},upper(sh_name_nc),': ESA & isardSAT (ISR)'))
if ~isempty(input_path_L2_STL)
    results=[SSH(:).ANALYTICAL_STL];
    plot([results.mean_error_fitting],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    title(strcat('Mean error on fitted',{' '},upper(sh_name_nc),': ESA & isardSAT (ISR) & Starlab (STL)'))
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SSH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean_error_fitting],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SSH(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.mean_error_fitting],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel_str='Track'; ylabel('\epsilon_{SSH} [m]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SSH.png'))
%% ---------------------------  SIGMA0 ------------------------------------
%$$$$$$$$$$$$$$$$$ Comparison retracker w.r.t ESA $$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SIGMA0(:).THRESHOLD];
    plot([results.RMSE_error_L2],linestyle_threshold)
    legend_text=[legend_text,{'L2 ESA-ISR (threshold-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0(:).OCOG];
    plot([results.RMSE_error_L2],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ESA-ISR (OCOG-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(input_path_L2_STL)
    results=[SIGMA0(:).ANALYTICAL_STL];
    plot([results.RMSE_error_L2],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    hold on;
    grid on;
    title('RMSE error on \sigma^0: ESA - isardSAT & ESA - Starlab','Interpreter','Tex')
else
    title('RMSE error on \sigma^0: ESA - isardSAT','Interpreter','Tex')
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.RMSE_error_L2],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.RMSE_error_L2],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end

legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('RMSE_{\sigma^0} [dB]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'RMSE_SIGMA0_ESA_ISD.png'))
%&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
if ~isempty(idx_int_thres)
    results=[SIGMA0(:).THRESHOLD];
    plot([results.mean_error_L2],linestyle_threshold)
    legend_text=[legend_text,{'L2 ESA-ISR (threshold-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0(:).OCOG];
    plot([results.mean_error_L2],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ESA-ISR (OCOG-retracker)'}];
    hold on;
    grid on;
end
if ~isempty(input_path_L2_STL)
    results=[SIGMA0(:).ANALYTICAL_STL];
    plot([results.mean_error_L2],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    hold on;
    grid on;
    title('Mean error on \sigma^0: ESA - isardSAT & ESA - Starlab','Interpreter','Tex')
else
    title('Mean error on \sigma^0: ESA - isardSAT','Interpreter','Tex')
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0(:).ANALYTICAL];
    plot([results.mean_error_L2],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean_error_L2],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ESA-ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end

legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('\epsilon_{\sigma^0} [dB]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'Mean_error_SIGMA0_ESA_ISD.png'))

%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
results=[SIGMA0(:).ESA_L2];
plot([results.rmse_fitting],linestyle_ESA)
legend_text={'L2 ESA'};
hold on;
grid on;
if ~isempty(idx_int_thres)
    results=[SIGMA0(:).THRESHOLD];
    plot([results.rmse_fitting],linestyle_threshold)
    legend_text=[legend_text,{'L2 ISR (threshold-retracker)'}];
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0(:).OCOG];
    plot([results.rmse_fitting],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ISR (OCOG-retracker)'}];
end
title('RMSE error on fitted \sigma^0: ESA & isardSAT (ISR)','Interpreter','Tex')
if ~isempty(input_path_L2_STL)
    results=[SIGMA0(:).ANALYTICAL_STL];
    plot([results.rmse_fitting],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    title('RMSE error on fitted \sigma^0: ESA & isardSAT (ISR) & Starlab (STL)','Interpreter','Tex')
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
end
legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('RMSE_{\sigma^0} [dB]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SIGMA0.png'))
%&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
results=[SIGMA0(:).ESA_L2];
plot([results.mean_error_fitting],linestyle_ESA)
legend_text={'L2 ESA'};
hold on;
grid on;
if ~isempty(idx_int_thres)
    results=[SIGMA0(:).THRESHOLD];
    plot([results.mean_error_fitting],linestyle_threshold)
    legend_text=[legend_text,{'L2 ISR (threshold-retracker)'}];
end
if ~isempty(idx_int_ocog)
    results=[SIGMA0(:).OCOG];
    plot([results.mean_error_fitting],linestyle_OCOG)
    legend_text=[legend_text,{'L2 ISR (OCOG-retracker)'}];
end
title('Mean error on fitted \sigma^0: ESA & isardSAT (ISR)','Interpreter','Tex')
if ~isempty(input_path_L2_STL)
    results=[SIGMA0(:).ANALYTICAL_STL];
    plot([results.mean_error_fitting],linestyle_STL)
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    title('Mean error on fitted \sigma^0: ESA & isardSAT (ISR) & Starlab (STL)','Interpreter','Tex')
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SIGMA0(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean_error_fitting],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SIGMA0(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.mean_error_fitting],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
end

legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('\epsilon_{\sigma^0} [dB]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SIGMA0.png'))

%% ---------------------------  SWH ---------------------------------------
%$$$$$$$$$$$$$$$$$$$$$$$$$ Error w.r.t STL $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if ~isempty(input_path_L2_STL) && ~isempty(idx_int_analytical_SWH_MSSfixed)
    res.SWH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2_STL
    figure;
    results=[SWH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean_error_L2_STL],linestyle_analytical_SWH_MSSfixed)
    %legend_text={strcat('L2 ISR (',title_name_SWH_MSSfixed,')')};
    title('Mean error on SWH: Starlab (STL) - isardSAT (ISR)','Interpreter','Tex')
    %legend(legend_text(~cellfun(@isempty,legend_text)));
    ylabel('\epsilon_{SWH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)
    print('-dpng',strcat(path_comparison_results,'Mean_error_ISD_STL_SWH.png'))
end


%$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
title('RMSE error on fitted SWH: ESA & isardSAT (ISR)','Interpreter','Tex')
if ~isempty(input_path_L2_STL)
    results=[SWH(:).ANALYTICAL_STL];
    plot([results.rmse_fitting],linestyle_STL)
    hold on;
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];
    title('RMSE error on fitted SWH: ESA & isardSAT (ISR) & Starlab (STL)','Interpreter','Tex')
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SWH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.rmse_fitting],linestyle_analytical_SWH_MSSfixed)    
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SWH(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.rmse_fitting],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end
legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('RMSE_{SWH} [m]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SWH.png'))
%&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure;
legend_text={''};
title('Mean error on fitted SWH: ESA & isardSAT (ISR)','Interpreter','Tex')
if ~isempty(input_path_L2_STL)
    results=[SWH(:).ANALYTICAL_STL];
    plot([results.mean_error_fitting],linestyle_STL)
    hold on;
    legend_text=[legend_text,{'L2 STL (a Sentinel-3)'}];    
    title('Mean error on fitted SWH: ESA & isardSAT (ISR) & Starlab (STL)','Interpreter','Tex')
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[SWH(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean_error_fitting],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[SWH(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.mean_error_fitting],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end
legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('\epsilon_{SWH} [m]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SWH.png'))

%% ----------------------  PEARSON CORR COEF. -----------------------------
%$$$$$$$$$$$$$$$$$$$$$$$$$ MEAN VALUE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure;
legend_text={''};
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    results=[COR(:).ANALYTICAL_SWH_MSSfixed];
    plot([results.mean],linestyle_analytical_SWH_MSSfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_SWH_MSSfixed,')')}];
    hold on;
    grid on;
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    results=[COR(:).ANALYTICAL_MSS_SWHfixed];
    plot([results.mean],linestyle_analytical_MSS_SWHfixed)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_MSS_SWHfixed,')')}];
    hold on;
    grid on;
end
title('Mean value pearson corr. coeff. \rho_{pearson}: ESA & isardSAT','Interpreter','Tex')
legend(legend_text(~cellfun(@isempty,legend_text)));
ylabel('\rho_{pearson} [%]','Interpreter','Tex');
xlabel(xlabel_str)
% %-------------- set the x-tick locations ----------------------------------
% % Reduce the size of the axis so that all the labels fit in the figure.
% pos = get(gca,'Position');
% set(gca,'Position',[pos(1), .2, pos(3) .65])
% Xt=1:1:filesBulk.nFilesNC_valid;
% Xl=[1 filesBulk.nFilesNC_valid];
% Yl=get(gca,'ylim');
% set(gca,'XTick',Xt,'XLim',Xl);
% ax = axis;    % Current axis limits
% axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
% % Place the text labels
% t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
% set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%       'Rotation',45);
% % Remove the default labels
% set(gca,'XTickLabel','')
% % Get the Extent of each text object.  This
% % loop is unavoidable.
% for i = 1:length(t)
%   ext(i,:) = get(t(i),'Extent');
% end
% % Determine the lowest point.  The X-label will be
% % placed so that the top is aligned with this point.
% LowYPoint = min(ext(:,2));
% % Place the axis label at this point
% XMidPoint = Xl(1)+abs(diff(Xl))/2;
% tl = text(XMidPoint,LowYPoint,xlabel_str, ...
%           'VerticalAlignment','top', ...
%           'HorizontalAlignment','center');
print('-dpng',strcat(path_comparison_results,'Mean_COR.png'))
close all
end
time_end=toc(time_init);
minutes_processing = floor(time_end/60);
secs_processing = time_end - minutes_processing*60;
disp(['Validation/processing time for ',input_path_L2_ISD,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);    
end

