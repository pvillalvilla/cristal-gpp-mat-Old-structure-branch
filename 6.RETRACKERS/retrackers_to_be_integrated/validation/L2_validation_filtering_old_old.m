%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VALIDATION L2/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this algorithm is to cross-check the L2 products:
% isardSAT and ESA for CryoSat-2 data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res]=L2_validation_filtering_old_old(filename_L2_ISD,filename_L2_ESA,path_results_comparison,varargin)

%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<3 || nargin>(3+4*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('filename_L2_STL',{''},@(x)ischar(x));
p.addParamValue('retrackers',{''},@(x)iscellstr(x)); %is a cell array with different names of different retrackers in L2 product
p.addParamValue('figures_visible',0);
p.addParamValue('flag_outliers_removal',0);
p.parse(varargin{:});
filename_L2_STL=char(p.Results.filename_L2_STL);
retrackers=p.Results.retrackers;
figures_visible=p.Results.figures_visible;
flag_outliers_removal=p.Results.flag_outliers_removal;
clear p;


close all;
set_default_plot;
%figures_visible=1;

%removal of outliers in ESA SSH retrievals L2
threshold_std=3.0;

%removal of outliers using percentile & IQR (Interquartile Range)
% outliers data<(percentil_low-IQR_times*IQR) | data>(percentil_high+IQR_times*IQR)
IQR_times=1.5; %number of IQR 
outlier_percentil_low=25.0;
outlier_percentil_high=75.0;

%sliding window definitions: smoothing function used as fitting
% number of samples at each side of current position
sliding_window_SSH=5;
sliding_window_SWH=5;
sliding_window_sigma0=5;
sliding_window_COR=5;

%Nbins
nbins=250;

%linestyle definition
linestyle_ESA='or';
linestyle_analytical='*b';
linestyle_threshold='+k';
linestyle_OCOG_ice='dc';
linestyle_STL='sy';
linestyle_ESA_smooth='-g';
linestyle_ESA_smooth_nocorr='-y';
linestyle_analytical_smooth='-g';
linestyle_threshold_smooth='-y';
linestyle_OCOG_ice_smooth='-m';


if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

%% -------------- Create the folder specific for the data -----------------------------
name_file_L2_ISD=strsplit(filename_L2_ISD,'/');
name_file_L2_ISD=name_file_L2_ISD(end);
name_file_L2_ESA=strsplit(filename_L2_ESA,'/');
name_file_L2_ESA=name_file_L2_ESA(end);
aux=strsplit(char(name_file_L2_ISD),'.'); aux=char(aux(1));
output_path=path_results_comparison;%strcat(path_results_comparison,aux,'\');
clear aux;
mkdir(output_path);
%mkdir([output_path 'subset\'])
%-------------- tEXT FILE COMPARISON --------------------------------------
fid = fopen(strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Evaluation.txt'), 'w');
fprintf(fid,'$---------------- Evaluation -----------------------------------$\n');
fprintf(fid,'ISD input file: '); fprintf(fid,'%s\n',char(name_file_L2_ISD));
fprintf(fid,'ESA input file: '); fprintf(fid,'%s\n',char(name_file_L2_ESA));



%% -------------- Constants definition ------------------------------------


%% --------------- Read isardSAT L2 product -------------------------------
%---------------- Geometry variables --------------------------------------
ISD_lat_surf=double(ncread(filename_L2_ISD,'lat_20_ku')).';
ISD_lon_surf=double(ncread(filename_L2_ISD,'lon_20_ku')).';
ISD_num_surfaces=length(ISD_lat_surf);
%---------------- Geophysical parameters ----------------------------------
%loop on different optional retrackers
for i_ret=1:length(retrackers)
    switch char(retrackers(i_ret))
        case {'ANALYTICAL','SAMOSA'}
            ISD_L2_SSH_analytical=double(ncread(filename_L2_ISD,'ssh_analytical_20_ku')).';
            ISD_L2_SWH_analytical=double(ncread(filename_L2_ISD,'swh_analytical_20_ku')).';
            ISD_L2_sigma0_analytical=double(ncread(filename_L2_ISD,'sig0_analytical_20_ku')).'; 
            ISD_L2_COR_analytical=double(ncread(filename_L2_ISD,'Pearson_corr_analytical_20_ku')).'; 
        case {'THRESHOLD'}
            %TBD
            ISD_L2_SSH_threshold=double(ncread(filename_L2_ISD,'ssh_threshold_20_ku')).';            
            ISD_L2_sigma0_threshold=double(ncread(filename_L2_ISD,'sig0_threshold_20_ku')).'; 
        case {'OCOG_ice'}
            %TBD
            ISD_L2_SSH_OCOG_ice=double(ncread(filename_L2_ISD,'ssh_ocog_ice_20_ku')).';            
            ISD_L2_sigma0_OCOG_ice=double(ncread(filename_L2_ISD,'sig0_ocog_ice_20_ku')).'; 
    end
end


%% ------------------- Read L2 ESA product --------------------------------
[~,CS2]=Cryo_L2_read(filename_L2_ESA);
s=size(CS2.MEA.surf_height_r1_20Hz);
records_db=s(2);
num_bursts_db=s(1);
ESA_num_surfaces=records_db*num_bursts_db;
%-------------- Geometry parameters ---------------------------------------
ESA_lat_surf=reshape(CS2.MEA.LAT_20Hz,[1,ESA_num_surfaces]);
ESA_lon_surf=reshape(CS2.MEA.LON_20Hz,[1,ESA_num_surfaces]);

%-------------- Geophysical parameters ------------------------------------
ESA_L2_SSH_r1=reshape(CS2.MEA.surf_height_r1_20Hz,[1,ESA_num_surfaces]);
ESA_L2_sigma0_r1=reshape(CS2.MEA.backsc_sig_r1_20Hz,[1,ESA_num_surfaces]);
%undo corrections from ESA
ESA_L2_SSH_r1_nocorr=ESA_L2_SSH_r1+reshape((ones(num_bursts_db,1))*(CS2.COR.total_ocean.'),[1,ESA_num_surfaces]);



%% -------------------- Read L2 STL product -------------------------------
if ~isempty(filename_L2_STL)
    STL_L2_SSH=double(ncread(filename_L2_STL,'ssh').');
    STL_L2_sigma0=double(ncread(filename_L2_STL,'Pu').');
    ncid = netcdf.open(filename_L2_STL,'NC_NOWRITE');
    varid = netcdf.inqVarID(ncid,'Lat');
    STL_L2_lat_surf=netcdf.getVar(ncid,varid);    
end




%% ------------------- FILTERING BY LATITUDE ------------------------------
% %Forcing the number of surfaces of the ISD and ESA product be the same
% %assuming first surface not contemplated in ESA product
% idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
% indices_lat_lon_zeros=find(idx_not_lat_lon_zeros==0);
% if ISD_num_surfaces <=ESA_num_surfaces
    % ISD_indexes_int=ones(1,ISD_num_surfaces);
    % ISD_indexes_int(1)=0;
    % ESA_indexes_int=zeros(1,ESA_num_surfaces);
    % ESA_indexes_int(1:ISD_num_surfaces-1)=1;        
    % for i_index=1:length(indices_lat_lon_zeros)
        % if (indices_lat_lon_zeros(i_index)+1)<=ISD_num_surfaces
            % ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
        % end
    % end
    % ISD_indexes_int=logical(ISD_indexes_int);
    % ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
% else
    % %the number of surfaces ESA limits
    % ISD_indexes_int=ones(1,ESA_num_surfaces);
    % ISD_indexes_int(1)=0;
    % ESA_indexes_int=ones(1,ESA_num_surfaces);    
    % for i_index=1:length(indices_lat_lon_zeros)
        % if (indices_lat_lon_zeros(i_index)+1)<=ISD_num_surfaces
            % ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
        % end
    % end
    % ISD_indexes_int=logical(ISD_indexes_int);
    % ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
% end

% Assuming L2 ISD processed filtering by geographical location
% Take the first closest surface of ESA as reference and force same number
% of surfaces as ISD
ISD_indexes_int=ones(1,ISD_num_surfaces);
[~,idx_min]=min(abs(ESA_lat_surf-ISD_lat_surf(1)));
ESA_indexes_int=zeros(1,ESA_num_surfaces);
ESA_indexes_int(idx_min(1):min(idx_min(1)+ISD_num_surfaces-1,ESA_num_surfaces))=1;
idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
indices_lat_lon_zeros=find(idx_not_lat_lon_zeros==0);
for i_index=1:length(indices_lat_lon_zeros)
    if (indices_lat_lon_zeros(i_index))<=ISD_num_surfaces
        ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
    end
end
ISD_indexes_int=logical(ISD_indexes_int);
ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;



ESA_num_surfaces_filtered=length(find(ESA_indexes_int)==1);
ISD_num_surfaces_filtered=length(find(ISD_indexes_int)==1);

idx_int_ESA=find(ESA_indexes_int);
idx_int_ISD=find(ISD_indexes_int);


%% --------------------------- OUTLIERS COMPUTATION -----------------------
idx_outliers_ESA_L2_SSH=zeros(1,ESA_num_surfaces_filtered);
idx_outliers_ESA_L2_sigma0=zeros(1,ESA_num_surfaces_filtered);
idx_int_analytical=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL')) | ~cellfun(@isempty,strfind(retrackers,'SAMOSA')), 1);
if ~isempty(idx_int_analytical)
    idx_outliers_ISD_L2_SSH_analytical      = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_analytical   = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_SWH_analytical      = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_COR_analytical      = zeros(1,ISD_num_surfaces_filtered);
    
end
idx_int_thres=find(~cellfun(@isempty,strfind(retrackers,'THRESHOLD')), 1);
if ~isempty(idx_int_thres)
    idx_outliers_ISD_L2_SSH_threshold       = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_threshold    = zeros(1,ISD_num_surfaces_filtered);
end
idx_int_ocog=find(~cellfun(@isempty,strfind(retrackers,'OCOG_ice')), 1);
if ~isempty(idx_int_ocog)
    idx_outliers_ISD_L2_SSH_OCOG_ice        = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_OCOG_ice        = zeros(1,ISD_num_surfaces_filtered);
end

if flag_outliers_removal
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ----------------------------- SSH ---------------------------------------
	IQR=iqr(ESA_L2_SSH_r1(ESA_indexes_int));
	percentiles = prctile(ESA_L2_SSH_r1(ESA_indexes_int),[outlier_percentil_low outlier_percentil_high]);
	idx_outliers_ESA_L2_SSH=idx_outliers_ESA_L2_SSH | (ESA_L2_SSH_r1(ESA_indexes_int)<(percentiles(1)-IQR_times*IQR) | ESA_L2_SSH_r1(ESA_indexes_int)>(percentiles(2)+IQR_times*IQR));
	%force the outliers as missing data or NaN
	ESA_L2_SSH_r1(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;
	ESA_L2_SSH_r1_nocorr(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;

	% ----------------------------- SIGMA0 ------------------------------------
	IQR=iqr(ESA_L2_sigma0_r1(ESA_indexes_int));
	percentiles = prctile(ESA_L2_sigma0_r1(ESA_indexes_int),[outlier_percentil_low outlier_percentil_high]);
	idx_outliers_ESA_L2_sigma0=idx_outliers_ESA_L2_sigma0 | (ESA_L2_sigma0_r1(ESA_indexes_int)<(percentiles(1)-IQR_times*IQR) | ESA_L2_sigma0_r1(ESA_indexes_int)>(percentiles(2)+IQR_times*IQR));
	%force the outliers as missing data or NaN
	ESA_L2_sigma0_r1(idx_int_ESA(idx_outliers_ESA_L2_sigma0))=NaN;
	

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ISD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i_ret=1:length(retrackers)
		switch char(retrackers(i_ret))
			case {'ANALYTICAL','SAMOSA'}
				% compute the errors w.r.t fitting on the data using a smooth
				% function
				%------------------- SSH ----------------------------------
				IQR=iqr(ISD_L2_SSH_analytical(ISD_indexes_int));
				percentiles = prctile(ISD_L2_SSH_analytical(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_SSH_analytical=idx_outliers_ISD_L2_SSH_analytical | (ISD_L2_SSH_analytical(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_SSH_analytical(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_SSH_analytical(idx_int_ISD(idx_outliers_ISD_L2_SSH_analytical))=NaN;

				%----------------- sigma0 ---------------------------------
				IQR=iqr(ISD_L2_sigma0_analytical(ISD_indexes_int));
				percentiles = prctile(ISD_L2_sigma0_analytical(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_sigma0_analytical=idx_outliers_ISD_L2_sigma0_analytical | (ISD_L2_sigma0_analytical(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_sigma0_analytical(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_sigma0_analytical(idx_int_ISD(idx_outliers_ISD_L2_sigma0_analytical))=NaN;
				
				%----------------- SWH ------------------------------------
				IQR=iqr(ISD_L2_SWH_analytical(ISD_indexes_int));
				percentiles = prctile(ISD_L2_SWH_analytical(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_SWH_analytical=idx_outliers_ISD_L2_SWH_analytical | (ISD_L2_SWH_analytical(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_SWH_analytical(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_SWH_analytical(idx_int_ISD(idx_outliers_ISD_L2_SWH_analytical))=NaN;
				
				% ---------- COR: pearson correlation coefficient ---------
				IQR=iqr(ISD_L2_COR_analytical(ISD_indexes_int));
				percentiles = prctile(ISD_L2_COR_analytical(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_COR_analytical=idx_outliers_ISD_L2_COR_analytical | (ISD_L2_COR_analytical(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_COR_analytical(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_COR_analytical(idx_int_ISD(idx_outliers_ISD_L2_COR_analytical))=NaN;
			  
			case {'THRESHOLD'}
				%-------------------- SSH ---------------------------------
				IQR=iqr(ISD_L2_SSH_threshold(ISD_indexes_int));
				percentiles = prctile(ISD_L2_SSH_threshold(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_SSH_threshold=idx_outliers_ISD_L2_SSH_threshold | (ISD_L2_SSH_threshold(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_SSH_threshold(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_SSH_threshold(idx_int_ISD(idx_outliers_ISD_L2_SSH_threshold))=NaN;
						  
				%------------------ sigma0 --------------------------------
				IQR=iqr(ISD_L2_sigma0_threshold(ISD_indexes_int));
				percentiles = prctile(ISD_L2_sigma0_threshold(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_sigma0_threshold=idx_outliers_ISD_L2_sigma0_threshold | (ISD_L2_sigma0_threshold(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_sigma0_threshold(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_sigma0_threshold(idx_int_ISD(idx_outliers_ISD_L2_sigma0_threshold))=NaN;
				
			case {'OCOG_ice'}
				%------------------------ SSH  ----------------------------
				IQR=iqr(ISD_L2_SSH_OCOG_ice(ISD_indexes_int));
				percentiles = prctile(ISD_L2_SSH_OCOG_ice(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_SSH_OCOG_ice=idx_outliers_ISD_L2_SSH_OCOG_ice | (ISD_L2_SSH_OCOG_ice(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_SSH_OCOG_ice(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));
				ISD_L2_SSH_OCOG_ice(idx_int_ISD(idx_outliers_ISD_L2_SSH_OCOG_ice))=NaN;
					 
				%-------------------- sigma0 ------------------------------
				IQR=iqr(ISD_L2_sigma0_OCOG_ice(ISD_indexes_int));
				percentiles = prctile(ISD_L2_sigma0_OCOG_ice(ISD_indexes_int),[outlier_percentil_low outlier_percentil_high]);
				idx_outliers_ISD_L2_sigma0_OCOG_ice=idx_outliers_ISD_L2_sigma0_OCOG_ice | (ISD_L2_sigma0_OCOG_ice(ISD_indexes_int)<(percentiles(1)-IQR_times*IQR) | ISD_L2_sigma0_OCOG_ice(ISD_indexes_int)>(percentiles(2)+IQR_times*IQR));            
				ISD_L2_sigma0_OCOG_ice(idx_int_ISD(idx_outliers_ISD_L2_sigma0_OCOG_ice))=NaN;
		end
	end
end



%% ------------------- COMPARISON -----------------------------------------         
%--------------------------------------------------------------------------
%------ Compute the errors around a fitting of the SSH & sigma0 -----------
%--------------------------------------------------------------------------
% A smoothing window is used
%--------------------------------------------------------------------------
%-----------------------------ESA------------------------------------------
for i_wfm=1:ESA_num_surfaces_filtered
    % compute the errors w.r.t fitting on the data using a smooth
    % function
    %------------------------------ SSH -----------------------------------
    start               = max(i_wfm-sliding_window_SSH,1);
    finish              = min(i_wfm+sliding_window_SSH,ESA_num_surfaces_filtered);
    % mean_win            = nanmean(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));
    % std_win             = nanstd(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));    
    % if ESA_L2_SSH_r1(idx_int_ESA(i_wfm))>(mean_win+threshold_std*std_win) || ESA_L2_SSH_r1(idx_int_ESA(i_wfm))<(mean_win-threshold_std*std_win)
        % idx_outliers_ESA_L2_SSH(i_wfm) = 1;
    % end    
    ESA_L2_SSH_r1_smoothed(i_wfm) = nanmean(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));
    ESA_L2_SSH_r1_nocorr_smoothed(i_wfm) = nanmean(ESA_L2_SSH_r1_nocorr(idx_int_ESA(start:finish)));
    
    %-------------------- sigma0 from L2 ESA product ----------------------
    start               = max(i_wfm-sliding_window_sigma0,1);
    finish              = min(i_wfm+sliding_window_sigma0,ESA_num_surfaces_filtered);
    % mean_win            = nanmean(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
    % std_win             = nanstd(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
    % if ESA_L2_sigma0_r1(idx_int_ESA(i_wfm))>(mean_win+threshold_std*std_win) || ESA_L2_sigma0_r1(idx_int_ESA(i_wfm))<(mean_win-threshold_std*std_win)
        % idx_outliers_ESA_L2_sigma0(i_wfm) = 1;
    % end
    ESA_L2_sigma0_r1_smoothed(i_wfm) = nanmean(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
end

%--------------------------------------------------------------------------
%--------------------------- ISD ------------------------------------------
idx_int_analytical=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL')) | ~cellfun(@isempty,strfind(retrackers,'SAMOSA')), 1);
if ~isempty(idx_int_analytical)
    ISD_L2_SSH_analytical_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_analytical_smoothed       = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_SWH_analytical_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_COR_analytical_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    
end
idx_int_thres=find(~cellfun(@isempty,strfind(retrackers,'THRESHOLD')), 1);
if ~isempty(idx_int_thres)
    ISD_L2_SSH_threshold_smoothed           = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_threshold_smoothed        = zeros(1,ISD_num_surfaces_filtered); 
end
idx_int_ocog=find(~cellfun(@isempty,strfind(retrackers,'OCOG_ice')), 1);
if ~isempty(idx_int_ocog)
    ISD_L2_SSH_OCOG_ice_smoothed            = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_OCOG_ice_smoothed         = zeros(1,ISD_num_surfaces_filtered);
end
for i_wfm=1:ISD_num_surfaces_filtered
    start_SSH               = max(i_wfm-sliding_window_SSH,1);
    finish_SSH              = min(i_wfm+sliding_window_SSH,ISD_num_surfaces_filtered);
    start_sigma0            = max(i_wfm-sliding_window_sigma0,1);
    finish_sigma0           = min(i_wfm+sliding_window_sigma0,ISD_num_surfaces_filtered);
    start_SWH               = max(i_wfm-sliding_window_SWH,1);
    finish_SWH              = min(i_wfm+sliding_window_SWH,ISD_num_surfaces_filtered);
    start_COR               = max(i_wfm-sliding_window_COR,1);
    finish_COR              = min(i_wfm+sliding_window_COR,ISD_num_surfaces_filtered);
    for i_ret=1:length(retrackers)
        switch char(retrackers(i_ret))
            case {'ANALYTICAL','SAMOSA'}
                % compute the errors w.r.t fitting on the data using a smooth
                % function
                %------------------- SSH ----------------------------------
                % mean_win            = nanmean(ISD_L2_SSH_analytical(idx_int_ISD(start_SSH:finish_SSH)));
                % std_win             = nanstd(ISD_L2_SSH_analytical(idx_int_ISD(start_SSH:finish_SSH)));
                % if ISD_L2_SSH_analytical(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SSH_analytical(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SSH_analytical(i_wfm) = 1;
                % end
                ISD_L2_SSH_analytical_smoothed(i_wfm)   = nanmean(ISD_L2_SSH_analytical(idx_int_ISD(start_SSH:finish_SSH)));
                
                %----------------- sigma0 ---------------------------------
                % mean_win            = nanmean(ISD_L2_sigma0_analytical(idx_int_ISD(start_sigma0:finish_sigma0)));
                % std_win             = nanstd(ISD_L2_sigma0_analytical(idx_int_ISD(start_sigma0:finish_sigma0)));
                % if ISD_L2_sigma0_analytical(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_sigma0_analytical(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_sigma0_analytical(i_wfm) = 1;
                % end
                ISD_L2_sigma0_analytical_smoothed(i_wfm) = nanmean(ISD_L2_sigma0_analytical(idx_int_ISD(start_sigma0:finish_sigma0)));
                %----------------- SWH ------------------------------------
                % mean_win            = nanmean(ISD_L2_SWH_analytical(idx_int_ISD(start_SWH:finish_SWH)));
                % std_win             = nanstd(ISD_L2_SWH_analytical(idx_int_ISD(start_SWH:finish_SWH)));
                % if ISD_L2_SWH_analytical(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SWH_analytical(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SWH_analytical(i_wfm) = 1;
                % end
                ISD_L2_SWH_analytical_smoothed(i_wfm) = nanmean(ISD_L2_SWH_analytical(idx_int_ISD(start_SWH:finish_SWH)));
                % ---------- COR: pearson correlation coefficient ---------
                % mean_win            = nanmean(ISD_L2_COR_analytical(idx_int_ISD(start_COR:finish_COR)));
                % std_win             = nanstd(ISD_L2_COR_analytical(idx_int_ISD(start_COR:finish_COR)));
                % if ISD_L2_COR_analytical(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_COR_analytical(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_COR_analytical(i_wfm) = 1;
                % end
                ISD_L2_COR_analytical_smoothed(i_wfm) = nanmean(ISD_L2_COR_analytical(idx_int_ISD(start_COR:finish_COR)));
            case {'THRESHOLD'}
                %-------------------- SSH ---------------------------------
                % mean_win            = nanmean(ISD_L2_SSH_threshold(idx_int_ISD(start_SSH:finish_SSH)));
                % std_win             = nanstd(ISD_L2_SSH_threshold(idx_int_ISD(start_SSH:finish_SSH)));
                % if ISD_L2_SSH_threshold(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SSH_threshold(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SSH_threshold(i_wfm) = 1;
                % end
                ISD_L2_SSH_threshold_smoothed(i_wfm)   = nanmean(ISD_L2_SSH_threshold(idx_int_ISD(start_SSH:finish_SSH)));
                %------------------ sigma0 --------------------------------
                % mean_win            = nanmean(ISD_L2_sigma0_threshold(idx_int_ISD(start_sigma0:finish_sigma0)));
                % std_win             = nanstd(ISD_L2_sigma0_threshold(idx_int_ISD(start_sigma0:finish_sigma0)));
                % if ISD_L2_sigma0_threshold(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_sigma0_threshold(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_sigma0_threshold(i_wfm) = 1;
                % end
                ISD_L2_sigma0_threshold_smoothed(i_wfm) = nanmean(ISD_L2_sigma0_threshold(idx_int_ISD(start_sigma0:finish_sigma0)));                
           case {'OCOG_ice'}
                %------------------------ SSH  ----------------------------
                % mean_win            = nanmean(ISD_L2_SSH_OCOG_ice(idx_int_ISD(start_SSH:finish_SSH)));
                % std_win             = nanstd(ISD_L2_SSH_OCOG_ice(idx_int_ISD(start_SSH:finish_SSH)));
                % if ISD_L2_SSH_OCOG_ice(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SSH_OCOG_ice(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SSH_OCOG_ice(i_wfm) = 1;
                % end
                ISD_L2_SSH_OCOG_ice_smoothed(i_wfm)   = nanmean(ISD_L2_SSH_OCOG_ice(idx_int_ISD(start_SSH:finish_SSH)));
                %-------------------- sigma0 ------------------------------
                % mean_win            = nanmean(ISD_L2_sigma0_OCOG_ice(idx_int_ISD(start_sigma0:finish_sigma0)));
                % std_win             = nanstd(ISD_L2_sigma0_OCOG_ice(idx_int_ISD(start_sigma0:finish_sigma0)));
                % if ISD_L2_sigma0_OCOG_ice(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_sigma0_OCOG_ice(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_sigma0_OCOG_ice(i_wfm) = 1;
                % end
                ISD_L2_sigma0_OCOG_ice_smoothed(i_wfm) = nanmean(ISD_L2_sigma0_OCOG_ice(idx_int_ISD(start_sigma0:finish_sigma0)));                                
        end
    end
end
   




%% --------------------------- SSH ----------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- --- --------------------------------------------$\n');
fprintf(fid,'$---------------- SSH --------------------------------------------$\n');
fprintf(fid,'$---------------- --- --------------------------------------------$\n');

%remove outliers:
ESA_L2_SSH_r1_filtered=ESA_L2_SSH_r1(ESA_indexes_int);
%idx_outliers_corr=ESA_L2_SSH_r1_filtered>(nanmean(ESA_L2_SSH_r1_filtered)+threshold_std*nanstd(ESA_L2_SSH_r1_filtered)) | ESA_L2_SSH_r1_filtered<(nanmean(ESA_L2_SSH_r1_filtered)-threshold_std*nanstd(ESA_L2_SSH_r1_filtered));

if ESA_num_surfaces_filtered==ISD_num_surfaces_filtered
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    fprintf(fid,'$---------------- RETRACKERS COMPARISON ESA -------------------------------------------$\n');
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    %----------------- analytical -------------------------------------
    if ~isempty(idx_int_analytical)
        fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_SSH_analytical_filtered=ISD_L2_SSH_analytical(ISD_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_SSH_analytical;
        res.SSH.ANALYTICAL.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_analytical_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISD (analytical-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL.RMSE_error_L2);
        res.SSH.ANALYTICAL.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_analytical_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISD (analytical-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL.mean_error_L2);        
    end
    if ~isempty(idx_int_thres)
        fprintf(fid,'$---------------- THRESHOLD  RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_SSH_threshold_filtered=ISD_L2_SSH_threshold(ISD_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_SSH_threshold;
        res.SSH.THRESHOLD.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_threshold_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISD (threshold-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.THRESHOLD.RMSE_error_L2);
        res.SSH.THRESHOLD.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_threshold_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISD (threshold-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.THRESHOLD.mean_error_L2);
    end
    if ~isempty(idx_int_ocog)
        fprintf(fid,'$---------------- OCOG_ice    RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_SSH_OCOG_ice_filtered=ISD_L2_SSH_OCOG_ice(ISD_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_SSH_OCOG_ice;
        res.SSH.OCOG_ice.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_OCOG_ice_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISD (OCOG_ice-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG_ice.RMSE_error_L2);
        res.SSH.OCOG_ice.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_OCOG_ice_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISD (OCOG_ice-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG_ice.mean_error_L2);
    end    
end
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical)
    fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.SSH.ANALYTICAL.mean_error_fitting=nanmean(ISD_L2_SSH_analytical_filtered(~idx_outliers_ISD_L2_SSH_analytical)-ISD_L2_SSH_analytical_smoothed(~idx_outliers_ISD_L2_SSH_analytical));
    fprintf(fid,'Mean error SSH fitting (analytical-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL.mean_error_fitting);
    res.SSH.ANALYTICAL.rmse_fitting=sqrt(nanmean((ISD_L2_SSH_analytical_filtered(~idx_outliers_ISD_L2_SSH_analytical)-ISD_L2_SSH_analytical_smoothed(~idx_outliers_ISD_L2_SSH_analytical)).^2));
    fprintf(fid,'RMSE SSH fitting (analytical-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL.rmse_fitting);
end
if ~isempty(idx_int_thres)
    fprintf(fid,'$---------------- THRESHOLD  RETRACKER -------------------------------------------$\n');
    res.SSH.THRESHOLD.mean_error_fitting=nanmean(ISD_L2_SSH_threshold_filtered(~idx_outliers_ISD_L2_SSH_threshold)-ISD_L2_SSH_threshold_smoothed(~idx_outliers_ISD_L2_SSH_threshold));
    fprintf(fid,'Mean error SSH fitting (threshold-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.THRESHOLD.mean_error_fitting);
    res.SSH.THRESHOLD.rmse_fitting=sqrt(nanmean((ISD_L2_SSH_threshold_filtered(~idx_outliers_ISD_L2_SSH_threshold)-ISD_L2_SSH_threshold_smoothed(~idx_outliers_ISD_L2_SSH_threshold)).^2));
    fprintf(fid,'RMSE SSH fitting (threshold-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.THRESHOLD.rmse_fitting);
end
if ~isempty(idx_int_ocog)
    fprintf(fid,'$---------------- OCOG_ICE   RETRACKER -------------------------------------------$\n');
    res.SSH.OCOG_ice.mean_error_fitting=nanmean(ISD_L2_SSH_OCOG_ice_filtered(~idx_outliers_ISD_L2_SSH_OCOG_ice)-ISD_L2_SSH_OCOG_ice_smoothed(~idx_outliers_ISD_L2_SSH_OCOG_ice));
    fprintf(fid,'Mean error SSH fitting (OCOG_ice-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG_ice.mean_error_fitting);
    res.SSH.OCOG_ice.rmse_fitting=sqrt(nanmean((ISD_L2_SSH_OCOG_ice_filtered(~idx_outliers_ISD_L2_SSH_OCOG_ice)-ISD_L2_SSH_OCOG_ice_smoothed(~idx_outliers_ISD_L2_SSH_OCOG_ice)).^2));
    fprintf(fid,'RMSE SSH fitting (OCOG_ice-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG_ice.rmse_fitting);
end

fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
res.SSH.ESA_L2.mean_error_fitting=nanmean(ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_SSH_r1_smoothed(~idx_outliers_ESA_L2_SSH));
fprintf(fid,'Mean error SSH fitting (L2-ESA) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ESA_L2.mean_error_fitting);
res.SSH.ESA_L2.rmse_fitting=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_SSH_r1_smoothed(~idx_outliers_ESA_L2_SSH)).^2));
fprintf(fid,'RMSE SSH fitting (L2-ESA) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ESA_L2.rmse_fitting);

%--------------------------------------------------------------------------
%-------------------- HISTOGRAM COMPUTATION -------------------------------
%--------------------------------------------------------------------------
if ~isempty(idx_int_analytical)
    %$---------------- ANALYTICAL RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_analytical]=hist(ISD_L2_SSH_analytical(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_analytical=counts./(sum(counts).*(centers_ISD_L2_SSH_analytical(2)-centers_ISD_L2_SSH_analytical(1)));
end
if ~isempty(idx_int_thres)
    %$---------------- THRESHOLD  RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_threshold]=hist(ISD_L2_SSH_threshold(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_threshold=counts./(sum(counts).*(centers_ISD_L2_SSH_threshold(2)-centers_ISD_L2_SSH_threshold(1)));
end
if ~isempty(idx_int_ocog)
    %$---------------- OCOG_ICE   RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_OCOG_ice]=hist(ISD_L2_SSH_OCOG_ice(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_OCOG_ice=counts./(sum(counts).*(centers_ISD_L2_SSH_OCOG_ice(2)-centers_ISD_L2_SSH_OCOG_ice(1)));
end
[counts,centers_ESA_L2_SSH_r1_filtered]=hist(ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH),nbins);
pdf_ESA_L2_SSH_r1_filtered=counts./(sum(counts).*(centers_ESA_L2_SSH_r1_filtered(2)-centers_ESA_L2_SSH_r1_filtered(1)));

%---------------------- ploting -------------------------------------------
%--------------------- Comparison smoothed vs non-smoothed ----------------
%------------------------ ESA ---------------------------------------------
figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_filtered,linestyle_ESA);
hold on;
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_smoothed,linestyle_ESA_smooth);
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr(ESA_indexes_int),'Marker','^','Color', [1 0.5 0.2]);
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr_smoothed,linestyle_ESA_smooth_nocorr);
text_in_textbox=strcat('Fitting ESA (corr.)--> RMSE: ',num2str(res.SSH.ESA_L2.rmse_fitting,'%.4g'),', Mean: ',num2str(res.SSH.ESA_L2.mean_error_fitting,'%.4g'));
legend('L2-ESA','L2-ESA - smoothed','L2-ESA (undoing corrections)','L2-ESA (undoing corrections) smooth');
xlabel('Latitude [deg.]'); ylabel('SSH [m]');
title(strcat('Comparison SSH: ESA (non-smoothed vs smoothed)'));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on','Units','points');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ESA_SSH_smoothed.png'))

%------------------------ ISD ---------------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};
if ~isempty(idx_int_thres)    
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold_filtered,linestyle_threshold);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold_smoothed,linestyle_threshold_smooth);
    legend_text=[legend_text,{'L2-ISD threshold','L2-ISD threshold (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting Threshold--> RMSE [m]: ',num2str(res.SSH.THRESHOLD.rmse_fitting,'%.4g'),', Mean [m]: ',num2str(res.SSH.THRESHOLD.mean_error_fitting,'%.4g'))}];
end
if ~isempty(idx_int_ocog)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG_ice_filtered,linestyle_OCOG_ice);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG_ice_smoothed,linestyle_OCOG_ice_smooth);
    legend_text=[legend_text,{'L2-ISD OCOG_ice','L2-ISD OCOG_ice (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting OCOG_ice--> RMSE [m]: ',num2str(res.SSH.OCOG_ice.rmse_fitting,'%.4g'),', Mean [m]: ',num2str(res.SSH.OCOG_ice.mean_error_fitting,'%.4g'))}];
end
if ~isempty(idx_int_analytical)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_filtered,linestyle_analytical);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_smoothed,linestyle_analytical_smooth);
    legend_text=[legend_text,{'L2-ISD analytical','L2-ISD analytical (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting Analytical--> RMSE [m]: ',num2str(res.SSH.ANALYTICAL.rmse_fitting,'%.4g'),', Mean [m]: ',num2str(res.SSH.ANALYTICAL.mean_error_fitting,'%.4g'))}];
end
legend(legend_text(~cellfun(@isempty,legend_text)));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
xlabel('Latitude [deg.]'); ylabel('SSH [m]');
title(strcat('Comparison SSH: isardSAT (non-smoothed vs smoothed)'));
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_SSH_smoothed.png'))

%---------------- Plot comparison ESA & ISD -------------------------------
figure;  
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),linestyle_ESA);
hold on;
legend_text={'L2-ESA'};
text_in_textbox={''};
if ~isempty(idx_int_thres)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold(ISD_indexes_int),linestyle_threshold);
    legend_text=[legend_text,{'L2-ISD threshold'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Comp. Threshold--> RMSE [m]: ',num2str(res.SSH.THRESHOLD.RMSE_error_L2,'%.4g'),', Mean [m]: ',num2str(res.SSH.THRESHOLD.mean_error_L2,'%.4g'))}];        
end
if ~isempty(idx_int_ocog)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG_ice(ISD_indexes_int),linestyle_OCOG_ice);
    legend_text=[legend_text,{'L2-ISD OCOG_ice'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Comp. OCOG_ice--> RMSE [m]: ',num2str(res.SSH.OCOG_ice.RMSE_error_L2,'%.4g'),', Mean [m]: ',num2str(res.SSH.OCOG_ice.mean_error_L2,'%.4g'))}];        
end
if ~isempty(idx_int_analytical)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical(ISD_indexes_int),linestyle_analytical);
    legend_text=[legend_text,{'L2-ISD analytical'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Comp. Analytical--> RMSE [m]: ',num2str(res.SSH.ANALYTICAL.RMSE_error_L2,'%.4g'),', Mean [m]: ',num2str(res.SSH.ANALYTICAL.mean_error_L2,'%.4g'))}];    
end

if ~isempty(filename_L2_STL)
    plot(STL_L2_lat_surf(ISD_indexes_int),STL_L2_SSH(ISD_indexes_int),linestyle_STL)
    legend_text=[legend_text,{'L2-STL'}];
    title(strcat('Comparison SSH: ESA & isardSAT & STARLAB'));
else
    title(strcat('Comparison SSH: ESA & isardSAT'));
end
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel('SSH [m]');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_SSH_ESA_ISD.png'))


%-------------------- HISTOGRAM PLOTING -----------------------------------
figure;  
plot(centers_ESA_L2_SSH_r1_filtered,pdf_ESA_L2_SSH_r1_filtered,strcat(linestyle_ESA,'-'));
hold on;
legend_text={'L2-ESA'};
if ~isempty(idx_int_thres)
    plot(centers_ISD_L2_SSH_threshold,pdf_ISD_L2_SSH_threshold,strcat(linestyle_threshold,'-'));
    legend_text=[legend_text,{'L2-ISD threshold'}];
end
if ~isempty(idx_int_ocog)
    plot(centers_ISD_L2_SSH_OCOG_ice,pdf_ISD_L2_SSH_OCOG_ice,strcat(linestyle_OCOG_ice,'-'));
    legend_text=[legend_text,{'L2-ISD OCOG_ice'}];
end
if ~isempty(idx_int_analytical)
    plot(centers_ISD_L2_SSH_analytical,pdf_ISD_L2_SSH_analytical,strcat(linestyle_analytical,'-'));
    legend_text=[legend_text,{'L2-ISD analytical'}];
end
title(strcat('Comparison SSH histograms: ESA & isardSAT'),'Interpreter','Tex');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('SSH [m]'); ylabel('PDF');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_PDF_SSH_ESA_ISD.png'))

    
%% ------------------------------ SIGMA0 ----------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$---------------- SIGMA0  --------------------------------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
ESA_L2_sigma0_r1_filtered=ESA_L2_sigma0_r1(ESA_indexes_int);
idx_outliers_ESA_L2_SSH=zeros(1,ESA_num_surfaces_filtered);
if ESA_num_surfaces_filtered==ISD_num_surfaces_filtered
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    fprintf(fid,'$---------------- RETRACKERS COMPARISON ESA -------------------------------------------$\n');
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    %----------------- analytical -------------------------------------
    if ~isempty(idx_int_analytical)
        fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_analytical_filtered=ISD_L2_sigma0_analytical(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_sigma0_analytical;
        res.SIGMA0.ANALYTICAL.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (analytical-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL.RMSE_error_L2);
        res.SIGMA0.ANALYTICAL.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (analytical-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL.mean_error_L2);        
    end
    if ~isempty(idx_int_thres)
        fprintf(fid,'$---------------- THRESHOLD  RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_threshold_filtered=ISD_L2_sigma0_threshold(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_sigma0_threshold;
        res.SIGMA0.THRESHOLD.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_threshold_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (threshold-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.THRESHOLD.RMSE_error_L2);
        res.SIGMA0.THRESHOLD.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_threshold_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (threshold-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.THRESHOLD.mean_error_L2);
    end
    if ~isempty(idx_int_ocog)
        fprintf(fid,'$---------------- OCOG_ice    RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_OCOG_ice_filtered=ISD_L2_sigma0_OCOG_ice(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_sigma0_OCOG_ice;
        res.SIGMA0.OCOG_ice.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_OCOG_ice_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (OCOG_ice-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG_ice.RMSE_error_L2);
        res.SIGMA0.OCOG_ice.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_OCOG_ice_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (OCOG_ice-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG_ice.mean_error_L2);
    end    
end
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical)
    fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.SIGMA0.ANALYTICAL.mean_error_fitting=nanmean(ISD_L2_sigma0_analytical(ISD_indexes_int)-ISD_L2_sigma0_analytical_smoothed);
    fprintf(fid,'Mean error sigma0 fitting (analytical-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL.mean_error_fitting);
    res.SIGMA0.ANALYTICAL.rmse_fitting=sqrt(nanmean((ISD_L2_sigma0_analytical(ISD_indexes_int)-ISD_L2_sigma0_analytical_smoothed).^2));
    fprintf(fid,'RMSE sigma0 fitting (analytical-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL.rmse_fitting);
end
if ~isempty(idx_int_thres)
    fprintf(fid,'$---------------- THRESHOLD  RETRACKER -------------------------------------------$\n');
    res.SIGMA0.THRESHOLD.mean_error_fitting=nanmean(ISD_L2_sigma0_threshold(ISD_indexes_int)-ISD_L2_sigma0_threshold_smoothed);
    fprintf(fid,'Mean error sigma0 fitting (threshold-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.THRESHOLD.mean_error_fitting);
    res.SIGMA0.THRESHOLD.rmse_fitting=sqrt(nanmean((ISD_L2_sigma0_threshold(ISD_indexes_int)-ISD_L2_sigma0_threshold_smoothed).^2));
    fprintf(fid,'RMSE sigma0 fitting (threshold-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.THRESHOLD.rmse_fitting);
end
if ~isempty(idx_int_ocog)
    fprintf(fid,'$---------------- OCOG_ICE   RETRACKER -------------------------------------------$\n');
    res.SIGMA0.OCOG_ice.mean_error_fitting=nanmean(ISD_L2_sigma0_OCOG_ice(ISD_indexes_int)-ISD_L2_sigma0_OCOG_ice_smoothed);
    fprintf(fid,'Mean error sigma0 fitting (OCOG_ice-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG_ice.mean_error_fitting);
    res.SIGMA0.OCOG_ice.rmse_fitting=sqrt(nanmean((ISD_L2_sigma0_OCOG_ice(ISD_indexes_int)-ISD_L2_sigma0_OCOG_ice_smoothed).^2));
    fprintf(fid,'RMSE sigma0 fitting (OCOG_ice-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG_ice.rmse_fitting);
end

fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
res.SIGMA0.ESA_L2.mean_error_fitting=nanmean(ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_sigma0_r1_smoothed(~idx_outliers_ESA_L2_SSH));
fprintf(fid,'Mean error sigma0 fitting (L2-ESA) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ESA_L2.mean_error_fitting);
res.SIGMA0.ESA_L2.rmse_fitting=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_sigma0_r1_smoothed(~idx_outliers_ESA_L2_SSH)).^2));
fprintf(fid,'RMSE sigma0 fitting (L2-ESA) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ESA_L2.rmse_fitting);

%--------------------------------------------------------------------------
%-------------------- HISTOGRAM COMPUTATION -------------------------------
%--------------------------------------------------------------------------
if ~isempty(idx_int_analytical)
    %$---------------- ANALYTICAL RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_analytical]=hist(ISD_L2_sigma0_analytical(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_analytical=counts./(sum(counts).*(centers_ISD_L2_sigma0_analytical(2)-centers_ISD_L2_sigma0_analytical(1)));
end
if ~isempty(idx_int_thres)
    %$---------------- THRESHOLD  RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_threshold]=hist(ISD_L2_sigma0_threshold(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_threshold=counts./(sum(counts).*(centers_ISD_L2_sigma0_threshold(2)-centers_ISD_L2_sigma0_threshold(1)));
end
if ~isempty(idx_int_ocog)
    %$---------------- OCOG_ICE   RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_OCOG_ice]=hist(ISD_L2_sigma0_OCOG_ice(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_OCOG_ice=counts./(sum(counts).*(centers_ISD_L2_sigma0_OCOG_ice(2)-centers_ISD_L2_sigma0_OCOG_ice(1)));
end
[counts,centers_ESA_L2_sigma0_r1_filtered]=hist(ESA_L2_sigma0_r1_filtered,nbins);
pdf_ESA_L2_sigma0_r1_filtered=counts./(sum(counts).*(centers_ESA_L2_sigma0_r1_filtered(2)-centers_ESA_L2_sigma0_r1_filtered(1)));


%---------------------- ploting -------------------------------------------
%--------------------- Comparison smoothed vs non-smoothed ----------------
%------------------------ ESA ---------------------------------------------
figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),linestyle_ESA);
hold on;
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1_smoothed,linestyle_ESA_smooth);
text_in_textbox={strcat('Fitting ESA (corr.)--> RMSE [dB]: ',num2str(res.SIGMA0.ESA_L2.rmse_fitting,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.ESA_L2.mean_error_fitting,'%.4g'))};
legend('L2-ESA','L2-ESA - smoothed');
xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
title(strcat('Comparison \sigma^0: ESA (non-smoothed vs smoothed)'),'Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ESA_SIGMA0_smoothed.png'))

%------------------------ ISD ---------------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};
if ~isempty(idx_int_analytical)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical(ISD_indexes_int),linestyle_analytical);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_smoothed,linestyle_analytical_smooth);
    legend_text=[legend_text,{'L2-ISD analytical','L2-ISD analytical (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting Analytical--> RMSE [dB]: ',num2str(res.SIGMA0.ANALYTICAL.rmse_fitting,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.ANALYTICAL.mean_error_fitting,'%.4g'))}];
end
if ~isempty(idx_int_thres)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold(ISD_indexes_int),linestyle_threshold);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold_smoothed,linestyle_threshold_smooth);
    legend_text=[legend_text,{'L2-ISD threshold','L2-ISD threshold (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting Threshold--> RMSE [dB]: ',num2str(res.SIGMA0.THRESHOLD.rmse_fitting,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.THRESHOLD.mean_error_fitting,'%.4g'))}];
    hold on;
end
if ~isempty(idx_int_ocog)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG_ice(ISD_indexes_int),linestyle_OCOG_ice);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG_ice_smoothed,linestyle_OCOG_ice_smooth);
    legend_text=[legend_text,{'L2-ISD OCOG_ice','L2-ISD OCOG_ice (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting OCOG_ice--> RMSE [dB]: ',num2str(res.SIGMA0.OCOG_ice.rmse_fitting,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.OCOG_ice.mean_error_fitting,'%.4g'))}];

end
legend(legend_text(~cellfun(@isempty,legend_text)));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
title(strcat('Comparison \sigma^0: isardSAT (non-smoothed vs smoothed)'),'Interpreter','Tex');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_SIGMA0_smoothed.png'))

%---------------- Plot comparison ESA & ISD -------------------------------
figure;  
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),linestyle_ESA);
hold on;
legend_text={'L2-ESA'};
text_in_textbox={''};
if ~isempty(idx_int_thres)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold(ISD_indexes_int),linestyle_threshold);
    legend_text=[legend_text,{'L2-ISD threshold'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Comp. Threshold--> RMSE [dB]: ',num2str(res.SIGMA0.THRESHOLD.RMSE_error_L2,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.THRESHOLD.mean_error_L2,'%.4g'))}];    
end
if ~isempty(idx_int_ocog)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG_ice(ISD_indexes_int),linestyle_OCOG_ice);
    legend_text=[legend_text,{'L2-ISD OCOG_ice'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Comp. OCOG_ice--> RMSE [dB]: ',num2str(res.SIGMA0.OCOG_ice.RMSE_error_L2,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.OCOG_ice.mean_error_L2,'%.4g'))}];        
end
if ~isempty(idx_int_analytical)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical(ISD_indexes_int),linestyle_analytical);
    legend_text=[legend_text,{'L2-ISD analytical'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Comp. Analytical--> RMSE [dB]: ',num2str(res.SIGMA0.ANALYTICAL.RMSE_error_L2,'%.4g'),', Mean [dB]: ',num2str(res.SIGMA0.ANALYTICAL.mean_error_L2,'%.4g'))}];
end

if ~isempty(filename_L2_STL)
    plot(STL_L2_lat_surf(ISD_indexes_int),STL_L2_sigma0(ISD_indexes_int),linestyle_STL)
    legend_text=[legend_text,{'L2-STL'}];
    title(strcat('Comparison \sigma^0: ESA & isardSAT & STARLAB'),'Interpreter','Tex');
else
    title(strcat('Comparison \sigma^0: ESA & isardSAT'),'Interpreter','Tex');
end
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_SIGMA0_ESA_ISD.png'))

%-------------------- HISTOGRAM PLOTING -----------------------------------
figure;  
plot(centers_ESA_L2_sigma0_r1_filtered,pdf_ESA_L2_sigma0_r1_filtered,strcat(linestyle_ESA,'-'));
hold on;
legend_text={'L2-ESA'};
if ~isempty(idx_int_thres)
    plot(centers_ISD_L2_sigma0_threshold,pdf_ISD_L2_sigma0_threshold,strcat(linestyle_threshold,'-'));
    legend_text=[legend_text,{'L2-ISD threshold'}];
end
if ~isempty(idx_int_ocog)
    plot(centers_ISD_L2_sigma0_OCOG_ice,pdf_ISD_L2_sigma0_OCOG_ice,strcat(linestyle_OCOG_ice,'-'));
    legend_text=[legend_text,{'L2-ISD OCOG_ice'}];
end
if ~isempty(idx_int_analytical)
    plot(centers_ISD_L2_sigma0_analytical,pdf_ISD_L2_sigma0_analytical,strcat(linestyle_analytical,'-'));
    legend_text=[legend_text,{'L2-ISD analytical'}];
end
title(strcat('Comparison \sigma^0 histograms: ESA & isardSAT'),'Interpreter','Tex');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('\sigma^0 [dB]','Interpreter','Tex'); ylabel('PDF');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_PDF_sigma0_ESA_ISD.png'))

%% ------------------------------ SWH -------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$----------------   SWH   --------------------------------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical)
    fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.SWH.ANALYTICAL.mean_error_fitting=nanmean(ISD_L2_SWH_analytical(ISD_indexes_int)-ISD_L2_SWH_analytical_smoothed);
    fprintf(fid,'Mean error SWH fitting (analytical-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL.mean_error_fitting);
    res.SWH.ANALYTICAL.rmse_fitting=sqrt(nanmean((ISD_L2_SWH_analytical(ISD_indexes_int)-ISD_L2_SWH_analytical_smoothed).^2));
    fprintf(fid,'RMSE SWH fitting (analytical-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL.rmse_fitting);
end
%--------------------------------------------------------------------------
%-------------------- HISTOGRAM COMPUTATION -------------------------------
%--------------------------------------------------------------------------
if ~isempty(idx_int_analytical)
    %$---------------- ANALYTICAL RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SWH_analytical]=hist(ISD_L2_SWH_analytical(ISD_indexes_int),nbins);
    pdf_ISD_L2_SWH_analytical=counts./(sum(counts).*(centers_ISD_L2_SWH_analytical(2)-centers_ISD_L2_SWH_analytical(1)));
end
%----------------------------- ploting ------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};
if ~isempty(idx_int_analytical)
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_SWH_analytical(ISD_indexes_int)),linestyle_analytical);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_SWH_analytical_smoothed),linestyle_analytical_smooth);
    legend_text=[legend_text,{'L2-ISD analytical','L2-ISD analytical (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting Analytical--> RMSE [m]: ',num2str(res.SWH.ANALYTICAL.rmse_fitting,'%.4g'),', Mean [m]: ',num2str(res.SWH.ANALYTICAL.mean_error_fitting,'%.4g'))}];
end
legend(legend_text(~cellfun(@isempty,legend_text)));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
xlabel('Latitude [deg.]'); ylabel('SWH [m]','Interpreter','Tex');
title(strcat('Comparison SWH: isardSAT (non-smoothed vs smoothed)'),'Interpreter','Tex');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_SWH_smoothed.png'))
%-------------------- HISTOGRAM PLOTING -----------------------------------
figure;  
legend_text={''};
if ~isempty(idx_int_analytical)
    plot(centers_ISD_L2_SWH_analytical,pdf_ISD_L2_SWH_analytical,strcat(linestyle_analytical,'-'));
    hold on;
    legend_text=[legend_text,{'L2-ISD analytical'}];
end
title(strcat('Histogram SWH: isardSAT'),'Interpreter','Tex');
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('SWH [m]'); ylabel('PDF');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_PDF_SWH_ISD.png'))

%% ---------- COR:PEARSON CORRELATION COEFFICIENT -------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$--------   COR:PEARSON CORRELATION COEFFICIENT   --------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical)
    fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.COR.ANALYTICAL.mean=nanmean(ISD_L2_COR_analytical(ISD_indexes_int));
    fprintf(fid,'Mean value Pearson Correlation Coefficient (analytical-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL.mean);
    res.COR.ANALYTICAL.mean_error_fitting=nanmean(ISD_L2_COR_analytical(ISD_indexes_int)-ISD_L2_COR_analytical_smoothed);
    fprintf(fid,'Mean error Pearson Correlation Coefficient fitting (analytical-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL.mean_error_fitting);
    res.COR.ANALYTICAL.rmse_fitting=sqrt(nanmean((ISD_L2_COR_analytical(ISD_indexes_int)-ISD_L2_COR_analytical_smoothed).^2));
    fprintf(fid,'RMSE Pearson Correlation Coefficient fitting (analytical-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL.rmse_fitting);
end
%----------------------------- ploting ------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};
if ~isempty(idx_int_analytical)
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_COR_analytical(ISD_indexes_int)),linestyle_analytical);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_COR_analytical_smoothed),linestyle_analytical_smooth);
    legend_text=[legend_text,{'L2-ISD analytical','L2-ISD analytical (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        {strcat('Fitting Analytical--> RMSE [%]: ',num2str(res.COR.ANALYTICAL.rmse_fitting,'%.4g'),', Mean [%]: ',num2str(res.SWH.ANALYTICAL.mean_error_fitting,'%.4g'))}];
end
legend(legend_text(~cellfun(@isempty,legend_text)));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1),0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
xlabel('Latitude [deg.]'); ylabel('\rho_{pearson} [%]','Interpreter','Tex');
title(strcat('Comparison \rho_{pearson}: isardSAT (non-smoothed vs smoothed)'),'Interpreter','Tex');
print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_COR_RHO_smoothed.png'))

%% --------------------------- Saving Information -------------------------
fclose(fid);
save(strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'L2_Evaluation.mat'),'res');

end