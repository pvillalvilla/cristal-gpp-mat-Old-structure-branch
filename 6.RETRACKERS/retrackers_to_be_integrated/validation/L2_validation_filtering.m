%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VALIDATION L2/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this algorithm is to cross-check the L2 products:
% isardSAT and ESA for CryoSat-2 data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res]=L2_validation_filtering(filename_L2_ISD,filename_L2_ESA,path_results_comparison,varargin)
global generate_kml
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<3 || nargin>(3+8*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParamValue('filename_L2_STL',{''},@(x)ischar(x));
p.addParamValue('retrackers',{''},@(x)iscellstr(x)); %is a cell array with different names of different retrackers in L2 product
p.addParamValue('figures_visible',0);
p.addParamValue('flag_outliers_removal',0);
p.addParamValue('type_outliers_removal','percentiles');
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('geo_mask',[]);
p.addParamValue('annotation_box_active',1);
filename_mask_KML=p.Results.filename_mask_KML;
p.parse(varargin{:});
filename_L2_STL=char(p.Results.filename_L2_STL);
retrackers=p.Results.retrackers;
figures_visible=p.Results.figures_visible;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=char(p.Results.type_outliers_removal);
sh_name_nc=p.Results.sh_name_nc;
geo_mask=p.Results.geo_mask;
annotation_box_active=p.Results.annotation_box_active;
clear p;


%close all;
set_default_plot;
%figures_visible=1;

%removal of outliers in ESA SSH retrievals L2
threshold_std=3.0;

%removal of outliers 
%-------using percentile & IQR (Interquartile Range)
% outliers data<(percentil_low-IQR_times*IQR) | data>(percentil_high+IQR_times*IQR)
IQR_times=1.5; %number of IQR 
outlier_percentil_low=25.0;
outlier_percentil_high=75.0;
%--------using hampel filter
hampel_wind=3;% size of window half size
hampel_sigma=3; %number of std deviations to which a sample differ from local median



%sliding window definitions: smoothing function used as fitting
% number of samples at each side of current position
sliding_window_SSH=10;
sliding_window_SWH=10;
sliding_window_sigma0=10;
sliding_window_COR=10;

%Nbins
nbins=250;

%linestyle definition
linestyle_ESA='or';
linestyle_analytical_SWH_MSSfixed='*b';
linestyle_analytical_MSS_SWHfixed='^g';
linestyle_threshold='+k';
linestyle_OCOG='dc';
linestyle_STL='sg';
linestyle_ESA_smooth='-g';
linestyle_ESA_smooth_nocorr='-y';
linestyle_analytical_SWH_MSSfixed_smooth='-c';
linestyle_analytical_MSS_SWHfixed_smooth='-r';
linestyle_threshold_smooth='-y';
linestyle_OCOG_smooth='-c';
linestyle_STL_smooth='-m';

%hardcoded to ensure annotation box aligned with axis for Fontsize 30
if ~isempty(filename_L2_STL)
    y_add_textbox=0.07;
else
    y_add_textbox=0.015;
end



if figures_visible
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end

%% -------------- Create the folder specific for the data -----------------------------
name_file_L2_ISD=strsplit(filename_L2_ISD,'/');
name_file_L2_ISD=char(name_file_L2_ISD(end));
name_file_L2_ISD=name_file_L2_ISD(17:17+30);
name_file_L2_ESA=strsplit(filename_L2_ESA,'/');
name_file_L2_ESA=name_file_L2_ESA(end);
%aux=strsplit(char(name_file_L2_ISD),'.'); aux=char(aux(1));
output_path=strcat(path_results_comparison,name_file_L2_ISD,'/');%path_results_comparison;
clear aux;
mkdir(output_path);
%mkdir([output_path 'subset\'])
%-------------- tEXT FILE COMPARISON --------------------------------------
%fid = fopen(strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Evaluation.txt'), 'w');
fid = fopen(strcat(output_path,name_file_L2_ISD,'_Evaluation.txt'), 'w');
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
        case {'ANALYTICAL_SWH'}
            ISD_L2_SSH_analytical_SWH_MSSfixed=double(ncread(filename_L2_ISD,strcat(sh_name_nc,'_analytical_SWH_MSSfixed_20_ku'))).';
            ISD_L2_SWH_analytical_SWH_MSSfixed=double(ncread(filename_L2_ISD,'swh_analytical_SWH_MSSfixed_20_ku')).';
            ISD_L2_sigma0_analytical_SWH_MSSfixed=double(ncread(filename_L2_ISD,'sig0_analytical_SWH_MSSfixed_20_ku')).'; 
            ISD_L2_COR_analytical_SWH_MSSfixed=double(ncread(filename_L2_ISD,'Pearson_corr_analytical_SWH_MSSfixed_20_ku')).'; 
            %title_name_SWH_MSSfixed='analytical-fit-SWH-and-MSS-fixed';
            title_name_SWH_MSSfixed='SWH fit';
        case {'ANALYTICAL_MSS'}
            ISD_L2_SSH_analytical_MSS_SWHfixed=double(ncread(filename_L2_ISD,strcat(sh_name_nc,'_analytical_MSS_SWHfixed_20_ku'))).';
%            ISD_L2_SWH_analytical_MSS_SWHfixed=double(ncread(filename_L2_ISD,'swh_analytical_MSS_SWHfixed_20_ku')).';
            ISD_L2_sigma0_analytical_MSS_SWHfixed=double(ncread(filename_L2_ISD,'sig0_analytical_MSS_SWHfixed_20_ku')).';
            ISD_L2_COR_analytical_MSS_SWHfixed=double(ncread(filename_L2_ISD,'Pearson_corr_analytical_MSS_SWHfixed_20_ku')).';
            title_name_MSS_SWHfixed='MSS fit';
        case {'THRESHOLD'}
            %TBD
            ISD_L2_SSH_threshold=double(ncread(filename_L2_ISD,strcat(sh_name_nc,'_threshold_20_ku'))).';            
            ISD_L2_sigma0_threshold=double(ncread(filename_L2_ISD,'sig0_threshold_20_ku')).'; 
        case {'OCOG'}
            %TBD
            ISD_L2_SSH_OCOG=double(ncread(filename_L2_ISD,strcat(sh_name_nc,'_OCOG_20_ku'))).';            
            ISD_L2_sigma0_OCOG=double(ncread(filename_L2_ISD,'sig0_OCOG_20_ku')).'; 
    end
end

%write the surfaces on a KML
if generate_kml
    ISD_alt_surf=ISD_lat_surf;
    ISD_alt_surf(:)=0;
    lla2kmlWaveforms_noimage(strcat(output_path,name_file_L2_ISD,'_Geolocated_track.kml'),name_file_L2_ISD,ISD_lat_surf,ISD_lon_surf,ISD_alt_surf, '.');
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
    STL_L2_SSH=double(ncread(filename_L2_STL,'ssh_corr').');
    STL_L2_sigma0=double(ncread(filename_L2_STL,'sigma0_L2').');
    STL_L2_SWH=double(ncread(filename_L2_STL,'swh').');
    STL_L2_lat_surf=double(ncread(filename_L2_STL,'Latitude').');   
    STL_L2_FIT=100.0-double(ncread(filename_L2_STL,'gof').'); %READING MISFIT COVNERT TO EQUIVALENT COR
end


%% ------------------- FILTERING BY LATITUDE ------------------------------
%--------------------------------------------------------------------------
%Forcing the number of surfaces of the ISD and ESA product be the same
%assuming first surface not contemplated in ESA product
idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
indices_lat_lon_zeros=find(idx_not_lat_lon_zeros==0);
if ISD_num_surfaces <=ESA_num_surfaces
    ISD_indexes_int=ones(1,ISD_num_surfaces);
    ISD_indexes_int(1)=0;
    ESA_indexes_int=zeros(1,ESA_num_surfaces);
    ESA_indexes_int(1:ISD_num_surfaces-1)=1;        
    for i_index=1:length(indices_lat_lon_zeros)
        if (indices_lat_lon_zeros(i_index)+1)<=ISD_num_surfaces
            ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
        end
    end
    ISD_indexes_int=logical(ISD_indexes_int);
    ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
else
    %the number of surfaces ESA limits
    ISD_indexes_int=zeros(1,ISD_num_surfaces);
    ISD_indexes_int(2:ESA_num_surfaces+1)=1;
    ISD_indexes_int(1)=0;
    ESA_indexes_int=ones(1,ESA_num_surfaces);    
    for i_index=1:length(indices_lat_lon_zeros)
        if (indices_lat_lon_zeros(i_index)+1)<=ISD_num_surfaces
            ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
        end
    end
    ISD_indexes_int=logical(ISD_indexes_int);
    ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
end

% %--------------------------------------------------------------------------
% 
% %--------------------------------------------------------------------------
% %---------------- geohraphical filtering ----------------------------------
% %reference is alawys ISD
% % ---- Checking whether the track within geomask if available ---------
% if ~isempty(geo_mask)
%     ISD_lon_surf_bis=ISD_lon_surf;
% %     idx_lt_0= ISD_lon_surf_bis<0;
% %     %longitudes +- values (-180,180)
% %     if any(idx_lt_0)
% %         ISD_lon_surf_bis(idx_lt_0)=ISD_lon_surf_bis(idx_lt_0)+360.0;
% %     end
% % 
% %     clear idx_lt_0;
%     
%     idx_int_geo=inpolygon(ISD_lon_surf_bis,ISD_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));    
%     if ~any(idx_int_geo)
%         disp(strcat('Track,',{' '},name_file_L2_ISD,{' '},'outside the limits of the geographical mask'))
%         return;
%     end
%     ESA_indexes_int(ESA_indexes_int==1)=idx_int_geo(ISD_indexes_int==1);
%     ISD_indexes_int=ISD_indexes_int & idx_int_geo;
% end
% %--------------------------------------------------------------------------

% %--------------------------------------------------------------------------
% %---------------- geohraphical filtering ----------------------------------
% %reference is alawys ISD
% % ---- Checking whether the track within geomask if available ---------
% ISD_indexes_int=ones(1,ISD_num_surfaces);
% ESA_indexes_int=ones(1,ESA_num_surfaces);
% if ~isempty(geo_mask)   
%     %ISD
%     idx_int_geo_ISD=inpolygon(ISD_lon_surf,ISD_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));    
%     if ~any(idx_int_geo_ISD)
%         disp(strcat('Track ISD,',{' '},name_file_L2_ISD,{' '},'outside the limits of the geographical mask'))
%         return;
%     end    
%     ISD_indexes_int=ISD_indexes_int & idx_int_geo_ISD;
%     
%     %ESA
%     idx_int_geo_ESA=inpolygon(ESA_lon_surf,ESA_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));    
%     if ~any(idx_int_geo_ESA)
%         disp(strcat('Track ESA,',{' '},name_file_L2_ESA,{' '},'outside the limits of the geographical mask'))
%         return;
%     end        
%     ESA_indexes_int=ESA_indexes_int & idx_int_geo_ESA;
% end
% %--------------------------------------------------------------------------

% %--------------------------------------------------------------------------
% % Assuming L2 ISD processed filtering by geographical location
% % Take the first closest surface of ESA as reference and force same number
% % of surfaces as ISD
% ISD_indexes_int=ones(1,ISD_num_surfaces);
% [~,idx_min]=min(abs(ESA_lat_surf-ISD_lat_surf(1)));
% ESA_indexes_int=zeros(1,ESA_num_surfaces);
% ESA_indexes_int(idx_min(1):min(idx_min(1)+ISD_num_surfaces-1,ESA_num_surfaces))=1;
% idx_not_lat_lon_zeros=~(ESA_lat_surf==0 & ESA_lon_surf==0);
% indices_lat_lon_zeros=find(idx_not_lat_lon_zeros==0);
% for i_index=1:length(indices_lat_lon_zeros)
%     if (indices_lat_lon_zeros(i_index))<=ISD_num_surfaces
%         ISD_indexes_int(indices_lat_lon_zeros(i_index)+1)=0;
%     end
% end
% ISD_indexes_int=logical(ISD_indexes_int);
% ESA_indexes_int=logical(ESA_indexes_int) & idx_not_lat_lon_zeros;
% %--------------------------------------------------------------------------



ESA_num_surfaces_filtered=length(find(ESA_indexes_int)==1);
ISD_num_surfaces_filtered=length(find(ISD_indexes_int)==1);

idx_int_ESA=find(ESA_indexes_int);
idx_int_ISD=find(ISD_indexes_int);

idx_int_analytical_SWH_MSSfixed=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL_SWH')), 1);
idx_int_analytical_MSS_SWHfixed=find(~cellfun(@isempty,strfind(retrackers,'ANALYTICAL_MSS')), 1);
idx_int_thres=find(~cellfun(@isempty,strfind(retrackers,'THRESHOLD')), 1);
idx_int_ocog=find(~cellfun(@isempty,strfind(retrackers,'OCOG')), 1);

%% --------------------------- APPLY THE GEO CORR ACCORING TO L2 ESA ------
if ISD_num_surfaces_filtered==ESA_num_surfaces_filtered
    %reads the geophysical correction for L2 ISD and ESA and applies the missing or removes the non-valid corrections as per L2 ESA
    script_geo_corr_a_ESA_L2
    %only valid for the same number of surfaces (forced)
end



%% --------------------------- OUTLIERS COMPUTATION -----------------------
idx_outliers_ESA_L2_SSH=zeros(1,ESA_num_surfaces_filtered);
idx_outliers_ESA_L2_sigma0=zeros(1,ESA_num_surfaces_filtered);

if ~isempty(idx_int_analytical_SWH_MSSfixed)
    idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed      = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed   = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed      = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_COR_analytical_SWH_MSSfixed      = zeros(1,ISD_num_surfaces_filtered);
    
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed      = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_analytical_MSS_SWHfixed   = zeros(1,ISD_num_surfaces_filtered);    
    idx_outliers_ISD_L2_COR_analytical_MSS_SWHfixed      = zeros(1,ISD_num_surfaces_filtered);
    
end

if ~isempty(idx_int_thres)
    idx_outliers_ISD_L2_SSH_threshold       = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_threshold    = zeros(1,ISD_num_surfaces_filtered);
end

if ~isempty(idx_int_ocog)
    idx_outliers_ISD_L2_SSH_OCOG        = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_ISD_L2_sigma0_OCOG        = zeros(1,ISD_num_surfaces_filtered);
end

if ~isempty(filename_L2_STL)
    %should have the same # surfaces as ISD generated from L1B ISD
    idx_outliers_STL_L2_SSH        = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_STL_L2_SWH        = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_STL_L2_sigma0     = zeros(1,ISD_num_surfaces_filtered);
    idx_outliers_STL_L2_FIT        = zeros(1,ISD_num_surfaces_filtered);
end

if flag_outliers_removal
	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ----------------------------- SSH ---------------------------------------
    switch type_outliers_removal
        case 'percentiles'
            [ESA_L2_SSH_r1(ESA_indexes_int),idx_outliers_ESA_L2_SSH]=outliers_by_percentiles(ESA_L2_SSH_r1(ESA_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
            ESA_L2_SSH_r1_nocorr(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;            
        case 'hampel'
            %[ESA_L2_SSH_r1(ESA_indexes_int)] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),hampel_wind,hampel_sigma);            
            [~,idx_outliers_ESA_L2_SSH] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),hampel_wind,hampel_sigma);            
            ESA_L2_SSH_r1(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;
            %[ESA_L2_SSH_r1_nocorr(ESA_indexes_int),indx] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr(ESA_indexes_int),hampel_wind,hampel_sigma);                        
            ESA_L2_SSH_r1_nocorr(idx_int_ESA(idx_outliers_ESA_L2_SSH))=NaN;
    end


	% ----------------------------- SIGMA0 ------------------------------------
    switch type_outliers_removal
        case 'percentiles'
            [ESA_L2_sigma0_r1(ESA_indexes_int),idx_outliers_ESA_L2_sigma0]=outliers_by_percentiles(ESA_L2_sigma0_r1(ESA_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);            
        case 'hampel'
            %[ESA_L2_sigma0_r1(ESA_indexes_int),indx] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),hampel_wind,hampel_sigma);
            [~,idx_outliers_ESA_L2_sigma0] = hampel(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),hampel_wind,hampel_sigma);
            ESA_L2_sigma0_r1(idx_int_ESA(idx_outliers_ESA_L2_sigma0))=NaN;
    end

	

	% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ISD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i_ret=1:length(retrackers)
		switch char(retrackers(i_ret))
			case {'ANALYTICAL_SWH'}
                switch type_outliers_removal
                    case 'percentiles'
                        % compute the errors w.r.t fitting on the data using a smooth
                        % function
                        %------------------- SSH ----------------------------------
                        [ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed]=outliers_by_percentiles(ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        
                        %----------------- sigma0 ---------------------------------
                        [ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed]=outliers_by_percentiles(ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        
                        %----------------- SWH ------------------------------------
                        [ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed]=outliers_by_percentiles(ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        
                        % ---------- COR: pearson correlation coefficient ---------
                        [ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int),idx_outliers_ISD_L2_COR_analytical_SWH_MSSfixed]=outliers_by_percentiles(ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                    case 'hampel'
                        %------------------ SSH ---------------------------
                        %[ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);
                        [~,idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);
                        ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed))=NaN;
                        %----------------- sigma0 ---------------------------------
                        %[ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                               
                        [~,idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                               
                        ISD_L2_sigma0_analytical_SWH_MSSfixed(idx_int_ISD(idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed))=NaN;
                        %----------------- SWH ------------------------------------
                        %[ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                                              
                        [~,idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                                               
                        ISD_L2_SWH_analytical_SWH_MSSfixed(idx_int_ISD(idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed))=NaN;
                        % ---------- COR: pearson correlation coefficient ---------
                        %[ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                     
                        [~,idx_outliers_ISD_L2_COR_analytical_SWH_MSSfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                     
                        ISD_L2_COR_analytical_SWH_MSSfixed(idx_int_ISD(idx_outliers_ISD_L2_COR_analytical_SWH_MSSfixed))=NaN;
                end

            case {'ANALYTICAL_MSS'}
                switch type_outliers_removal
                    case 'percentiles'
                        %------------------- SSH ----------------------------------
                        [ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed]=outliers_by_percentiles(ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        %----------------- sigma0 ---------------------------------
                        [ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),idx_outliers_ISD_L2_sigma0_analytical_MSS_SWHfixed]=outliers_by_percentiles(ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        % ---------- COR: pearson correlation coefficient ---------
                        [ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int),idx_outliers_ISD_L2_COR_analytical_MSS_SWHfixed]=outliers_by_percentiles(ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                    case 'hampel'
                        %------------------- SSH ----------------------------------
                        %[ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                        
                        [~,idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                        
                        ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed))=NaN;
                        %----------------- sigma0 ---------------------------------
                        %[ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                
                        [~,idx_outliers_ISD_L2_sigma0_analytical_MSS_SWHfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                
                        ISD_L2_sigma0_analytical_MSS_SWHfixed(idx_int_ISD(idx_outliers_ISD_L2_sigma0_analytical_MSS_SWHfixed))=NaN;
                        % ---------- COR: pearson correlation coefficient ---------
                        %[ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                        
                        [~,idx_outliers_ISD_L2_COR_analytical_MSS_SWHfixed] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                        
                        ISD_L2_COR_analytical_MSS_SWHfixed(idx_int_ISD(idx_outliers_ISD_L2_COR_analytical_MSS_SWHfixed))=NaN;
                end
                % compute the errors w.r.t fitting on the data using a smooth
                % function
                
            case {'THRESHOLD'}
                switch type_outliers_removal
                    case 'percentiles'
                        %-------------------- SSH ---------------------------------
                        [ISD_L2_SSH_threshold(ISD_indexes_int),idx_outliers_ISD_L2_SSH_threshold]=outliers_by_percentiles(ISD_L2_SSH_threshold(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        %------------------ sigma0 --------------------------------
                        [ISD_L2_sigma0_threshold(ISD_indexes_int),idx_outliers_ISD_L2_sigma0_threshold]=outliers_by_percentiles(ISD_L2_sigma0_threshold(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                    case 'hampel'
                        %-------------------- SSH ---------------------------------
                        %[ISD_L2_SSH_threshold(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold(ISD_indexes_int),hampel_wind,hampel_sigma);                                                
                        [~,idx_outliers_ISD_L2_SSH_threshold] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold(ISD_indexes_int),hampel_wind,hampel_sigma);                                                
                        ISD_L2_SSH_threshold(idx_int_ISD(idx_outliers_ISD_L2_SSH_threshold))=NaN;
                        %------------------ sigma0 --------------------------------
                        %[ISD_L2_sigma0_threshold(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold(ISD_indexes_int),hampel_wind,hampel_sigma);                        
                        [~,idx_outliers_ISD_L2_sigma0_threshold] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold(ISD_indexes_int),hampel_wind,hampel_sigma);                        
                        ISD_L2_sigma0_threshold(idx_int_ISD(idx_outliers_ISD_L2_sigma0_threshold))=NaN;
                end

				
			case {'OCOG'}
                switch type_outliers_removal
                    case 'percentiles'
                        %------------------------ SSH  ----------------------------
                        [ISD_L2_SSH_OCOG(ISD_indexes_int),idx_outliers_ISD_L2_SSH_OCOG]=outliers_by_percentiles(ISD_L2_SSH_OCOG(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        %-------------------- sigma0 ------------------------------
                        [ISD_L2_sigma0_OCOG(ISD_indexes_int),idx_outliers_ISD_L2_sigma0_OCOG]=outliers_by_percentiles(ISD_L2_sigma0_OCOG(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                    case 'hampel'
                        %------------------------ SSH  ----------------------------
                        %[ISD_L2_SSH_OCOG(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG(ISD_indexes_int),hampel_wind,hampel_sigma);                                                 
                        [~,idx_outliers_ISD_L2_SSH_OCOG] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG(ISD_indexes_int),hampel_wind,hampel_sigma);                                                 
                        ISD_L2_SSH_OCOG(idx_int_ISD(idx_outliers_ISD_L2_SSH_OCOG))=NaN;
                        %-------------------- sigma0 ------------------------------
                        %[ISD_L2_sigma0_OCOG(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                         
                        [~,idx_outliers_ISD_L2_sigma0_OCOG] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                         
                        ISD_L2_sigma0_OCOG(idx_int_ISD(idx_outliers_ISD_L2_sigma0_OCOG))=NaN;
                end

		end
    end
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&& STARLAB &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if ~isempty(filename_L2_STL)        
        switch type_outliers_removal
                 case 'percentiles'
                        % compute the errors w.r.t fitting on the data using a smooth
                        % function
                        %------------------- SSH ----------------------------------
                        [STL_L2_SSH(ISD_indexes_int),idx_outliers_STL_L2_SSH]=outliers_by_percentiles(STL_L2_SSH(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        
                        %----------------- sigma0 ---------------------------------
                        [STL_L2_sigma0(ISD_indexes_int),idx_outliers_STL_L2_sigma0]=outliers_by_percentiles(STL_L2_sigma0(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        
                        %----------------- SWH ------------------------------------
                        [STL_L2_SWH(ISD_indexes_int),idx_outliers_STL_L2_SWH]=outliers_by_percentiles(STL_L2_SWH(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);

                        %----------------- FIT ------------------------------------
                        [STL_L2_FIT(ISD_indexes_int),idx_outliers_STL_L2_FIT]=outliers_by_percentiles(STL_L2_FIT(ISD_indexes_int),outlier_percentil_low,outlier_percentil_high,IQR_times);

                    case 'hampel'
                        %------------------ SSH ---------------------------
                        %[ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);
                        [~,idx_outliers_STL_L2_SSH] = hampel(STL_L2_lat_surf(ISD_indexes_int),STL_L2_SSH(ISD_indexes_int),hampel_wind,hampel_sigma);
                        STL_L2_SSH(idx_int_ISD(idx_outliers_STL_L2_SSH))=NaN;
                        %----------------- sigma0 ---------------------------------
                        %[ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                               
                        [~,idx_outliers_STL_L2_sigma0] = hampel(STL_L2_lat_surf(ISD_indexes_int),STL_L2_sigma0(ISD_indexes_int),hampel_wind,hampel_sigma);                                               
                        STL_L2_sigma0(idx_int_ISD(idx_outliers_STL_L2_sigma0))=NaN;
                        %----------------- SWH ------------------------------------
                        %[ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                                              
                        [~,idx_outliers_STL_L2_SWH] = hampel(STL_L2_lat_surf(ISD_indexes_int),STL_L2_SWH(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                                               
                        STL_L2_SWH(idx_int_ISD(idx_outliers_STL_L2_SWH))=NaN; 
                        %----------------- FIT ------------------------------------
                        %[ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),indx] = hampel(ISD_lat_surf(ISD_indexes_int),ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                                              
                        [~,idx_outliers_STL_L2_FIT] = hampel(STL_L2_lat_surf(ISD_indexes_int),STL_L2_FIT(ISD_indexes_int),hampel_wind,hampel_sigma);                                                                                               
                        STL_L2_FIT(idx_int_ISD(idx_outliers_STL_L2_FIT))=NaN;                         
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
    finish              = min(i_wfm+sliding_window_SSH-1,ESA_num_surfaces_filtered);
    % mean_win            = nanmean(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));
    % std_win             = nanstd(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));    
    % if ESA_L2_SSH_r1(idx_int_ESA(i_wfm))>(mean_win+threshold_std*std_win) || ESA_L2_SSH_r1(idx_int_ESA(i_wfm))<(mean_win-threshold_std*std_win)
        % idx_outliers_ESA_L2_SSH(i_wfm) = 1;
    % end    
    ESA_L2_SSH_r1_smoothed(i_wfm) = nanmean(ESA_L2_SSH_r1(idx_int_ESA(start:finish)));
    ESA_L2_SSH_r1_nocorr_smoothed(i_wfm) = nanmean(ESA_L2_SSH_r1_nocorr(idx_int_ESA(start:finish)));
    
    %-------------------- sigma0 from L2 ESA product ----------------------
    start               = max(i_wfm-sliding_window_sigma0,1);
    finish              = min(i_wfm+sliding_window_sigma0-1,ESA_num_surfaces_filtered);
    % mean_win            = nanmean(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
    % std_win             = nanstd(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
    % if ESA_L2_sigma0_r1(idx_int_ESA(i_wfm))>(mean_win+threshold_std*std_win) || ESA_L2_sigma0_r1(idx_int_ESA(i_wfm))<(mean_win-threshold_std*std_win)
        % idx_outliers_ESA_L2_sigma0(i_wfm) = 1;
    % end
    ESA_L2_sigma0_r1_smoothed(i_wfm) = nanmean(ESA_L2_sigma0_r1(idx_int_ESA(start:finish)));
end

%--------------------------------------------------------------------------
%--------------------------- ISD ------------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    ISD_L2_SSH_analytical_SWH_MSSfixed_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_analytical_SWH_MSSfixed_smoothed       = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_SWH_analytical_SWH_MSSfixed_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_COR_analytical_SWH_MSSfixed_smoothed          = zeros(1,ISD_num_surfaces_filtered);   
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    ISD_L2_SSH_analytical_MSS_SWHfixed_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_analytical_MSS_SWHfixed_smoothed       = zeros(1,ISD_num_surfaces_filtered);    
    ISD_L2_COR_analytical_MSS_SWHfixed_smoothed          = zeros(1,ISD_num_surfaces_filtered);    
end

if ~isempty(idx_int_thres)
    ISD_L2_SSH_threshold_smoothed           = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_threshold_smoothed        = zeros(1,ISD_num_surfaces_filtered); 
end

if ~isempty(idx_int_ocog)
    ISD_L2_SSH_OCOG_smoothed            = zeros(1,ISD_num_surfaces_filtered);
    ISD_L2_sigma0_OCOG_smoothed         = zeros(1,ISD_num_surfaces_filtered);
end
if ~isempty(filename_L2_STL)
    %should have the same # surfaces as ISD generated from L1B ISD
    STL_L2_SSH_smoothed          = zeros(1,ISD_num_surfaces_filtered);
    STL_L2_sigma0_smoothed       = zeros(1,ISD_num_surfaces_filtered);
    STL_L2_SWH_smoothed          = zeros(1,ISD_num_surfaces_filtered); 
    STL_L2_FIT_smoothed          = zeros(1,ISD_num_surfaces_filtered); 
end
for i_wfm=1:ISD_num_surfaces_filtered
    start_SSH               = max(i_wfm-sliding_window_SSH,1);
    finish_SSH              = min(i_wfm+sliding_window_SSH-1,ISD_num_surfaces_filtered);
    start_sigma0            = max(i_wfm-sliding_window_sigma0,1);
    finish_sigma0           = min(i_wfm+sliding_window_sigma0-1,ISD_num_surfaces_filtered);
    start_SWH               = max(i_wfm-sliding_window_SWH,1);
    finish_SWH              = min(i_wfm+sliding_window_SWH-1,ISD_num_surfaces_filtered);
    start_COR               = max(i_wfm-sliding_window_COR,1);
    finish_COR              = min(i_wfm+sliding_window_COR-1,ISD_num_surfaces_filtered);
    for i_ret=1:length(retrackers)
        switch char(retrackers(i_ret))
            case {'ANALYTICAL_SWH'}
                % compute the errors w.r.t fitting on the data using a smooth
                % function
                %------------------- SSH ----------------------------------
                % mean_win            = nanmean(ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(start_SSH:finish_SSH)));
                % std_win             = nanstd(ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(start_SSH:finish_SSH)));
                % if ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed(i_wfm) = 1;
                % end
                ISD_L2_SSH_analytical_SWH_MSSfixed_smoothed(i_wfm)   = nanmean(ISD_L2_SSH_analytical_SWH_MSSfixed(idx_int_ISD(start_SSH:finish_SSH)));
                
                %----------------- sigma0 ---------------------------------
                % mean_win            = nanmean(ISD_L2_sigma0_analytical_SWH_MSSfixed(idx_int_ISD(start_sigma0:finish_sigma0)));
                % std_win             = nanstd(ISD_L2_sigma0_analytical_SWH_MSSfixed(idx_int_ISD(start_sigma0:finish_sigma0)));
                % if ISD_L2_sigma0_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_sigma0_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed(i_wfm) = 1;
                % end
                ISD_L2_sigma0_analytical_SWH_MSSfixed_smoothed(i_wfm) = nanmean(ISD_L2_sigma0_analytical_SWH_MSSfixed(idx_int_ISD(start_sigma0:finish_sigma0)));
                %----------------- SWH ------------------------------------
                % mean_win            = nanmean(ISD_L2_SWH_analytical_SWH_MSSfixed(idx_int_ISD(start_SWH:finish_SWH)));
                % std_win             = nanstd(ISD_L2_SWH_analytical_SWH_MSSfixed(idx_int_ISD(start_SWH:finish_SWH)));
                % if ISD_L2_SWH_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SWH_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed(i_wfm) = 1;
                % end
                ISD_L2_SWH_analytical_SWH_MSSfixed_smoothed(i_wfm) = nanmean(ISD_L2_SWH_analytical_SWH_MSSfixed(idx_int_ISD(start_SWH:finish_SWH)));
                % ---------- COR: pearson correlation coefficient ---------
                % mean_win            = nanmean(ISD_L2_COR_analytical_SWH_MSSfixed(idx_int_ISD(start_COR:finish_COR)));
                % std_win             = nanstd(ISD_L2_COR_analytical_SWH_MSSfixed(idx_int_ISD(start_COR:finish_COR)));
                % if ISD_L2_COR_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_COR_analytical_SWH_MSSfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_COR_analytical_SWH_MSSfixed(i_wfm) = 1;
                % end
                ISD_L2_COR_analytical_SWH_MSSfixed_smoothed(i_wfm) = nanmean(ISD_L2_COR_analytical_SWH_MSSfixed(idx_int_ISD(start_COR:finish_COR)));
            case {'ANALYTICAL_MSS'}
                % compute the errors w.r.t fitting on the data using a smooth
                % function
                %------------------- SSH ----------------------------------
                % mean_win            = nanmean(ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(start_SSH:finish_SSH)));
                % std_win             = nanstd(ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(start_SSH:finish_SSH)));
                % if ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed(i_wfm) = 1;
                % end
                ISD_L2_SSH_analytical_MSS_SWHfixed_smoothed(i_wfm)   = nanmean(ISD_L2_SSH_analytical_MSS_SWHfixed(idx_int_ISD(start_SSH:finish_SSH)));
                
                %----------------- sigma0 ---------------------------------
                % mean_win            = nanmean(ISD_L2_sigma0_analytical_MSS_SWHfixed(idx_int_ISD(start_sigma0:finish_sigma0)));
                % std_win             = nanstd(ISD_L2_sigma0_analytical_MSS_SWHfixed(idx_int_ISD(start_sigma0:finish_sigma0)));
                % if ISD_L2_sigma0_analytical_MSS_SWHfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_sigma0_analytical_MSS_SWHfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_sigma0_analytical_MSS_SWHfixed(i_wfm) = 1;
                % end
                ISD_L2_sigma0_analytical_MSS_SWHfixed_smoothed(i_wfm) = nanmean(ISD_L2_sigma0_analytical_MSS_SWHfixed(idx_int_ISD(start_sigma0:finish_sigma0)));
                % ---------- COR: pearson correlation coefficient ---------
                % mean_win            = nanmean(ISD_L2_COR_analytical_MSS_SWHfixed(idx_int_ISD(start_COR:finish_COR)));
                % std_win             = nanstd(ISD_L2_COR_analytical_MSS_SWHfixed(idx_int_ISD(start_COR:finish_COR)));
                % if ISD_L2_COR_analytical_MSS_SWHfixed(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_COR_analytical_MSS_SWHfixed(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_COR_analytical_MSS_SWHfixed(i_wfm) = 1;
                % end
                ISD_L2_COR_analytical_MSS_SWHfixed_smoothed(i_wfm) = nanmean(ISD_L2_COR_analytical_MSS_SWHfixed(idx_int_ISD(start_COR:finish_COR)));                
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
           case {'OCOG'}
                %------------------------ SSH  ----------------------------
                % mean_win            = nanmean(ISD_L2_SSH_OCOG(idx_int_ISD(start_SSH:finish_SSH)));
                % std_win             = nanstd(ISD_L2_SSH_OCOG(idx_int_ISD(start_SSH:finish_SSH)));
                % if ISD_L2_SSH_OCOG(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_SSH_OCOG(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_SSH_OCOG(i_wfm) = 1;
                % end
                ISD_L2_SSH_OCOG_smoothed(i_wfm)   = nanmean(ISD_L2_SSH_OCOG(idx_int_ISD(start_SSH:finish_SSH)));
                %-------------------- sigma0 ------------------------------
                % mean_win            = nanmean(ISD_L2_sigma0_OCOG(idx_int_ISD(start_sigma0:finish_sigma0)));
                % std_win             = nanstd(ISD_L2_sigma0_OCOG(idx_int_ISD(start_sigma0:finish_sigma0)));
                % if ISD_L2_sigma0_OCOG(idx_int_ISD(i_wfm))>(mean_win+threshold_std*std_win) || ISD_L2_sigma0_OCOG(idx_int_ISD(i_wfm))<(mean_win-threshold_std*std_win)
                    % idx_outliers_ISD_L2_sigma0_OCOG(i_wfm) = 1;
                % end
                ISD_L2_sigma0_OCOG_smoothed(i_wfm) = nanmean(ISD_L2_sigma0_OCOG(idx_int_ISD(start_sigma0:finish_sigma0)));                                
        end
    end
    if ~isempty(filename_L2_STL)
        STL_L2_SSH_smoothed(i_wfm)   = nanmean(STL_L2_SSH(idx_int_ISD(start_SSH:finish_SSH)));
        STL_L2_SWH_smoothed(i_wfm)   = nanmean(STL_L2_SWH(idx_int_ISD(start_SSH:finish_SSH)));
        STL_L2_sigma0_smoothed(i_wfm)= nanmean(STL_L2_sigma0(idx_int_ISD(start_SSH:finish_SSH)));
        STL_L2_FIT_smoothed(i_wfm)= nanmean(STL_L2_FIT(idx_int_ISD(start_SSH:finish_SSH)));
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
    %----------------- analytical_SWH_MSSfixed -------------------------------------
    if ~isempty(idx_int_analytical_SWH_MSSfixed)
        fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_SSH_analytical_SWH_MSSfixed_filtered=ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed;
        res.SSH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_analytical_SWH_MSSfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISD (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2);
        res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_analytical_SWH_MSSfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISD (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_L2);        
    end
    if ~isempty(idx_int_analytical_MSS_SWHfixed)
        fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_SSH_analytical_MSS_SWHfixed_filtered=ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed;
        res.SSH.ANALYTICAL_MSS_SWHfixed.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_analytical_MSS_SWHfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISD (analytical_MSS_SWHfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_MSS_SWHfixed.RMSE_error_L2);
        res.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_analytical_MSS_SWHfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISD (analytical_MSS_SWHfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_L2);
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
        fprintf(fid,'$---------------- OCOG    RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_SSH_OCOG_filtered=ISD_L2_SSH_OCOG(ISD_indexes_int);
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_ISD_L2_SSH_OCOG;
        res.SSH.OCOG.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_OCOG_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-ISD (OCOG-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG.RMSE_error_L2);
        res.SSH.OCOG.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-ISD_L2_SSH_OCOG_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-ISD (OCOG-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG.mean_error_L2);
    end    
end
if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    fprintf(fid,'$---------------- RETRACKERS COMPARISON STL -------------------------------------------$\n');
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    %----------------- analytical_SWH_MSSfixed -------------------------------------
    if ~isempty(idx_int_analytical_SWH_MSSfixed)
        fprintf(fid,'$---------------- ANALYTICAL RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ISD  ---------------------------
        STL_L2_SSH_filtered=STL_L2_SSH(ISD_indexes_int);
        idx_outliers=idx_outliers_STL_L2_SSH | idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed;
        res.SSH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2_STL=sqrt(nanmean((STL_L2_SSH_filtered(~idx_outliers)-ISD_L2_SSH_analytical_SWH_MSSfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH STL-ISD (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2_STL);
        res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_L2_STL=(nanmean((STL_L2_SSH_filtered(~idx_outliers)-ISD_L2_SSH_analytical_SWH_MSSfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH STL-ISD (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_L2_STL);
        
        %------------ Comparison with ESA ---------------------------------
        %----------- Comparison with ESA L2 ---------------------------
        fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
        idx_outliers=idx_outliers_ESA_L2_SSH | idx_outliers_STL_L2_SSH;
        res.SSH.ANALYTICAL_STL.RMSE_error_L2=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-STL_L2_SSH_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SSH ESA-STL (ANALYTICAL_STL) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_STL.RMSE_error_L2);
        res.SSH.ANALYTICAL_STL.mean_error_L2=(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers)-STL_L2_SSH_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SSH ESA-STL (ANALYTICAL_STL) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_STL.mean_error_L2);        
    end
end

fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
    res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting=nanmean(ISD_L2_SSH_analytical_SWH_MSSfixed_filtered(~idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed)-ISD_L2_SSH_analytical_SWH_MSSfixed_smoothed(~idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed));
    fprintf(fid,'Mean error SSH fitting (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting);
    res.SSH.ANALYTICAL_SWH_MSSfixed.rmse_fitting=sqrt(nanmean((ISD_L2_SSH_analytical_SWH_MSSfixed_filtered(~idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed)-ISD_L2_SSH_analytical_SWH_MSSfixed_smoothed(~idx_outliers_ISD_L2_SSH_analytical_SWH_MSSfixed)).^2));
    fprintf(fid,'RMSE SSH fitting (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_SWH_MSSfixed.rmse_fitting);
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    fprintf(fid,'$---------------- ANALYTICAL_MSS_SWHfixed RETRACKER -------------------------------------------$\n');
    res.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_fitting=nanmean(ISD_L2_SSH_analytical_MSS_SWHfixed_filtered(~idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed)-ISD_L2_SSH_analytical_MSS_SWHfixed_smoothed(~idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed));
    fprintf(fid,'Mean error SSH fitting (analytical_MSS_SWHfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_fitting);
    res.SSH.ANALYTICAL_MSS_SWHfixed.rmse_fitting=sqrt(nanmean((ISD_L2_SSH_analytical_MSS_SWHfixed_filtered(~idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed)-ISD_L2_SSH_analytical_MSS_SWHfixed_smoothed(~idx_outliers_ISD_L2_SSH_analytical_MSS_SWHfixed)).^2));
    fprintf(fid,'RMSE SSH fitting (analytical_MSS_SWHfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_MSS_SWHfixed.rmse_fitting);
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
    res.SSH.OCOG.mean_error_fitting=nanmean(ISD_L2_SSH_OCOG_filtered(~idx_outliers_ISD_L2_SSH_OCOG)-ISD_L2_SSH_OCOG_smoothed(~idx_outliers_ISD_L2_SSH_OCOG));
    fprintf(fid,'Mean error SSH fitting (OCOG-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG.mean_error_fitting);
    res.SSH.OCOG.rmse_fitting=sqrt(nanmean((ISD_L2_SSH_OCOG_filtered(~idx_outliers_ISD_L2_SSH_OCOG)-ISD_L2_SSH_OCOG_smoothed(~idx_outliers_ISD_L2_SSH_OCOG)).^2));
    fprintf(fid,'RMSE SSH fitting (OCOG-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SSH.OCOG.rmse_fitting);
end

if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- STARLAB ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.SSH.ANALYTICAL_STL.mean_error_fitting=nanmean(STL_L2_SSH_filtered(~idx_outliers_STL_L2_SSH)-STL_L2_SSH_smoothed(~idx_outliers_STL_L2_SSH));
    fprintf(fid,'Mean error SSH fitting (analytical STL) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_STL.mean_error_fitting);
    res.SSH.ANALYTICAL_STL.rmse_fitting=sqrt(nanmean((STL_L2_SSH_filtered(~idx_outliers_STL_L2_SSH)-STL_L2_SSH_smoothed(~idx_outliers_STL_L2_SSH)).^2));
    fprintf(fid,'RMSE SSH fitting (analytical STL) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ANALYTICAL_STL.rmse_fitting);    
end

fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
res.SSH.ESA_L2.mean_error_fitting=nanmean(ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_SSH_r1_smoothed(~idx_outliers_ESA_L2_SSH));
fprintf(fid,'Mean error SSH fitting (L2-ESA) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ESA_L2.mean_error_fitting);
res.SSH.ESA_L2.rmse_fitting=sqrt(nanmean((ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH)-ESA_L2_SSH_r1_smoothed(~idx_outliers_ESA_L2_SSH)).^2));
fprintf(fid,'RMSE SSH fitting (L2-ESA) [m]: '); fprintf(fid,'%.18g\n',res.SSH.ESA_L2.rmse_fitting);

%--------------------------------------------------------------------------
%-------------------- HISTOGRAM COMPUTATION -------------------------------
%--------------------------------------------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_analytical_SWH_MSSfixed]=hist(ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_analytical_SWH_MSSfixed=counts./(sum(counts).*(centers_ISD_L2_SSH_analytical_SWH_MSSfixed(2)-centers_ISD_L2_SSH_analytical_SWH_MSSfixed(1)));
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    %$---------------- ANALYTICAL_MSS_SWHfixed RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_analytical_MSS_SWHfixed]=hist(ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_analytical_MSS_SWHfixed=counts./(sum(counts).*(centers_ISD_L2_SSH_analytical_MSS_SWHfixed(2)-centers_ISD_L2_SSH_analytical_MSS_SWHfixed(1)));
end
if ~isempty(idx_int_thres)
    %$---------------- THRESHOLD  RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_threshold]=hist(ISD_L2_SSH_threshold(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_threshold=counts./(sum(counts).*(centers_ISD_L2_SSH_threshold(2)-centers_ISD_L2_SSH_threshold(1)));
end
if ~isempty(idx_int_ocog)
    %$---------------- OCOG   RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SSH_OCOG]=hist(ISD_L2_SSH_OCOG(ISD_indexes_int),nbins);
    pdf_ISD_L2_SSH_OCOG=counts./(sum(counts).*(centers_ISD_L2_SSH_OCOG(2)-centers_ISD_L2_SSH_OCOG(1)));
end
[counts,centers_ESA_L2_SSH_r1_filtered]=hist(ESA_L2_SSH_r1_filtered(~idx_outliers_ESA_L2_SSH),nbins);
pdf_ESA_L2_SSH_r1_filtered=counts./(sum(counts).*(centers_ESA_L2_SSH_r1_filtered(2)-centers_ESA_L2_SSH_r1_filtered(1)));

%---------------------- ploting -------------------------------------------
%--------------------- Comparison smoothed vs non-smoothed ----------------
% %------------------------ ESA ---------------------------------------------
% figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_filtered,linestyle_ESA);
% hold on;
% plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_smoothed,linestyle_ESA_smooth);
% %plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr(ESA_indexes_int),'Marker','^','Color', [1 0.5 0.2]);
% %plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1_nocorr_smoothed,linestyle_ESA_smooth_nocorr);
% text_in_textbox=strcat('Fitting ESA (corr.)--> RMSE [m]:',{' '},num2str(res.SSH.ESA_L2.rmse_fitting,'%.4g'),', Mean [m]:',{' '},num2str(res.SSH.ESA_L2.mean_error_fitting,'%.4g'));
% %legend('L2-ESA','L2-ESA smoothed','L2-ESA (undoing corrections)','L2-ESA (undoing corrections) smoothed');
% legend('L2-ESA','L2-ESA (smoothed)');
% xlabel('Latitude [deg.]'); ylabel(strcat(upper(sh_name_nc),' [m]'));
% title(strcat('Comparison ',upper(sh_name_nc),': ESA (non-smoothed vs smoothed)'));
% xlim=get(gca,'xlim');
% ylim=get(gca,'ylim');
% [figx,figy] = dsxy2figxy(gca, xlim, ylim);
% dim = [figx(1),figy(1),0.1,0.1];
% if annotation_box_active
% annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on','Units','points');
% end
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ESA_SSH_smoothed.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_ESA_SSH.png'))

% %------------------------ ISD ---------------------------------------------
% figure; 
% legend_text={''};
% text_in_textbox={''};
% if ~isempty(idx_int_thres)    
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold_filtered,linestyle_threshold);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold_smoothed,linestyle_threshold_smooth);
%     legend_text=[legend_text,{'L2-ISD threshold','L2-ISD threshold (smoothed)'}];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting Threshold--> RMSE [m]:',{' '},num2str(res.SSH.THRESHOLD.rmse_fitting,'%.4g'),', Mean [m]:',{' '},num2str(res.SSH.THRESHOLD.mean_error_fitting,'%.4g'))];
% end
% if ~isempty(idx_int_ocog)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG_filtered,linestyle_OCOG);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG_smoothed,linestyle_OCOG_smooth);
%     legend_text=[legend_text,{'L2-ISD OCOG','L2-ISD OCOG (smoothed)'}];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting OCOG--> RMSE [m]:',{' '},num2str(res.SSH.OCOG.rmse_fitting,'%.4g'),', Mean [m]:',{' '},num2str(res.SSH.OCOG.mean_error_fitting,'%.4g'))];
% end
% if ~isempty(idx_int_analytical_SWH_MSSfixed)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_SWH_MSSfixed_filtered,linestyle_analytical_SWH_MSSfixed);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_SWH_MSSfixed_smoothed,linestyle_analytical_SWH_MSSfixed_smooth);
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_SWH_MSSfixed),strcat('L2-ISD',{' '},title_name_SWH_MSSfixed,' (smoothed)')];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting',{' '},title_name_SWH_MSSfixed,' --> RMSE [m]:',{' '},num2str(res.SSH.ANALYTICAL_SWH_MSSfixed.rmse_fitting,'%.4g'),', Mean [m]: ',num2str(res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting,'%.4g'))];
% end
% if ~isempty(idx_int_analytical_MSS_SWHfixed)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_MSS_SWHfixed_filtered,linestyle_analytical_MSS_SWHfixed);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_MSS_SWHfixed_smoothed,linestyle_analytical_MSS_SWHfixed_smooth);
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_MSS_SWHfixed),strcat('L2-ISD',{' '},title_name_MSS_SWHfixed,' (smoothed)')];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting',{' '},title_name_MSS_SWHfixed,' --> RMSE [m]:',{' '},num2str(res.SSH.ANALYTICAL_MSS_SWHfixed.rmse_fitting,'%.4g'),', Mean [m]: ',num2str(res.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_fitting,'%.4g'))];
% end
% legend(legend_text(~cellfun(@isempty,legend_text)));
% xlim=get(gca,'xlim');
% ylim=get(gca,'ylim');
% [figx,figy] = dsxy2figxy(gca, xlim, ylim);
% dim = [figx(1),figy(1),0.1,0.1];
% text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
% if annotation_box_active
% annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
% end
% xlabel('Latitude [deg.]'); ylabel(strcat(upper(sh_name_nc),' [m]'));
% title(strcat('Comparison ',upper(sh_name_nc),': isardSAT (non-smoothed vs smoothed)'));
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_SSH_smoothed.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_ISD_SSH.png'))

%---------------- Plot comparison ESA & ISD -------------------------------
figure;  
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_SSH_r1(ESA_indexes_int),linestyle_ESA);
hold on;
legend_text={'L2-ESA'};
text_in_textbox={''};
text_in_textbox=[text_in_textbox,...
        strcat('ESA --> RMSE [m]:',{' '},num2str(res.SSH.ESA_L2.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.ESA_L2.mean_error_fitting,'%.4g'))];        
if ~isempty(idx_int_thres)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_threshold(ISD_indexes_int),linestyle_threshold);
    legend_text=[legend_text,{'L2-ISR threshold'}];
    text_in_textbox=[text_in_textbox,...
        strcat('Threshold --> RMSE [m]:',{' '},num2str(res.SSH.THRESHOLD.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.THRESHOLD.mean_error_fitting,'%.4g'))];        
end
if ~isempty(idx_int_ocog)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_OCOG(ISD_indexes_int),linestyle_OCOG);
    legend_text=[legend_text,{'L2-ISR OCOG'}];
    text_in_textbox=[text_in_textbox,...
        strcat('OCOG --> RMSE [m]:',{' '},num2str(res.SSH.OCOG.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.OCOG.mean_error_fitting,'%.4g'))];        
end
if ~isempty(filename_L2_STL)
    plot(STL_L2_lat_surf(ISD_indexes_int),STL_L2_SSH(ISD_indexes_int),linestyle_STL)
    legend_text=[legend_text,{'L2-STL SAMOSA-3 (S-3)'}];
    text_in_textbox=[text_in_textbox,...
        strcat({'STL SAMOSA-3 (S-3) '},'--> RMSE [m]:',{' '},num2str(res.SSH.ANALYTICAL_STL.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.ANALYTICAL_STL.mean_error_fitting,'%.4g'))];    
    title(strcat('Comparison',{' '},upper(sh_name_nc),': ESA & isardSAT & Starlab'));
else
    title(strcat('Comparison',{' '},upper(sh_name_nc),': ESA & isardSAT'));
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_SWH_MSSfixed(ISD_indexes_int),linestyle_analytical_SWH_MSSfixed);
    legend_text=[legend_text,strcat('L2-ISR',{' '},title_name_SWH_MSSfixed)];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_SWH_MSSfixed,' --> RMSE [m]:',{' '},num2str(res.SSH.ANALYTICAL_SWH_MSSfixed.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting,'%.4g'))];    
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_SSH_analytical_MSS_SWHfixed(ISD_indexes_int),linestyle_analytical_MSS_SWHfixed);
    legend_text=[legend_text,strcat('L2-ISR ',{' '},title_name_MSS_SWHfixed)];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_MSS_SWHfixed,' --> RMSE [m]:',{' '},num2str(res.SSH.ANALYTICAL_MSS_SWHfixed.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SSH.ANALYTICAL_MSS_SWHfixed.mean_error_fitting,'%.4g'))];    
end
%title(strcat('Comparison',{' '},upper(sh_name_nc),': ESA & isardSAT'));
legend(legend_text(~cellfun(@isempty,legend_text)));
xlabel('Latitude [deg.]'); ylabel(strcat(upper(sh_name_nc),' [m]'));
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+y_add_textbox,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
if annotation_box_active
    annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
end
%print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_SSH_ESA_ISD.png'))
print('-dpng',strcat(output_path,name_file_L2_ISD,'_cmp_SSH.png'))


% %-------------------- HISTOGRAM PLOTING -----------------------------------
% figure;  
% plot(centers_ESA_L2_SSH_r1_filtered,pdf_ESA_L2_SSH_r1_filtered,strcat(linestyle_ESA,'-'));
% hold on;
% legend_text={'L2-ESA'};
% if ~isempty(idx_int_thres)
%     plot(centers_ISD_L2_SSH_threshold,pdf_ISD_L2_SSH_threshold,strcat(linestyle_threshold,'-'));
%     legend_text=[legend_text,{'L2-ISD threshold'}];
% end
% if ~isempty(idx_int_ocog)
%     plot(centers_ISD_L2_SSH_OCOG,pdf_ISD_L2_SSH_OCOG,strcat(linestyle_OCOG,'-'));
%     legend_text=[legend_text,{'L2-ISD OCOG'}];
% end
% if ~isempty(idx_int_analytical_SWH_MSSfixed)
%     plot(centers_ISD_L2_SSH_analytical_SWH_MSSfixed,pdf_ISD_L2_SSH_analytical_SWH_MSSfixed,strcat(linestyle_analytical_SWH_MSSfixed,'-'));
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_SWH_MSSfixed)];
% end
% if ~isempty(idx_int_analytical_MSS_SWHfixed)
%     plot(centers_ISD_L2_SSH_analytical_MSS_SWHfixed,pdf_ISD_L2_SSH_analytical_MSS_SWHfixed,strcat(linestyle_analytical_MSS_SWHfixed,'-'));
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_MSS_SWHfixed)];
% end
% title(strcat('Comparison SSH histograms: ESA & isardSAT'),'Interpreter','Tex');
% legend(legend_text(~cellfun(@isempty,legend_text)));
% xlabel('SSH [m]'); ylabel('PDF');
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_PDF_SSH_ESA_ISD.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_cmp_SSH_PDF.png'))

    
%% ------------------------------ SIGMA0 ----------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$---------------- SIGMA0  --------------------------------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
ESA_L2_sigma0_r1_filtered=ESA_L2_sigma0_r1(ESA_indexes_int);
%idx_outliers_ESA_L2_SSH=zeros(1,ESA_num_surfaces_filtered);
if ESA_num_surfaces_filtered==ISD_num_surfaces_filtered
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    fprintf(fid,'$---------------- RETRACKERS COMPARISON ESA -------------------------------------------$\n');
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    %----------------- analytical_SWH_MSSfixed -------------------------------------
    if ~isempty(idx_int_analytical_SWH_MSSfixed)
        fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_analytical_SWH_MSSfixed_filtered=ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed;
        res.SIGMA0.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_SWH_MSSfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (analytical_SWH_MSSfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2);
        res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_SWH_MSSfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (analytical_SWH_MSSfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_L2);        
    end
    %----------------- analytical_MSS_SWHfixed -------------------------------------
    if ~isempty(idx_int_analytical_MSS_SWHfixed)
        fprintf(fid,'$---------------- ANALYTICAL_MSS_SWHfixed RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_analytical_MSS_SWHfixed_filtered=ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_ISD_L2_sigma0_analytical_MSS_SWHfixed;
        res.SIGMA0.ANALYTICAL_MSS_SWHfixed.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_MSS_SWHfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (analytical_MSS_SWHfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_MSS_SWHfixed.RMSE_error_L2);
        res.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_MSS_SWHfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (analytical_MSS_SWHfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_L2);        
    end
    if ~isempty(idx_int_thres)
        fprintf(fid,'$---------------- THRESHOLD  RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_threshold_filtered=ISD_L2_sigma0_threshold(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_ISD_L2_sigma0_threshold;
        res.SIGMA0.THRESHOLD.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_threshold_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (threshold-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.THRESHOLD.RMSE_error_L2);
        res.SIGMA0.THRESHOLD.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_threshold_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (threshold-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.THRESHOLD.mean_error_L2);
    end
    if ~isempty(idx_int_ocog)
        fprintf(fid,'$---------------- OCOG    RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        ISD_L2_sigma0_OCOG_filtered=ISD_L2_sigma0_OCOG(ISD_indexes_int);
		idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_ISD_L2_sigma0_OCOG;
        res.SIGMA0.OCOG.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_OCOG_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-ISD (OCOG-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG.RMSE_error_L2);
        res.SIGMA0.OCOG.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-ISD_L2_sigma0_OCOG_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-ISD (OCOG-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG.mean_error_L2);
    end    
end
if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    fprintf(fid,'$---------------- RETRACKERS COMPARISON STL -------------------------------------------$\n');
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    %----------------- analytical_SWH_MSSfixed -------------------------------------
    if ~isempty(idx_int_analytical_SWH_MSSfixed)
        fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------
        STL_L2_sigma0_filtered=STL_L2_sigma0(ISD_indexes_int);
		idx_outliers=idx_outliers_STL_L2_sigma0 | idx_outliers_ISD_L2_sigma0_analytical_SWH_MSSfixed;
        res.SIGMA0.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2_STL=sqrt(nanmean((STL_L2_sigma0_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_SWH_MSSfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 STL-ISD (analytical_SWH_MSSfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2);
        res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_L2_STL=(nanmean((STL_L2_sigma0_filtered(~idx_outliers)-ISD_L2_sigma0_analytical_SWH_MSSfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 STL-ISD (analytical_SWH_MSSfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_L2);        
        
        %----------- Comparison with ESA L2 ---------------------------        
		idx_outliers=idx_outliers_ESA_L2_sigma0 | idx_outliers_STL_L2_sigma0;
        res.SIGMA0.ANALYTICAL_STL.RMSE_error_L2=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-STL_L2_sigma0_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error sigma0 ESA-STL (ANALYTICAL_STL-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_STL.RMSE_error_L2);
        res.SIGMA0.ANALYTICAL_STL.mean_error_L2=(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers)-STL_L2_sigma0_filtered(~idx_outliers))));
        fprintf(fid,'Mean error sigma0 ESA-STL (ANALYTICAL_STL-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_STL.mean_error_L2);        
    end
end
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
    res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_fitting=nanmean(ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int)-ISD_L2_sigma0_analytical_SWH_MSSfixed_smoothed);
    fprintf(fid,'Mean error sigma0 fitting (analytical_SWH_MSSfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_fitting);
    res.SIGMA0.ANALYTICAL_SWH_MSSfixed.rmse_fitting=sqrt(nanmean((ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int)-ISD_L2_sigma0_analytical_SWH_MSSfixed_smoothed).^2));
    fprintf(fid,'RMSE sigma0 fitting (analytical_SWH_MSSfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_SWH_MSSfixed.rmse_fitting);
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    fprintf(fid,'$---------------- ANALYTICAL_MSS_SWHfixed RETRACKER -------------------------------------------$\n');
    res.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_fitting=nanmean(ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int)-ISD_L2_sigma0_analytical_MSS_SWHfixed_smoothed);
    fprintf(fid,'Mean error sigma0 fitting (analytical_MSS_SWHfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_fitting);
    res.SIGMA0.ANALYTICAL_MSS_SWHfixed.rmse_fitting=sqrt(nanmean((ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int)-ISD_L2_sigma0_analytical_MSS_SWHfixed_smoothed).^2));
    fprintf(fid,'RMSE sigma0 fitting (analytical_MSS_SWHfixed-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_MSS_SWHfixed.rmse_fitting);
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
    res.SIGMA0.OCOG.mean_error_fitting=nanmean(ISD_L2_sigma0_OCOG(ISD_indexes_int)-ISD_L2_sigma0_OCOG_smoothed);
    fprintf(fid,'Mean error sigma0 fitting (OCOG-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG.mean_error_fitting);
    res.SIGMA0.OCOG.rmse_fitting=sqrt(nanmean((ISD_L2_sigma0_OCOG(ISD_indexes_int)-ISD_L2_sigma0_OCOG_smoothed).^2));
    fprintf(fid,'RMSE sigma0 fitting (OCOG-retracker) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.OCOG.rmse_fitting);
end

if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- STARLAB ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.SIGMA0.ANALYTICAL_STL.mean_error_fitting=nanmean(STL_L2_sigma0_filtered(~idx_outliers_STL_L2_sigma0)-STL_L2_sigma0_smoothed(~idx_outliers_STL_L2_sigma0));
    fprintf(fid,'Mean error sigma0 fitting (analytical STL) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_STL.mean_error_fitting);
    res.SIGMA0.ANALYTICAL_STL.rmse_fitting=sqrt(nanmean((STL_L2_sigma0_filtered(~idx_outliers_STL_L2_sigma0)-STL_L2_sigma0_smoothed(~idx_outliers_STL_L2_sigma0)).^2));
    fprintf(fid,'RMSE sigma0 fitting (analytical STL) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ANALYTICAL_STL.rmse_fitting);    
end

fprintf(fid,'$---------------- ESA RETRACKER -------------------------------------------$\n');
res.SIGMA0.ESA_L2.mean_error_fitting=nanmean(ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_sigma0)-ESA_L2_sigma0_r1_smoothed(~idx_outliers_ESA_L2_sigma0));
fprintf(fid,'Mean error sigma0 fitting (L2-ESA) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ESA_L2.mean_error_fitting);
res.SIGMA0.ESA_L2.rmse_fitting=sqrt(nanmean((ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_sigma0)-ESA_L2_sigma0_r1_smoothed(~idx_outliers_ESA_L2_sigma0)).^2));
fprintf(fid,'RMSE sigma0 fitting (L2-ESA) [dB]: '); fprintf(fid,'%.18g\n',res.SIGMA0.ESA_L2.rmse_fitting);

%--------------------------------------------------------------------------
%-------------------- HISTOGRAM COMPUTATION -------------------------------
%--------------------------------------------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_analytical_SWH_MSSfixed]=hist(ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_analytical_SWH_MSSfixed=counts./(sum(counts).*(centers_ISD_L2_sigma0_analytical_SWH_MSSfixed(2)-centers_ISD_L2_sigma0_analytical_SWH_MSSfixed(1)));
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    %$---------------- ANALYTICAL_MSS_SWHfixed RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_analytical_MSS_SWHfixed]=hist(ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_analytical_MSS_SWHfixed=counts./(sum(counts).*(centers_ISD_L2_sigma0_analytical_MSS_SWHfixed(2)-centers_ISD_L2_sigma0_analytical_MSS_SWHfixed(1)));
end
if ~isempty(idx_int_thres)
    %$---------------- THRESHOLD  RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_threshold]=hist(ISD_L2_sigma0_threshold(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_threshold=counts./(sum(counts).*(centers_ISD_L2_sigma0_threshold(2)-centers_ISD_L2_sigma0_threshold(1)));
end
if ~isempty(idx_int_ocog)
    %$---------------- OCOG_ICE   RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_sigma0_OCOG]=hist(ISD_L2_sigma0_OCOG(ISD_indexes_int),nbins);
    pdf_ISD_L2_sigma0_OCOG=counts./(sum(counts).*(centers_ISD_L2_sigma0_OCOG(2)-centers_ISD_L2_sigma0_OCOG(1)));
end
[counts,centers_ESA_L2_sigma0_r1_filtered]=hist(ESA_L2_sigma0_r1_filtered(~idx_outliers_ESA_L2_sigma0),nbins);
pdf_ESA_L2_sigma0_r1_filtered=counts./(sum(counts).*(centers_ESA_L2_sigma0_r1_filtered(2)-centers_ESA_L2_sigma0_r1_filtered(1)));


%---------------------- ploting -------------------------------------------
%--------------------- Comparison smoothed vs non-smoothed ----------------
% %------------------------ ESA ---------------------------------------------
% figure; plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),linestyle_ESA);
% hold on;
% plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1_smoothed,linestyle_ESA_smooth);
% text_in_textbox=strcat('Fitting ESA (corr.)--> RMSE [dB]:',{' '},num2str(res.SIGMA0.ESA_L2.rmse_fitting,'%.4g'),', Mean [dB]:',{' '},num2str(res.SIGMA0.ESA_L2.mean_error_fitting,'%.4g'));
% legend('L2-ESA','L2-ESA (smoothed)');
% xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
% title(strcat('Comparison \sigma^0: ESA (non-smoothed vs smoothed)'),'Interpreter','Tex');
% xlim=get(gca,'xlim');
% ylim=get(gca,'ylim');
% [figx,figy] = dsxy2figxy(gca, xlim, ylim);
% dim = [figx(1),figy(1),0.1,0.1];
% text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
% if annotation_box_active
% annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
% end
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ESA_SIGMA0_smoothed.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_ESA_SIG0.png'))

% %------------------------ ISD ---------------------------------------------
% figure; 
% legend_text={''};
% text_in_textbox={''};
% if ~isempty(idx_int_thres)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold(ISD_indexes_int),linestyle_threshold);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold_smoothed,linestyle_threshold_smooth);
%     legend_text=[legend_text,{'L2-ISD threshold','L2-ISD threshold (smoothed)'}];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting Threshold--> RMSE [dB]:',{' '},num2str(res.SIGMA0.THRESHOLD.rmse_fitting,'%.4g'),', Mean [dB]:',{' '},num2str(res.SIGMA0.THRESHOLD.mean_error_fitting,'%.4g'))];
%     hold on;
% end
% if ~isempty(idx_int_ocog)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG(ISD_indexes_int),linestyle_OCOG);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG_smoothed,linestyle_OCOG_smooth);
%     legend_text=[legend_text,{'L2-ISD OCOG','L2-ISD OCOG (smoothed)'}];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting OCOG--> RMSE [dB]:',{' '},num2str(res.SIGMA0.OCOG.rmse_fitting,'%.4g'),', Mean [dB]:',{' '},num2str(res.SIGMA0.OCOG.mean_error_fitting,'%.4g'))];
% 
% end
% if ~isempty(idx_int_analytical_SWH_MSSfixed)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),linestyle_analytical_SWH_MSSfixed);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_SWH_MSSfixed_smoothed,linestyle_analytical_SWH_MSSfixed_smooth);
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_SWH_MSSfixed),strcat('L2-ISD',{' '},title_name_SWH_MSSfixed,' (smoothed)')];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting',{' '},title_name_SWH_MSSfixed,' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_SWH_MSSfixed.rmse_fitting,'%.4g'),', Mean [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_fitting,'%.4g'))];
% end
% if ~isempty(idx_int_analytical_MSS_SWHfixed)
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),linestyle_analytical_MSS_SWHfixed);
%     hold on;
%     plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_MSS_SWHfixed_smoothed,linestyle_analytical_MSS_SWHfixed_smooth);
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_MSS_SWHfixed),strcat('L2-ISD',{' '},title_name_MSS_SWHfixed,' (smoothed)')];
%     text_in_textbox=[text_in_textbox,...
%         strcat('Fitting',{' '},title_name_MSS_SWHfixed,' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_MSS_SWHfixed.rmse_fitting,'%.4g'),', Mean [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_fitting,'%.4g'))];
% end
% legend(legend_text(~cellfun(@isempty,legend_text)));
% xlim=get(gca,'xlim');
% ylim=get(gca,'ylim');
% [figx,figy] = dsxy2figxy(gca, xlim, ylim);
% dim = [figx(1),figy(1),0.1,0.1];
% text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
% if annotation_box_active
% annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
% end
% xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
% title(strcat('Comparison \sigma^0: isardSAT (non-smoothed vs smoothed)'),'Interpreter','Tex');
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_SIGMA0_smoothed.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_ISD_SIG0.png'))

%---------------- Plot comparison ESA & ISD -------------------------------
figure;  
plot(ESA_lat_surf(ESA_indexes_int),ESA_L2_sigma0_r1(ESA_indexes_int),linestyle_ESA);
hold on;
legend_text={'L2-ESA'};
text_in_textbox={''};
text_in_textbox=[text_in_textbox,...
    strcat('ESA --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ESA_L2.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.ESA_L2.mean_error_fitting,'%.4g'))];
if ~isempty(idx_int_thres)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_threshold(ISD_indexes_int),linestyle_threshold);
    legend_text=[legend_text,{'L2-ISR threshold'}];
    text_in_textbox=[text_in_textbox,...
        strcat('Threshold --> RMSE [dB]:',{' '},num2str(res.SIGMA0.THRESHOLD.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.THRESHOLD.mean_error_fitting,'%.4g'))];    
end
if ~isempty(idx_int_ocog)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_OCOG(ISD_indexes_int),linestyle_OCOG);
    legend_text=[legend_text,{'L2-ISR OCOG'}];
    text_in_textbox=[text_in_textbox,...
        strcat('OCOG --> RMSE [dB]:',{' '},num2str(res.SIGMA0.OCOG.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.OCOG.mean_error_fitting,'%.4g'))];        
end
if ~isempty(filename_L2_STL)
    plot(STL_L2_lat_surf(ISD_indexes_int),STL_L2_sigma0(ISD_indexes_int),linestyle_STL)
    legend_text=[legend_text,{'L2-STL SAMOSA-3 (S-3)'}];
    text_in_textbox=[text_in_textbox,...
        strcat({'STL SAMOSA-3 (S-3) '},' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_STL.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_STL.mean_error_fitting,'%.4g'))];
    title(strcat('Comparison \sigma^0: ESA & isardSAT & Starlab'),'Interpreter','Tex');
else
    title(strcat('Comparison \sigma^0: ESA & isardSAT'),'Interpreter','Tex');
end

if ~isempty(idx_int_analytical_SWH_MSSfixed)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_SWH_MSSfixed(ISD_indexes_int),linestyle_analytical_SWH_MSSfixed);
    legend_text=[legend_text,strcat('L2-ISR',{' '},title_name_SWH_MSSfixed)];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_SWH_MSSfixed,' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_SWH_MSSfixed.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_SWH_MSSfixed.mean_error_fitting,'%.4g'))];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    plot(ISD_lat_surf(ISD_indexes_int),ISD_L2_sigma0_analytical_MSS_SWHfixed(ISD_indexes_int),linestyle_analytical_MSS_SWHfixed);
    legend_text=[legend_text,strcat('L2-ISR',{' '},title_name_MSS_SWHfixed)];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_MSS_SWHfixed,' --> RMSE [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_MSS_SWHfixed.rmse_fitting,'%.4g'))];%,', Mean [dB]:',{' '},num2str(res.SIGMA0.ANALYTICAL_MSS_SWHfixed.mean_error_fitting,'%.4g'))];
end
%title(strcat('Comparison \sigma^0: ESA & isardSAT'),'Interpreter','Tex');
legend(legend_text(~cellfun(@isempty,legend_text))); %,'Location','southeast'
xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]','Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+y_add_textbox,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
if annotation_box_active
    annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
end
print('-dpng',strcat(output_path,name_file_L2_ISD,'_cmp_SIG0.png'))

% %-------------------- HISTOGRAM PLOTING -----------------------------------
% figure;  
% plot(centers_ESA_L2_sigma0_r1_filtered,pdf_ESA_L2_sigma0_r1_filtered,strcat(linestyle_ESA,'-'));
% hold on;
% legend_text={'L2-ESA'};
% if ~isempty(idx_int_thres)
%     plot(centers_ISD_L2_sigma0_threshold,pdf_ISD_L2_sigma0_threshold,strcat(linestyle_threshold,'-'));
%     legend_text=[legend_text,{'L2-ISD threshold'}];
% end
% if ~isempty(idx_int_ocog)
%     plot(centers_ISD_L2_sigma0_OCOG,pdf_ISD_L2_sigma0_OCOG,strcat(linestyle_OCOG,'-'));
%     legend_text=[legend_text,{'L2-ISD OCOG'}];
% end
% if ~isempty(idx_int_analytical_SWH_MSSfixed)
%     plot(centers_ISD_L2_sigma0_analytical_SWH_MSSfixed,pdf_ISD_L2_sigma0_analytical_SWH_MSSfixed,strcat(linestyle_analytical_SWH_MSSfixed,'-'));
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_SWH_MSSfixed)];
% end
% if ~isempty(idx_int_analytical_MSS_SWHfixed)
%     plot(centers_ISD_L2_sigma0_analytical_MSS_SWHfixed,pdf_ISD_L2_sigma0_analytical_MSS_SWHfixed,strcat(linestyle_analytical_MSS_SWHfixed,'-'));
%     legend_text=[legend_text,strcat('L2-ISD ',{' '},title_name_MSS_SWHfixed)];
% end
% title(strcat('Comparison \sigma^0 histograms: ESA & isardSAT'),'Interpreter','Tex');
% legend(legend_text(~cellfun(@isempty,legend_text)));
% xlabel('\sigma^0 [dB]','Interpreter','Tex'); ylabel('PDF');
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_PDF_sigma0_ESA_ISD.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_cmp_SIG0_PDF.png'))

%% ------------------------------ SWH -------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
fprintf(fid,'$----------------   SWH   --------------------------------------------$\n');
fprintf(fid,'$---------------- ------  --------------------------------------------$\n');
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    ISD_L2_SWH_analytical_SWH_MSSfixed_filtered=ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int);
end
if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    fprintf(fid,'$---------------- RETRACKERS COMPARISON STL -------------------------------------------$\n');
    fprintf(fid,'$---------------- ---------- ---------- --- -------------------------------------------$\n');
    STL_L2_SWH_filtered=STL_L2_SWH(ISD_indexes_int);
    %----------------- analytical_SWH_MSSfixed -------------------------------------
    if ~isempty(idx_int_analytical_SWH_MSSfixed)
        fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
        %----------- Comparison with ESA L2 ---------------------------        
		idx_outliers=idx_outliers_STL_L2_SSH | idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed;
        res.SWH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2_STL=sqrt(nanmean((STL_L2_SWH_filtered(~idx_outliers)-ISD_L2_SWH_analytical_SWH_MSSfixed_filtered(~idx_outliers)).^2));
        fprintf(fid,'RMSE error SWH STL-ISD (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL_SWH_MSSfixed.RMSE_error_L2_STL);
        res.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_L2_STL=(nanmean((STL_L2_SWH_filtered(~idx_outliers)-ISD_L2_SWH_analytical_SWH_MSSfixed_filtered(~idx_outliers))));
        fprintf(fid,'Mean error SWH STL-ISD (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_L2_STL);        
    end
end
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
fprintf(fid,'$---------------- RETRACKERS FITTING ERRORS -------------------------------------------$\n');
fprintf(fid,'$---------------- ---------- ------- ------ -------------------------------------------$\n');
%----------- Fitting errors ------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
    res.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting=nanmean(ISD_L2_SWH_analytical_SWH_MSSfixed_filtered(~idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed)-ISD_L2_SWH_analytical_SWH_MSSfixed_smoothed(~idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed));
    fprintf(fid,'Mean error SWH fitting (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting);
    res.SWH.ANALYTICAL_SWH_MSSfixed.rmse_fitting=sqrt(nanmean((ISD_L2_SWH_analytical_SWH_MSSfixed_filtered(~idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed)-ISD_L2_SWH_analytical_SWH_MSSfixed_smoothed(~idx_outliers_ISD_L2_SWH_analytical_SWH_MSSfixed)).^2));
    fprintf(fid,'RMSE SWH fitting (analytical_SWH_MSSfixed-retracker) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL_SWH_MSSfixed.rmse_fitting);
end
if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- STARLAB ANALYTICAL RETRACKER -------------------------------------------$\n');
    res.SWH.ANALYTICAL_STL.mean_error_fitting=nanmean(STL_L2_SWH_filtered(~idx_outliers_STL_L2_SWH)-STL_L2_SWH_smoothed(~idx_outliers_STL_L2_SWH));
    fprintf(fid,'Mean error SWH fitting (analytical STL) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL_STL.mean_error_fitting);
    res.SWH.ANALYTICAL_STL.rmse_fitting=sqrt(nanmean((STL_L2_SWH_filtered(~idx_outliers_STL_L2_SWH)-STL_L2_SWH_smoothed(~idx_outliers_STL_L2_SWH)).^2));
    fprintf(fid,'RMSE SWH fitting (analytical STL) [m]: '); fprintf(fid,'%.18g\n',res.SWH.ANALYTICAL_STL.rmse_fitting);
    
end
%--------------------------------------------------------------------------
%-------------------- HISTOGRAM COMPUTATION -------------------------------
%--------------------------------------------------------------------------
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    %$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$
    [counts,centers_ISD_L2_SWH_analytical_SWH_MSSfixed]=hist(ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int),nbins);
    pdf_ISD_L2_SWH_analytical_SWH_MSSfixed=counts./(sum(counts).*(centers_ISD_L2_SWH_analytical_SWH_MSSfixed(2)-centers_ISD_L2_SWH_analytical_SWH_MSSfixed(1)));
end
%----------------------------- ploting ------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};
if ~isempty(filename_L2_STL)
    plot(STL_L2_lat_surf(ISD_indexes_int),STL_L2_SWH(ISD_indexes_int),linestyle_STL)
    hold on;
    plot(STL_L2_lat_surf(ISD_indexes_int),abs(STL_L2_SWH_smoothed),linestyle_STL_smooth);
    legend_text=[legend_text,{'L2-STL SAMOSA-3 (S-3)','L2-STL SAMOSA-3 (S-3) (smoothed)'}];
    text_in_textbox=[text_in_textbox,...
        strcat({'STL SAMOSA-3 (S-3) '},' --> RMSE [m]:',{' '},num2str(res.SWH.ANALYTICAL_STL.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SWH.ANALYTICAL_STL.mean_error_fitting,'%.4g'))];
    title(strcat('Comparison SWH: isardSAT & Starlab'),'Interpreter','Tex');
else
    title(strcat('Comparison SWH: isardSAT'),'Interpreter','Tex');
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_SWH_analytical_SWH_MSSfixed(ISD_indexes_int)),linestyle_analytical_SWH_MSSfixed);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_SWH_analytical_SWH_MSSfixed_smoothed),linestyle_analytical_SWH_MSSfixed_smooth);
    legend_text=[legend_text,strcat('L2-ISR',{' '},title_name_SWH_MSSfixed),strcat('L2-ISR',{' '},title_name_SWH_MSSfixed,' (smoothed)')];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_SWH_MSSfixed,' --> RMSE [m]:',{' '},num2str(res.SWH.ANALYTICAL_SWH_MSSfixed.rmse_fitting,'%.4g'))];%,', Mean [m]:',{' '},num2str(res.SWH.ANALYTICAL_SWH_MSSfixed.mean_error_fitting,'%.4g'))];
end
%axis([2 18 0 4]);
%axis([14 26 0 5]);
%axis([-42 -26 1 6]);
legend(legend_text(~cellfun(@isempty,legend_text)));%,'Location','southeast'
xlabel('Latitude [deg.]'); ylabel('SWH [m]','Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+0.015,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
if annotation_box_active
    annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
end
%print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_SWH_smoothed.png'))
print('-dpng',strcat(output_path,name_file_L2_ISD,'_ISD_SWH.png'))
% %-------------------- HISTOGRAM PLOTING -----------------------------------
% figure;  
% legend_text={''};
% if ~isempty(idx_int_analytical_SWH_MSSfixed)
%     plot(centers_ISD_L2_SWH_analytical_SWH_MSSfixed,pdf_ISD_L2_SWH_analytical_SWH_MSSfixed,strcat(linestyle_analytical_SWH_MSSfixed,'-'));
%     hold on;
%     legend_text=[legend_text,strcat('L2-ISD',{' '},title_name_SWH_MSSfixed)];
% end
% title(strcat('Histogram SWH: isardSAT'),'Interpreter','Tex');
% legend(legend_text(~cellfun(@isempty,legend_text)));
% xlabel('SWH [m]'); ylabel('PDF');
% %print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_PDF_SWH_ISD.png'))
% print('-dpng',strcat(output_path,name_file_L2_ISD,'_cmp_SWH_PDF.png'))

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
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    fprintf(fid,'$---------------- ANALYTICAL_SWH_MSSfixed RETRACKER -------------------------------------------$\n');
    res.COR.ANALYTICAL_SWH_MSSfixed.mean=nanmean(ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int));
    fprintf(fid,'Mean value Pearson Correlation Coefficient (analytical_SWH_MSSfixed-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_SWH_MSSfixed.mean);
    res.COR.ANALYTICAL_SWH_MSSfixed.mean_error_fitting=nanmean(ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int)-ISD_L2_COR_analytical_SWH_MSSfixed_smoothed);
    fprintf(fid,'Mean error Pearson Correlation Coefficient fitting (analytical_SWH_MSSfixed-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_SWH_MSSfixed.mean_error_fitting);
    res.COR.ANALYTICAL_SWH_MSSfixed.rmse_fitting=sqrt(nanmean((ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int)-ISD_L2_COR_analytical_SWH_MSSfixed_smoothed).^2));
    fprintf(fid,'RMSE Pearson Correlation Coefficient fitting (analytical_SWH_MSSfixed-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_SWH_MSSfixed.rmse_fitting);
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    fprintf(fid,'$---------------- ANALYTICAL_MSS_SWHfixed RETRACKER -------------------------------------------$\n');
    res.COR.ANALYTICAL_MSS_SWHfixed.mean=nanmean(ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int));
    fprintf(fid,'Mean value Pearson Correlation Coefficient (analytical_MSS_SWHfixed-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_MSS_SWHfixed.mean);
    res.COR.ANALYTICAL_MSS_SWHfixed.mean_error_fitting=nanmean(ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int)-ISD_L2_COR_analytical_MSS_SWHfixed_smoothed);
    fprintf(fid,'Mean error Pearson Correlation Coefficient fitting (analytical_MSS_SWHfixed-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_MSS_SWHfixed.mean_error_fitting);
    res.COR.ANALYTICAL_MSS_SWHfixed.rmse_fitting=sqrt(nanmean((ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int)-ISD_L2_COR_analytical_MSS_SWHfixed_smoothed).^2));
    fprintf(fid,'RMSE Pearson Correlation Coefficient fitting (analytical_MSS_SWHfixed-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_MSS_SWHfixed.rmse_fitting);
end
if ~isempty(filename_L2_STL)
    fprintf(fid,'$---------------- ANALYTICAL_STL RETRACKER -------------------------------------------$\n');
    res.COR.ANALYTICAL_STL.mean=nanmean(STL_L2_FIT(ISD_indexes_int));
    fprintf(fid,'Mean value Pearson Correlation Coefficient (ANALYTICAL_STL-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_STL.mean);
    res.COR.ANALYTICAL_STL.mean_error_fitting=nanmean(STL_L2_FIT(ISD_indexes_int)-STL_L2_FIT_smoothed);
    fprintf(fid,'Mean error Pearson Correlation Coefficient fitting (ANALYTICAL_STL-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_STL.mean_error_fitting);
    res.COR.ANALYTICAL_STL.rmse_fitting=sqrt(nanmean((STL_L2_FIT(ISD_indexes_int)-STL_L2_FIT_smoothed).^2));
    fprintf(fid,'RMSE Pearson Correlation Coefficient fitting (ANALYTICAL_STL-retracker) [percentage]: '); fprintf(fid,'%.18g\n',res.COR.ANALYTICAL_STL.rmse_fitting);
end
%----------------------------- ploting ------------------------------------
figure; 
legend_text={''};
text_in_textbox={''};
if ~isempty(filename_L2_STL)
    plot(STL_L2_lat_surf(ISD_indexes_int),abs(STL_L2_FIT(ISD_indexes_int)),linestyle_STL);
    hold on;
    plot(STL_L2_lat_surf(ISD_indexes_int),abs(STL_L2_FIT_smoothed),linestyle_STL_smooth);
    legend_text=[legend_text,strcat('L2-STL',{' (100-misfit) '},' SAMOSA-3 (S-3)'),strcat('L2-STL',{' (100-misfit) '},' SAMOSA-3 (S-3)',' (smoothed)')];
    text_in_textbox=[text_in_textbox,...
        strcat('STL SAMOSA-3 (S-3) --> RMSE [%]:',{' '},num2str(res.COR.ANALYTICAL_STL.rmse_fitting,'%.4g'))];%,', Mean [%]:',{' '},num2str(res.COR.ANALYTICAL_STL.mean_error_fitting,'%.4g'))];
    title(strcat('Comparison goodness of fitting (gof): isardSAT & Starlab'),'Interpreter','Tex');
else
    title(strcat('Comparison goodness of fitting (gof): isardSAT'),'Interpreter','Tex');
end
if ~isempty(idx_int_analytical_SWH_MSSfixed)
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_COR_analytical_SWH_MSSfixed(ISD_indexes_int)),linestyle_analytical_SWH_MSSfixed);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_COR_analytical_SWH_MSSfixed_smoothed),linestyle_analytical_SWH_MSSfixed_smooth);
    legend_text=[legend_text,strcat('L2-ISR',{' \rho_{pearson} '},title_name_SWH_MSSfixed),strcat('L2-ISR',{' \rho_{pearson} '},title_name_SWH_MSSfixed,' (smoothed)')];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_SWH_MSSfixed,' --> RMSE [%]:',{' '},num2str(res.COR.ANALYTICAL_SWH_MSSfixed.rmse_fitting,'%.4g'))];%,', Mean [%]:',{' '},num2str(res.COR.ANALYTICAL_SWH_MSSfixed.mean_error_fitting,'%.4g'))];
end
if ~isempty(idx_int_analytical_MSS_SWHfixed)
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_COR_analytical_MSS_SWHfixed(ISD_indexes_int)),linestyle_analytical_MSS_SWHfixed);
    hold on;
    plot(ISD_lat_surf(ISD_indexes_int),abs(ISD_L2_COR_analytical_MSS_SWHfixed_smoothed),linestyle_analytical_MSS_SWHfixed_smooth);
    legend_text=[legend_text,strcat('L2-ISR',{' \rho_{pearson} '},title_name_MSS_SWHfixed),strcat('L2-ISR',{' \rho_{pearson} '},title_name_MSS_SWHfixed,' (smoothed)')];
    text_in_textbox=[text_in_textbox,...
        strcat({'ISR '},title_name_MSS_SWHfixed,' --> RMSE [%]:',{' '},num2str(res.COR.ANALYTICAL_MSS_SWHfixed.rmse_fitting,'%.4g'))];%,', Mean [%]:',{' '},num2str(res.COR.ANALYTICAL_MSS_SWHfixed.mean_error_fitting,'%.4g'))];
end
%axis([2 18 95 100]);
%axis([14 26 95 100]);
%axis([-42 -26 95 100]);
legend(legend_text(~cellfun(@isempty,legend_text))); %,'Location','southeast'
xlabel('Latitude [deg.]'); ylabel('gof [%]','Interpreter','Tex');
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[figx,figy] = dsxy2figxy(gca, xlim, ylim);
dim = [figx(1),figy(1)+0.015,0.1,0.1];
text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
if annotation_box_active
    annotation('textbox',dim,'String',text_in_textbox,'FitBoxToText','on');
end
%print('-dpng',strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'Comparison_ISD_COR_RHO_smoothed.png'))
print('-dpng',strcat(output_path,name_file_L2_ISD,'_ISD_COR.png'))

%% --------------------------- Saving Information -------------------------
fclose(fid);
%close all;
%save(strcat(output_path,strrep(char(name_file_L2_ISD),'.nc','_'),'L2_Evaluation.mat'),'res');
save(strcat(output_path,name_file_L2_ISD,'_L2_Evaluation.mat'),'res');

end