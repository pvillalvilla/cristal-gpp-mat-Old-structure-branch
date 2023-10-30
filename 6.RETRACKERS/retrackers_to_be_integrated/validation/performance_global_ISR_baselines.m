function performance_global_ISR_baselines(input_path_L2_ISR_bs,name_bs,path_comparison_results,varargin)
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
version_matlab=version;
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<3 || nargin>(3+17*2))
    error('Wrong number of input parameters');
end
p = inputParser;
p.addParamValue('figure_format','',@(x)ischar(x));
p.addParamValue('res_fig','',@(x)ischar(x));
p.addParamValue('win_size_detrending',20);
p.addParamValue('step_SWH',0.2);
p.addParamValue('step_SSH',0.5);
p.addParamValue('step_hr',0.5);
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('flag_outliers_removal',0);
p.addParamValue('type_outliers_removal','',@(x)ischar(x));
p.addParamValue('smooth_param',0);
p.addParamValue('define_min_max_SWH',0);
p.addParamValue('min_SWH',1.5);
p.addParamValue('max_SWH',4.5);
p.addParamValue('input_path_L1_ISR_bs',{''},@(x)iscellstr(x));
p.addParamValue('input_path_L2_ESA','',@(x)ischar(x));
p.addParamValue('input_path_L2_STL','',@(x)ischar(x));
p.addParamValue('input_path_L1_ESA','',@(x)ischar(x));
p.addParamValue('annotation_box_active',1);
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.parse(varargin{:});
figure_format=p.Results.figure_format;
res_fig=p.Results.res_fig;
win_size_detrending=p.Results.win_size_detrending;
step_SWH=p.Results.step_SWH;
step_SSH=p.Results.step_SSH;
step_hr=p.Results.step_hr;
sh_name_nc=p.Results.sh_name_nc;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=p.Results.type_outliers_removal;
smooth_param=p.Results.smooth_param;
define_min_max_SWH=p.Results.define_min_max_SWH;
min_SWH=p.Results.min_SWH;
max_SWH=p.Results.max_SWH;
input_path_L1_ISR_bs=p.Results.input_path_L1_ISR_bs;
input_path_L2_ESA=p.Results.input_path_L2_ESA;
input_path_L1_ESA=p.Results.input_path_L1_ESA;
input_path_L2_STL=p.Results.input_path_L2_STL;
annotation_box_active=p.Results.annotation_box_active;
filename_mask_KML=p.Results.filename_mask_KML;

y_add_textbox=0.01;
nanclr=[1 1 1];
min_std_SSH=1; %cm
max_std_SSH=10;%cm

min_std_SWH=10; %cm
max_std_SWH=50;%cm

min_std_sigma0=0; %cm
max_std_sigma0=0.7;%cm




%% ----------------- Hard coded variables definition ----------------------
%----------------- Define linstyles for bulk comparison -------------------
%linestyle_ESA='+k';
marker_ESA='+';
color_ESA=rgb('Black');
%linestyle_STL='h';
marker_STL='h';
color_STL=rgb('Lime');
%linestyle_bs={'^r','*b','og','xc','sy','.m'};
marker_bs={'^','*','o','x','s','.'};
color_bs=[rgb('red'); rgb('blue'); rgb('Green'); rgb('Cyan'); rgb('gold'); rgb('Magenta')];
fontsize_xlabel_tracks=8;
size_marker=12;
thick_marker=3.0;
generate_plots=1;
plot_fits_downsampling=50;
textbox_fontsize=10;
legend_fontsize=15;


switch lower(figure_format)
    case 'eps'
        file_ext='.eps';
        print_file='-depsc';
    case 'png'
        file_ext='.png';
        print_file='-dpng';
end
clear p;

%--------------- Define the outliers filtering parameters -----------------
%removal of outliers 
%-------using percentile & IQR (Interquartile Range)
% outliers data<(percentil_low-IQR_times*IQR) | data>(percentil_high+IQR_times*IQR)
IQR_times=1.5; %number of IQR 
outlier_percentil_low=25.0;
outlier_percentil_high=75.0;
%--------using hampel filter
hampel_wind=3;% size of window half size
hampel_sigma=3; %number of std deviations to which a sample differ from local median




%% --------- CREATE COMPARISON RESULTS PATH ---------------------------------
mkdir(path_comparison_results);

%% ------------ LOAD THE GEOGRAPHICAL MASK --------------------------------
if isempty(filename_mask_KML)
    within_geo_mask=1;
    geo_mask  = [];
else
    geo_mask  = kml2lla(filename_mask_KML);
end


%% ------------ REFERENCE BASELINE-1 FOR FILES SEARCH ---------------------
filesBulk.inputPath       =   char(input_path_L2_ISR_bs(1));
filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,'.nc')));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

i_files_valid=0;
N_baselines=length(input_path_L2_ISR_bs);

%N_samples_max=1;

%--------------------------------------------------------------------------
for i_file=1:filesBulk.nFilesNC
    empty_flag=ones(1,N_baselines-1);
    
    %taking first baseline as reference
    filename_L2_ISR=char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);
    data_string=filename_L2_ISR(17:17+30);
    filename_L2_ISR=strcat(char(input_path_L2_ISR_bs(1)),filename_L2_ISR);
    
    % ---- Checking whether the track within geomask is available ---------
    if ~isempty(filename_mask_KML)        
        if ~isempty(strfind(lower(char(name_bs(1))),'acdc'))
            %------------------- ACDC product -------------------------------------
            %---------------- Geometry variables --------------------------------------
            ISR_lat_surf=double(ncread(filename_L2_ISR,'lat_l1b_echo_sar_ku')).';
            ISR_lon_surf=double(ncread(filename_L2_ISR,'lon_l1b_echo_sar_ku')).';
        else
            ISR_lat_surf=double(ncread(filename_L2_ISR,'lat_20_ku')).';
            ISR_lon_surf=double(ncread(filename_L2_ISR,'lon_20_ku')).';            
        end
        
%         if i_file==1
%             figure; geoshow('landareas.shp','DisplayType','texturemap');
%             hold on;
%             geoshow(geo_mask.coord(:,2),geo_mask.coord(:,1)); hold on;            
%         end
%         geoshow(ISR_lat_surf,ISR_lon_surf,'DisplayType','Point','Marker','+','MarkerEdgeColor','red'); hold on;
                                        
        idx_int=inpolygon(ISR_lon_surf,ISR_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));
        if ~any(idx_int)
            disp(strcat('Track,',{' '},filename_L2_ISR,{' '},'outside the limits of the geographical mask'))
            continue;               
        end
    end
    
    for b=2:N_baselines
        %----------- checking whether the L2 ESA is available -------------
        if ~isempty(input_path_L2_ESA)
            if b==2
                input_L2_ESA_Files   = dir(fullfile(char(input_path_L2_ESA),['*' data_string(1:15) '*.DBL']));                
                if isempty(input_L2_ESA_Files)                    
                    %add one second to initial time acquisition
                    init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                    input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
                    if isempty(input_L2_ESA_Files)
                        %add one second to initial time acquisition
                        end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                        input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
                        if isempty(input_L2_ESA_Files)
                            input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
                            if isempty(input_L2_ESA_Files)
                                disp(strcat('File not available in L2 ESA data set: ',data_string));
                                continue;
                            end
                        end
                    end
                end
                %check whether files form Starlab are available
                if ~isempty(input_path_L2_STL)
                    input_L2_STL_Files   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
                    if isempty(input_L2_ESA_Files)
                        disp(char(strcat('File not available in ',{' '},name_bs(b),{' baseline data set: '},data_string)));
                        continue;
                    end
                end
            end
        end
        %------------ CHECKING ISR L2 AVAILABLE PRODUCT for other baselines------------        
        input_L2_ISR_Files   = dir(fullfile(char(input_path_L2_ISR_bs(b)),['*' data_string(1:15) '*.nc']));
        if ~isempty(input_L2_ISR_Files)
            empty_flag(b-1)=0;
        end
    end
    if ~any(empty_flag)
        i_files_valid=i_files_valid+1;
        filesBulk.indexFilesNC_valid(i_files_valid)=filesBulk.indexFilesNC(i_file);
    end
end

filesBulk.nFilesNC_valid=length(filesBulk.indexFilesNC_valid);
disp(strcat('Total number of input b1 (',name_bs(1),') files: ',num2str(filesBulk.nFilesNC)));
disp(strcat('Total number of valid files for evaluation (',strjoin(name_bs,' vs '),'): ',num2str(filesBulk.nFilesNC_valid)));


%% ------------------ RUN the performance ---------------------------------
%----------- Definition of variables --------------------------------------
SSH_std_whole=[];
SWH_std_whole=[];
sigma0_std_whole=[];
amp_fit_std_whole=[];
max_wvfm_std_whole=[];
if ~isempty(input_path_L2_ESA)
    ESA_SSH_std_whole=[];
    ESA_sigma0_std_whole=[];
end
if ~isempty(input_path_L2_STL)
   STL_SSH_std_whole=[];
   STL_sigma0_std_whole=[];
   STL_SWH_std_whole=[];
end
%mean values of stds over track
SSH_std_mean_whole=[];
SWH_std_mean_whole=[];
sigma0_std_mean_whole=[];
COR_std_mean_whole=[];
if ~isempty(input_path_L2_ESA)
    ESA_SSH_std_mean_whole=[];
    ESA_sigma0_std_mean_whole=[];
end
if ~isempty(input_path_L2_STL)
    STL_SSH_std_mean_whole=[];
    STL_sigma0_std_mean_whole=[];
    STL_SWH_std_mean_whole=[];
end


SSH_mean_whole=[];
SWH_mean_whole=[];
sigma0_mean_whole=[];
amp_fit_mean_whole=[];
COR_mean_whole=[];
max_wvfm_mean_whole=[];
if ~cellfun(@isempty,input_path_L1_ISR_bs)
    hr_mean_whole=[];
    nb_mean_whole=[];
end
if ~isempty(input_path_L2_ESA)
    ESA_SSH_mean_whole=[];
    ESA_sigma0_mean_whole=[];
    if ~isempty(input_path_L1_ESA)
        ESA_hr_mean_whole=[]; %height rate (radial veloctiy)
        ESA_nb_mean_whole=[];
    end
end
if ~isempty(input_path_L2_STL)
    STL_SSH_mean_whole=[];
    STL_sigma0_mean_whole=[];
    STL_SWH_mean_whole=[];
end

LAT_median_whole=[];
LON_median_whole=[];

SSH_whole=[];
SWH_whole=[];
sigma0_whole=[];
COR_whole=[];
amp_fit_whole=[];
max_wvfm_whole=[];
if ~cellfun(@isempty,input_path_L1_ISR_bs)
    hr_whole=[];
    nb_whole=[];
end
if ~isempty(input_path_L2_ESA)
    ESA_SSH_whole=[];
    ESA_sigma0_whole=[];
    if ~isempty(input_path_L1_ESA)
        ESA_hr_whole=[];
        ESA_nb_whole=[];
    end
    SSH_mean_error_ESA_ISR_whole=[];
    sigma0_mean_error_ESA_ISR_whole=[];
end
if ~isempty(input_path_L2_STL)
    STL_SSH_whole=[];
    STL_sigma0_whole=[];
    STL_SWH_whole=[];
    SSH_mean_error_STL_ISR_whole=[];
    SWH_mean_error_STL_ISR_whole=[];
    sigma0_mean_error_STL_ISR_whole=[];
end


if smooth_param
    SSH_smooth_whole=[]; %use the same window as the one for std and detrending
    SWH_smooth_whole=[];
    sigma0_smooth_whole=[];
    COR_smooth_whole=[];
    SSH_RMSE_whole=[];
    SWH_RMSE_whole=[];
    sigma0_RMSE_whole=[];
    COR_RMSE_whole=[];   
    if ~isempty(input_path_L2_ESA)
        ESA_SSH_smooth_whole=[];
        ESA_sigma0_smooth_whole=[];
        ESA_SSH_RMSE_whole=[];
        ESA_sigma0_RMSE_whole=[];
    end  
    if ~isempty(input_path_L2_STL)
        STL_SSH_smooth_whole=[];
        STL_sigma0_smooth_whole=[];
        STL_SWH_smooth_whole=[];
        STL_SWH_RMSE_whole=[];
        STL_SSH_RMSE_whole=[];
        STL_sigma0_RMSE_whole=[];
    end
end


if define_min_max_SWH
    max_SWH=-Inf;
    min_SWH=Inf;
end


num_surfaces_whole=NaN(1,filesBulk.nFilesNC_valid);

for i_files_valid=1:filesBulk.nFilesNC_valid
    
    try
        for b=1:N_baselines
            
            %% ----------- Taking the files names -----------------------------------
            if b==1
                filename_L2_ISR=char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name);
                data_string=filename_L2_ISR(17:17+30);
                disp(num2str(i_files_valid));
                disp(data_string);
                
                %% --------------- READ ESA product -----------------------------------
                %try
                if ~isempty(input_path_L2_ESA)
                    input_L2_ESA_Files   = dir(fullfile(char(input_path_L2_ESA),['*' data_string(1:15) '*.DBL']));
                    if isempty(input_L2_ESA_Files)
                        %add one second to initial time acquisition
                        init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                        input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
                        if isempty(input_L2_ESA_Files)
                            %add one second to initial time acquisition
                            end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                            input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
                            if isempty(input_L2_ESA_Files)
                                input_L2_ESA_Files   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
                            end
                        end
                    end
                    filename_L2_ESA=strcat(char(input_path_L2_ESA),input_L2_ESA_Files.name);
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
                    clear CS2;
                end
                %---------- read the radial velocity information L1------------
                if ~isempty(input_path_L1_ESA)
                    input_L1_ESA_Files   = dir(fullfile(char(input_path_L1_ESA),['*' data_string(1:15) '*.DBL']));
                    if isempty(input_L1_ESA_Files)
                        %add one second to initial time acquisition
                        init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                        input_L1_ESA_Files   = dir(fullfile(input_path_L1_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
                        if isempty(input_L1_ESA_Files)
                            %add one second to initial time acquisition
                            end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                            input_L1_ESA_Files   = dir(fullfile(input_path_L1_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
                            if isempty(input_L1_ESA_Files)
                                input_L1_ESA_Files   = dir(fullfile(input_path_L1_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
                            end
                        end
                    end
                    filename_L1_ESA=strcat(char(input_path_L1_ESA),input_L1_ESA_Files.name);
                    [~,CS2]=Cryo_L1b_read(filename_L1_ESA);
                    ESA_hr=reshape(CS2.GEO.H_rate,[1,ESA_num_surfaces]);
                    ESA_nb=reshape(CS2.SAR.N_averaged_echoes,[1,ESA_num_surfaces]);
                    
                    
                    clear CS;
                end
                
                %------------- Read Starlab file ------------------------------
                if ~isempty(input_path_L2_STL)
                    input_L2_STL_Files   = dir(fullfile(input_path_L2_STL,['*' data_string(1:15) '*.nc']));
                    filename_L2_STL=char(input_L2_STL_Files.name);
                    filename_L2_STL=strcat(char(input_path_L2_STL),filename_L2_STL);
                    STL_L2_SSH=double(ncread(filename_L2_STL,'ssh_corr').');
                    STL_L2_sigma0=double(ncread(filename_L2_STL,'sigma0_L2').');
                    STL_L2_SWH=double(ncread(filename_L2_STL,'swh').');
                    
                end
                
            else
                input_L2_ISR_Files   = dir(fullfile(char(input_path_L2_ISR_bs(b)),['*' data_string(1:15) '*.nc']));
                filename_L2_ISR=char(input_L2_ISR_Files.name);
                clear input_L2_ISR_b2_Files;
            end
            filename_L2_ISR=strcat(char(input_path_L2_ISR_bs(b)),filename_L2_ISR);
            
            if ~cellfun(@isempty,input_path_L1_ISR_bs)
                input_L1_ISR_Files   = dir(fullfile(char(input_path_L1_ISR_bs(b)),['*' data_string(1:15) '*.nc']));
                filename_L1_ISR=char(input_L1_ISR_Files.name);
                filename_L1_ISR=strcat(char(input_path_L1_ISR_bs(b)),filename_L1_ISR);
            end
            
            %% ------------- Reading the information for each file ----------------
            % ********************* BASELINE -1 ***************************************
            if ~isempty(strfind(lower(char(name_bs(b))),'acdc'))
                %------------------- ACDC product -------------------------------------
                %---------------- Geometry variables --------------------------------------
                %Assuming they should have the same surfaces b1 and b2
                if b==1
                    ISR_lat_surf=double(ncread(filename_L2_ISR,'lat_l1b_echo_sar_ku')).';
                    ISR_lon_surf=double(ncread(filename_L2_ISR,'lon_l1b_echo_sar_ku')).';
                end
                
                %------------- Geophysical parameters ---------------------------------
                SSH_dumm=double(ncread(filename_L2_ISR,'ssh_ACDC_20_ku')).';
                SWH_dumm=double(ncread(filename_L2_ISR,'swh_ACDC_20_ku')).';
                sigma0_dumm=double(ncread(filename_L2_ISR,'sig0_ACDC_20_ku')).';
                COR_dumm=double(ncread(filename_L2_ISR,'Pearson_corr_ACDC_20_ku')).';
                
                %------------ read waveforms information ----------------------
                if ~cellfun(@isempty,input_path_L1_ISR_bs)
                    ISD_i2q2_meas=double(ncread(filename_L2_ISR,'i2q2_meas_ku_ACDC_echo_sar_ku').');
                    ISD_waveform_scale=ncread(filename_L2_ISR,'waveform_scale_factor_ACDC_echo_sar_ku');
                    ISD_N_samples=length(ISD_i2q2_meas(1,:));
                    %apply scaling factor to waveforms
                    ISD_i2q2_meas=ISD_i2q2_meas.*repmat(ISD_waveform_scale,1,ISD_N_samples);
                    max_wvfm_dumm=(max(ISD_i2q2_meas,[],2)).';
                    clear ISD_i2q2_meas ISD_waveform_scale;
                    Pu_dumm=10.^(double(ncread(filename_L2_ISR,'Pu_analytical_ACDC_20_ku')).'./10);
                    amp_fit_dumm=Pu_dumm./max_wvfm_dumm;
                    if b==1
                        hr=double(ncread(filename_L1_ISR,'orb_alt_rate_l1b_echo_sar_ku')).'; %radial velocity
                    end
                    nb_dumm=double(ncread(filename_L1_ISR,'nb_stack_l1b_echo_sar_ku')).'; %number of beams
                end
                
            else
                %--------------------- L2 product -------------------------------------
                %---------------- Geometry variables ----------------------------------
                %Assuming they should have the same surfaces b1 and b2
                if b==1
                    ISR_lat_surf=double(ncread(filename_L2_ISR,'lat_20_ku')).';
                    ISR_lon_surf=double(ncread(filename_L2_ISR,'lon_20_ku')).';
                end
                
                %------------- Geophysical parameters ---------------------------------
                SSH_dumm=double(ncread(filename_L2_ISR,strcat(sh_name_nc,'_analytical_SWH_MSSfixed_20_ku'))).';
                SWH_dumm=double(ncread(filename_L2_ISR,'swh_analytical_SWH_MSSfixed_20_ku')).';
                sigma0_dumm=double(ncread(filename_L2_ISR,'sig0_analytical_SWH_MSSfixed_20_ku')).';
                COR_dumm=double(ncread(filename_L2_ISR,'Pearson_corr_analytical_SWH_MSSfixed_20_ku')).';
                
                
                %------------ read waveforms information ----------------------
                if ~cellfun(@isempty,input_path_L1_ISR_bs)
                    ISD_i2q2_meas=double(ncread(filename_L1_ISR,'i2q2_meas_ku_l1b_echo_sar_ku').');
                    ISD_waveform_scale=ncread(filename_L1_ISR,'waveform_scale_factor_l1b_echo_sar_ku');
                    ISD_N_samples=length(ISD_i2q2_meas(1,:));
                    ISD_i2q2_meas=ISD_i2q2_meas.*repmat(ISD_waveform_scale,1,ISD_N_samples);
                    max_wvfm_dumm=(max(ISD_i2q2_meas,[],2)).';
                    clear ISD_i2q2_meas ISD_waveform_scale;
                    Pu_dumm=10.^(double(ncread(filename_L2_ISR,'Pu_analytical_SWH_MSSfixed_20_ku')).'./10);
                    amp_fit_dumm=Pu_dumm./max_wvfm_dumm;
                    if b==1
                        hr=double(ncread(filename_L1_ISR,'orb_alt_rate_l1b_echo_sar_ku')).'; %radial velocity
                    end
                    nb_dumm=double(ncread(filename_L1_ISR,'nb_stack_l1b_echo_sar_ku')).'; %number of beams
                end
            end
            ISD_num_surfaces=length(ISR_lat_surf);
            
            
            %% ----------------- Algin ESA and ISR arrays of parameters -------
            if b==1
                %---------- Filtering by mask ---------------------------------
                if ~isempty(filename_mask_KML)
                    ISD_indexes_int=inpolygon(ISR_lon_surf,ISR_lat_surf,geo_mask.coord(:,1),geo_mask.coord(:,2));
                else
                    ISD_indexes_int=logical(ones(1,ISD_num_surfaces));
                end
                
                if ~isempty(input_path_L2_ESA)
                    if isempty(filename_mask_KML)
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
                    else
                        %Forcing same number of surfaces ESA/ISD for comparison: taking
                        %the first closest surface ESA to first ISD within the mask
                        [~,idx_min]=min(abs(ESA_lat_surf-ISR_lat_surf(find(ISD_indexes_int,1,'first'))));
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
                        
                    end
                    
                    ISD_num_surfaces=length(find(ISD_indexes_int)==1);
                    ESA_num_surfaces=length(find(ESA_indexes_int)==1);
                    
                    idx_int_ESA=find(ESA_indexes_int);
                    idx_int_ISD=find(ISD_indexes_int);
                    
                    % ESA replacing values to meet same number surfaces as ISR
                    %-------------- Geometry parameters ---------------------------------------
                    ESA_lat_surf=ESA_lat_surf(ESA_indexes_int);
                    ESA_lon_surf=ESA_lon_surf(ESA_indexes_int);
                    if ~cellfun(@isempty,input_path_L1_ISR_bs)
                        ESA_hr=ESA_hr(ESA_indexes_int);
                        ESA_nb=ESA_nb(ESA_indexes_int);
                    end
                    
                    
                    %-------------- Geophysical parameters ------------------------------------
                    ESA_L2_SSH_r1=ESA_L2_SSH_r1(ESA_indexes_int);
                    ESA_L2_sigma0_r1=ESA_L2_sigma0_r1(ESA_indexes_int);
                end
                
                
                ISR_lat_surf=ISR_lat_surf(ISD_indexes_int);
                ISR_lon_surf=ISR_lon_surf(ISD_indexes_int);
                
                
                num_surfaces_whole(i_files_valid)=ISD_num_surfaces;
                
                %check STL data avialbel
                if ~isempty(input_path_L2_STL)
                    STL_L2_SSH=STL_L2_SSH(ISD_indexes_int);
                    STL_L2_sigma0=STL_L2_sigma0(ISD_indexes_int);
                    STL_L2_SWH=STL_L2_SWH(ISD_indexes_int);
                end
                
            end
            
            %load the common surfaces of interest for ISD
            %------------- Geophysical parameters ---------------------------------
            SSH(b,:)=SSH_dumm(ISD_indexes_int);
            SWH(b,:)=SWH_dumm(ISD_indexes_int);
            sigma0(b,:)=sigma0_dumm(ISD_indexes_int);
            COR(b,:)=COR_dumm(ISD_indexes_int);
            if ~cellfun(@isempty,input_path_L1_ISR_bs)
                max_wvfm(b,:)=max_wvfm_dumm(ISD_indexes_int);
                Pu(b,:)=Pu_dumm(ISD_indexes_int);
                amp_fit(b,:)=amp_fit_dumm(ISD_indexes_int);
                
                if b==1
                    hr=hr(ISD_indexes_int); %radial velocity
                end
                nb(b,:)=nb_dumm(ISD_indexes_int); %number of beams
                
                clear max_wvfm_dumm Pu_dumm amp_fit_dumm;
            end
            clear SSH_dumm SWH_dumm sigma0_dumm COR_dumm;
            
            
            
            
            %% ------------------- Outliers filtering -----------------------------
            %----------------------------------------------------------------------
            if flag_outliers_removal
                switch type_outliers_removal
                    case 'percentiles'
                        if ~isempty(input_path_L2_ESA) && b==1
                            %------------------- SSH ----------------------------------
                            [ESA_L2_SSH_r1,idx_outliers_SSH_ESA]=outliers_by_percentiles(ESA_L2_SSH_r1,outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_SSH_ESA=find(~idx_outliers_SSH_ESA);
                            %----------------- sigma0 ---------------------------------
                            [ESA_L2_sigma0_r1,idx_outliers_sigma0_ESA]=outliers_by_percentiles(ESA_L2_sigma0_r1,outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_sigma0_ESA=find(~idx_outliers_sigma0_ESA);
                        end
                        
                        if ~isempty(input_path_L2_STL) && b==1
                            %------------------- SSH ----------------------------------
                            [STL_L2_SSH,idx_outliers_SSH_STL]=outliers_by_percentiles(STL_L2_SSH,outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_SSH_STL=find(~idx_outliers_SSH_STL);
                            %----------------- sigma0 ---------------------------------
                            [STL_L2_sigma0,idx_outliers_sigma0_STL]=outliers_by_percentiles(STL_L2_sigma0,outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_sigma0_STL=find(~idx_outliers_sigma0_STL);
                            %----------------- SWH ------------------------------------
                            [STL_L2_SWH,idx_outliers_SWH_STL]=outliers_by_percentiles(STL_L2_SWH,outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_SWH_STL=find(~idx_outliers_SWH_STL);
                        end
                        
                        % compute the errors w.r.t fitting on the data using a smooth
                        % function
                        %------------------- SSH ----------------------------------
                        [SSH(b,:),idx_outliers_SSH(b,:)]=outliers_by_percentiles(SSH(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        idx_nooutliers_SSH=find(~idx_outliers_SSH(b,:));
                        %----------------- sigma0 ---------------------------------
                        [sigma0(b,:),idx_outliers_sigma0(b,:)]=outliers_by_percentiles(sigma0(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        idx_nooutliers_sigma0=find(~idx_outliers_sigma0(b,:));
                        %----------------- SWH ------------------------------------
                        [SWH(b,:),idx_outliers_SWH(b,:)]=outliers_by_percentiles(SWH(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        idx_nooutliers_SWH=find(~idx_outliers_SWH(b,:));
                        %----------------- COR ------------------------------------
                        [COR(b,:),idx_outliers_COR(b,:)]=outliers_by_percentiles(COR(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                        idx_nooutliers_COR=find(~idx_outliers_COR(b,:));
                        
                        if ~cellfun(@isempty,input_path_L1_ISR_bs)
                            %----------------- amplitude fit ----------------------
                            [amp_fit(b,:),idx_outliers_amp_fit(b,:)]=outliers_by_percentiles(amp_fit(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_amp_fit=find(~idx_outliers_amp_fit(b,:));
                            %---------------- max waveforms -----------------------
                            [max_wvfm(b,:),idx_outliers_max_wvfm(b,:)]=outliers_by_percentiles(max_wvfm(b,:),outlier_percentil_low,outlier_percentil_high,IQR_times);
                            idx_nooutliers_max_wvfm=find(~idx_outliers_max_wvfm(b,:));
                        end
                        
                    case 'hampel'
                        if ~isempty(input_path_L2_ESA) && b==1
                            %------------------ SSH ---------------------------
                            [~,idx_outliers_SSH_ESA] = hampel(ESA_lat_surf,ESA_L2_SSH_r1,hampel_wind,hampel_sigma);
                            ESA_L2_SSH_r1(idx_outliers_SSH_ESA)=NaN;
                            idx_nooutliers_SSH_ESA=find(~idx_outliers_SSH_ESA);
                            %----------------- sigma0 ---------------------------------
                            [~,idx_outliers_sigma0_ESA] = hampel(ESA_lat_surf,ESA_L2_sigma0_r1,hampel_wind,hampel_sigma);
                            ESA_L2_sigma0_r1(idx_outliers_sigma0_ESA)=NaN;
                            idx_nooutliers_sigma0_ESA=find(~idx_outliers_sigma0_ESA);
                        end
                        
                        if ~isempty(input_path_L2_STL) && b==1
                            %------------------ SSH ---------------------------
                            [~,idx_outliers_SSH_STL] = hampel(ISR_lat_surf,STL_L2_SSH,hampel_wind,hampel_sigma);
                            STL_L2_SSH(idx_outliers_SSH_STL)=NaN;
                            idx_nooutliers_SSH_STL=find(~idx_outliers_SSH_STL);
                            %----------------- sigma0 ---------------------------------
                            [~,idx_outliers_sigma0_STL] = hampel(ISR_lat_surf,STL_L2_sigma0,hampel_wind,hampel_sigma);
                            STL_L2_sigma0(idx_outliers_sigma0_STL)=NaN;
                            idx_nooutliers_sigma0_STL=find(~idx_outliers_sigma0_STL);
                            %----------------- SWH ------------------------------------
                            [~,idx_outliers_SWH_STL] = hampel(ISR_lat_surf,STL_L2_SWH,hampel_wind,hampel_sigma);
                            STL_L2_SWH(idx_outliers_SWH_STL)=NaN;
                            idx_nooutliers_SWH_STL=find(~idx_outliers_SWH_STL);
                        end
                        
                        %------------------ SSH ---------------------------
                        [~,idx_outliers_SSH(b,:)] = hampel(ISR_lat_surf,SSH(b,:),hampel_wind,hampel_sigma);
                        SSH(b,idx_outliers_SSH(b,:))=NaN;
                        idx_nooutliers_SSH=find(~idx_outliers_SSH(b,:));
                        %----------------- sigma0 ---------------------------------
                        [~,idx_outliers_sigma0(b,:)] = hampel(ISR_lat_surf,sigma0(b,:),hampel_wind,hampel_sigma);
                        sigma0(b,idx_outliers_sigma0(b,:))=NaN;
                        idx_nooutliers_sigma0=find(~idx_outliers_sigma0(b,:));
                        %----------------- SWH ------------------------------------
                        [~,idx_outliers_SWH(b,:)] = hampel(ISR_lat_surf,SWH(b,:),hampel_wind,hampel_sigma);
                        SWH(b,idx_outliers_SWH(b,:))=NaN;
                        idx_nooutliers_SWH=find(~idx_outliers_SWH(b,:));
                        %----------------- COR ------------------------------------
                        [~,idx_outliers_COR(b,:)] = hampel(ISR_lat_surf,COR(b,:),hampel_wind,hampel_sigma);
                        COR(b,idx_outliers_COR(b,:))=NaN;
                        idx_nooutliers_COR=find(~idx_outliers_COR(b,:));
                        
                        if ~cellfun(@isempty,input_path_L1_ISR_bs)
                            %----------------- amplitude fit ----------------------
                            [~,idx_outliers_amp_fit(b,:)] = hampel(ISR_lat_surf,amp_fit(b,:),hampel_wind,hampel_sigma);
                            amp_fit(b,idx_outliers_amp_fit(b,:))=NaN;
                            idx_nooutliers_amp_fit=find(~idx_outliers_amp_fit(b,:));
                            
                            %----------------- max wvfm ---------------------------
                            [~,idx_outliers_max_wvfm(b,:)] = hampel(ISR_lat_surf,max_wvfm(b,:),hampel_wind,hampel_sigma);
                            max_wvfm(b,idx_outliers_max_wvfm(b,:))=NaN;
                            idx_nooutliers_max_wvfm=find(~idx_outliers_max_wvfm(b,:));
                            
                            %----------------- Pu ---------------------------
                            [~,idx_outliers_Pu(b,:)] = hampel(ISR_lat_surf,Pu(b,:),hampel_wind,hampel_sigma);
                            Pu(b,idx_outliers_Pu(b,:))=NaN;
                            idx_nooutliers_Pu=find(~idx_outliers_Pu(b,:));
                        end
                end
            end
            
            %% --------------- Smoothing geophysical retrievals ---------------
            if smooth_param
                
                if ~isempty(input_path_L2_ESA) && b==1
                    %----------------- SSH ------------------------------------
                    ESA_L2_SSH_r1_smooth=smooth(ESA_L2_SSH_r1,win_size_detrending).';
                    %----------------- SIGMA0 ---------------------------------
                    ESA_L2_sigma0_r1_smooth=smooth(ESA_L2_sigma0_r1,win_size_detrending).';
                end
                
                if ~isempty(input_path_L2_STL) && b==1
                    %----------------- SSH ------------------------------------
                    STL_L2_SSH_smooth=smooth(STL_L2_SSH,win_size_detrending).';
                    %----------------- SIGMA0 ---------------------------------
                    STL_L2_sigma0_smooth=smooth(STL_L2_sigma0,win_size_detrending).';
                    %----------------- SWH ---------------------------------
                    STL_L2_SWH_smooth=smooth(STL_L2_SWH,win_size_detrending).';
                end
                
                %------------------ SSH ---------------------------
                SSH_smooth(b,:) = smooth(SSH(b,:),win_size_detrending).';
                %----------------- sigma0 ---------------------------------
                sigma0_smooth(b,:) = smooth(sigma0(b,:),win_size_detrending).';
                %----------------- SWH ------------------------------------
                SWH_smooth(b,:) = smooth(SWH(b,:),win_size_detrending).';
                %----------------- COR ------------------------------------
                COR_smooth(b,:) = smooth(COR(b,:),win_size_detrending).';
            end
            
            
            %% ------------------- Compute the std BLOCK-WISE ---------------------
            %----------------------------------------------------------------------
            num_boxes=floor(ISD_num_surfaces/win_size_detrending);
            for i_box=1:(num_boxes+1)
                init_sample=max([(i_box-1)*win_size_detrending+1,1]);
                last_sample=min([(i_box-1)*win_size_detrending+win_size_detrending,ISD_num_surfaces]);
                SSH_std(b,i_box)=nanstd(detrend(SSH(b,idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
                SWH_std(b,i_box)=nanstd(detrend(SWH(b,idx_nooutliers_SWH(idx_nooutliers_SWH>=init_sample & idx_nooutliers_SWH<=last_sample))));
                sigma0_std(b,i_box)=nanstd(detrend(sigma0(b,idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));
                COR_std(b,i_box)=nanstd(detrend(COR(b,idx_nooutliers_COR(idx_nooutliers_COR>=init_sample & idx_nooutliers_COR<=last_sample))));
                
                
                
                SSH_mean(b,i_box)=nanmean((SSH(b,idx_nooutliers_SSH(idx_nooutliers_SSH>=init_sample & idx_nooutliers_SSH<=last_sample))));
                SWH_mean(b,i_box)=nanmean((SWH(b,idx_nooutliers_SWH(idx_nooutliers_SWH>=init_sample & idx_nooutliers_SWH<=last_sample))));
                sigma0_mean(b,i_box)=nanmean((sigma0(b,idx_nooutliers_sigma0(idx_nooutliers_sigma0>=init_sample & idx_nooutliers_sigma0<=last_sample))));
                COR_mean(b,i_box)=nanmean((COR(b,idx_nooutliers_COR(idx_nooutliers_COR>=init_sample & idx_nooutliers_COR<=last_sample))));
                
                if ~cellfun(@isempty,input_path_L1_ISR_bs)
                    amp_fit_std(b,i_box)=nanstd(detrend(amp_fit(b,idx_nooutliers_amp_fit(idx_nooutliers_amp_fit>=init_sample & idx_nooutliers_amp_fit<=last_sample))));
                    max_wvfm_std(b,i_box)=nanstd(detrend(max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample))));
                    max_wvfm_dB_std(b,i_box)=nanstd(detrend(10*log10(max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample)))));
                    Pu_std(b,i_box)=nanstd(detrend(Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample))));
                    Pu_dB_std(b,i_box)=nanstd(detrend(10*log10(Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample)))));
                    Pu_mean(b,i_box)=nanmean((Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample))));
                    Pu_dB_mean(b,i_box)=nanmean(10*log10(Pu(b,idx_nooutliers_Pu(idx_nooutliers_Pu>=init_sample & idx_nooutliers_Pu<=last_sample))));
                    
                    
                    
                    amp_fit_mean(b,i_box)=nanmean((amp_fit(b,idx_nooutliers_amp_fit(idx_nooutliers_amp_fit>=init_sample & idx_nooutliers_amp_fit<=last_sample))));
                    max_wvfm_mean(b,i_box)=nanmean((max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample))));
                    max_wvfm_dB_mean(b,i_box)=nanmean(10*log10(max_wvfm(b,idx_nooutliers_max_wvfm(idx_nooutliers_max_wvfm>=init_sample & idx_nooutliers_max_wvfm<=last_sample))));
                    
                    if b==1
                        hr_mean(i_box)=nanmean((hr(init_sample:last_sample)));
                    end
                    nb_mean(b,i_box)=nanmean((nb(b,init_sample:last_sample)));
                end
                
                if b==1
                    LAT_median(i_box)=nanmedian((ISR_lat_surf(init_sample:last_sample)));
                    LON_median(i_box)=nanmedian((ISR_lon_surf(init_sample:last_sample)));
                end
                
                
                if ~isempty(input_path_L2_ESA) && b==1
                    ESA_SSH_std(i_box)=nanstd(detrend(ESA_L2_SSH_r1(idx_nooutliers_SSH_ESA(idx_nooutliers_SSH_ESA>=init_sample & idx_nooutliers_SSH_ESA<=last_sample))));
                    ESA_sigma0_std(i_box)=nanstd(detrend(ESA_L2_sigma0_r1(idx_nooutliers_sigma0_ESA(idx_nooutliers_sigma0_ESA>=init_sample & idx_nooutliers_sigma0_ESA<=last_sample))));
                    ESA_SSH_mean(i_box)=nanmean((ESA_L2_SSH_r1(idx_nooutliers_SSH_ESA(idx_nooutliers_SSH_ESA>=init_sample & idx_nooutliers_SSH_ESA<=last_sample))));
                    ESA_sigma0_mean(i_box)=nanmean((ESA_L2_sigma0_r1(idx_nooutliers_sigma0_ESA(idx_nooutliers_sigma0_ESA>=init_sample & idx_nooutliers_sigma0_ESA<=last_sample))));
                end
                
                
                
                if ~isempty(input_path_L2_STL) && b==1
                    STL_SSH_std(i_box)=nanstd(detrend(STL_L2_SSH(idx_nooutliers_SSH_STL(idx_nooutliers_SSH_STL>=init_sample & idx_nooutliers_SSH_STL<=last_sample))));
                    STL_sigma0_std(i_box)=nanstd(detrend(STL_L2_sigma0(idx_nooutliers_sigma0_STL(idx_nooutliers_sigma0_STL>=init_sample & idx_nooutliers_sigma0_STL<=last_sample))));
                    STL_SWH_std(i_box)=nanstd(detrend(STL_L2_SWH(idx_nooutliers_SWH_STL(idx_nooutliers_SWH_STL>=init_sample & idx_nooutliers_SWH_STL<=last_sample))));
                    
                    STL_SSH_mean(i_box)=nanmean((STL_L2_SSH(idx_nooutliers_SSH_STL(idx_nooutliers_SSH_STL>=init_sample & idx_nooutliers_SSH_STL<=last_sample))));
                    STL_sigma0_mean(i_box)=nanmean((STL_L2_sigma0(idx_nooutliers_sigma0_STL(idx_nooutliers_sigma0_STL>=init_sample & idx_nooutliers_sigma0_STL<=last_sample))));
                    STL_SWH_mean(i_box)=nanmean((STL_L2_SWH(idx_nooutliers_SWH_STL(idx_nooutliers_SWH_STL>=init_sample & idx_nooutliers_SWH_STL<=last_sample))));
                end
                
                if ~isempty(input_path_L1_ESA) && b==1
                    ESA_hr_mean(i_box)=nanmean((ESA_hr(init_sample:last_sample)));
                    ESA_nb_mean(i_box)=nanmean((ESA_nb(init_sample:last_sample)));
                end
                
                
                
            end
            if define_min_max_SWH
                max_SWH=max([max_SWH,nanmax(SWH(:))]);
                min_SWH=min([min_SWH,nanmin(SWH(:))]);
            end
            
            %% --------- Compute the RMSE FITTING ---------------------------------
            if smooth_param
                
                if ~isempty(input_path_L2_ESA) && b==1
                    ESA_SSH_RMSE=sqrt(nanmean((ESA_L2_SSH_r1-ESA_L2_SSH_r1_smooth).^2));
                    ESA_sigma0_RMSE=sqrt(nanmean((ESA_L2_sigma0_r1-ESA_L2_sigma0_r1_smooth).^2));
                end
                
                if ~isempty(input_path_L2_STL) && b==1
                    STL_SSH_RMSE=sqrt(nanmean((STL_L2_SSH-STL_L2_SSH_smooth).^2));
                    STL_sigma0_RMSE=sqrt(nanmean((STL_L2_sigma0-STL_L2_sigma0_smooth).^2));
                    STL_SWH_RMSE=sqrt(nanmean((STL_L2_SWH-STL_L2_SWH_smooth).^2));
                end
                
                SSH_RMSE(b,1)=sqrt(nanmean((SSH(b,:)-SSH_smooth(b,:)).^2));
                SWH_RMSE(b,1)=sqrt(nanmean((SWH(b,:)-SWH_smooth(b,:)).^2));
                sigma0_RMSE(b,1)=sqrt(nanmean((sigma0(b,:)-sigma0_smooth(b,:)).^2));
                COR_RMSE(b,1)=sqrt(nanmean((COR(b,:)-COR_smooth(b,:)).^2));
            end
            
            
            %% --------- Compute the mean std over track equivalent to RMSE--------
            if ~isempty(input_path_L2_ESA) && b==1
                ESA_SSH_std_mean=nanmean(ESA_SSH_std);
                ESA_sigma0_std_mean=nanmean(ESA_sigma0_std);
            end
            
            if ~isempty(input_path_L2_STL) && b==1
                STL_SSH_std_mean=nanmean(STL_SSH_std);
                STL_sigma0_std_mean=nanmean(STL_sigma0_std);
                STL_SWH_std_mean=nanmean(STL_SWH_std);
            end
            
            SSH_std_mean(b,1)=nanmean(SSH_std(b,:));
            SWH_std_mean(b,1)=nanmean(SWH_std(b,:));
            sigma0_std_mean(b,1)=nanmean(sigma0_std(b,:));
            COR_std_mean(b,1)=nanmean(COR_std(b,:));
            
            
            
            
            %% ------------------ Compute the mean errors ESA/SLT-ISD--------------
            if ~isempty(input_path_L2_ESA)
                SSH_mean_error_ESA_ISR(b,1)=nanmean(ESA_SSH_mean-SSH_mean(b,:));
                sigma0_mean_error_ESA_ISR(b,1)=nanmean(ESA_sigma0_mean-sigma0_mean(b,:));
            end
            
            if ~isempty(input_path_L2_STL)
                SSH_mean_error_STL_ISR(b,1)=nanmean(STL_SSH_mean-SSH_mean(b,:));
                SWH_mean_error_STL_ISR(b,1)=nanmean(STL_SWH_mean-SWH_mean(b,:));
                sigma0_mean_error_STL_ISR(b,1)=nanmean(STL_sigma0_mean-sigma0_mean(b,:));
            end
            clear idx_nooutliers_SSH idx_nooutliers_SWH idx_nooutliers_sigma0 idx_nooutliers_COR;
            clear idx_nooutliers_SSH_ESA idx_nooutliers_sigma0_ESA;
            clear idx_nooutliers_SSH_STL idx_nooutliers_sigma0_STL idx_nooutliers_SWH_STL;
            clear idx_nooutliers_amp_fit idx_nooutliers_max_wvfm idx_nooutliers_Pu;
        end %loop over the channels
        max_SWH=max([max_SWH,nanmax(SWH(:))]);
        min_SWH=min([min_SWH,nanmin(SWH(:))]);
        
        
        
        %% ---------------- PLOTING RESULTS -----------------------------------
        %----------------------------------------------------------------------
        %--------------------- SSH --------------------------------------------
        if generate_plots==1
            if  mod(i_files_valid,plot_fits_downsampling)==0 || i_files_valid==1
                f1=figure;
                legend_text={''};
                text_in_textbox={''};
                text_errors_ESA_ISR={''};
                text_errors_STL_ISR={''};
                if ~isempty(input_path_L2_ESA)
                    plot(ESA_lat_surf,ESA_L2_SSH_r1,'Marker',marker_ESA,'Color',color_ESA);
                    hold on;
                    legend_text=[legend_text,'ESA'];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'ESA: RMSE [m]: '},num2str(ESA_SSH_RMSE,'%.4g'),{' , std [m]: '},num2str(ESA_SSH_std_mean,'%.4g'))];
                end
                
                if ~isempty(input_path_L2_STL)
                    
                    plot(ISR_lat_surf,STL_L2_SSH,'Marker',marker_STL,'Color',color_STL);
                    hold on;
                    legend_text=[legend_text,'STL'];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'STL: RMSE [m]: '},num2str(STL_SSH_RMSE,'%.4g'),{' , std [m]: '},num2str(STL_SSH_std_mean,'%.4g'))];
                end
                
                for b=1:N_baselines
                    plot(ISR_lat_surf,SSH(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
                    hold on;
                    legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'ISR '},name_bs(b),{': RMSE [m]: '},num2str(SSH_RMSE(b,1),'%.4g'),{' , std [m]: '},num2str(SSH_std_mean(b,1),'%.4g'))];
                    if ~isempty(input_path_L2_ESA)
                        text_errors_ESA_ISR=[text_errors_ESA_ISR,...
                            strcat({'Mean error ESA-ISR '},name_bs(b),{' [m]: '},num2str(SSH_mean_error_ESA_ISR(b,1),'%.4g'))];
                    end
                    if ~isempty(input_path_L2_STL)
                        text_errors_STL_ISR=[text_errors_STL_ISR,...
                            strcat({'Mean error STL-ISR '},name_bs(b),{' [m]: '},num2str(SSH_mean_error_STL_ISR(b,1),'%.4g'))];
                    end
                end
                title(strcat(upper(sh_name_nc),': Baselines comparison'));
                leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
                pos_leg=get(leg,'Position');
                xlabel('Latitude [deg.]'); ylabel(strcat(upper(sh_name_nc),' [m]'));
                text_in_textbox=[text_in_textbox,text_errors_ESA_ISR,text_errors_STL_ISR];
                text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
                
                if annotation_box_active
                    pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),pos_leg(4)],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
                end
                print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SSH',file_ext))
                save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SSH','.mat'),...
                    'input_path_L2_ESA','ESA_lat_surf','ESA_L2_SSH_r1','marker_ESA','color_ESA',...
                    'input_path_L2_STL','ISR_lat_surf','STL_L2_SSH','marker_STL','color_STL',...
                    'N_baselines','SSH','marker_bs','color_bs',...
                    'name_bs','sh_name_nc','legend_text','text_in_textbox');
                close(f1)
                %--------------------- SWH --------------------------------------------
                f1=figure;
                legend_text={''};
                text_in_textbox={''};
                text_errors_STL_ISR={''};
                if ~isempty(input_path_L2_STL)
                    plot(ISR_lat_surf,STL_L2_SWH,'Marker',marker_STL,'Color',color_STL);
                    hold on;
                    legend_text=[legend_text,'STL'];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'STL: RMSE [m]: '},num2str(STL_SWH_RMSE,'%.4g'),{' , std [m]: '},num2str(STL_SWH_std_mean,'%.4g'))];
                end
                for b=1:N_baselines
                    plot(ISR_lat_surf,SWH(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
                    hold on;
                    legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'ISR '},name_bs(b),{': RMSE [m]: '},num2str(SWH_RMSE(b,1),'%.4g'),{' , std [m]: '},num2str(SWH_std_mean(b,1),'%.4g'))];
                    if ~isempty(input_path_L2_STL)
                        text_errors_STL_ISR=[text_errors_STL_ISR,...
                            strcat({'Mean error STL-ISR '},name_bs(b),{' [m]: '},num2str(SWH_mean_error_STL_ISR(b,1),'%.4g'))];
                    end
                end
                title(strcat('SWH (H_s): Baselines comparison'));
                leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
                pos_leg=get(leg,'Position');
                xlabel('Latitude [deg.]'); ylabel('SWH [m]');
                text_in_textbox=[text_in_textbox,text_errors_STL_ISR];
                text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
                if annotation_box_active
                    pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
                end
                print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SWH',file_ext))
                save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_SWH','.mat'),...
                    'input_path_L2_STL','ISR_lat_surf','STL_L2_SWH','marker_STL','color_STL',...
                    'N_baselines','SWH','marker_bs','color_bs',...
                    'name_bs','sh_name_nc','legend_text','text_in_textbox',...
                    'path_comparison_results','data_string');
                close(f1);
                
                %--------------------- sigma0 --------------------------------------------
                f1=figure;
                legend_text={''};
                text_in_textbox={''};
                text_errors_ESA_ISR={''};
                text_errors_STL_ISR={''};
                if ~isempty(input_path_L2_ESA)
                    
                    plot(ESA_lat_surf,ESA_L2_sigma0_r1,'Marker',marker_ESA,'Color',color_ESA);
                    hold on;
                    legend_text=[legend_text,'ESA'];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'ESA: RMSE [dB]: '},num2str(ESA_sigma0_RMSE,'%.4g'),{' , std [dB]: '},num2str(ESA_sigma0_std_mean,'%.4g'))];
                end
                
                if ~isempty(input_path_L2_STL)
                    
                    plot(ISR_lat_surf,STL_L2_sigma0,'Marker',marker_STL,'Color',color_STL);
                    hold on;
                    legend_text=[legend_text,'STL'];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'STL: RMSE [dB]: '},num2str(STL_sigma0_RMSE,'%.4g'),{' , std [dB]: '},num2str(STL_sigma0_std_mean,'%.4g'))];
                end
                
                for b=1:N_baselines
                    plot(ISR_lat_surf,sigma0(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
                    hold on;
                    legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'ISR '},name_bs(b),{': RMSE [dB]: '},num2str(sigma0_RMSE(b,1),'%.4g'),{' , std [dB]: '},num2str(sigma0_std_mean(b,1),'%.4g'))];
                    
                    if ~isempty(input_path_L2_ESA)
                        text_errors_ESA_ISR=[text_errors_ESA_ISR,...
                            strcat({'Mean error ESA-ISR '},name_bs(b),{' [dB]: '},num2str(sigma0_mean_error_ESA_ISR(b,1),'%.4g'))];
                    end
                    if ~isempty(input_path_L2_STL)
                        text_errors_STL_ISR=[text_errors_STL_ISR,...
                            strcat({'Mean error STL-ISR '},name_bs(b),{' [dB]: '},num2str(sigma0_mean_error_STL_ISR(b,1),'%.4g'))];
                    end
                    
                end
                title(strcat('\sigma^0: Baselines comparison'));
                leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
                pos_leg=get(leg,'Position');
                xlabel('Latitude [deg.]'); ylabel('\sigma^0 [dB]');
                text_in_textbox=[text_in_textbox,text_errors_ESA_ISR,text_errors_STL_ISR];
                text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
                
                if annotation_box_active
                    pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
                end
                print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_sigma0',file_ext))
                save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_sigma0','.mat'),...
                    'input_path_L2_ESA','ESA_lat_surf','ESA_L2_sigma0_r1','marker_ESA','color_ESA',...
                    'input_path_L2_STL','ISR_lat_surf','STL_L2_sigma0','marker_STL','color_STL',...
                    'N_baselines','sigma0','marker_bs','color_bs',...
                    'name_bs','sh_name_nc','legend_text','text_in_textbox',...
                    'path_comparison_results','data_string');
                close(f1)
                
                
                %--------------------- nb beams --------------------------------------------
                f1=figure;
                legend_text={''};
                text_in_textbox={''};
                text_errors_ESA_ISR={''};
                text_errors_STL_ISR={''};
                if ~isempty(input_path_L2_ESA)
                    plot(ESA_lat_surf,ESA_nb,'Marker',marker_ESA,'Color',color_ESA);
                    hold on;
                    legend_text=[legend_text,'ESA'];
                end
                
                for b=1:N_baselines
                    plot(ISR_lat_surf,nb(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
                    hold on;
                    legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
                    
                    if ~isempty(input_path_L2_ESA)
                        text_errors_ESA_ISR=[text_errors_ESA_ISR,...
                            strcat({'Mean error ESA-ISR '},name_bs(b),{' : '},num2str(nanmean(ESA_nb-nb(b,:)),'%.4g'))];
                    end
                    
                end
                title(strcat('Contributing Beams: Baselines comparison'));
                leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
                pos_leg=get(leg,'Position');
                xlabel('Latitude [deg.]'); ylabel('n_b');
                text_in_textbox=[text_in_textbox,text_errors_ESA_ISR];
                text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
                
                if annotation_box_active
                    pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
                end
                print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_nb',file_ext))
                save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_nb','.mat'),...
                    'input_path_L2_ESA','ESA_lat_surf','ESA_nb','marker_ESA','color_ESA',...
                    'N_baselines','nb','marker_bs','color_bs',...
                    'name_bs','sh_name_nc','legend_text','text_in_textbox',...
                    'path_comparison_results','data_string');
                close(f1)
                
                %--------------------- COR --------------------------------------------
                f1=figure;
                legend_text={''};
                text_in_textbox={''};
                
                for b=1:N_baselines
                    plot(ISR_lat_surf,COR(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
                    hold on;
                    legend_text=[legend_text,strcat({'ISR '},name_bs(b))];
                    text_in_textbox=[text_in_textbox,...
                        strcat({'ISR '},name_bs(b),{' --> RMSE [%]: '},num2str(COR_RMSE(b,1),'%.4g'),{' , std [%]: '},num2str(COR_std_mean(b,1),'%.4g'))];
                end
                title(strcat('\rho (Pearson correlation coeff.): Baselines comparison'));
                leg=legend(legend_text(~cellfun(@isempty,legend_text)),'Location','northeastoutside','Fontsize',legend_fontsize);
                pos_leg=get(leg,'Position');
                xlabel('Latitude [deg.]'); ylabel('\rho [%]');
                text_in_textbox=text_in_textbox(~cellfun(@isempty,text_in_textbox));
                if annotation_box_active
                    pos_a1=annotation('textbox',[pos_leg(1),pos_leg(2)-pos_leg(4)-0.15,pos_leg(3),0.1],'String',text_in_textbox,'FitBoxToText','on','FontSize',textbox_fontsize);
                end
                print(print_file,res_fig,strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_COR',file_ext))
                save(strcat(path_comparison_results,data_string,'_cmp_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),'_COR','.mat'),...
                    'N_baselines','COR','marker_bs','color_bs',...
                    'name_bs','sh_name_nc','legend_text','text_in_textbox',...
                    'path_comparison_results','data_string');
                close(f1);
            end
        end
        %% -------------- Combine the information -----------------------------
        %----------------------------------------------------------------------
        %------------ Combine the stds for each track -------------------------
        SSH_std_whole=[SSH_std_whole,SSH_std];
        
        SWH_std_whole=[SWH_std_whole,SWH_std];
        
        sigma0_std_whole=[sigma0_std_whole,sigma0_std];
        
        if ~cellfun(@isempty,input_path_L1_ISR_bs)
            amp_fit_std_whole=[amp_fit_std_whole,amp_fit_std];
            max_wvfm_std_whole=[max_wvfm_std_whole,max_wvfm_std];
        end
        
        if ~isempty(input_path_L2_ESA)
            ESA_SSH_std_whole=[ESA_SSH_std_whole,ESA_SSH_std];
            ESA_sigma0_std_whole=[ESA_sigma0_std_whole,ESA_sigma0_std];
        end
        
        if ~isempty(input_path_L2_STL)
            STL_SSH_std_whole=[STL_SSH_std_whole,STL_SSH_std];
            STL_sigma0_std_whole=[STL_sigma0_std_whole,STL_sigma0_std];
            STL_SWH_std_whole=[STL_SWH_std_whole,STL_SWH_std];
        end
        
        
        %------------ Combine the mean stds for each track --------------------
        SSH_std_mean_whole=[SSH_std_mean_whole,SSH_std_mean];
        
        SWH_std_mean_whole=[SWH_std_mean_whole,SWH_std_mean];
        
        sigma0_std_mean_whole=[sigma0_std_mean_whole,sigma0_std_mean];
        
        COR_std_mean_whole=[COR_std_mean_whole,COR_std_mean];
        
        if ~isempty(input_path_L2_ESA)
            ESA_SSH_std_mean_whole=[ESA_SSH_std_mean_whole,ESA_SSH_std_mean];
            ESA_sigma0_std_mean_whole=[ESA_sigma0_std_mean_whole,ESA_sigma0_std_mean];
        end
        
        if ~isempty(input_path_L2_STL)
            STL_SSH_std_mean_whole=[STL_SSH_std_mean_whole,STL_SSH_std_mean];
            STL_sigma0_std_mean_whole=[STL_sigma0_std_mean_whole,STL_sigma0_std_mean];
            STL_SWH_std_mean_whole=[STL_SWH_std_mean_whole,STL_SWH_std_mean];
        end
        
        
        
        %------------ Combine the mean for each track ---------------------
        SSH_mean_whole=[SSH_mean_whole,SSH_mean];
        
        SWH_mean_whole=[SWH_mean_whole,SWH_mean];
        
        sigma0_mean_whole=[sigma0_mean_whole,sigma0_mean];
        
        if ~cellfun(@isempty,input_path_L1_ISR_bs)
            amp_fit_mean_whole=[amp_fit_mean_whole,amp_fit_mean];
            max_wvfm_mean_whole=[max_wvfm_mean_whole,max_wvfm_mean];
            hr_mean_whole=[hr_mean_whole,hr_mean];
            nb_mean_whole=[nb_mean_whole,nb_mean];
        end
        
        COR_mean_whole=[COR_mean_whole,COR_mean];
        
        LAT_median_whole=[LAT_median_whole,LAT_median];
        LON_median_whole=[LON_median_whole,LON_median];
        
        if ~isempty(input_path_L2_ESA)
            ESA_SSH_mean_whole=[ESA_SSH_mean_whole,ESA_SSH_mean];
            ESA_sigma0_mean_whole=[ESA_sigma0_mean_whole,ESA_sigma0_mean];
        end
        
        if ~isempty(input_path_L2_STL)
            STL_SSH_mean_whole=[STL_SSH_mean_whole,STL_SSH_mean];
            STL_sigma0_mean_whole=[STL_sigma0_mean_whole,STL_sigma0_mean];
            STL_SWH_mean_whole=[STL_SWH_mean_whole,STL_SWH_mean];
        end
        
        %------------ Combine the retrieved param for each track --------------
        SSH_whole=[SSH_whole,SSH];
        SWH_whole=[SWH_whole,SWH];
        sigma0_whole=[sigma0_whole,sigma0];
        COR_whole=[COR_whole,COR];
        
        if ~isempty(input_path_L2_ESA)
            ESA_SSH_whole=[ESA_SSH_whole,ESA_L2_SSH_r1];
            ESA_sigma0_whole=[ESA_sigma0_whole,ESA_L2_sigma0_r1];
        end
        
        if ~isempty(input_path_L2_STL)
            STL_SSH_whole=[STL_SSH_whole,STL_L2_SSH];
            STL_sigma0_whole=[STL_sigma0_whole,STL_L2_sigma0];
            STL_SWH_whole=[STL_SWH_whole,STL_L2_SWH];
        end
        
        if ~cellfun(@isempty,input_path_L1_ISR_bs)
            amp_fit_whole=[amp_fit_whole,amp_fit];
            max_wvfm_whole=[max_wvfm_whole,max_wvfm];
            hr_whole=[hr_whole,hr];
            nb_whole=[nb_whole,nb];
        end
        
        if ~isempty(input_path_L1_ESA)
            ESA_hr_mean_whole=[ESA_hr_mean_whole,ESA_hr_mean];
            ESA_hr_whole=[ESA_hr_whole,ESA_hr];
            
            ESA_nb_mean_whole=[ESA_nb_mean_whole,ESA_nb_mean];
            ESA_nb_whole=[ESA_nb_whole,ESA_nb];
        end
        
        %smoothed version
        if smooth_param
            SSH_smooth_whole=[SSH_smooth_whole,SSH_smooth];
            SWH_smooth_whole=[SWH_smooth_whole,SWH_smooth];
            sigma0_smooth_whole=[sigma0_smooth_whole,sigma0_smooth];
            COR_smooth_whole=[COR_smooth_whole,COR_smooth];
            
            SSH_RMSE_whole=[SSH_RMSE_whole,SSH_RMSE];
            SWH_RMSE_whole=[SWH_RMSE_whole,SWH_RMSE];
            sigma0_RMSE_whole=[sigma0_RMSE_whole,sigma0_RMSE];
            COR_RMSE_whole=[COR_RMSE_whole,COR_RMSE];
            
            if ~isempty(input_path_L2_ESA)
                ESA_SSH_smooth_whole=[ESA_SSH_smooth_whole,ESA_L2_SSH_r1_smooth];
                ESA_sigma0_smooth_whole=[ESA_sigma0_smooth_whole,ESA_L2_sigma0_r1_smooth];
                
                ESA_SSH_RMSE_whole=[ESA_SSH_RMSE_whole,ESA_SSH_RMSE];
                ESA_sigma0_RMSE_whole=[ESA_sigma0_RMSE_whole,ESA_sigma0_RMSE];
            end
            if ~isempty(input_path_L2_STL)
                STL_SSH_smooth_whole=[STL_SSH_smooth_whole,STL_L2_SSH_smooth];
                STL_SWH_smooth_whole=[STL_SWH_smooth_whole,STL_L2_SWH_smooth];
                STL_sigma0_smooth_whole=[STL_sigma0_smooth_whole,STL_L2_sigma0_smooth];
                
                STL_SSH_RMSE_whole=[STL_SSH_RMSE_whole,STL_SSH_RMSE];
                STL_sigma0_RMSE_whole=[STL_sigma0_RMSE_whole,STL_sigma0_RMSE];
                STL_SWH_RMSE_whole=[STL_SWH_RMSE_whole,STL_SWH_RMSE];
            end
        end
        
        %----------------- Combine the mean errors ----------------------------
        %% ------------------ Compute the mean errors ESA/SLT-ISD--------------
        if ~isempty(input_path_L2_ESA)
            SSH_mean_error_ESA_ISR_whole=[SSH_mean_error_ESA_ISR_whole,SSH_mean_error_ESA_ISR];
            sigma0_mean_error_ESA_ISR_whole=[sigma0_mean_error_ESA_ISR_whole,sigma0_mean_error_ESA_ISR];
        end
        
        if ~isempty(input_path_L2_STL)
            SSH_mean_error_STL_ISR_whole=[SSH_mean_error_STL_ISR_whole,SSH_mean_error_STL_ISR];
            SWH_mean_error_STL_ISR_whole=[SWH_mean_error_STL_ISR_whole,SWH_mean_error_STL_ISR];
            sigma0_mean_error_STL_ISR_whole=[sigma0_mean_error_STL_ISR_whole,sigma0_mean_error_STL_ISR];
        end
        
        
        clear SSH SWH sigma0 SSH_smooth SWH_smooth sigma0_smooth
        clear ESA_L2_SSH_r1 ESA_L2_sigma0_r1 ESA_L2_sigma0_r1_smooth ESA_L2_sigma0_r1_smooth
        clear hr hr_mean nb nb_mean ESA_hr ESA_hr_mean ESA_nb ESA_nb_mean;
        clear STL_L2_SSH STL_L2_SWH STL_L2_sigma0 STL_L2_SSH_smooth STL_L2_SWH_smooth STL_L2_sigma0_smooth
        clear Pu amp_fit max_wvfm COR COR_smooth
        clear idx_outliers_SSH idx_outliers_SWH idx_outliers_sigma0 idx_outliers_COR idx_outliers_amp_fit idx_outliers_max_wvfm idx_outliers_Pu;
        clear idx_outliers_SSH_ESA idx_outliers_sigma0_ESA;
        clear idx_outliers_SSH_STL idx_outliers_sigma0_STL idx_outliers_SWH_STL;
        
        clear SSH_mean_error_ESA_ISR sigma0_mean_error_ESA_ISR SSH_mean_error_STL_ISR SWH_mean_error_STL_ISR sigma0_mean_error_STL_ISR;
        clear SSH_std_mean SWH_std_mean sigma0_std_mean COR_std_mean;
        clear ESA_SSH_std_mean ESA_sigma0_std_mean STL_SSH_std_mean STL_sigma0_std_mean STL_SWH_std_mean;
        clear ESA_SSH_RMSE ESA_sigma0_RMSE STL_SSH_RMSE STL_sigma0_RMSE STL_SWH_RMSE;
        clear SSH_RMSE SWH_RMSE sigma0_RMSE COR_RMSE;
        clear STL_SSH_std STL_sigma0_std STL_SWH_std STL_SSH_mean STL_sigma0_mean STL_SWH_mean ESA_hr_mean;
        clear ESA_SSH_std ESA_sigma0_std ESA_SSH_mean ESA_sigma0_mean;
        clear hr_mean LAT_median LON_median;
        clear amp_fit_mean max_wvfm_mean max_wvfm_dB_mean;
        clear amp_fit_std max_wvfm_std max_wvfm_dB_std Pu_std Pu_dB_std Pu_mean Pu_dB_mean;
        clear SSH_mean SWH_mean sigma0_mean COR_mean;
        clear SSH_std SWH_std sigma0_std COR_std;
        
        clear SSH_std SWH_std sigma0_std amp_fit_std max_wvfm_std SSH_mean SWH_mean sigma0_mean amp_fit_mean max_wvfm_mean COR_mean;
        clear ESA_SSH_std ESA_sigma0_std ESA_SSH_mean ESA_sigma0_mean;
        clear STL_SSH_std STL_sigma0_std STL_SWH_std STL_SSH_mean STL_sigma0_mean STL_SWH_mean;
        clear ISR_lat_surf ISR_lon_surf LAT_mean LON_mean;
        clear SSH_RMSE SWH_RMSE sigma0_RMSE COR_RMSE;
        clear ESA_SSH_RMSE ESA_sigma0_RMSE;
        clear STL_SSH_RMSE STL_SWH_RMSE STL_sigma0_RMSE;
        clear SSH_std_mean SWH_std_mean sigma0_std_mean COR_std_mean;
        clear ESA_SSH_std_mean ESA_sigma0_std_mean;
        clear STL_SSH_std_mean STL_SWH_std_mean STL_sigma0_std_mean;
        
        clear SSH_mean_error_ESA_ISR sigma0_mean_error_ESA_ISR;
        clear SSH_mean_error_STL_ISR sigma0_mean_error_STL_ISR SWH_mean_error_STL_ISR;
        
    catch
        disp(strcat('Some errors ocurred for:',data_string))
        clear SSH SWH sigma0 SSH_smooth SWH_smooth sigma0_smooth
        clear ESA_L2_SSH_r1 ESA_L2_sigma0_r1 ESA_L2_sigma0_r1_smooth ESA_L2_sigma0_r1_smooth
        clear hr hr_mean nb nb_mean ESA_hr ESA_hr_mean ESA_nb ESA_nb_mean;
        clear STL_L2_SSH STL_L2_SWH STL_L2_sigma0 STL_L2_SSH_smooth STL_L2_SWH_smooth STL_L2_sigma0_smooth
        clear Pu amp_fit max_wvfm COR COR_smooth
        clear idx_outliers_SSH idx_outliers_SWH idx_outliers_sigma0 idx_outliers_COR idx_outliers_amp_fit idx_outliers_max_wvfm idx_outliers_Pu;
        clear idx_outliers_SSH_ESA idx_outliers_sigma0_ESA;
        clear idx_outliers_SSH_STL idx_outliers_sigma0_STL idx_outliers_SWH_STL;
        
        clear SSH_mean_error_ESA_ISR sigma0_mean_error_ESA_ISR SSH_mean_error_STL_ISR SWH_mean_error_STL_ISR sigma0_mean_error_STL_ISR;
        clear SSH_std_mean SWH_std_mean sigma0_std_mean COR_std_mean;
        clear ESA_SSH_std_mean ESA_sigma0_std_mean STL_SSH_std_mean STL_sigma0_std_mean STL_SWH_std_mean;
        clear ESA_SSH_RMSE ESA_sigma0_RMSE STL_SSH_RMSE STL_sigma0_RMSE STL_SWH_RMSE;
        clear SSH_RMSE SWH_RMSE sigma0_RMSE COR_RMSE;
        clear STL_SSH_std STL_sigma0_std STL_SWH_std STL_SSH_mean STL_sigma0_mean STL_SWH_mean ESA_hr_mean;
        clear ESA_SSH_std ESA_sigma0_std ESA_SSH_mean ESA_sigma0_mean;
        clear hr_mean LAT_median LON_median;
        clear amp_fit_mean max_wvfm_mean max_wvfm_dB_mean;
        clear amp_fit_std max_wvfm_std max_wvfm_dB_std Pu_std Pu_dB_std Pu_mean Pu_dB_mean;
        clear SSH_mean SWH_mean sigma0_mean COR_mean;
        clear SSH_std SWH_std sigma0_std COR_std;
        
        clear SSH_std SWH_std sigma0_std amp_fit_std max_wvfm_std SSH_mean SWH_mean sigma0_mean amp_fit_mean max_wvfm_mean COR_mean;
        clear ESA_SSH_std ESA_sigma0_std ESA_SSH_mean ESA_sigma0_mean;
        clear STL_SSH_std STL_sigma0_std STL_SWH_std STL_SSH_mean STL_sigma0_mean STL_SWH_mean;
        clear ISR_lat_surf ISR_lon_surf LAT_mean LON_mean;
        clear SSH_RMSE SWH_RMSE sigma0_RMSE COR_RMSE;
        clear ESA_SSH_RMSE ESA_sigma0_RMSE;
        clear STL_SSH_RMSE STL_SWH_RMSE STL_sigma0_RMSE;
        clear SSH_std_mean SWH_std_mean sigma0_std_mean COR_std_mean;
        clear ESA_SSH_std_mean ESA_sigma0_std_mean;
        clear STL_SSH_std_mean STL_SWH_std_mean STL_sigma0_std_mean;
        
        clear SSH_mean_error_ESA_ISR sigma0_mean_error_ESA_ISR;
        clear SSH_mean_error_STL_ISR sigma0_mean_error_STL_ISR SWH_mean_error_STL_ISR;
        continue;
    end
        
end %loop over tracks
save(strcat(path_comparison_results,'Workspace_performance_comparison_baselines_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_')));



load('C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\results\Phase_2\test\2013\performance_comparison\agulhas\comparison_S3_old_S3_new\Workspace_performance_comparison_baselines_S_3_old_vs_S_3_new.mat');

%% ------------------- Noise performance vs SWH ---------------------------
%--------------------------------------------------------------------------
%-------------------- Compute edges of binning ----------------------------
edges_SWH=min_SWH:step_SWH:max_SWH;
n_total_bins=length(edges_SWH);

% Mean of std using specific binning on SWH basis of the stds
SSH_mean_std_binned_SWH=NaN(N_baselines,n_total_bins);
SWH_mean_std_binned_SWH=NaN(N_baselines,n_total_bins);
sigma0_mean_std_binned_SWH=NaN(N_baselines,n_total_bins);

if ~cellfun(@isempty,input_path_L1_ISR_bs)
    amp_fit_std_std_binned_SWH=NaN(N_baselines,n_total_bins);
    max_wvfm_std_std_binned_SWH=NaN(N_baselines,n_total_bins);
end

COR_mean_mean_binned_SWH=NaN(N_baselines,n_total_bins);


SSH_std_std_binned_SWH=NaN(N_baselines,n_total_bins);
SWH_std_std_binned_SWH=NaN(N_baselines,n_total_bins);
sigma0_std_std_binned_SWH=NaN(N_baselines,n_total_bins);


% if ~isempty(input_path_L2_ESA)
%     ESA_SSH_mean_std_binned_SWH=NaN(N_baselines,n_total_bins);
%     ESA_sigma0_mean_std_binned_SWH=NaN(N_baselines,n_total_bins); 
%     ESA_SSH_std_std_binned_SWH=NaN(N_baselines,n_total_bins);
%     ESA_sigma0_std_std_binned_SWH=NaN(N_baselines,n_total_bins);
% end

if ~isempty(input_path_L2_STL)
    STL_SSH_mean_std_binned_SWH=NaN(1,n_total_bins);
    STL_sigma0_mean_std_binned_SWH=NaN(1,n_total_bins); 
    STL_SWH_mean_std_binned_SWH=NaN(1,n_total_bins); 
    STL_SSH_std_std_binned_SWH=NaN(1,n_total_bins);
    STL_sigma0_std_std_binned_SWH=NaN(1,n_total_bins);
    STL_SWH_std_std_binned_SWH=NaN(1,n_total_bins);
end

%--------- Sorting mean SWH for the diferent boxes-------------------------
%taking as reference the conventional SWH
sorting_SWH=discretize(SWH_mean_whole(1,:),edges_SWH);
for b=1:N_baselines
    %sorting_SWH=discretize(SWH_mean_whole(b,:),edges_SWH);
    if ~isempty(input_path_L2_STL) && b==1
        %STL_sorting_SWH=discretize(STL_SWH_mean_whole,edges_SWH);
        STL_sorting_SWH=sorting_SWH;
    end
    for i_bin_SWH=1:n_total_bins        
        idx_int=sorting_SWH == i_bin_SWH;
        if ~isempty(idx_int)
            SSH_mean_std_binned_SWH(b,i_bin_SWH)=nanmean(SSH_std_whole(b,idx_int));
            SSH_std_std_binned_SWH(b,i_bin_SWH)=nanstd(SSH_std_whole(b,idx_int));
            
            SWH_mean_std_binned_SWH(b,i_bin_SWH)=nanmean(SWH_std_whole(b,idx_int));
            SWH_std_std_binned_SWH(b,i_bin_SWH)=nanstd(SWH_std_whole(b,idx_int));
            
            sigma0_mean_std_binned_SWH(b,i_bin_SWH)=nanmean(sigma0_std_whole(b,idx_int));
            sigma0_std_std_binned_SWH(b,i_bin_SWH)=nanstd(sigma0_std_whole(b,idx_int));
            
            if ~cellfun(@isempty,input_path_L1_ISR_bs)
                amp_fit_mean_std_binned_SWH(b,i_bin_SWH)=nanmean(amp_fit_std_whole(b,idx_int));
                amp_fit_std_std_binned_SWH(b,i_bin_SWH)=nanstd(amp_fit_std_whole(b,idx_int));
                
                max_wvfm_mean_std_binned_SWH(b,i_bin_SWH)=nanmean(max_wvfm_std_whole(b,idx_int));
                max_wvfm_std_std_binned_SWH(b,i_bin_SWH)=nanstd(max_wvfm_std_whole(b,idx_int));
                
            end
            
            COR_mean_mean_binned_SWH(b,i_bin_SWH)=nanmean(COR_mean_whole(b,idx_int));            
            
        end
        %------------------- starlab --------------------------------------
        if ~isempty(input_path_L2_STL) && b==1
            idx_int=STL_sorting_SWH == i_bin_SWH;
            if ~isempty(idx_int)
                STL_SSH_mean_std_binned_SWH(1,i_bin_SWH)=nanmean(STL_SSH_std_whole(b,idx_int));
                STL_SSH_std_std_binned_SWH(1,i_bin_SWH)=nanstd(STL_SSH_std_whole(b,idx_int));
                
                STL_SWH_mean_std_binned_SWH(1,i_bin_SWH)=nanmean(STL_SWH_std_whole(b,idx_int));
                STL_SWH_std_std_binned_SWH(1,i_bin_SWH)=nanstd(STL_SWH_std_whole(b,idx_int));
                
                STL_sigma0_mean_std_binned_SWH(1,i_bin_SWH)=nanmean(STL_sigma0_std_whole(b,idx_int));
                STL_sigma0_std_std_binned_SWH(1,i_bin_SWH)=nanstd(STL_sigma0_std_whole(b,idx_int));
            end
        end
        
    end
    
end


%% ---------------- Noise performance/mean values vs radial velocity ------
%--------------------------------------------------------------------------
%-------------------- Compute edges of binning ----------------------------
if ~cellfun(@isempty,input_path_L1_ISR_bs)
    
    idx_minus_hr=hr_whole<0.0;
    dumm=minmax(hr_whole(idx_minus_hr));
    min_minus_hr=dumm(1);
    max_minus_hr=dumm(2);
    edges_minus_hr=min_minus_hr:step_hr:max_minus_hr;
    n_total_bins_minus_hr=length(edges_minus_hr);
    dumm=minmax(hr_whole(~idx_minus_hr));
    min_plus_hr=dumm(1);
    max_plus_hr=dumm(2);
    edges_plus_hr=min_plus_hr:step_hr:max_plus_hr;
    n_total_bins_plus_hr=length(edges_plus_hr);
    
    
    % Mean of std using specific binning on radial velocity basis of the stds
    %positive radial velocities
    SSH_mean_std_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    SWH_mean_std_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    sigma0_mean_std_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    
    SSH_mean_mean_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    SWH_mean_mean_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    sigma0_mean_mean_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    COR_mean_mean_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    nb_mean_mean_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    
    SSH_std_std_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    SWH_std_std_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    sigma0_std_std_binned_plus_hr=NaN(N_baselines,n_total_bins_plus_hr);
    
    %negative radial velocities
    SSH_mean_std_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    SWH_mean_std_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    sigma0_mean_std_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    
    SSH_mean_mean_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    SWH_mean_mean_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    sigma0_mean_mean_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    COR_mean_mean_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    nb_mean_mean_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    
    SSH_std_std_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    SWH_std_std_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    sigma0_std_std_binned_minus_hr=NaN(N_baselines,n_total_bins_minus_hr);
    
    
    if ~isempty(input_path_L2_STL)
        % Mean of std using specific binning on radial velocity basis of the stds
        %positive radial velocities
        STL_SSH_mean_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_SWH_mean_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_sigma0_mean_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        STL_SSH_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_SWH_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_sigma0_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_COR_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        STL_SSH_std_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_SWH_std_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        STL_sigma0_std_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        %negative radial velocities
        STL_SSH_mean_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        STL_SWH_mean_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        STL_sigma0_mean_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        
        STL_SSH_mean_mean_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        STL_SWH_mean_mean_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        STL_sigma0_mean_mean_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        
        
        STL_SSH_std_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        STL_SWH_std_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        STL_sigma0_std_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
    end
    
    if ~isempty(input_path_L2_ESA)
        % Mean of std using specific binning on radial velocity basis of the stds
        %positive radial velocities
        ESA_SSH_mean_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        ESA_sigma0_mean_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        ESA_SSH_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        ESA_sigma0_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        ESA_SSH_std_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        ESA_sigma0_std_std_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        ESA_nb_mean_mean_binned_plus_hr=NaN(1,n_total_bins_plus_hr);
        
        %negative radial velocities
        ESA_SSH_mean_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        ESA_sigma0_mean_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        
        ESA_SSH_mean_mean_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        ESA_sigma0_mean_mean_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        
        ESA_SSH_std_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        ESA_sigma0_std_std_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
        
        ESA_nb_mean_mean_binned_minus_hr=NaN(1,n_total_bins_minus_hr);
    end
    
    
    %--------- Sorting mean SWH for the diferent boxes-------------------------
    %taking as reference the conventional SWH
    sorting_plus_hr=discretize(hr_mean_whole,edges_plus_hr);
    sorting_minus_hr=discretize(hr_mean_whole,edges_minus_hr);
    
    if ~isempty(input_path_L2_ESA)
        ESA_sorting_minus_hr=discretize(ESA_hr_mean_whole,edges_minus_hr);
        ESA_sorting_plus_hr=discretize(ESA_hr_mean_whole,edges_plus_hr);
    end
    for b=1:N_baselines
        % -------------- plus radial velocity -----------------------------
        for i_bin_hr=1:n_total_bins_plus_hr
            idx_int=sorting_plus_hr == i_bin_hr;
            if ~isempty(idx_int)
                SSH_mean_std_binned_plus_hr(b,i_bin_hr)=nanmean(SSH_std_whole(b,idx_int));
                SSH_std_std_binned_plus_hr(b,i_bin_hr)=nanstd(SSH_std_whole(b,idx_int));
                
                SWH_mean_std_binned_plus_hr(b,i_bin_hr)=nanmean(SWH_std_whole(b,idx_int));
                SWH_std_std_binned_plus_hr(b,i_bin_hr)=nanstd(SWH_std_whole(b,idx_int));
                
                sigma0_mean_std_binned_plus_hr(b,i_bin_hr)=nanmean(sigma0_std_whole(b,idx_int));
                sigma0_std_std_binned_plus_hr(b,i_bin_hr)=nanstd(sigma0_std_whole(b,idx_int));
                
                SSH_mean_mean_binned_plus_hr(b,i_bin_hr)=nanmean(SSH_mean_whole(b,idx_int));
                SWH_mean_mean_binned_plus_hr(b,i_bin_hr)=nanmean(SWH_mean_whole(b,idx_int));
                sigma0_mean_mean_binned_plus_hr(b,i_bin_hr)=nanmean(sigma0_mean_whole(b,idx_int));
                COR_mean_mean_binned_plus_hr(b,i_bin_hr)=nanmean(COR_mean_whole(b,idx_int));
                
                nb_mean_mean_binned_plus_hr(b,i_bin_hr)=nanmean(nb_mean_whole(b,idx_int));
                
                %----------------------- STL ------------------------------
                if ~isempty(input_path_L2_STL)
                    STL_SSH_mean_std_binned_plus_hr(1,i_bin_hr)=nanmean(STL_SSH_std_whole(idx_int));
                    STL_SSH_std_std_binned_plus_hr(1,i_bin_hr)=nanstd(STL_SSH_std_whole(idx_int));
                    
                    STL_SWH_mean_std_binned_plus_hr(1,i_bin_hr)=nanmean(STL_SWH_std_whole(idx_int));
                    STL_SWH_std_std_binned_plus_hr(1,i_bin_hr)=nanstd(STL_SWH_std_whole(idx_int));
                    
                    STL_sigma0_mean_std_binned_plus_hr(1,i_bin_hr)=nanmean(STL_sigma0_std_whole(idx_int));
                    STL_sigma0_std_std_binned_plus_hr(1,i_bin_hr)=nanstd(STL_sigma0_std_whole(idx_int));
                    
                    STL_SSH_mean_mean_binned_plus_hr(1,i_bin_hr)=nanmean(STL_SSH_mean_whole(idx_int));
                    STL_SWH_mean_mean_binned_plus_hr(1,i_bin_hr)=nanmean(STL_SWH_mean_whole(idx_int));
                    STL_sigma0_mean_mean_binned_plus_hr(1,i_bin_hr)=nanmean(STL_sigma0_mean_whole(idx_int));
                    
                end
                
                %------------------- ESA ---------------------------------
                if ~isempty(input_path_L2_ESA)
                    idx_int=ESA_sorting_plus_hr == i_bin_hr;
                    
                    if ~isempty(idx_int)
                        ESA_SSH_mean_std_binned_plus_hr(1,i_bin_hr)=nanmean(ESA_SSH_std_whole(idx_int));
                        ESA_SSH_std_std_binned_plus_hr(1,i_bin_hr)=nanstd(ESA_SSH_std_whole(idx_int));
                        
                        ESA_sigma0_mean_std_binned_plus_hr(1,i_bin_hr)=nanmean(ESA_sigma0_std_whole(idx_int));
                        ESA_sigma0_std_std_binned_plus_hr(1,i_bin_hr)=nanstd(ESA_sigma0_std_whole(idx_int));
                        
                        ESA_SSH_mean_mean_binned_plus_hr(1,i_bin_hr)=nanmean(ESA_SSH_mean_whole(idx_int));
                        ESA_sigma0_mean_mean_binned_plus_hr(1,i_bin_hr)=nanmean(ESA_sigma0_mean_whole(idx_int));
                        
                        ESA_nb_mean_mean_binned_plus_hr(1,i_bin_hr)=nanmean(ESA_nb_mean_whole(idx_int));
                        
                    end
                    
                end
                
            end
             clear idx_int;
        end
       
        % -------------- minus radial velocity -----------------------------
        for i_bin_hr=1:n_total_bins_minus_hr
            idx_int=sorting_minus_hr == i_bin_hr;
            if ~isempty(idx_int)
                SSH_mean_std_binned_minus_hr(b,i_bin_hr)=nanmean(SSH_std_whole(b,idx_int));
                SSH_std_std_binned_minus_hr(b,i_bin_hr)=nanstd(SSH_std_whole(b,idx_int));
                
                SWH_mean_std_binned_minus_hr(b,i_bin_hr)=nanmean(SWH_std_whole(b,idx_int));
                SWH_std_std_binned_minus_hr(b,i_bin_hr)=nanstd(SWH_std_whole(b,idx_int));
                
                sigma0_mean_std_binned_minus_hr(b,i_bin_hr)=nanmean(sigma0_std_whole(b,idx_int));
                sigma0_std_std_binned_minus_hr(b,i_bin_hr)=nanstd(sigma0_std_whole(b,idx_int));
                
                SSH_mean_mean_binned_minus_hr(b,i_bin_hr)=nanmean(SSH_mean_whole(b,idx_int));
                SWH_mean_mean_binned_minus_hr(b,i_bin_hr)=nanmean(SWH_mean_whole(b,idx_int));
                sigma0_mean_mean_binned_minus_hr(b,i_bin_hr)=nanmean(sigma0_mean_whole(b,idx_int));
                COR_mean_mean_binned_minus_hr(b,i_bin_hr)=nanmean(COR_mean_whole(b,idx_int));
                
                nb_mean_mean_binned_minus_hr(b,i_bin_hr)=nanmean(nb_mean_whole(b,idx_int));
                %---------------------- STL -------------------------------
                if ~isempty(input_path_L2_STL)
                    STL_SSH_mean_std_binned_minus_hr(1,i_bin_hr)=nanmean(STL_SSH_std_whole(idx_int));
                    STL_SSH_std_std_binned_minus_hr(1,i_bin_hr)=nanstd(STL_SSH_std_whole(idx_int));
                    
                    STL_SWH_mean_std_binned_minus_hr(1,i_bin_hr)=nanmean(STL_SWH_std_whole(idx_int));
                    STL_SWH_std_std_binned_minus_hr(1,i_bin_hr)=nanstd(STL_SWH_std_whole(idx_int));
                    
                    STL_sigma0_mean_std_binned_minus_hr(1,i_bin_hr)=nanmean(STL_sigma0_std_whole(idx_int));
                    STL_sigma0_std_std_binned_minus_hr(1,i_bin_hr)=nanstd(STL_sigma0_std_whole(idx_int));
                    
                    STL_SSH_mean_mean_binned_minus_hr(1,i_bin_hr)=nanmean(STL_SSH_mean_whole(idx_int));
                    STL_SWH_mean_mean_binned_minus_hr(1,i_bin_hr)=nanmean(STL_SWH_mean_whole(idx_int));
                    STL_sigma0_mean_mean_binned_minus_hr(1,i_bin_hr)=nanmean(STL_sigma0_mean_whole(idx_int));
                    
                end
                
                %------------------- ESA ---------------------------------
                if ~isempty(input_path_L2_ESA)
                    idx_int=ESA_sorting_minus_hr == i_bin_hr;
                    
                    if ~isempty(idx_int)
                        ESA_SSH_mean_std_binned_minus_hr(1,i_bin_hr)=nanmean(ESA_SSH_std_whole(idx_int));
                        ESA_SSH_std_std_binned_minus_hr(1,i_bin_hr)=nanstd(ESA_SSH_std_whole(idx_int));
                        
                        ESA_sigma0_mean_std_binned_minus_hr(1,i_bin_hr)=nanmean(ESA_sigma0_std_whole(idx_int));
                        ESA_sigma0_std_std_binned_minus_hr(1,i_bin_hr)=nanstd(ESA_sigma0_std_whole(idx_int));
                        
                        ESA_SSH_mean_mean_binned_minus_hr(1,i_bin_hr)=nanmean(ESA_SSH_mean_whole(idx_int));
                        ESA_sigma0_mean_mean_binned_minus_hr(1,i_bin_hr)=nanmean(ESA_sigma0_mean_whole(idx_int));
                        
                        ESA_nb_mean_mean_binned_minus_hr(1,i_bin_hr)=nanmean(ESA_nb_mean_whole(idx_int));
                    end
                    
                end
                
            end
             clear idx_int;
        end
        
    end
end
save(strcat(path_comparison_results,'Workspace_performance_comparison_baselines_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_')));

min_std_SSH=4;
max_std_SSH=14;

min_std_SWH=20;
max_std_SWH=100;

%% ----------------- Ploting noise performance vs SWH ---------------------
%-------------- SSH vs SWH ------------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_SWH,STL_SSH_mean_std_binned_SWH*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_SWH,SSH_mean_std_binned_SWH(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([min_SWH max_SWH min_std_SSH max_std_SSH]);
title(strcat('Mean $\sigma_{SSH}$ versus $\hat{H}_{s}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('$\hat{H}_{s}$ [m]','Interpreter','latex'); ylabel('$\hat{\sigma}_{SSH}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SSH_vs_SWH_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%-------------- SWH vs SWH ------------------------------------------------
figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_SWH,STL_SWH_mean_std_binned_SWH*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_SWH,SWH_mean_std_binned_SWH(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_SWH max_SWH min_std_SWH max_std_SWH]);
title(strcat('Mean $\sigma_{H_{s}}$ versus $\hat{H}_{s}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('$\hat{H}_{s}$ [m]','Interpreter','latex'); ylabel('$\hat{\sigma}_{H_{s}}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SWH_vs_SWH_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
%-------------- sigma0 vs SWH ------------------------------------------------
figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_SWH,STL_sigma0_mean_std_binned_SWH,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_SWH,sigma0_mean_std_binned_SWH(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_SWH max_SWH min_std_sigma0 max_std_sigma0]);
title(strcat('Mean $\sigma_{\sigma^{0}}$ versus $\hat{H}_{s}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('$\hat{H}_{s}$ [m]','Interpreter','latex'); ylabel('$\hat{\sigma}_{\sigma^{0}}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_sigma0_vs_SWH_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))


%% ----------------- Ploting noise performance vs hr ---------------------
%-------------- SSH vs v_radial ------------------------------------------------
min_std_SSH=4;
max_std_SSH=30;
%minus
legend_text={''};
f1=figure;
if ~isempty(input_path_L2_ESA)
    plot(edges_minus_hr,ESA_SSH_mean_std_binned_minus_hr*100.0,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(edges_minus_hr,STL_SSH_mean_std_binned_minus_hr*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_minus_hr,SSH_mean_std_binned_minus_hr(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_minus_hr max_minus_hr min_std_SSH max_std_SSH]);
title(strcat('Mean $\sigma_{SSH}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}_{SSH}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SSH_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_plus_hr,ESA_SSH_mean_std_binned_plus_hr*100.0,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(edges_plus_hr,STL_SSH_mean_std_binned_plus_hr*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,SSH_mean_std_binned_plus_hr(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_std_SSH max_std_SSH]);
title(strcat('Mean $\sigma_{SSH}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}_{SSH}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SSH_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%-------------- SWH vs vrad ------------------------------------------------
%minus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_minus_hr,STL_SWH_mean_std_binned_minus_hr*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_minus_hr,SWH_mean_std_binned_minus_hr(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_minus_hr max_minus_hr min_std_SWH max_std_SWH]);
title(strcat('Mean $\sigma_{H_{s}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}_{H_{s}}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SWH_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_plus_hr,STL_SWH_mean_std_binned_plus_hr*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,SWH_mean_std_binned_plus_hr(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_std_SWH max_std_SWH]);
title(strcat('Mean $\sigma_{H_{s}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}_{H_{s}}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SWH_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%-------------- sigma0 vs vrad ------------------------------------------------
%minus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_minus_hr,ESA_sigma0_mean_std_binned_minus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(edges_minus_hr,STL_sigma0_mean_std_binned_minus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_minus_hr,sigma0_mean_std_binned_minus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_minus_hr max_minus_hr min_std_sigma0 max_std_sigma0]);
title(strcat('Mean $\sigma_{\sigma^{0}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}_{\sigma^{0}}$ [dB]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_sigma0_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_plus_hr,ESA_sigma0_mean_std_binned_plus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(edges_plus_hr,STL_sigma0_mean_std_binned_plus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,sigma0_mean_std_binned_plus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_std_sigma0 max_std_sigma0]);
title(strcat('Mean $\sigma_{\sigma^{0}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}_{\sigma^{0}}$ [dB]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_sigma0_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);


%% ----------------- Ploting mean values  vs hr ---------------------------
%-------------- SSH vs v_radial ------------------------------------------------
min_SSH=Inf;
max_SSH=-Inf;
%minus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_minus_hr,ESA_SSH_mean_mean_binned_minus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
    min_SSH=min([min_SSH,min(ESA_SSH_mean_mean_binned_minus_hr),min(ESA_SSH_mean_mean_binned_plus_hr)]);
    max_SSH=max([max_SSH,max(ESA_SSH_mean_mean_binned_minus_hr),max(ESA_SSH_mean_mean_binned_plus_hr)]);
end
if ~isempty(input_path_L2_STL)
    plot(edges_minus_hr,STL_SSH_mean_mean_binned_minus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
    min_SSH=min([min_SSH,min(STL_SSH_mean_mean_binned_minus_hr),min(STL_SSH_mean_mean_binned_plus_hr)]);
    max_SSH=max([max_SSH,max(STL_SSH_mean_mean_binned_minus_hr),max(STL_SSH_mean_mean_binned_plus_hr)]);
end
for b=1:N_baselines
    plot(edges_minus_hr,SSH_mean_mean_binned_minus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on; 
    min_SSH=min([min_SSH,min(SSH_mean_mean_binned_minus_hr(b,:)),min(SSH_mean_mean_binned_plus_hr(b,:))]);
    max_SSH=max([max_SSH,max(SSH_mean_mean_binned_minus_hr(b,:)),max(SSH_mean_mean_binned_plus_hr(b,:))]);
end
axis([min_minus_hr max_minus_hr min_SSH max_SSH]);
title(strcat('Mean SSH versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{SSH}$ [m]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_SSH_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_plus_hr,ESA_SSH_mean_mean_binned_plus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(edges_plus_hr,STL_SSH_mean_mean_binned_plus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,SSH_mean_mean_binned_plus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_SSH max_SSH]);
title(strcat('Mean SSH versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{SSH}$ [m]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_SSH_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%-------------- SWH vs vrad ------------------------------------------------
%minus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_minus_hr,STL_SWH_mean_mean_binned_minus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_minus_hr,SWH_mean_mean_binned_minus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_minus_hr max_minus_hr min_SWH max_SWH]);
title(strcat('Mean ${H_{s}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{H}_{s}$ [m]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_SWH_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(edges_plus_hr,STL_SWH_mean_mean_binned_plus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,SWH_mean_mean_binned_plus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_SWH max_SWH]);
title(strcat('Mean ${H_{s}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{H}_{s}$ [m]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_SWH_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%-------------- sigma0 vs vrad ------------------------------------------------
min_sigma0=Inf;
max_sigma0=-Inf;
%minus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_minus_hr,ESA_sigma0_mean_mean_binned_minus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
    min_sigma0=min([min_sigma0,min(ESA_sigma0_mean_mean_binned_minus_hr),min(ESA_sigma0_mean_mean_binned_plus_hr)]);
    max_sigma0=max([max_sigma0,max(ESA_sigma0_mean_mean_binned_minus_hr),max(ESA_sigma0_mean_mean_binned_plus_hr)]);
end
if ~isempty(input_path_L2_STL)
    plot(edges_minus_hr,STL_sigma0_mean_mean_binned_minus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
    min_sigma0=min([min_sigma0,min(STL_sigma0_mean_mean_binned_minus_hr),min(STL_sigma0_mean_mean_binned_plus_hr)]);
    max_sigma0=max([max_sigma0,max(STL_sigma0_mean_mean_binned_minus_hr),max(STL_sigma0_mean_mean_binned_plus_hr)]);
end
for b=1:N_baselines
    plot(edges_minus_hr,sigma0_mean_mean_binned_minus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
    min_sigma0=min([min_sigma0,min(sigma0_mean_mean_binned_minus_hr(b,:)),min(sigma0_mean_mean_binned_plus_hr(b,:))]);
    max_sigma0=max([max_sigma0,max(sigma0_mean_mean_binned_minus_hr(b,:)),max(sigma0_mean_mean_binned_plus_hr(b,:))]);
end
axis([min_minus_hr max_minus_hr min_sigma0 max_sigma0]);
title(strcat('Mean $\sigma^{0}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}^{0}$ [dB]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_sigma0_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_plus_hr,ESA_sigma0_mean_mean_binned_plus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(edges_plus_hr,STL_sigma0_mean_mean_binned_plus_hr,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,sigma0_mean_mean_binned_plus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_sigma0 max_sigma0]);
title(strcat('Mean $\sigma^{0}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{\sigma}^{0}$ [dB]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_sigma0_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);

min_nb=150;
max_nb=250;
%-------------- nb beams vs vrad ------------------------------------------------
%minus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_minus_hr,ESA_nb_mean_mean_binned_minus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
for b=1:N_baselines
    plot(edges_minus_hr,nb_mean_mean_binned_minus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_minus_hr max_minus_hr min_nb max_nb]);
title(strcat('Mean Number of contributing beams ${n_{b}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{n}_{b}$','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_nb_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%plus
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(edges_plus_hr,ESA_nb_mean_mean_binned_plus_hr,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
for b=1:N_baselines
    plot(edges_plus_hr,nb_mean_mean_binned_plus_hr(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;  
end
axis([min_plus_hr max_plus_hr min_nb max_nb]);
title(strcat('Mean Number of contributing beams ${n_{b}}$ versus $v_{r}$'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthWest');
xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\hat{n}_{b}$','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'mean_nb_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);

end

