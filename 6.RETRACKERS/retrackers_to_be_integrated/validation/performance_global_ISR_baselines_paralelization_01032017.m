function performance_global_ISR_baselines_paralelization_01032017(input_path_L2_ISR_bs,name_bs,path_comparison_results,varargin)
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
version_matlab=version;
ttotal=tic;
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<3 || nargin>(3+25*2))
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
p.addParamValue('input_path_L2_ISD_L1B_ESA','',@(x)ischar(x));
p.addParamValue('input_path_L2_STL','',@(x)ischar(x));
p.addParamValue('input_path_L1_ESA','',@(x)ischar(x));
p.addParamValue('annotation_box_active',1);
p.addParamValue('filename_mask_KML','',@(x)ischar(x));
p.addParamValue('num_pools',1);
p.addParamValue('plot_downsampling',50);
p.addParamValue('generate_plots',0);
p.addParamValue('generate_kml',0);
p.addParamValue('filter_land_surf_type',0);
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
input_path_L2_ISD_L1B_ESA=p.Results.input_path_L2_ISD_L1B_ESA;
input_path_L1_ESA=p.Results.input_path_L1_ESA;
input_path_L2_STL=p.Results.input_path_L2_STL;
annotation_box_active=p.Results.annotation_box_active;
filename_mask_KML=p.Results.filename_mask_KML;
num_pools=p.Results.num_pools;
plot_downsampling=p.Results.plot_downsampling;
generate_plots=p.Results.generate_plots;
generate_kml=p.Results.generate_kml;
filter_land_surf_type=p.Results.filter_land_surf_type;


%% ----------------- Hard coded variables definition ----------------------
y_add_textbox=0.01;
nanclr=[1 1 1];

min_std_SSH=4; %cm
max_std_SSH=50;%cm

min_std_SWH=10; %cm
max_std_SWH=50;%cm

min_std_sigma0=0; %cm
max_std_sigma0=0.7;%cm

min_mean_error_SSH=-100; %cm
max_mean_error_SSH=100; %cm

min_mean_error_SWH=-150; %cm
max_mean_error_SWH=150; %cm

min_mean_error_sigma0=-2; %dB
max_mean_error_sigma0=10; %dB

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
textbox_fontsize=8;
legend_fontsize=14;


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
    %disp(i_file)
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
SSH_std_whole={};
SWH_std_whole={};
sigma0_std_whole={};
amp_fit_std_whole={};
max_wvfm_std_whole={};
SSH_std_mean_whole={};
SWH_std_mean_whole={};
sigma0_std_mean_whole={};
COR_std_mean_whole={};
SSH_mean_whole={};
SWH_mean_whole={};
sigma0_mean_whole={};
COR_mean_whole={};
LAT_median_whole={};
LON_median_whole={};
SSH_whole={};
SWH_whole={};
sigma0_whole={};
COR_whole={};


if ~cellfun(@isempty,input_path_L1_ISR_bs)
    amp_fit_mean_whole={};
    max_wvfm_mean_whole={};
    hr_mean_whole={};
    nb_mean_whole={};
    amp_fit_whole={};
    max_wvfm_whole={};
    hr_whole={};
    nb_whole={};
else
    amp_fit_mean_whole={[]};
    max_wvfm_mean_whole={[]};
    hr_mean_whole={[]};
    nb_mean_whole={[]};
    amp_fit_whole={[]};
    max_wvfm_whole={[]};
    hr_whole={[]};
    nb_whole={[]};
end



if ~isempty(input_path_L2_ESA)
    ESA_SSH_std_whole={};
    ESA_sigma0_std_whole={};
    ESA_SSH_std_mean_whole={};
    ESA_sigma0_std_mean_whole={};
    ESA_SSH_mean_whole={};
    ESA_sigma0_mean_whole={};
    ESA_SSH_whole={};
    ESA_sigma0_whole={};
    SSH_mean_error_ESA_ISR_whole={};
    sigma0_mean_error_ESA_ISR_whole={};
    
    if ~isempty(input_path_L1_ESA)
        ESA_hr_mean_whole={}; %height rate (radial veloctiy)
        ESA_nb_mean_whole={};
        ESA_hr_whole={};
        ESA_nb_whole={};
    else
        ESA_hr_mean_whole={[]}; %height rate (radial veloctiy)
        ESA_nb_mean_whole={[]};
        ESA_hr_whole={[]};
        ESA_nb_whole={[]};
    end
    
else
    ESA_SSH_std_whole={[]};
    ESA_sigma0_std_whole={[]};
    ESA_SSH_std_mean_whole={[]};
    ESA_sigma0_std_mean_whole={[]};
    ESA_SSH_mean_whole={[]};
    ESA_sigma0_mean_whole={[]};
    ESA_hr_mean_whole={[]}; %height rate (radial veloctiy)
    ESA_nb_mean_whole={[]};
end

if ~isempty(input_path_L2_STL)
    STL_SSH_std_whole={};
    STL_sigma0_std_whole={};
    STL_SWH_std_whole={};
    STL_SSH_std_mean_whole={};
    STL_sigma0_std_mean_whole={};
    STL_SWH_std_mean_whole={};
    STL_SSH_mean_whole={};
    STL_sigma0_mean_whole={};
    STL_SWH_mean_whole={};
    STL_SSH_whole={};
    STL_sigma0_whole={};
    STL_SWH_whole={};
    SSH_mean_error_STL_ISR_whole={};
    SWH_mean_error_STL_ISR_whole={};
    sigma0_mean_error_STL_ISR_whole={};
else
    STL_SSH_std_whole={[]};
    STL_sigma0_std_whole={[]};
    STL_SWH_std_whole={[]};
    STL_SSH_std_mean_whole={[]};
    STL_sigma0_std_mean_whole={[]};
    STL_SWH_std_mean_whole={[]};
    STL_SSH_mean_whole={[]};
    STL_sigma0_mean_whole={[]};
    STL_SWH_mean_whole={[]};
    STL_SSH_whole={[]};
    STL_sigma0_whole={[]};
    STL_SWH_whole={[]};
    SSH_mean_error_STL_ISR_whole={[]};
    SWH_mean_error_STL_ISR_whole={[]};
    sigma0_mean_error_STL_ISR_whole={[]};
end

if smooth_param
    SSH_smooth_whole={}; %use the same window as the one for std and detrending
    SWH_smooth_whole={};
    sigma0_smooth_whole={};
    COR_smooth_whole={};
    SSH_RMSE_whole={};
    SWH_RMSE_whole={};
    sigma0_RMSE_whole={};
    COR_RMSE_whole={};   
    if ~isempty(input_path_L2_ESA)
        ESA_SSH_smooth_whole={};
        ESA_sigma0_smooth_whole={};
        ESA_SSH_RMSE_whole={};
        ESA_sigma0_RMSE_whole={};
    else
        ESA_SSH_smooth_whole={[]};
        ESA_sigma0_smooth_whole={[]};
        ESA_SSH_RMSE_whole={[]};
        ESA_sigma0_RMSE_whole={[]};        
    end
    if ~isempty(input_path_L2_STL)
        STL_SSH_smooth_whole={};
        STL_sigma0_smooth_whole={};
        STL_SWH_smooth_whole={};
        STL_SWH_RMSE_whole={};
        STL_SSH_RMSE_whole={};
        STL_sigma0_RMSE_whole={};
    else
        STL_SSH_smooth_whole={[]};
        STL_sigma0_smooth_whole={[]};
        STL_SWH_smooth_whole={[]};
        STL_SWH_RMSE_whole={[]};
        STL_SSH_RMSE_whole={[]};
        STL_sigma0_RMSE_whole={[]};
    end
else
    SSH_smooth_whole={[]}; %use the same window as the one for std and detrending
    SWH_smooth_whole={[]};
    sigma0_smooth_whole={[]};
    COR_smooth_whole={[]};
    SSH_RMSE_whole={[]};
    SWH_RMSE_whole={[]};
    sigma0_RMSE_whole={[]};
    COR_RMSE_whole={[]};
    
    ESA_SSH_smooth_whole={[]};
    ESA_sigma0_smooth_whole={[]};
    ESA_SSH_RMSE_whole={[]};
    ESA_sigma0_RMSE_whole={[]};
    STL_SSH_smooth_whole={[]};
    STL_sigma0_smooth_whole={[]};
    STL_SWH_smooth_whole={[]};
    STL_SWH_RMSE_whole={[]};
    STL_SSH_RMSE_whole=[];
    STL_sigma0_RMSE_whole=[];
    
end


if ~define_min_max_SWH
    max_SWH=-Inf;
    min_SWH=Inf;
end


num_surfaces_whole={};
if num_pools~=1
    %create pools
    if str2double(version_matlab(end-5:end-2))>2013
        parpool(num_pools);
    else
        matlabpool('open',num_pools);
    end
    parfor i_files_valid=1:filesBulk.nFilesNC_valid
        
         try
            
            [SSH_std_whole{i_files_valid},SWH_std_whole{i_files_valid},sigma0_std_whole{i_files_valid},amp_fit_std_whole{i_files_valid},...
                max_wvfm_std_whole{i_files_valid},ESA_SSH_std_whole{i_files_valid},ESA_sigma0_std_whole{i_files_valid},...
                STL_SSH_std_whole{i_files_valid},STL_sigma0_std_whole{i_files_valid},STL_SWH_std_whole{i_files_valid},...
                SSH_std_mean_whole{i_files_valid},SWH_std_mean_whole{i_files_valid},sigma0_std_mean_whole{i_files_valid},COR_std_mean_whole{i_files_valid},...
                ESA_SSH_std_mean_whole{i_files_valid},ESA_sigma0_std_mean_whole{i_files_valid},...
                STL_SSH_std_mean_whole{i_files_valid},STL_sigma0_std_mean_whole{i_files_valid},STL_SWH_std_mean_whole{i_files_valid},...
                SSH_mean_whole{i_files_valid},SWH_mean_whole{i_files_valid},sigma0_mean_whole{i_files_valid},...
                amp_fit_mean_whole{i_files_valid},max_wvfm_mean_whole{i_files_valid},hr_mean_whole{i_files_valid},nb_mean_whole{i_files_valid},...
                COR_mean_whole{i_files_valid},LAT_median_whole{i_files_valid},LON_median_whole{i_files_valid},...
                ESA_SSH_mean_whole{i_files_valid},ESA_sigma0_mean_whole{i_files_valid},...
                STL_SSH_mean_whole{i_files_valid},STL_sigma0_mean_whole{i_files_valid},STL_SWH_mean_whole{i_files_valid},...
                SSH_whole{i_files_valid},SWH_whole{i_files_valid},sigma0_whole{i_files_valid},COR_whole{i_files_valid},...
                ESA_SSH_whole{i_files_valid},ESA_sigma0_whole{i_files_valid},...
                STL_SSH_whole{i_files_valid},STL_sigma0_whole{i_files_valid},STL_SWH_whole{i_files_valid},...
                amp_fit_whole{i_files_valid},max_wvfm_whole{i_files_valid},hr_whole{i_files_valid},nb_whole{i_files_valid},...
                ESA_hr_mean_whole{i_files_valid},ESA_hr_whole{i_files_valid},ESA_nb_mean_whole{i_files_valid},ESA_nb_whole{i_files_valid},...
                SSH_smooth_whole{i_files_valid},SWH_smooth_whole{i_files_valid},sigma0_smooth_whole{i_files_valid},COR_smooth_whole{i_files_valid},...
                SSH_RMSE_whole{i_files_valid},SWH_RMSE_whole{i_files_valid},sigma0_RMSE_whole{i_files_valid},COR_RMSE_whole{i_files_valid},...
                ESA_SSH_smooth_whole{i_files_valid},ESA_sigma0_smooth_whole{i_files_valid},ESA_SSH_RMSE_whole{i_files_valid},ESA_sigma0_RMSE_whole{i_files_valid},...
                STL_SSH_smooth_whole{i_files_valid},STL_SWH_smooth_whole{i_files_valid},STL_sigma0_smooth_whole{i_files_valid},...
                STL_SSH_RMSE_whole{i_files_valid},STL_sigma0_RMSE_whole{i_files_valid},STL_SWH_RMSE_whole{i_files_valid},...
                SSH_mean_error_ESA_ISR_whole{i_files_valid},sigma0_mean_error_ESA_ISR_whole{i_files_valid},...
                SSH_mean_error_STL_ISR_whole{i_files_valid},SWH_mean_error_STL_ISR_whole{i_files_valid},sigma0_mean_error_STL_ISR_whole{i_files_valid},...
                num_surfaces_whole{i_files_valid}]=run_performance_ISR_baselines(N_baselines,filesBulk,input_path_L2_ISR_bs,input_path_L2_ESA,input_path_L1_ESA,input_path_L2_STL,input_path_L1_ISR_bs,...
                                                    name_bs,filename_mask_KML,...
                                                    flag_outliers_removal,type_outliers_removal,outlier_percentil_low,outlier_percentil_high,IQR_times,hampel_wind,hampel_sigma,...
                                                    smooth_param,win_size_detrending,...
                                                    i_files_valid,plot_downsampling,marker_ESA,color_ESA,marker_STL,color_STL,marker_bs,color_bs,sh_name_nc,...
                                                    path_comparison_results,generate_plots,print_file,res_fig,file_ext,legend_fontsize,textbox_fontsize,annotation_box_active,...
                                                    geo_mask,generate_kml,filter_land_surf_type);
            
            
        catch
            continue;
        end

        
    end %loop over tracks
    
    %close pools
    if str2double(version_matlab(end-5:end-2))>2013
        poolobj = gcp('nocreate');
        delete(poolobj);
    else
        matlabpool('close');
    end
else
    for i_files_valid=1:filesBulk.nFilesNC_valid        
%         try
            [SSH_std_whole{i_files_valid},SWH_std_whole{i_files_valid},sigma0_std_whole{i_files_valid},amp_fit_std_whole{i_files_valid},...
                max_wvfm_std_whole{i_files_valid},ESA_SSH_std_whole{i_files_valid},ESA_sigma0_std_whole{i_files_valid},...
                STL_SSH_std_whole{i_files_valid},STL_sigma0_std_whole{i_files_valid},STL_SWH_std_whole{i_files_valid},...
                SSH_std_mean_whole{i_files_valid},SWH_std_mean_whole{i_files_valid},sigma0_std_mean_whole{i_files_valid},COR_std_mean_whole{i_files_valid},...
                ESA_SSH_std_mean_whole{i_files_valid},ESA_sigma0_std_mean_whole{i_files_valid},...
                STL_SSH_std_mean_whole{i_files_valid},STL_sigma0_std_mean_whole{i_files_valid},STL_SWH_std_mean_whole{i_files_valid},...
                SSH_mean_whole{i_files_valid},SWH_mean_whole{i_files_valid},sigma0_mean_whole{i_files_valid},...
                amp_fit_mean_whole{i_files_valid},max_wvfm_mean_whole{i_files_valid},hr_mean_whole{i_files_valid},nb_mean_whole{i_files_valid},...
                COR_mean_whole{i_files_valid},LAT_median_whole{i_files_valid},LON_median_whole{i_files_valid},...
                ESA_SSH_mean_whole{i_files_valid},ESA_sigma0_mean_whole{i_files_valid},...
                STL_SSH_mean_whole{i_files_valid},STL_sigma0_mean_whole{i_files_valid},STL_SWH_mean_whole{i_files_valid},...
                SSH_whole{i_files_valid},SWH_whole{i_files_valid},sigma0_whole{i_files_valid},COR_whole{i_files_valid},...
                ESA_SSH_whole{i_files_valid},ESA_sigma0_whole{i_files_valid},...
                STL_SSH_whole{i_files_valid},STL_sigma0_whole{i_files_valid},STL_SWH_whole{i_files_valid},...
                amp_fit_whole{i_files_valid},max_wvfm_whole{i_files_valid},hr_whole{i_files_valid},nb_whole{i_files_valid},...
                ESA_hr_mean_whole{i_files_valid},ESA_hr_whole{i_files_valid},ESA_nb_mean_whole{i_files_valid},ESA_nb_whole{i_files_valid},...
                SSH_smooth_whole{i_files_valid},SWH_smooth_whole{i_files_valid},sigma0_smooth_whole{i_files_valid},COR_smooth_whole{i_files_valid},...
                SSH_RMSE_whole{i_files_valid},SWH_RMSE_whole{i_files_valid},sigma0_RMSE_whole{i_files_valid},COR_RMSE_whole{i_files_valid},...
                ESA_SSH_smooth_whole{i_files_valid},ESA_sigma0_smooth_whole{i_files_valid},ESA_SSH_RMSE_whole{i_files_valid},ESA_sigma0_RMSE_whole{i_files_valid},...
                STL_SSH_smooth_whole{i_files_valid},STL_SWH_smooth_whole{i_files_valid},STL_sigma0_smooth_whole{i_files_valid},...
                STL_SSH_RMSE_whole{i_files_valid},STL_sigma0_RMSE_whole{i_files_valid},STL_SWH_RMSE_whole{i_files_valid},...
                SSH_mean_error_ESA_ISR_whole{i_files_valid},sigma0_mean_error_ESA_ISR_whole{i_files_valid},...
                SSH_mean_error_STL_ISR_whole{i_files_valid},SWH_mean_error_STL_ISR_whole{i_files_valid},sigma0_mean_error_STL_ISR_whole{i_files_valid},...
                num_surfaces_whole{i_files_valid}]=run_performance_ISR_baselines(N_baselines,filesBulk,input_path_L2_ISR_bs,input_path_L2_ESA,input_path_L1_ESA,input_path_L2_STL,input_path_L1_ISR_bs,...
                name_bs,filename_mask_KML,...
                flag_outliers_removal,type_outliers_removal,outlier_percentil_low,outlier_percentil_high,IQR_times,hampel_wind,hampel_sigma,...
                smooth_param,win_size_detrending,...
                i_files_valid,plot_downsampling,marker_ESA,color_ESA,marker_STL,color_STL,marker_bs,color_bs,sh_name_nc,...
                path_comparison_results,generate_plots,print_file,res_fig,file_ext,legend_fontsize,textbox_fontsize,annotation_box_active,...
                geo_mask,generate_kml,filter_land_surf_type);
            
            
%         catch
%             continue;
%         end
          
    end %loop over tracks
    
end
%------------- re-organize the info from cell arrays to arrays ------------
    SSH_std_whole=[SSH_std_whole{:}];
    SWH_std_whole=[SWH_std_whole{:}];
    sigma0_std_whole=[sigma0_std_whole{:}];
    amp_fit_std_whole=[amp_fit_std_whole{:}];
    max_wvfm_std_whole=[max_wvfm_std_whole{:}];
    ESA_SSH_std_whole=[ESA_SSH_std_whole{:}];
    ESA_sigma0_std_whole=[ESA_sigma0_std_whole{:}];
    STL_SSH_std_whole=[STL_SSH_std_whole{:}];
    STL_sigma0_std_whole=[STL_sigma0_std_whole{:}];
    STL_SWH_std_whole=[STL_SWH_std_whole{:}];
    SSH_std_mean_whole=[SSH_std_mean_whole{:}];
    SWH_std_mean_whole=[SWH_std_mean_whole{:}];
    sigma0_std_mean_whole=[sigma0_std_mean_whole{:}];
    COR_std_mean_whole=[COR_std_mean_whole{:}];
    ESA_SSH_std_mean_whole=[ESA_SSH_std_mean_whole{:}];
    ESA_sigma0_std_mean_whole=[ESA_sigma0_std_mean_whole{:}];
    STL_SSH_std_mean_whole=[STL_SSH_std_mean_whole{:}];
    STL_sigma0_std_mean_whole=[STL_sigma0_std_mean_whole{:}];
    STL_SWH_std_mean_whole=[STL_SWH_std_mean_whole{:}];
    SSH_mean_whole=[SSH_mean_whole{:}];
    SWH_mean_whole=[SWH_mean_whole{:}];
    sigma0_mean_whole=[sigma0_mean_whole{:}];
    amp_fit_mean_whole=[amp_fit_mean_whole{:}];
    max_wvfm_mean_whole=[max_wvfm_mean_whole{:}];
    hr_mean_whole=[hr_mean_whole{:}];
    nb_mean_whole=[nb_mean_whole{:}];
    COR_mean_whole=[COR_mean_whole{:}];
    LAT_median_whole=[LAT_median_whole{:}];
    LON_median_whole=[LON_median_whole{:}];
    ESA_SSH_mean_whole=[ESA_SSH_mean_whole{:}];
    ESA_sigma0_mean_whole=[ESA_sigma0_mean_whole{:}];
    STL_SSH_mean_whole=[STL_SSH_mean_whole{:}];
    STL_sigma0_mean_whole=[STL_sigma0_mean_whole{:}];
    STL_SWH_mean_whole=[STL_SWH_mean_whole{:}];
    SSH_whole=[SSH_whole{:}];
    SWH_whole=[SWH_whole{:}];
    sigma0_whole=[sigma0_whole{:}];
    COR_whole=[COR_whole{:}];
    ESA_SSH_whole=[ESA_SSH_whole{:}];
    ESA_sigma0_whole=[ESA_sigma0_whole{:}];
    STL_SSH_whole=[STL_SSH_whole{:}];
    STL_sigma0_whole=[STL_sigma0_whole{:}];
    STL_SWH_whole=[STL_SWH_whole{:}];
    amp_fit_whole=[amp_fit_whole{:}];
    max_wvfm_whole=[max_wvfm_whole{:}];
    hr_whole=[hr_whole{:}];
    nb_whole=[nb_whole{:}];
    ESA_hr_mean_whole=[ESA_hr_mean_whole{:}];
    ESA_hr_whole=[ESA_hr_whole{:}];
    ESA_nb_mean_whole=[ESA_nb_mean_whole{:}];
    ESA_nb_whole=[ESA_nb_whole{:}];
    SSH_smooth_whole=[SSH_smooth_whole{:}];
    SWH_smooth_whole=[SWH_smooth_whole{:}];
    sigma0_smooth_whole=[sigma0_smooth_whole{:}];
    COR_smooth_whole=[COR_smooth_whole{:}];
    SSH_RMSE_whole=[SSH_RMSE_whole{:}];
    SWH_RMSE_whole=[SWH_RMSE_whole{:}];
    sigma0_RMSE_whole=[sigma0_RMSE_whole{:}];
    COR_RMSE_whole=[COR_RMSE_whole{:}];
    ESA_SSH_smooth_whole=[ESA_SSH_smooth_whole{:}];
    ESA_sigma0_smooth_whole=[ESA_sigma0_smooth_whole{:}];
    ESA_SSH_RMSE_whole=[ESA_SSH_RMSE_whole{:}];
    ESA_sigma0_RMSE_whole=[ESA_sigma0_RMSE_whole{:}];
    STL_SSH_smooth_whole=[STL_SSH_smooth_whole{:}];
    STL_SWH_smooth_whole=[STL_SWH_smooth_whole{:}];
    STL_sigma0_smooth_whole=[STL_sigma0_smooth_whole{:}];
    STL_SSH_RMSE_whole=[STL_SSH_RMSE_whole{:}];
    STL_sigma0_RMSE_whole=[STL_sigma0_RMSE_whole{:}];
    STL_SWH_RMSE_whole=[STL_SWH_RMSE_whole{:}];
    SSH_mean_error_ESA_ISR_whole=[SSH_mean_error_ESA_ISR_whole{:}];
    sigma0_mean_error_ESA_ISR_whole=[sigma0_mean_error_ESA_ISR_whole{:}];
    SSH_mean_error_STL_ISR_whole=[SSH_mean_error_STL_ISR_whole{:}];
    SWH_mean_error_STL_ISR_whole=[SWH_mean_error_STL_ISR_whole{:}];
    sigma0_mean_error_STL_ISR_whole=[sigma0_mean_error_STL_ISR_whole{:}];
    num_surfaces_whole=[num_surfaces_whole{:}];

save(strcat(path_comparison_results,'Workspace_performance_comparison_baselines_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_')));

%load('C:\Users\eduard.makhoul\isardSAT\projects\SCOOP\processing\results\Phase_2\test\2013\performance_comparison\agulhas\comparison_S3_old_S3_new_bis\Workspace_performance_comparison_baselines_S_3_old_vs_S_3_new.mat');



% min_SWH=1.5;
% max_SWH=5.5;


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
%sorting_SWH=discretize(SWH_mean_whole(2,:),edges_SWH);
for b=1:N_baselines
    sorting_SWH=discretize(SWH_mean_whole(b,:),edges_SWH);
    if ~isempty(input_path_L2_STL) && b==1
        STL_sorting_SWH=discretize(STL_SWH_mean_whole,edges_SWH);
        %STL_sorting_SWH=sorting_SWH;
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

%save(strcat(path_comparison_results,'Workspace_performance_comparison_baselines_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_')));







init_track=1;
last_track=length(ESA_SSH_std_mean_whole);

%% ---------------------  Values along-tracks -----------------------------
%--------------------------------------------------------------------------
%--------------------------NOISE PERFORMANCE ------------------------------
%--------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&& SSH &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%---------------------- std -----------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(ESA_SSH_std_mean_whole*100.0,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(STL_SSH_std_mean_whole*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(SSH_std_mean_whole(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([init_track last_track min_std_SSH max_std_SSH]);
title(strcat('Mean $\sigma_{SSH}$ over the different tracks'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('Track','Interpreter','latex'); ylabel('$\hat{\sigma}_{SSH}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SSH_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%---------------------- RMSE ----------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(ESA_SSH_RMSE_whole*100.0,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(STL_SSH_RMSE_whole*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(SSH_RMSE_whole(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([init_track last_track min_std_SSH max_std_SSH]);
title(strcat('RMSE SSH over the different tracks'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('Track','Interpreter','latex'); ylabel('$RMSE_{SSH}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'RMSE_SSH_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);

% &&&&&&&&&&&&&&&&&&&&&&&&&&& SWH &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%---------------------- std -----------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(STL_SWH_std_mean_whole*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(SWH_std_mean_whole(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([init_track last_track min_std_SWH max_std_SWH]);
title(strcat('Mean $\sigma_{H_{s}}$ over the different tracks'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('Track','Interpreter','latex'); ylabel('$\hat{\sigma}_{H_{s}}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_SWH_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%---------------------- RMSE ----------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_STL)
    plot(STL_SWH_RMSE_whole*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(SWH_RMSE_whole(b,:)*100.0,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([init_track last_track min_std_SWH max_std_SWH]);
title(strcat('RMSE $H_{s}$ over the different tracks'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('Track','Interpreter','latex'); ylabel('$RMSE_{H_{s}}$ [cm]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'RMSE_SWH_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);

% &&&&&&&&&&&&&&&&&&&&&&&&&&& Sigma0 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%---------------------- std -----------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(ESA_sigma0_std_mean_whole,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(STL_sigma0_std_mean_whole,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(sigma0_std_mean_whole(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([init_track last_track min_std_sigma0 max_std_sigma0]);
title(strcat('Mean $\sigma_{\sigma^0}$ over the different tracks'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('Track','Interpreter','latex'); ylabel('$\hat{\sigma}_{\sigma^0}$ [dB]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'std_sigma0_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);
%---------------------- RMSE ----------------------------------------------
f1=figure;
legend_text={''};
if ~isempty(input_path_L2_ESA)
    plot(ESA_sigma0_RMSE_whole,'Marker',char(marker_ESA),'Color',color_ESA);
    legend_text=[legend_text,'ESA'];
    hold on;
end
if ~isempty(input_path_L2_STL)
    plot(STL_SSH_RMSE_whole*100.0,'Marker',char(marker_STL),'Color',color_STL);
    legend_text=[legend_text,'STL'];
    hold on;
end
for b=1:N_baselines
    plot(sigma0_RMSE_whole(b,:),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
    hold on;
end
axis([init_track last_track min_std_sigma0 max_std_sigma0]);
title(strcat('RMSE $\sigma^0$ over the different tracks'),'Interpreter','latex');
legend_text=[legend_text,name_bs];
legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
xlabel('Track','Interpreter','latex'); ylabel('$RMSE_{\sigma^0}$ [dB]','Interpreter','latex');
print(print_file,res_fig,strcat(path_comparison_results,'RMSE_sigma0_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
close(f1);

%--------------------------------------------------------------------------
%--------------------------Mean errors ------------------------------------
%--------------------------------------------------------------------------
% &&&&&&&&&&&&&&&&&&&&&&&&&&& SSH &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if ~isempty(input_path_L2_ESA) || ~isempty(input_path_L2_STL)
    f1=figure;
    legend_text={''};
    for b=1:N_baselines
        if ~isempty(input_path_L2_ESA)
            plot(SSH_mean_error_ESA_ISR_whole(b,:)*100.0,'LineStyle','none','Marker',char(marker_bs(b)),'Color',color_ESA);
            legend_text=[legend_text,strcat('ESA-ISR',{' '},name_bs(b))];
            hold on;
        end
        if ~isempty(input_path_L2_STL)
            plot(SSH_mean_error_STL_ISR_whole(b,:)*100.0,'LineStyle','none','Marker',char(marker_bs(b)),'Color',color_STL);
            legend_text=[legend_text,strcat('STL-ISR',{' '},name_bs(b))];
            hold on;
        end        
    end
    axis([init_track last_track min_mean_error_SSH max_mean_error_SSH]);
    title(strcat('Mean error $\epsilon_{SSH}$ over the different tracks'),'Interpreter','latex');
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','SouthEast');
    xlabel('Track','Interpreter','latex'); ylabel('$\hat{\epsilon}_{SSH}$ [cm]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'mean_error_SSH_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);
end
% &&&&&&&&&&&&&&&&&&&&&&&&&&&& SWH &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if ~isempty(input_path_L2_STL)
    f1=figure;
    legend_text={''};
    for b=1:N_baselines
        plot(SWH_mean_error_STL_ISR_whole(b,:)*100.0,'LineStyle','none','Marker',char(marker_bs(b)),'Color',color_STL);
        legend_text=[legend_text,strcat('STL-ISR',{' '},name_bs(b))];
        hold on;
    end
    axis([init_track last_track min_mean_error_SWH max_mean_error_SWH]);
    title(strcat('Mean error $\epsilon_{H_{s}}$ over the different tracks'),'Interpreter','latex');
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
    xlabel('Track','Interpreter','latex'); ylabel('$\hat{\epsilon}_{H_{s}}$ [cm]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'mean_error_SWH_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);
end
% &&&&&&&&&&&&&&&&&&&&&&&&&&& sigma0 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
if ~isempty(input_path_L2_ESA) || ~isempty(input_path_L2_STL)
    f1=figure;
    legend_text={''};
    for b=1:N_baselines
        if ~isempty(input_path_L2_ESA)
            plot(sigma0_mean_error_ESA_ISR_whole(b,:),'LineStyle','none','Marker',char(marker_bs(b)),'Color',color_ESA);
            legend_text=[legend_text,strcat('ESA-ISR',{' '},name_bs(b))];
            hold on;
        end
        if ~isempty(input_path_L2_STL)
            plot(sigma0_mean_error_STL_ISR_whole(b,:),'LineStyle','none','Marker',char(marker_bs(b)),'Color',color_STL);
            legend_text=[legend_text,strcat('STL-ISR',{' '},name_bs(b))];
            hold on;
        end        
    end
    axis([init_track last_track min_mean_error_sigma0 max_mean_error_sigma0]);
    title(strcat('Mean error $\epsilon_{\sigma^0}$ over the different tracks'),'Interpreter','latex');
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','NorthEast');
    xlabel('Track','Interpreter','latex'); ylabel('$\hat{\epsilon}_{\sigma^0}$ [dB]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'mean_error_sigma0_vs_track_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);
end


min_std_SSH=4;
max_std_SSH=30;
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

if ~isempty(input_path_L2_ESA)
    %% ----------------- Ploting difference mean values ESA vs hr ----------------
    %-------------- SSH vs v_radial ------------------------------------------------
    min_error_SSH=Inf;
    max_error_SSH=-Inf;
    %minus
    f1=figure;
    legend_text={''};
    if ~isempty(input_path_L2_STL)
        plot(edges_minus_hr,(ESA_SSH_mean_mean_binned_minus_hr-STL_SSH_mean_mean_binned_minus_hr)*100,'Marker',char(marker_STL),'Color',color_STL);
        legend_text=[legend_text,'ESA - STL'];
        hold on;
        min_error_SSH=min([min_error_SSH,min(ESA_SSH_mean_mean_binned_minus_hr-STL_SSH_mean_mean_binned_minus_hr).*100,min(ESA_SSH_mean_mean_binned_plus_hr-STL_SSH_mean_mean_binned_plus_hr).*100]);
        max_error_SSH=max([max_error_SSH,max(ESA_SSH_mean_mean_binned_minus_hr-STL_SSH_mean_mean_binned_minus_hr).*100,max(ESA_SSH_mean_mean_binned_plus_hr-STL_SSH_mean_mean_binned_plus_hr).*100]);
    end
    for b=1:N_baselines
        plot(edges_minus_hr,(ESA_SSH_mean_mean_binned_minus_hr-SSH_mean_mean_binned_minus_hr(b,:)).*100,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
        hold on;
        min_error_SSH=min([min_error_SSH,min(ESA_SSH_mean_mean_binned_minus_hr-SSH_mean_mean_binned_minus_hr(b,:)).*100,min(ESA_SSH_mean_mean_binned_plus_hr-SSH_mean_mean_binned_plus_hr(b,:)).*100]);
        max_error_SSH=max([max_error_SSH,max(ESA_SSH_mean_mean_binned_minus_hr-SSH_mean_mean_binned_minus_hr(b,:)).*100,max(ESA_SSH_mean_mean_binned_plus_hr-SSH_mean_mean_binned_plus_hr(b,:)).*100]);
    end
    axis([min_minus_hr max_minus_hr min_error_SSH max_error_SSH]);
    title(strcat('Difference mean SSH versus $v_{r}$'),'Interpreter','latex');
    legend_text=[legend_text,strcat({'ESA - '},name_bs)];
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','SouthEast');
    xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\epsilon_{SSH}$ [cm]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'difference_mean_SSH_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);
    %plus
        f1=figure;
    legend_text={''};
    if ~isempty(input_path_L2_STL)
        plot(edges_plus_hr,(ESA_SSH_mean_mean_binned_plus_hr-STL_SSH_mean_mean_binned_plus_hr).*100,'Marker',char(marker_STL),'Color',color_STL);
        legend_text=[legend_text,'ESA - STL'];
        hold on;
    end
    for b=1:N_baselines
        plot(edges_plus_hr,(ESA_SSH_mean_mean_binned_plus_hr-SSH_mean_mean_binned_plus_hr(b,:)).*100,'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
        hold on;
    end
    axis([min_plus_hr max_plus_hr min_error_SSH max_error_SSH]);
    title(strcat('Difference mean SSH versus $v_{r}$'),'Interpreter','latex');
    legend_text=[legend_text,strcat({'ESA - '},name_bs)];
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','SouthEast');
    xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\epsilon_{SSH}$ [cm]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'difference_mean_SSH_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);

    %-------------- sigma0 vs vrad ------------------------------------------------
    min_sigma0=Inf;
    max_sigma0=-Inf;
    %minus
    f1=figure;
    legend_text={''};
    if ~isempty(input_path_L2_STL)
        plot(edges_minus_hr,(ESA_sigma0_mean_mean_binned_minus_hr-STL_sigma0_mean_mean_binned_minus_hr),'Marker',char(marker_STL),'Color',color_STL);
        legend_text=[legend_text,'ESA - STL'];
        hold on;
        min_sigma0=min([min_sigma0,min(ESA_sigma0_mean_mean_binned_minus_hr-STL_sigma0_mean_mean_binned_minus_hr),min(ESA_sigma0_mean_mean_binned_plus_hr-STL_sigma0_mean_mean_binned_plus_hr)]);
        max_sigma0=max([max_sigma0,max(ESA_sigma0_mean_mean_binned_minus_hr-STL_sigma0_mean_mean_binned_minus_hr),max(ESA_sigma0_mean_mean_binned_plus_hr-STL_sigma0_mean_mean_binned_plus_hr)]);
    end
    for b=1:N_baselines
        plot(edges_minus_hr,(ESA_sigma0_mean_mean_binned_minus_hr-sigma0_mean_mean_binned_minus_hr(b,:)),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
        hold on;
        min_sigma0=min([min_sigma0,min(ESA_sigma0_mean_mean_binned_minus_hr-sigma0_mean_mean_binned_minus_hr(b,:)),min(ESA_sigma0_mean_mean_binned_plus_hr-sigma0_mean_mean_binned_plus_hr(b,:))]);
        max_sigma0=max([max_sigma0,max(ESA_sigma0_mean_mean_binned_minus_hr-sigma0_mean_mean_binned_minus_hr(b,:)),max(ESA_sigma0_mean_mean_binned_plus_hr-sigma0_mean_mean_binned_plus_hr(b,:))]);
    end
    axis([min_minus_hr max_minus_hr min_sigma0 max_sigma0]);
    title(strcat('Difference mean $\sigma^0$ versus $v_{r}$'),'Interpreter','latex');
    legend_text=[legend_text,strcat({'ESA - '},name_bs)];
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','Southeast');
    xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\epsilon_{\sigma^0}$ [dB]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'difference_mean_sigma0_vs_minus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);
    %plus
        f1=figure;
    legend_text={''};
    if ~isempty(input_path_L2_STL)
        plot(edges_plus_hr,(ESA_sigma0_mean_mean_binned_plus_hr-STL_sigma0_mean_mean_binned_plus_hr),'Marker',char(marker_STL),'Color',color_STL);
        legend_text=[legend_text,'ESA - STL'];
        hold on;
    end
    for b=1:N_baselines
        plot(edges_plus_hr,(ESA_sigma0_mean_mean_binned_plus_hr-sigma0_mean_mean_binned_plus_hr(b,:)),'Marker',char(marker_bs(b)),'Color',color_bs(b,:));
        hold on;
    end
    axis([min_plus_hr max_plus_hr min_sigma0 max_sigma0]);
    title(strcat('Difference mean $\sigma^0$ versus $v_{r}$'),'Interpreter','latex');
    legend_text=[legend_text,strcat({'ESA - '},name_bs)];
    legend(legend_text(~cellfun(@isempty,legend_text)),'Location','Southeast');
    xlabel('$v_{r}$ [m/s]','Interpreter','latex'); ylabel('$\epsilon_{\sigma^0}$ [dB]','Interpreter','latex');
    print(print_file,res_fig,strcat(path_comparison_results,'difference_mean_sigma0_vs_plus_hr_',strjoin(strrep(strrep(name_bs,' ','_'),'-','_'),'_vs_'),file_ext))
    close(f1);

end

time = toc(ttotal);
minutes_processing = floor(time/60);
secs_processing = time - minutes_processing*60;
disp(['Total processing time: ', num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);

end

