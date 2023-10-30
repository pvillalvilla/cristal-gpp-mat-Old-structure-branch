function [icebergs_properties_L1B_DDP, icebergs_properties_L1BS, iceberg_L1BS_flag, icebergs_properties_FF_ML] = ...
          iceberg_detection(filesBulk,L1B,cnf,cnf_L2,chd,cst)
      
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code allows to run an iceberg detection algorithm over open ocean to
% detect signals on the thermal noise part of the waveforms caused by the
% presence of icebergs. It can run on SAR-L1B, SAR-L1Bs and FF-ML-L1B data
% for CRISTAL.
% -------------------------------------------------------------------------
% 
% Author:           Juan Pedro López-Zaragoza / isardSAT
%
% Reviewer:         ---- / isardSAT
%
% Last revision:    Juan Pedro López-Zaragoza / isardSAT V1 03/06/2021
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       - filesbulk = variable containing paths and datafiles variables.
%       - L1B = struct containing the output data of the L1A to L1B processor.
%       - L1Bs = 
%       - FF-ML = 
%       - cnf = configuration file with several L1 configuration parameters.
%       - cnf_L2 = onfiguration file with several L2 configuration parameters. Several iceberg algorithm configuration parameters are found here.
%       - chd = file containing several characterization parameters.
%       - cst = file containing values of constants.
% OUTPUT:
%       - icebergs_properties_L1B_DDP = struct array containing the following properties of the iceberg candidates detected in SAR-L1B data:
%           -- alongtrack_index: along track indexes of the SAR-L1B waveform in which the candidate is contained.
%           -- region: number of the area after the connected components image creation in which the candidate is located.
%           -- rangebin_ind: range indexes of the L1B waveform in which the candidate is contained.
%           -- power: mean power (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- power: mean power (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- stdpower: standard deviation of the power (in dBs).
%           -- sigma0: mean sigma0 (backscatter) (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- stdsigma0: standard deviation of sigma0 (in dBs).
%           -- freeboard: corresponding freeboard of the iceberg candidate assuming that the iceberg is located just below the satellite track.
%       - icebergs_properties_L1BS = struct array containing the following properties of the iceberg candidates detected in SAR-L1BS data:
%           -- region: number of the area after the connected components image creation in which the candidate is located.
%           -- alongtrack_index: along track indexes of the L1B waveform in which the candidate is contained.
%           -- rangebin_ind: range indexes of the L1B waveform in which the candidate is contained.
%           -- power: mean power (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- power: mean power (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- stdpower: standard deviation of the power (in dBs).
%           -- sigma0: mean sigma0 (backscatter) (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- stdsigma0: standard deviation of sigma0 (in dBs).
%           -- freeboard: corresponding freeboard of the iceberg candidate assuming that the iceberg is located just below the satellite track.
%       - iceberg_L1BS_flag = cell array containg L2 file name for each baseline
%       - icebergs_properties_FF_ML = cell array of size N_baselines containing baselines names
%           -- region: number of the area after the connected components image creation in which the candidate is located.
%           -- alongtrack_index: along track indexes of the L1B waveform in which the candidate is contained.
%           -- rangebin_ind: range indexes of the L1B waveform in which the candidate is contained.
%           -- power: mean power (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- power: mean power (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- stdpower: standard deviation of the power (in dBs).
%           -- sigma0: mean sigma0 (backscatter) (in dBs) of the pixels corresponding to 'alongtrack_index' & 'rangebin_ind' which contain the iceberg candidates.
%           -- stdsigma0: standard deviation of sigma0 (in dBs).
%           -- freeboard: corresponding freeboard of the iceberg candidate assuming that the iceberg is located just below the satellite track.       
% -------------------------------------------------------------------------
% Versions control:
% V1.0: First algorithm version. It works with CRISTAL L1B-SAR data only,
% for one input file which can be divided in segments
%--------------------------------------------------------------------------
% Some variables initialization
icebergs_properties_L1B_DDP=[];
iceberg_L1BS_flag = [];
icebergs_properties_L1BS_DDP=[];
icebergs_properties_L1B_FF_ML=[];

%% L1B DDP ICEBERG DETECTION ALGORITHM

if cnf_L2.iceberg_L1B_DDP_flag
    
    %% Divide the data in segments of a given length to look for icebergs
    L1B_DDP_data_length=size(L1B.scaled_waveforms,1);
    
    jj=0;
    for ii=1:cnf_L2.iceberg_L1B_DDP_segments_length:L1B_DDP_data_length
        jj=jj+1;
        if jj==ceil(L1B_DDP_data_length/cnf_L2.iceberg_L1B_DDP_segments_length) % for last segment of data
            L1B_DDP_segments(jj).lat(1:length(ii:numel(L1B.lat)))=L1B.lat(ii:numel(L1B.lat))';
            L1B_DDP_segments(jj).lon(1:length(ii:numel(L1B.lon)))=L1B.lat(ii:numel(L1B.lon))';
            L1B_DDP_segments(jj).alt(1:length(ii:numel(L1B.alt)))=L1B.lat(ii:numel(L1B.alt))';
            L1B_DDP_segments(jj).range_ku_l1b_echo(1:length(ii:numel(L1B.range_ku_l1b_echo)))=L1B.range_ku_l1b_echo(ii:numel(L1B.range_ku_l1b_echo))';
            L1B_DDP_segments(jj).sigma0(1:length(ii:numel(L1B.sigma0)))=L1B.sigma0(ii:numel(L1B.sigma0))';
            
            L1B_DDP_segments(jj).scaled_waveforms(1:length(ii:size(L1B.scaled_waveforms,1)),1:size(L1B.scaled_waveforms,2))=L1B.scaled_waveforms(ii:size(L1B.scaled_waveforms,1),1:size(L1B.scaled_waveforms,2));
            L1B_DDP_segments(jj).scaled_waveforms=L1B_DDP_segments(jj).scaled_waveforms';
        else % filling with the data segments
            L1B_DDP_segments(jj).lat(1:cnf_L2.iceberg_L1B_DDP_segments_length)=L1B.lat(ii:ii+cnf_L2.iceberg_L1B_DDP_segments_length-1);
            L1B_DDP_segments(jj).lon(1:cnf_L2.iceberg_L1B_DDP_segments_length)=L1B.lon(ii:ii+cnf_L2.iceberg_L1B_DDP_segments_length-1)';
            L1B_DDP_segments(jj).alt(1:cnf_L2.iceberg_L1B_DDP_segments_length)=L1B.alt(ii:ii+cnf_L2.iceberg_L1B_DDP_segments_length-1)';
            L1B_DDP_segments(jj).range_ku_l1b_echo(1:cnf_L2.iceberg_L1B_DDP_segments_length)=L1B.range_ku_l1b_echo(ii:ii+cnf_L2.iceberg_L1B_DDP_segments_length-1)';
            L1B_DDP_segments(jj).sigma0(1:cnf_L2.iceberg_L1B_DDP_segments_length)=L1B.sigma0(ii:ii+cnf_L2.iceberg_L1B_DDP_segments_length-1)';
            
            L1B_DDP_segments(jj).scaled_waveforms(1:cnf_L2.iceberg_L1B_DDP_segments_length,1:size(L1B.scaled_waveforms,2))=L1B.scaled_waveforms(ii:ii+cnf_L2.iceberg_L1B_DDP_segments_length-1,1:size(L1B.scaled_waveforms,2));
            L1B_DDP_segments(jj).scaled_waveforms=L1B_DDP_segments(jj).scaled_waveforms';
        end
    end
    
    disp(['------------Data divided in ', num2str(jj), ' segments------------'])
    
    %% Latitude range selection
    z=0;
    % Starting loop for the different data segments
    for kk=1:ceil(L1B_DDP_data_length/cnf_L2.iceberg_L1B_DDP_segments_length)
        disp(['------------Processing segment ', num2str(kk),'/', num2str(ceil(L1B_DDP_data_length/cnf_L2.iceberg_L1B_DDP_segments_length)), ' ------------'])
        
        if any(L1B_DDP_segments(kk).lat > 20) || any(L1B_DDP_segments(kk).lat < -50)
            % Some latitude within the processing range, process the segment
            latitude_check(kk)=1;
            
            %% Align Waveforms and build window
            
            Bw=cnf.range_BW_to_process; % in Hz
            % dont we have zp factor? 'range_oversampling_factor' in S6 data. We
            % need it to build the window of the tracker right?
            zp_factor_DDP=cnf.zp_fact_range; % Is this?
            c_cst=cst.c;
            
            clear geoid_correction_DDP elev_DDP top_elev_DDP bottom_elev_DDP elev_axis_DDP wvf_raw_shifted_DDP
            geoid_correction_DDP=geoidheight(L1B_DDP_segments(kk).lat, L1B_DDP_segments(kk).lon);
            elev_DDP = (L1B_DDP_segments(kk).alt-L1B_DDP_segments(kk).range_ku_l1b_echo)-geoid_correction_DDP;
            top_elev_DDP  = elev_DDP(1)+(size(L1B_DDP_segments(kk).scaled_waveforms,1)/2/Bw/zp_factor_DDP*c_cst/2)*3;
            bottom_elev_DDP    = elev_DDP(1)-((size(L1B_DDP_segments(kk).scaled_waveforms,1)/2)/Bw/zp_factor_DDP*c_cst/2)*3;
            elev_axis_DDP = (top_elev_DDP:-1/Bw/zp_factor_DDP*c_cst/2:bottom_elev_DDP+1/Bw/zp_factor_DDP*c_cst/2);
            
            wvf_raw_shifted_DDP = cat(1,0.*L1B_DDP_segments(kk).scaled_waveforms,0.*L1B_DDP_segments(kk).scaled_waveforms,0.*L1B_DDP_segments(kk).scaled_waveforms);
            clear shift_mat
            
            [~,reference_elevation_surf]=max(elev_DDP);
            for i_surf=1:length(L1B_DDP_segments(kk).range_ku_l1b_echo)
                shift_mat(i_surf)=round((elev_DDP(reference_elevation_surf)-elev_DDP(i_surf))/c_cst*zp_factor_DDP*2*Bw); %This extra 2 is for zeropadding
                %shift_mat(i_surf)=round((elev_DDP(1)-elev_DDP(i_surf))/c_cst*2*Bw);
                %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
                wvf_raw_shifted_DDP(size(L1B_DDP_segments(kk).scaled_waveforms,1) + 1 +shift_mat(i_surf):2*size(L1B_DDP_segments(kk).scaled_waveforms,1)+shift_mat(i_surf),i_surf)= L1B_DDP_segments(kk).scaled_waveforms(:,i_surf);
            end
            
            %% L1B power waveform plot creation
            figure;
            imagesc(L1B_DDP_segments(kk).lat, elev_axis_DDP, 10.*log10(wvf_raw_shifted_DDP));
            title('DDP')
            set(gca,'YDir', 'normal');
            hcb=colorbar;
            hcb.Label.String = 'Rx Power [dB]';
            % xlim([-59.5033  -59.4054])
            xlabel('Latitude [degrees]')
            ylabel('Window Elevation [m]')
            load('colormap_blues.mat');
            colormap(colormap_blues);
            
            %% Iceberg detection L1B algorithm
            % Initialize some variables
            i=0;
            j=0;
            clear wvf_raw_shifted_DDP_norm wvf_raw_shifted_DDP_norm_xrecord i2_DDP closeBW_DDP X_DDP
            
            % WVs NORMALIZATION. Here we normalize the waveforms
            % Compute maximum value in each column
            L1B_max_values_xrow= max(wvf_raw_shifted_DDP);
            % Compute total maximum value and its wvform index
            [L1B_max_total,L1B_ind_max_total]=max(L1B_max_values_xrow);
            % Normalize it (we use this for the threshold definition)
            wvf_raw_shifted_DDP_norm =  wvf_raw_shifted_DDP ./ L1B_max_total;
            
            %Normalize by record along-track (we don't use this for the threshold definition)
            for k=1:size(wvf_raw_shifted_DDP,2)
                wvf_raw_shifted_DDP_norm_xrecord(:,k) = wvf_raw_shifted_DDP(:,k) ./ L1B_max_values_xrow(k);
            end
            % We manually select a treshold of 4% (Tournardre presentation) in
            % power wrt the maximum to define a positive from an iceberg
            iceberg_threshold=0.04;
            
            % BINARY IMAGE CREATION. We only keep the values above the selected threshold to create a binary image
            %i2= reduced_fp_0928_DDP_p_norm>iceberg_threshold;
            i2_DDP= wvf_raw_shifted_DDP_norm>iceberg_threshold;
            figure;
            imagesc(i2_DDP);
            colormap(gray);
            %freezeColors;
            title ('Binary image');
            
            % BINARY IMAGE FILLING. Here we fill the empty sections of the binary image to define connected areas later on
            % Here we define vertical lines to fill the areas connected in less than the selected number of pixels
            se = strel('line',4,90);
            closeBW_DDP = imclose(i2_DDP,se);
            figure;
            imagesc(closeBW_DDP);
            colormap(gray);
            %freezeColors;
            title ('Binary image, holes filled');
            
            % CONNECTED COMPONENTS IMAGE CREATION. Here we define the different areas which are connected
            X_DDP=bwlabel(closeBW_DDP);
            figure;
            imagesc(X_DDP);
            colormap(jet);
            %freezeColors;
            title ('Connected components');
            
            % Here we create a histogram with all the coordinates of the connected
            % components map, and save the position of the biggest area (excluding
            % the zeros, ie background)
            clear rangebin_indexes alongtrack_index rangebinocean rangebin100ocean rangebin87ocean
            [yy,xx]=hist(reshape (X_DDP,1,size(X_DDP,1)*size(X_DDP,2)),size(X_DDP,1)*size(X_DDP,2));
            [maxim,posmax]=max(yy(2:end));
            % Numeric code of ocean/big area
            oceanWV_areacode=round(xx(posmax));
            % Coordinates of ocean/big area
            [rangebin_indexes, alongtrack_index]=find(X_DDP==round(xx(posmax)));
            for j=1:size(X_DDP,2)
                [rangebin_indexes, alongtrack_index]=find(X_DDP(:,j)==round(xx(posmax)));
                % Find the 1st point in range of the ocean/big area
                rangebinocean(j)=min(rangebin_indexes);
                % Peak of ocean
                rangebin100ocean(j)=find(wvf_raw_shifted_DDP_norm_xrecord(:,j)==1);
                % 87% of the peak of ocean
                a=wvf_raw_shifted_DDP_norm_xrecord(1:rangebin100ocean(j),j);
                [val,idx]=min(abs(a-0.87));
                rangebin87ocean(j)=idx;
            end
            
            % Loop through all the along track records in the data segment
            for j=1:size(X_DDP,2)
                % Find the different regions ID present within the same along track record
                unique_regions=unique(X_DDP(:,j));
                % Loop through the different regions (in the same along track record) to select the ones with iceberg candidates
                for jj=2:length(unique_regions) % We skip 1st region, since it is the zeros region (values below threshold/non-positive values)
                    region=unique_regions(jj);
                    % Coordinates of the current region under analysis
                    [rangebin_indexes, alongtrack_index]=find(X_DDP(:,j)==region);
                    % We only select the regions which are before the ocean leading edge and with region code different from the ocean/big area
                    if min(rangebin_indexes) < min(rangebinocean(alongtrack_index)) && region~=oceanWV_areacode
                        z=z+1;
                        abs_alon_track_ind=(cnf_L2.iceberg_L1B_DDP_segments_length*(kk-1))+j;
                        % Save properties of positives in a struct
                        icebergs_properties_L1B_DDP(z).alongtrack_index=abs_alon_track_ind;
                        icebergs_properties_L1B_DDP(z).regionID=region;
                        icebergs_properties_L1B_DDP(z).rangebin_ind=rangebin_indexes;
                        clear power_iceberg;
                        clear sigma0;
                        clear nadir_freeboard_iceberg;
                        for ii=1:length(rangebin_indexes)
                            power_iceberg(ii)=10.*log10(wvf_raw_shifted_DDP(rangebin_indexes(ii),j));
                            sigma0(ii)=power_iceberg(ii)+L1B_DDP_segments(kk).sigma0(j);
                            % freeboard here is the distance from the positive rangebin until the 87% of the peak of the waveform edge
                            nadir_freeboard_iceberg(ii)=abs(-elev_axis_DDP(rangebin87ocean(j)) - -elev_axis_DDP(rangebin_indexes(ii)));
                        end
                        icebergs_properties_L1B_DDP(z).power=mean(power_iceberg);
                        icebergs_properties_L1B_DDP(z).stdpower=std(power_iceberg);
                        icebergs_properties_L1B_DDP(z).sigma0=mean(sigma0);
                        icebergs_properties_L1B_DDP(z).stdsigma0=std(sigma0);
                        icebergs_properties_L1B_DDP(z).nadir_freeboard=max(nadir_freeboard_iceberg);
                    else
                        % When no iceberg is found, set properties to default
                        z=z+1;
                        abs_alon_track_ind=(cnf_L2.iceberg_L1B_DDP_segments_length*(kk-1))+j;
                        % Save properties of positives in a struct
                        icebergs_properties_L1B_DDP(z).alongtrack_index=abs_alon_track_ind;
                        icebergs_properties_L1B_DDP(z).regionID=region;
                        icebergs_properties_L1B_DDP(z).rangebin_ind=[];
                        icebergs_properties_L1B_DDP(z).power=[];
                        icebergs_properties_L1B_DDP(z).stdpower=[];
                        icebergs_properties_L1B_DDP(z).sigma0=[];
                        icebergs_properties_L1B_DDP(z).stdsigma0=[];
                        icebergs_properties_L1B_DDP(z).nadir_freeboard=[];
                    end
                end
            end
        else
            % No latitude within the processing range, jump to next data segment
            latitude_check(kk)=0;
        end % end of loop to check latitudes within segment
    end % end of loop over data segments
end

%% L1BS ICEBERG DETECTION ALGORITHM
if L1BS_flag
    
    icebergs_properties_L1BS=[];
    
    %% Reading netCDF L1B data (in case the user didn't run the L1B chain, we have to read the data again here since we use some L1B parameters)
    DDP_file=dir([ input_path '\*HR__1B*.nc']);
    
    % Load DDP waveforms
    fp_0928_DDP = [input_path '\' DDP_file.name];
    fp_0928_DDP_lat = ncread(fp_0928_DDP, 'data_20/ku/latitude');
    fp_0928_DDP_lon = ncread(fp_0928_DDP, 'data_20/ku/longitude');
    fp_0928_DDP_alt = ncread(fp_0928_DDP, 'data_20/ku/altitude');
    fp_0928_DDP_trk_range = ncread(fp_0928_DDP, 'data_20/ku/tracker_range_calibrated');
    fp_0928_DDP_p_raw = ncread(fp_0928_DDP, 'data_20/ku/power_waveform');
    fp_0928_DDP_wsc = ncread(fp_0928_DDP, 'data_20/ku/waveform_scale_factor');
    %fp_0928_DDP_p   = 10*log10(fp_0928_DDP_p_raw.*((ones(512,1)*fp_0928_DDP_wsc.')));
    fp_0928_DDP_p   = (fp_0928_DDP_p_raw.*((ones(size(fp_0928_DDP_p_raw,1),1)*fp_0928_DDP_wsc.')));
    zp_factor_DDP = ncread(fp_0928_DDP, '/global/ku/range_oversampling_factor');
    
    geoid_correction_DDP=geoidheight(fp_0928_DDP_lat.',fp_0928_DDP_lon.');
    lla2kml([input_path '\' DDP_file.name],fp_0928_DDP_lat(1:100:end),fp_0928_DDP_lon(1:100:end),100+zeros(length(fp_0928_DDP_lon(1:100:end)),1),'.');
    
    
    %% Read netCDF L1Bs data
    disp('------------Reading L1BS data------------')
    %initialize some variables before the loop through all the stacks
    i=0;
    j=0;
    % *****Define the alongtrack indexes to be analysed (THIS IS HARDCODED, HAS TO BE MODIFIED)
    index_L1B_start=10101;
    index_L1B_end=11100;
    for index_L1B=index_L1B_start:index_L1B_end
        % Some variables initialization
        j=j+1;
        iceberg_flag=0;
        
        % Reading L1BS file and detecting if there is more than 1 input file
        L1BS_file=dir([ input_path '\*HR__1S*.nc']);
        if(numel(L1BS_file)>1)
            error('Please provide just one L1BS input file at a time');
        end
            
        fp_0928_L1BS_refa_1 = [input_path L1BS_file.name];
        fp_0928_L1BS_refa_1_lat = ncread(fp_0928_L1BS_refa_1, 'data_20/ku/latitude');
        index_L1BS=find(fp_0928_L1BS_refa_1_lat==fp_0928_DDP_lat(index_L1B));
             
        try
            samples_L1BS = length(ncread(fp_0928_L1BS_refa_1,'samples'));
        catch
            samples_L1BS = length(ncread(fp_0928_L1BS_refa_1,'samples_ov'));
        end
        looks_L1BS = length(ncread(fp_0928_L1BS_refa_1,'looks'));
        %reading L1BS vars
       
        i_look   = ncread(fp_0928_L1BS_refa_1,'data_20/ku/look_i_samples',[1 1 index_L1BS],[samples_L1BS looks_L1BS 1]);
        q_look   = ncread(fp_0928_L1BS_refa_1,'data_20/ku/look_q_samples',[1 1 index_L1BS],[samples_L1BS looks_L1BS 1]);
        
        iq_scale_L1BS    = ncread(fp_0928_L1BS_refa_1,'data_20/ku/iq_scale_factor',[1 index_L1BS],[looks_L1BS 1]);
        
        % Create stack matrix in different flavours
        stack_L1BS   = ((i_look+1i*q_look).*(ones(samples_L1BS,1)*iq_scale_L1BS.')); % complex stack
        stack_L1BS_dB = 10*log10(abs(stack_L1BS.').^2); % in dBs
        stack_L1BS_real = (abs(stack_L1BS.').^2); % real part, not dBs
        
        % Plot of stack in dBs
        if L1BS_checking_plots
            figure;
            imagesc(1:samples_L1BS,1:looks_L1BS,  10*log10(abs(stack_L1BS.').^2))
            set(gca,'XLim',[1 samples_L1BS])
            set(gca,'YLim',[1 looks_L1BS])
            colorbar;
            xlab = get(gca,'XLabel'); set(xlab,'String','Range index')
            ylab = get(gca,'YLabel'); set(ylab,'String','Doppler beam index')
            title(['RMC stack #' num2str(index_L1BS)  ' record, latitude: ' num2str(fp_0928_L1BS_refa_1_lat(index_L1BS)) ])
            hcb=colorbar; colormap(jet)
            hcb.Label.String = 'Stack Power [dB]';
            addlogoisardSAT('plot');
        end
        
        %% Iceberg detection L1Bs chain
        
        if index_L1B==10440
        %if index_L1B==10123
        %if index_L1B==10255
            stop=1;
        end
        % WVs NORMALIZATION. Here we normalize the waveforms
        % Compute maximum value in each column
        L1BS_max_values_xrow= max(stack_L1BS_real);
        % Compute total maximum value and its wvform index
        [L1BS_max_total(index_L1B),L1BS_ind_max_total(index_L1B)]=max(L1BS_max_values_xrow);
        % Normalize it
        stack_L1BS_dB_norm =  stack_L1BS_real ./ L1BS_max_total(index_L1B);
        
        % Plot of normalized stack
        if L1BS_checking_plots
            figure;
            imagesc(1:samples_L1BS,1:looks_L1BS,  stack_L1BS_dB_norm)
            set(gca,'XLim',[1 samples_L1BS])
            set(gca,'YLim',[1 looks_L1BS])
            colorbar;
            xlab = get(gca,'XLabel'); set(xlab,'String','Range index')
            ylab = get(gca,'YLabel'); set(ylab,'String','Doppler beam index')
            title(['RMC normalized stack #' num2str(index_L1BS)  ' record, latitude: ' num2str(fp_0928_L1BS_refa_1_lat(index_L1BS)) ])
            hcb=colorbar; colormap(jet)
            hcb.Label.String = 'Normalized stack';
            addlogoisardSAT('plot');
        end
        
        % We manually select a treshold of 4% (Tournardre presentation) in
        % power wrt the maximum to define a positive from an iceberg
        iceberg_threshold=0.04;
        %iceberg_threshold=0.05;
        % BINARY IMAGE CREATION. We only keep the values above the selected threshold to create a binary image
        i2= stack_L1BS_dB_norm>iceberg_threshold;
        if L1BS_checking_plots
            figure;
            imagesc(i2);
            colormap(gray);
            freezeColors;
            title ('Binary image');
            axis off;
        end
        % BINARY IMAGE FILLING. Here we fill the empty sections of the binary image to define connected areas later on
        % Here we define vertical lines to fill the areas connected in less than the selected number of pixels
        se = strel('cube',8);
        se = strel('line',10,90);
        se = strel('rectangle',[10 5]);
        closeBW = imclose(i2,se);
        if L1BS_checking_plots
            figure;
            imagesc(closeBW);title ('Binary image, holes filled');
            axis off;
        end
        % CONNECTED COMPONENTS IMAGE CREATION. Here we define the different areas which are connected
        X=bwlabel(closeBW, 8);
        if L1BS_checking_plots
            figure;
            imagesc(X);
            colormap(jet);
            freezeColors;
            title ('Connected components');
            axis off;
        end
    
        % Here we create a histogram with all the coordinates of the connected
        % components map, and save the position of the biggest area (excluding
        % the zeros, ie background)
        [yy,xx]=hist(reshape (X,1,size(X,1)*size(X,2)),size(X,1)*size(X,2));
        [maxim,posmax]=max(yy(2:end));
        % Numeric code of ocean/big area
        oceanWV_areacode=round(xx(posmax));
        % Coordinates of ocean/big area
        [dopplerbeam_index,rangebin_indexes]=find(X==round(xx(posmax)));
        % Find the 1st point in range of the ocean/big area
        rangebinocean=min(rangebin_indexes);
        
        % Loop through all the different regions in the stack to select the regions with iceberg candidates
        for region=(1+min(min(X))):max(max(X))
            % Coordinates of the current region under analysis
            [dopplerbeam_index,rangebin_indexes]=find(X==region);
            % We compute region area to exclude 1 pixel positives in doppler beam indexes that are not related to icebergs, but to high noise
            uniquevalues_doppler_beam=unique(dopplerbeam_index);
            doppler_region_area=numel(uniquevalues_doppler_beam);
            % We only select the regions which are before the ocean leading
            % edge and with region code different from the ocean/big area and with area bigger than 1 pixel (which are typically noise)
            if rangebinocean > min(rangebin_indexes) && region~=oceanWV_areacode && doppler_region_area > 1
                % There is an iceberg region candidate within this stack
                iceberg_flag=1;
                
%                 %This flag saves a 1 if the region we analyse is left of the oceanWV region code
%                 %It saves a code for each region inside each stack (region 0, ie background, is not considered)
%                 flag_iceberg_L1BS(index_L1B,region).region=region; % THINK HOW TO YIELD ONLY A BINARY OUTPUT FOR EACK STACK
%                 flag_iceberg_L1BS(index_L1B,region).flag=1;
                
                % save properties of positives in a struct
                i=i+1;
                % THINK WHICH INDEX/RECORD NUMBER TO SAVE, ONLY 1
                icebergs_properties_L1BS(i).alongtrack_record=j;
                icebergs_properties_L1BS(i).index_L1B=index_L1B;
                icebergs_properties_L1BS(i).stack_L1BS=index_L1BS;
                icebergs_properties_L1BS(i).region=region;
                icebergs_properties_L1BS(i).dopplerbeam_ind=dopplerbeam_index;
                icebergs_properties_L1BS(i).rangebin_ind=rangebin_indexes;   
                for ii=1:length(dopplerbeam_index)
                    power_iceberg(ii)=stack_L1BS_dB(icebergs_properties_L1BS(i).dopplerbeam_ind(ii),icebergs_properties_L1BS(i).rangebin_ind(ii));
                end     
                icebergs_properties_L1BS(i).power=mean(power_iceberg(~isinf(power_iceberg)));
                icebergs_properties_L1BS(i).sigma=std(power_iceberg(~isinf(power_iceberg)));
            else
                %                 flag_iceberg_L1BS(index_L1B,region).region=region;
                %                 flag_iceberg_L1BS(index_L1B,region).flag=0;
            end
            
        end
        % Save a single iceberg output flag for each along track sample
        iceberg_L1BS_flag(j).alongtrack_record=j;
        iceberg_L1BS_flag(j).iceberg_flag=iceberg_flag;
        
        disp(['Processed stack ', num2str(index_L1B),'/',num2str(index_L1B_end)]);
    end
    
    for i=1:numel(iceberg_L1BS_flag)
        icebergs_array(i)=iceberg_L1BS_flag(i).iceberg_flag;
    end
    numberstacks_icebergs=nnz(icebergs_array);
    
end

%% L1B FF ICEBERG DETECTION ALGORITHM
if FF_flag
    %% Reading netCDF data
    disp('------------Reading FF L1B data------------')
    
    SL_file=dir([ input_path '\*L1B_FF_SL_*.nc']);
    % Reading L1B-FF-ML file and detecting if there is more than 1 input file 
    ML_file=dir([ input_path '\*L1B_FF_ML_*.nc']);        
    if(numel(ML_file)>1)
        error('Please provide just one FF-ML-L1B input file at a time');
    end
    
    % Load ML waveforms
    fp_0928_ML_refa_1 = [input_path '\' ML_file.name];
    fp_0928_ML_refa_1_lat = ncread(fp_0928_ML_refa_1, 'data/ku/latitude');
    fp_0928_ML_refa_1_lon = ncread(fp_0928_ML_refa_1, 'data/ku/longitude');
    fp_0928_ML_refa_1_alt = ncread(fp_0928_ML_refa_1, 'data/ku/altitude');
    fp_0928_ML_refa_1_trk_range = ncread(fp_0928_ML_refa_1, 'data/ku/tracker_range_calibrated');
    fp_0928_ML_refa_1_p_raw = ncread(fp_0928_ML_refa_1, 'data/ku/power_waveform');
    fp_0928_ML_refa_1_wsc = ncread(fp_0928_ML_refa_1, 'data/ku/waveform_scale_factor');
    fp_0928_ML_refa_1_p = (fp_0928_ML_refa_1_p_raw).*(ones(size(fp_0928_ML_refa_1_p_raw,1),1)*fp_0928_ML_refa_1_wsc');
    %fp_0928_ML_refa_1_p = 10*log10(fp_0928_ML_refa_1_p_raw).*(ones(256,1)*fp_0928_ML_refa_1_wsc');
    zp_factor_ML = ncread(fp_0928_ML_refa_1, '/global/ku/range_oversampling_factor');
    sig0_scaling_factor_FF_ML = ncread(fp_0928_ML_refa_1, 'data/ku/sig0_scaling_factor');
      
%     % Load SL waveforms
%     fp_0928_SL_refa_1 = [input_path '\' SL_file.name];
%     fp_0928_SL_refa_1_lat = ncread(fp_0928_SL_refa_1, 'data/ku/latitude');
%     fp_0928_SL_refa_1_alt = ncread(fp_0928_SL_refa_1, 'data/ku/altitude');
%     fp_0928_SL_refa_1_trk_range = ncread(fp_0928_SL_refa_1, 'data/ku/tracker_range_calibrated');
%     fp_0928_SL_refa_1_i_raw = ncread(fp_0928_SL_refa_1, 'data/ku/i_samples');
%     fp_0928_SL_refa_1_q_raw = ncread(fp_0928_SL_refa_1, 'data/ku/q_samples');
%     fp_0928_SL_refa_1_iqscf = ncread(fp_0928_SL_refa_1, 'data/ku/iq_scale_factor');
%     zp_factor_SL = ncread(fp_0928_SL_refa_1, '/global/ku/range_oversampling_factor');
%     
%     fp_0928_SL_refa_1_i = fp_0928_SL_refa_1_i_raw.*(ones(size(fp_0928_SL_refa_1_q_raw,1),1)*fp_0928_SL_refa_1_iqscf');
%     fp_0928_SL_refa_1_q = fp_0928_SL_refa_1_q_raw.*(ones(size(fp_0928_SL_refa_1_q_raw,1),1)*fp_0928_SL_refa_1_iqscf');
%     clear fp_0928_SL_refa_1_i_raw fp_0928_SL_refa_1_q_raw
%     fp_0928_SL_refa_1_p_raw = sqrt(real((fp_0928_SL_refa_1_i + j*fp_0928_SL_refa_1_q).*(fp_0928_SL_refa_1_i - j*fp_0928_SL_refa_1_q)));
%     clear fp_0928_SL_refa_1_i fp_0928_SL_refa_1_q
%     % Estimate power from SL waveforms
%     fp_0928_SL_refa_1_p = (fp_0928_SL_refa_1_p_raw);
%     %fp_0928_SL_refa_1_p = 10*log10(fp_0928_SL_refa_1_p_raw);
    
    
        
    % fp_0928_ML_refa_1_p = 10*log10(fp_0928_ML_refa_1_p/max(max(fp_0928_ML_refa_1_p)));
    % fp_0928_SL_refa_1_p = 10*log10(fp_0928_SL_refa_1_p/max(max(fp_0928_SL_refa_1_p)));
    % fp_0928_DDP_p = 10*log10(fp_0928_DDP_p/max(max(fp_0928_DDP_p)));
    
    % fp_0928_ML_refa_1_p = 10*log10(fp_0928_ML_refa_1_p);
    % fp_0928_SL_refa_1_p = 10*log10(fp_0928_SL_refa_1_p);
    % fp_0928_DDP_p = 10*log10(fp_0928_DDP_p);
    % align power between DDP and FF
    %ratio_ML = max(max(fp_0928_DDP_p))/max(max(fp_0928_ML_refa_1_p));
    %ratio_SL = max(max(fp_0928_DDP_p))/max(max(fp_0928_SL_refa_1_p));
    %fp_0928_ML_refa_1_p=fp_0928_ML_refa_1_p.*ratio_ML;
    %fp_0928_SL_refa_1_p=fp_0928_SL_refa_1_p.*ratio_ML;
    
    
    %% Select the region to analyse, only open water samples
    % WE HAVE TO MODIFY HOW TO WORK IN REDUCED SEGMENTS
   % reduced_fp_0928_DDP_p=fp_0928_DDP_p(1:512,10101:11100);
    reduced_fp_0928_ML_refa_1_p_raw=fp_0928_ML_refa_1_p_raw(1:512,1175:9978);
    reduced_fp_0928_ML_refa_1_p=fp_0928_ML_refa_1_p(1:512,1175:9978);
    
 
%     reduced_fp_0928_SL_refa_1_p=abs(fp_0928_SL_refa_1_p(1:512,36397:309305).^2);
%     reduced_fp_0928_SL_refa_1_p_raw=abs(fp_0928_SL_refa_1_p_raw(1:512,36397:309305).^2);
%     

    % raw data plots
    figure;
    imagesc(reduced_fp_0928_ML_refa_1_p_raw);
    xlabel('Along track index')
    ylabel('Range index')
    title('FF ML waveforms raw');
%     figure;
%     imagesc(fp_0928_SL_refa_1_p_raw(1:512,36397:309305));
%     xlabel('Along track index')
%     ylabel('Range index')
%     title('FF SL waveforms');
    
%     %corrected? data plots
%     figure;
%     imagesc(reduced_fp_0928_ML_refa_1_p);
%     xlabel('Along track index')
%     ylabel('Range index')
%     title('FF ML waveforms *ratioML');
    
%     figure;
%     imagesc(reduced_fp_0928_SL_refa_1_p);
%     xlabel('Along track index')
%     ylabel('Range index')

    
    %% Correct for window delay
    % SL waveforms
    
%     elev_SL = (fp_0928_SL_refa_1_alt-fp_0928_SL_refa_1_trk_range);
%     top_elev_SL  = elev_SL(1)+size(fp_0928_SL_refa_1_p,1)/2/Bw/zp_factor_SL*c_cst/2;
%     bottom_elev_SL    = elev_SL(1)-(size(fp_0928_SL_refa_1_p,1)/2)/Bw/zp_factor_SL*c_cst/2;
%     elev_axis_SL = (top_elev_SL-1/Bw/zp_factor_SL*c_cst/2:-1/Bw/zp_factor_SL*c_cst/2:bottom_elev_SL);
%     
%     clear wvf_raw_shifted_SL
%     wvf_raw_shifted_SL=cat(1,0.*fp_0928_SL_refa_1_p,0.*fp_0928_SL_refa_1_p,0.*fp_0928_SL_refa_1_p);
%     clear shift_mat
%     for i_surf=1:length(fp_0928_SL_refa_1_trk_range)
%         %shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/c_cst*2*2*Bw); %This extra 2 is for zeropadding
%         shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/c_cst*2*Bw*zp_factor_SL);
%         %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
%         wvf_raw_shifted_SL(size(fp_0928_SL_refa_1_p,1) + 1 +shift_mat(i_surf):(size(fp_0928_SL_refa_1_p,1)*2+shift_mat(i_surf)),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
%     end
%     
    % ML waveforms
    elev_ML = (fp_0928_ML_refa_1_alt-fp_0928_ML_refa_1_trk_range);
    top_elev_ML  = elev_ML(1)+size(fp_0928_ML_refa_1_p,1)/2/Bw/zp_factor_ML*c_cst/2;
    bottom_elev_ML    = elev_ML(1)-(size(fp_0928_ML_refa_1_p,1)/2)/Bw/zp_factor_ML*c_cst/2;
    elev_axis_ML = (top_elev_ML:-1/Bw/zp_factor_ML*c_cst/2:bottom_elev_ML-1/Bw/zp_factor_ML*c_cst/2);
    
    clear wvf_raw_shifted_ML
    wvf_raw_shifted_ML=cat(1,0.*fp_0928_ML_refa_1_p,0.*fp_0928_ML_refa_1_p,0.*fp_0928_ML_refa_1_p);
    clear shift_mat
    for i_surf=1:length(fp_0928_ML_refa_1_trk_range)
        %shift_mat(i_surf)=round((elev_SL(1)-elev_SL(i_surf))/c_cst*2*2*Bw); %This extra 2 is for zeropadding
        shift_mat(i_surf)=round((elev_ML(1)-elev_ML(i_surf))/c_cst*2*Bw*zp_factor_ML);
        %wvf_raw_shifted_SL(513+shift_mat(i_surf):1024+shift_mat(i_surf),i_surf)= fp_0928_SL_refa_1_p(:,i_surf) ;
        wvf_raw_shifted_ML(  size(fp_0928_ML_refa_1_p,1) + 1+shift_mat(i_surf):(2*size(fp_0928_ML_refa_1_p,1)+shift_mat(i_surf)),i_surf)= fp_0928_ML_refa_1_p(:,i_surf) ;
    end
    %% L1B power waveform plot creation
    % Reduced segment
    reduced_wvf_raw_shifted_ML=wvf_raw_shifted_ML(:,1175:9978);
    
    figure;
    imagesc(fp_0928_ML_refa_1_lat, elev_axis_ML, 10.*log10(reduced_wvf_raw_shifted_ML));
    title('FF ML reduced')
    set(gca,'YDir', 'normal');
    hcb=colorbar;
    hcb.Label.String = 'Rx Power [dB]';
    % xlim([-59.5033  -59.4054])
    xlabel('Latitude [degrees]')
    ylabel('Window Elevation [m]')
    load('colormap_blues.mat');
    colormap(colormap_blues);
    
    %% Iceberg detection FF L1B chain
    % Initialize some variables
    i=0;
    j=0;
    % WVs NORMALIZATION. Here we normalize the waveforms  
    % Compute maximum value in each column
    FF_ML_max_values_xrow= max(reduced_wvf_raw_shifted_ML);
    % Compute total maximum value and its wvform index
    [FF_ML_max_total,FF_ML_ind_max_total]=max(FF_ML_max_values_xrow);
    % Normalize it
    reduced_wvf_raw_shifted_ML_norm =  reduced_wvf_raw_shifted_ML ./ FF_ML_max_total;
    
    %Normalize by record along-track
    for k=1:size(reduced_wvf_raw_shifted_ML,2)
        reduced_wvf_raw_shifted_ML_norm_xrecord(:,k) = reduced_wvf_raw_shifted_ML(:,k) ./ FF_ML_max_values_xrow(k);      
    end
%    % Some plots (DELETE IN FINAL VERSION)
%    figure;
%    mesh(reduced_fp_0928_DDP_p_norm);
%    % Create a histogram to select threshold (DELETE IN FINAL VERSION IF WE MANUALLY SELECT THE TRESHOLD TO 4%)
%     [f,xi] = ksdensity(reshape (reduced_fp_0928_DDP_p_norm,1,size(reduced_fp_0928_DDP_p_norm,1)*size(reduced_fp_0928_DDP_p_norm,2)));
%     [histy,histx]=hist(reshape (reduced_fp_0928_DDP_p_norm,1,size(reduced_fp_0928_DDP_p_norm,1)*size(reduced_fp_0928_DDP_p_norm,2)),50);
%     figure; %histogram and threshold selection plot
%     plot(xi,f*max(histy)/max(f),'r'); hold all; bar(histx,histy);
%     [minim,pos]=findpeaks(-1*f);
    % We manually select a treshold of 4% (Tournardre presentation) in
    % power wrt the maximum to define a positive from an iceberg
    iceberg_threshold=0.04; 
    % BINARY IMAGE CREATION. We only keep the values above the selected threshold to create a binary image
    i2= reduced_wvf_raw_shifted_ML_norm_xrecord>iceberg_threshold;
    figure;
    imagesc(i2); 
    colormap(gray); 
    freezeColors; 
    title ('Binary image'); 
    % BINARY IMAGE FILLING. Here we fill the empty sections of the binary image to define connected areas later on
    % Here we define vertical lines to fill the areas connected in less than the selected number of pixels
    se = strel('line',3,90);
    closeBW = imclose(i2,se);
    figure;
    imagesc(closeBW);
    title ('Binary image, holes filled');
    % CONNECTED COMPONENTS IMAGE CREATION. Here we define the different areas which are connected
    X=bwlabel(closeBW);
    figure; 
    imagesc(X);
    colormap(jet); 
    freezeColors; 
    title ('Connected components'); 
    
    % Here we create a histogram with all the coordinates of the connected
    % components map, and save the position of the biggest area (excluding
    % the zeros, ie background)
    [yy,xx]=hist(reshape (X,1,size(X,1)*size(X,2)),size(X,1)*size(X,2));
    [maxim,posmax]=max(yy(2:end));
    % Numeric code of ocean/big area
    clear oceanWV_areacode;
    oceanWV_areacode=round(xx(posmax));
    % Coordinates of ocean/big area
    clear rangebin_indexes;
    clear alongtrack_index;
    clear rangebinocean;
    clear rangebin100ocean;
    clear rangebin87ocean;
    
    [rangebin_indexes, alongtrack_index]=find(X==round(xx(posmax)));
    for j=1:size(X,2)
        [rangebin_indexes, alongtrack_index]=find(X(:,j)==round(xx(posmax)));
        % Find the 1st point in range of the ocean/big area
        rangebinocean(j)=min(rangebin_indexes);
        % Peak of ocean
        rangebin100ocean(j)=find(reduced_wvf_raw_shifted_ML_norm_xrecord(:,j)==1,1,'first');
        % 87% of the peak of ocean
        a=reduced_wvf_raw_shifted_ML_norm_xrecord(1:rangebin100ocean(j),j);
        [val2,idx2]=min(abs(a-0.87));
        rangebin87ocean(j)=idx2;
    end
   
    %TESTING: BINARY IMAGE STARTING FROM THE 87% THRESHOLD ONCE WE CREATED THE PREVIOUS IMAGE
    oceanwv_threshold=0.87;
    i2new=i2;
    i3= reduced_wvf_raw_shifted_ML_norm_xrecord>oceanwv_threshold;
    for k=1:size(i3,2)
        % Find the 1st point in the binary image of the 87% threshold
        rangebinocean87threshold=find(i3(:,k),1,'first');
        for kk=rangebinocean(k):rangebinocean87threshold
            % Move the ocean WV starting area to the 87% threshold
            i2new(kk,k)=i2(rangebinocean(k),k)==0;
        end
    end
    figure;
    imagesc(i2new); 
    colormap(gray); 
    freezeColors; 
    title ('Binary image OCEANWV threshold FF ML'); 
    
    icebergs_properties_FF_ML=[];
    % Loop through all the different regions to select the ones with iceberg candidates
    for region=(1+min(min(X))):max(max(X))
        % Coordinates of the current region under analysis
        [rangebin_indexes, alongtrack_index]=find(X==region);
        % We only select the regions which are before the ocean leading edge and with region code different from the ocean/big area
        if min(rangebin_indexes) < min(rangebinocean) && region~=oceanWV_areacode
            i=i+1;
            % Save properties of positives in a struct
            % Save properties of positives in a struct
            icebergs_properties_FF_ML(i).region=region;
            icebergs_properties_FF_ML(i).alongtrack_index=alongtrack_index;
            icebergs_properties_FF_ML(i).rangebin_ind=rangebin_indexes;
            clear power_iceberg;
            clear sigma0;
            clear freeboard_iceberg;
            for ii=1:length(rangebin_indexes)
                power_iceberg(ii)=10.*log10(reduced_wvf_raw_shifted_ML(icebergs_properties_FF_ML(i).rangebin_ind(ii),icebergs_properties_FF_ML(i).alongtrack_index(ii)));
                sigma0(ii)=power_iceberg(ii)+sig0_scaling_factor_FF_ML(alongtrack_index(ii));
                %freeboard_iceberg(ii)=abs(-elev_axis_ML(rangebin87ocean(alongtrack_index(ii))) - -elev_axis_ML(rangebin_indexes(ii)));
            end          
            icebergs_properties_FF_ML(i).power=mean(power_iceberg);
            icebergs_properties_FF_ML(i).stdpower=std(power_iceberg);
            icebergs_properties_FF_ML(i).sigma0=mean(sigma0);
            icebergs_properties_FF_ML(i).stdsigma0=std(sigma0);
            %icebergs_properties_FF_ML(i).freeboard=max(freeboard_iceberg);
        else

        end
    end
    %%plot of waveforms with markers
    figure;imagesc(10.*log10(reduced_wvf_raw_shifted_DDP));
    clear i;
    for i=1:length(icebergs_properties_L1B)
        hold on; scatter([icebergs_properties_L1B(i).alongtrack_index],[icebergs_properties_L1B(i).rangebin_ind],20,'red')
    end  
end

end