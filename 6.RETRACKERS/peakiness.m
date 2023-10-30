function PP = peakiness(scaled_waveforms,cnf,chd)


N_samples	= size(scaled_waveforms,1);
N_records   = size(scaled_waveforms,2);
[mean_noise,~]     = compute_noise_floor(scaled_waveforms,cnf);
[max_value,end_LE_pos] = peak_retracker (scaled_waveforms);
noise_threshold= 1.5;

for i_record=1:N_records
    valid_samples = scaled_waveforms(:,i_record)>mean_noise(i_record)*noise_threshold;
    
%     if cnf.limitlengthLE
%         % JPLZ: create a filter to not select the start of the LE too soon if
%         % thermal noise part is too noisy. We set a max distance of
%         % samples for the length of the LE
%         %valid_samples_indxs=find(valid_samples);
%         for j=1:length(valid_samples)
%             if valid_samples(j)==1
%                valid_samples_filter(j)=j;
%             else
%                valid_samples_filter(j)=0; 
%             end
%             
%             if valid_samples_filter(j)~=0 && valid_samples_filter(j) < end_LE_pos(i_record) &&  end_LE_pos(i_record)-valid_samples_filter(j) > cnf.maxdistendLE
%                 valid_samples_filtered(j)=0;
%             else
%                 valid_samples_filtered(j)=valid_samples_filter(j);
%             end
%         end
%         valid_samples=ismember(1:numel(valid_samples),valid_samples_filtered);
%     end
    
    PP(i_record)    =  max(scaled_waveforms(valid_samples,i_record))/mean(scaled_waveforms(valid_samples,i_record));
    %max(scaled_waveforms(:,i_record))/sum(scaled_waveforms(:,i_record)); % JPLZ: PP based on DiBella2021 paper
   
end
% if cnf.filterPP
%     % Filter waveforms by PP, according to  Lawrence, I. 2018. The snow depth corrections taken only apply for PP between 3.5-4.6 (Ka band)
%     % and PP 4.5-7.5 (Ku) band. Waveforms with greater PP are considered ambiguous regarding snow depth
%     switch chd.band
%         case 'Ku'
%             PP_filtered_log = PP > 4.5 & PP < 7.5;
%             indxs_PP_filtered=find(PP_filtered_log);
%             PP_filtered=PP.*PP_filtered_log;
%             PP(PP_filtered==0)=NaN;
%         case 'Ka'
%             PP_filtered_log = PP > 3.5 & PP < 4.6;
%             indxs_PP_filtered=find(PP_filtered_log);
%             PP_filtered=PP.*PP_filtered_log;
%             PP(PP_filtered==0)=NaN;
%         otherwise
%             disp('Nor Band defined')
%     end
% end
end