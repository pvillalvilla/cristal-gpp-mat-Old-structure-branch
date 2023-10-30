% v1.0 First Version
% v1.1 2019/06/10 Added noise threshold and case when everything is below noise
function [value,pos] = threshold_first_maximum_restracker(scaled_waveforms,threshold,cnf)

% Find the leading edge point with %THRESHOLD 

% 1 Find peak
% 2 Compute noise and smooth
% 3 Find start leadind edge
% 4 Linear Interpolate leading edge
interpolation_step = 5;
% 5 Find sample of the Threshold.
[~,records]     = size(scaled_waveforms);
value           = zeros(1,records);
pos             = zeros(1,records);

[max_value,end_LE_pos] = peak_retracker (scaled_waveforms);
[mean_noise,~]     = compute_noise_floor(scaled_waveforms,cnf);
noise_threshold= 1.2;
for i=1:records
    
    above_noise= find(scaled_waveforms(:,i) > mean(mean_noise)*noise_threshold);
    % JPLZ: create a filter to not select the start of the LE too soon if
    % thermal noise part is too noisy. We set a max distance of cnf.maxdistendLE
    % samples for the length of the LE
%     if cnf.limitlengthLE
%         for j=1:length(above_noise)
%             if above_noise(j) < end_LE_pos(i) &&  end_LE_pos(i)-above_noise(j) > cnf.maxdistendLE
%                 above_noise_filtered(j)=0;
%             else
%                 above_noise_filtered(j)=above_noise(j);
%             end
%         end
%         filt_indxs=find(above_noise_filtered);
%         above_noise=above_noise_filtered(filt_indxs);
%     end
    
    if (isempty(above_noise))
        start_LE_pos(i) = end_LE_pos(i)-3;
    else
        start_LE_pos(i)   = min(above_noise);
        if(start_LE_pos(i) == end_LE_pos(i))
            start_LE_pos(i) = start_LE_pos(i)-1;          
        end
    end
    length_LE           = length(start_LE_pos(i):end_LE_pos(i));
    LE_axis             = start_LE_pos(i):end_LE_pos(i);
    %JPLZ: added condition to avoid 0 or negative values in the waveform index for badly behaving waveforms
    if any(LE_axis==0 | LE_axis<0)
      indxs_dump=find(LE_axis<1);
      LE_axis=LE_axis(length(indxs_dump)+1:length(LE_axis));
      if length(LE_axis)<2
          LE_axis(2)=2;
      end
    end
    inperpolated_LE_axis = start_LE_pos(i):1/interpolation_step:length_LE+start_LE_pos(i)-1;
    LE_interpolated     = interp1(LE_axis,scaled_waveforms(LE_axis,i),inperpolated_LE_axis,'linear');
    LE_interpolated_smoothed     = smooth(LE_interpolated,15);
    
    [~,interp_pos(i)]   = min(abs(LE_interpolated_smoothed-(threshold*max_value(i))));
    pos(i)              = inperpolated_LE_axis(interp_pos(i));
    value(i)            = LE_interpolated(interp_pos(i));
    clear LE_interpolated inperpolated_LE_axis LE_interpolated_smoothed;
end


end