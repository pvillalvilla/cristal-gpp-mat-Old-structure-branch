function [L2, validity_tfmra] = TFMRA_retracker (L1B, chd, cst, cnf_L2, L2)
[~,records]     = size(L2.scaled_waveform);
value           = zeros(1,records);
pos             = zeros(1,records);

for i=1:records
    %checking waveform valid due to L1B processing related to noise
    %filtering outer beams
    if ~any(~isnan(L2.scaled_waveform(:,i)))
        L2.range_tfmra(i) = NaN;
        L2.epoch_tfmra(i) = NaN;
        L2.sig0_tfmra(i) = NaN;
        L2.amplitude_tfmra(i) = NaN;
        validity_tfmra(i) = 0; %invalid due to L1B processing
        continue;
    end
    % Find the leading edge point with %THRESHOLD
    % 1 Find peak
    % 2 Compute noise and smooth
    % 3 Find start leading edge
    % 4 Linear Interpolate leading edge
    % 5 Find sample of the Threshold.
    
    [max_value(i),end_LE_pos(i)] = max(L2.scaled_waveform(:,i));
    above_noise= find(L2.scaled_waveform(:,i) > mean(L2.noise_floor)* cnf_L2.noise_threshold_tfmra);
    % JPLZ: create a filter to not select the start of the LE too soon if
    % thermal noise part is too noisy. We set a max distance of cnf.maxdistendLE
    % samples for the length of the LE
    if cnf_L2.limitlengthLE
        for j=1:length(above_noise)
            if above_noise(j) < end_LE_pos(i) &&  end_LE_pos(i)-above_noise(j) > cnf_L2.maxdistendLE
                above_noise_filtered(j)=0;
            else
                above_noise_filtered(j)=above_noise(j);
            end
        end
        filt_indxs=find(above_noise_filtered);
        above_noise=above_noise_filtered(filt_indxs);
    end
    
    if (isempty(above_noise))
        start_LE_pos(i) = end_LE_pos(i)-3;
    else
        start_LE_pos(i) = min(above_noise);
        if (start_LE_pos(i) == end_LE_pos(i))
            start_LE_pos(i) = start_LE_pos(i)-1;
        end
    end
    
    length_LE = length(start_LE_pos(i):end_LE_pos(i));
    LE_axis = start_LE_pos(i):end_LE_pos(i);
    
    %JPLZ: added condition to avoid 0 or negative values in the waveform index for badly behaving waveforms
    if any(LE_axis==0 | LE_axis<0)
      indxs_dump=find(LE_axis<1);
      LE_axis=LE_axis(length(indxs_dump)+1:length(LE_axis));
      if length(LE_axis)<2
          LE_axis(2)=2;
      end
    end
    
    interpolated_LE_axis = start_LE_pos(i):1/cnf_L2.oversampling_factor_tfmra:length_LE+start_LE_pos(i)-1;
    LE_interpolated = interp1(LE_axis,L2.power_waveform(i,LE_axis),interpolated_LE_axis,'linear');
    LE_interpolated_smoothed = smooth(LE_interpolated,cnf_L2.boxcar_width_tfmra);
    [~,interp_pos(i)] = min(abs(LE_interpolated_smoothed-(cnf_L2.threshold_level_tfmra*max_value(i))));
    pos(i) = interpolated_LE_axis(interp_pos(i));
    max_value(i) = LE_interpolated(interp_pos(i));
    
    T0 = 1/(chd.uso_freq_nom*chd.alt_freq_multiplier); %[s]
    L2.epoch_tfmra(i) = (pos(i)-cnf_L2.tracker_range_ref_sample)*T0; %[s]
    L2.range_tfmra(i) = L2.epoch_tfmra(i)*cst.c/2.+L1B.tracker_range_calibrated(i); %[m]
    L2.amplitude_tfmra(i) = max_value(i);
    L2.sig0_tfmra(i) = 10*log(L2.amplitude_tfmra(i)*L2.waveform_scale_factor(i))+L1B.sig0_scaling_factor(i); %[dB]
    validity_tfmra(i) = 1;

end
end
