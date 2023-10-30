function [mean_noise,std_floor] = compute_noise_floor (scaled_waveforms, cnf)


mean_noise = mean(scaled_waveforms(cnf.range_index1:cnf.range_index2,:),1);
std_floor  = std(scaled_waveforms(cnf.range_index1:cnf.range_index2,:),1);

end


