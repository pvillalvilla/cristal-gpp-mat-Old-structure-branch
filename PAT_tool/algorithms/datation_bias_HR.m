function [datation_error_value, datation_error_met] = datation_bias_HR(tracker_range_calibrated, ...
    slant_range_correction_applied, point_target_idx, chd, cst)

theo_range = chd.expected_range;
range_without_slant = tracker_range_calibrated(point_target_idx) - ...
    (slant_range_correction_applied(point_target_idx,:).')*chd.T0_nom*cst.c/2;


for i_beam=1:length(range_without_slant)
    time_beam(i_beam)= L1A(L1BS(point_target_idx).burst_index(i_beam)).time_sar_ku;
end

% time interpolation:
interp_fact=0.001;
delta_time=mean(diff(time_beam));
delta_time_interp= delta_time*interp_fact;
times_STACK = min(time_beam):delta_time_interp:max(time_beam);

[range_theo_coef,SS,MUMU] = polyfit(time_beam,theo_range,2);
range_theo_interp = polyval(range_theo_coef, times_STACK,SS,MUMU);
[range_meas_coef,SS,MUMU] = polyfit(time_beam(valid_beams),range_without_slant(valid_beams),2);
range_meas_interp = polyval(range_meas_coef, times_STACK,SS,MUMU); 

if cnf.RMSE_minimisation_flag % datation_error_A1_v2
    deltaR = fitfunctions_yaxis (range_meas_interp, range_theo_interp);
    deltaT = -1*fitfunctions_xaxis ((range_meas_interp-deltaR), range_theo_interp); %--20110711-- reduce Y diff
    if deltaT==0
      deltaT = fitfunctions_xaxis (range_theo_interp, range_meas_interp-deltaR);
    end   
    datation_error_value = deltaT*delta_time_interp *1e6;%[micros]

else  % datation_error_A1_v1
    [~,pos_theo]=min(range_theo_interp);
    [~,pos_meas]=min(range_meas_interp);
    datation_error_value = (times_STACK(pos_meas)-times_STACK(pos_theo))*1e6; %[micros]
end

datation_error_met = datation_error_value < chd.datation_bias;

end