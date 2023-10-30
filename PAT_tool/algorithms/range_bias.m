function [range_bias, range_bias_req_met] = range_bias(tracker_range_calibrated, ...
    slant_range_correction_applied, power_waveform, Gap_flag, range_PT_expected, required_range_bias, left_beams, right_beams, ...
    point_target_idx, chd, cst)
    range_without_slant = tracker_range_calibrated(point_target_idx) - ...
        (slant_range_correction_applied(point_target_idx,:).')*chd.T0_nom*cst.c/2;

    pos_max = max(power_waveform(point_target_idx,:));
    if(pos_max<50)
            left_beams=pos_max-1;
    else
    valid_beams   = find(Gap_flag(pos_max-left_beams:pos_max+right_beams));
    range_bias = abs(mean(range_without_slant(valid_beams))-range_PT_expected);
    range_bias_req_met = range_bias < required_range_bias;
end