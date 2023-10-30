%Compute the difference between sidelobe levles (left and right) of an 1D Power IRF or antenna
%pattern
function [side_lobes_ratio, ratio_req_met] = power_side_lobes_ratio(wfm_AC_interp, pos_max_along_focused, chd)

cut_IRF_across = (abs(wfm_AC_interp(pos_max_along_focused,:)))^2;

[pos_value,pos_max] = max(cut_IRF_across);

[value_peak,loc_peak]=findpeaks(cut_IRF_across);
if isempty(find(loc_peak==pos_max))
    %findpeaks not able to retrieve the max peak position when
    %having backfolding at edges of array
    loc_peak =sort([loc_peak;pos_max]);
    idx_pos_max = find(loc_peak==pos_max);
    idx_rest = find(loc_peak~=pos_max);
    dumm=NaN(length(loc_peak),1);
    dumm(idx_pos_max)=max_val(i_beam,i_surf);
    dumm(idx_rest)=value_peak;
    value_peak=dumm;
    clear dumm
end
idx_max = find(value_peak==max(value_peak));
if idx_max>1 && idx_max<length(value_peak)
    side_lobes_ratio = [10*log10(value_peak(idx_max-1))-10*log10(value_peak(idx_max+1))];
else
    %shift the array to avoid backfolding and compute the peaks
    %locations again
    AF_shifted=fftshift(cut_IRF_across,1);
    [value_peak,loc_peak]=findpeaks(AF_shifted);
    [~,max_pos_az_bis]=max(AF_shifted);
    if isempty(find(loc_peak==max_pos_az_bis))
        %findpeaks not able to retrieve the max peak position when
        %having backfolding at edges of array
        loc_peak =sort([loc_peak;max_pos_az_bis]);
        idx_pos_max = find(loc_peak==max_pos_az_bis);
        idx_rest = find(loc_peak~=max_pos_az_bis);
        dumm=NaN(length(loc_peak),1);
        dumm(idx_pos_max)=max_val(i_beam);
        dumm(idx_rest)=value_peak;
        value_peak=dumm;
        clear dumm
    end
    idx_max = find(value_peak==max(value_peak));
    side_lobes_ratio = [10*log10(value_peak(idx_max-1))-10*log10(value_peak(idx_max+1))];
    
end

ratio_req_met = abs(side_lobes_ratio) < chd.side_lobes_ratio;

end