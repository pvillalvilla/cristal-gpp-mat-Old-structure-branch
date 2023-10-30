function [burst_air_pslr, burst_air_pslr_met] = burst_air_PSLR(wfm_AC_interp, pos_max_along_focused, chd)
    cut_IRF_across = abs(wfm_AC_interp(pos_max_along_focused,:));
    [burst_air_pslr] = psl(cut_IRF_across.^2);
    burst_air_pslr_met = burst_air_pslr < chd.burst_air_PSLR;
end