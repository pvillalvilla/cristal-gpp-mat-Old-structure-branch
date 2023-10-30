% Author: Eduard Makhoul
%to compute the beamwidth of antenna pattern or PTR or IRF with a single
%main lobe
function [hpbw,angle_l,angle_h]=cc_ldb(AF,phi,varargin)
p = inputParser;
p.addParamValue('l',3.0); %cut at 3dB
p.addParamValue('interp_method','linear'); %cut at 3dB
p.parse(varargin{:});
interp_method = p.Results.interp_method;
l=(p.Results.l);
clear p;

[max_value,indx_max]=max(AF);
indx_max=indx_max(1);

indx_h =find(AF >=(max_value-l), 1, 'last');
indx_l =find(AF >=(max_value-l), 1, 'first');

if indx_l>1 && indx_h<length(AF)
    indx_l=indx_l-1;
    indx_h=indx_h+1;
    
    angle_l=interp1(AF(indx_l:indx_max),phi(indx_l:indx_max),max_value-l,interp_method);
    angle_h=interp1(AF(indx_max:indx_h),phi(indx_max:indx_h),max_value-l,interp_method);

%     angle_l = indx_l;
%     angle_h = indx_h;
    hpbw=angle_h-angle_l;
    
else
    AF_shifted = fftshift(AF,1);
    [max_value,indx_max]=max(AF_shifted);
    indx_max=indx_max(1);
    indx_h =find(AF_shifted >=(max_value-l), 1, 'last');
    indx_l =find(AF_shifted >=(max_value-l), 1, 'first');
    
    if indx_l >1
        indx_l=indx_l-1;        
    end
    
    if indx_h<length(AF_shifted)
        indx_h=indx_h+1;
    end
    
    angle_l=interp1(AF_shifted(indx_l:indx_max),phi(indx_l:indx_max),max_value-l,interp_method);
    angle_h=interp1(AF_shifted(indx_max:indx_h),phi(indx_max:indx_h),max_value-l,interp_method);
    hpbw=angle_h-angle_l;
    
%     %detect potential backfolding of beamwidth at begining end array
%     [value_peaks,loc_peaks]=findpeaks(AF);
%     if isempty(find(loc_peaks==indx_max))
%         %findpeaks not able to retrieve the max peak position when
%         %having backfolding at edges of array
%         loc_peaks =sort([loc_peaks;indx_max]);
%         idx_pos_max = find(loc_peaks==indx_max);
%         idx_rest = find(loc_peaks~=indx_max);
%         dumm=NaN(length(loc_peaks),1);
%         dumm(idx_pos_max)=max_value;
%         dumm(idx_rest)=value_peaks;
%         value_peaks=dumm;
%         clear dumm
%     else
%         idx_pos_max = find(loc_peaks==indx_max);
%     end
%     
%     if idx_pos_max==1
%         %part of the main-lobe is folded in the end of the array
%         %left side-lobe at the end of the array
%         peak_max_idx = find(value_peaks==max_value);
%         indx_h =find(AF(1:loc_peaks(peak_max_idx+1)) >=(max_value-l), 1, 'last');
%         indx_l =loc_peaks(end)-1+find(AF(loc_peaks(end):end)>=(max_value-l), 1, 'first');
%         
%         indx_l=indx_l-1;
%         indx_h=indx_h+1;
%         
%         angle_l=interp1(AF(indx_l:end),phi(indx_l:end)-phi(end),max_value-l,'spline');
%         angle_h=interp1(AF(indx_max:indx_h),phi(indx_max:indx_h),max_value-l,'spline');
%         hpbw=angle_h-angle_l;
%     elseif idx_pos_max==length(loc_peaks)
%         %part of the main-lobe is folded in the beginning of the array
%         %right sidelobe is folded back at the beginning of the array
%         peak_max_idx = find(value_peaks==max_value);
%         indx_l =loc_peaks(peak_max_idx-1)+...
%             find(AF(loc_peaks(peak_max_idx-1):end) >=(max_value-l), 1, 'first');
%         indx_h =find(AF(1:loc_peaks(1))>=(max_value-l), 1, 'last');
%         
%         indx_l=indx_l-1;
%         indx_h=indx_h+1;
%         
%         angle_l=interp1(AF(indx_l:indx_max),phi(indx_l:indx_max)-phi(indx_max),max_value-l,'spline');
%         angle_h=interp1(AF(1:indx_h),phi(1:indx_h),max_value-l,'spline');
%         hpbw=angle_h-angle_l;
%         
%     end
    
end





end