function stack_phase_diff = compute_phase_difference (L1BS, beams_zp_fft,beams_zp_fft_2, maxpos_A1, maxpos_A2)
interpolation_margin= 1:3;
interpolation_step= 0.01;
interpolation_length=length(interpolation_margin(1):interpolation_step:interpolation_margin(3));
diff_phase_12=zeros(1,L1BS.N_beams_stack);
% angle_meas=zeros(1,L1BS.N_beams_stack);

A1_Re   = real(beams_zp_fft);
A2_Re   = real(beams_zp_fft_2);
A1_Img  = imag(beams_zp_fft);
A2_Img  = imag(beams_zp_fft_2);

Re1     = zeros(L1BS.N_beams_stack,interpolation_length);
Re11    = zeros(1,L1BS.N_beams_stack);
Img1    = zeros(L1BS.N_beams_stack,interpolation_length);
Img11   = zeros(1,L1BS.N_beams_stack);
Re2     = zeros(L1BS.N_beams_stack,interpolation_length);
Re22    = zeros(1,L1BS.N_beams_stack);
Img2    = zeros(L1BS.N_beams_stack,interpolation_length);
Img22   = zeros(1,L1BS.N_beams_stack);

Phase1  = zeros(1,L1BS.N_beams_stack);
Phase2  = zeros(1,L1BS.N_beams_stack);
stack_phase_diff = zeros(1,L1BS.N_beams_stack);

for nwvf = 1: L1BS.N_beams_stack
    Re1(nwvf,:)  = interp1(interpolation_margin, A1_Re(nwvf, round(maxpos_A1(nwvf))-1:round(maxpos_A1(nwvf))+1),interpolation_margin(1):interpolation_step:interpolation_margin(3) , 'spline');
    Re11(nwvf)   =  Re1(nwvf,round((maxpos_A1(nwvf)-(round(maxpos_A1(nwvf))-1))/interpolation_step+1));
    Img1(nwvf,:) = interp1(interpolation_margin, A1_Img(nwvf, round(maxpos_A1(nwvf))-1:round(maxpos_A1(nwvf))+1),interpolation_margin(1):interpolation_step:interpolation_margin(3) , 'spline');
    Img11(nwvf)  =  Img1(nwvf,round((maxpos_A1(nwvf)-(round(maxpos_A1(nwvf))-1))/interpolation_step+1));
    Re2(nwvf,:)  = interp1(interpolation_margin, A2_Re(nwvf, round(maxpos_A2(nwvf))-1:round(maxpos_A2(nwvf))+1),interpolation_margin(1):interpolation_step:interpolation_margin(3) , 'spline');
    Re22(nwvf)   =  Re2(nwvf,round((maxpos_A2(nwvf)-(round(maxpos_A2(nwvf))-1))/interpolation_step +1 ));
    Img2(nwvf,:) = interp1(interpolation_margin, A2_Img(nwvf, round(maxpos_A2(nwvf))-1:round(maxpos_A2(nwvf))+1),interpolation_margin(1):interpolation_step:interpolation_margin(3) , 'spline');
    Img22(nwvf)  =  Img2(nwvf,round((maxpos_A2(nwvf)-(round(maxpos_A2(nwvf))-1))/interpolation_step+1));
    
    %Phase difference version 2
    Phase1(nwvf)=atan(Img11(nwvf)/Re11(nwvf));
    Phase2(nwvf)=atan(Img22(nwvf)/Re22(nwvf));
    
    %Unwrapp the tangent wrapps on [pi/2:-pi/2]
%     if((Re11(nwvf)<0))
%         Phase1(nwvf)=Phase1(nwvf)+pi;
%         %                 Phase1_full(nwvf)=Phase1_full(nwvf)+pi;
%     end
    
%     if(Re22(nwvf)<0)
%         Phase2(nwvf)=Phase2(nwvf)+pi;
%         %                 Phase2_full(nwvf)=Phase2_full(nwvf)+pi;
%     end
    stack_phase_diff(nwvf)= Phase1(nwvf)-Phase2(nwvf);
    if(diff_phase_12(nwvf)>pi)
        stack_phase_diff(nwvf)=stack_phase_diff(nwvf)-2*pi;
    elseif(diff_phase_12(nwvf)<-pi)
        stack_phase_diff(nwvf)=stack_phase_diff(nwvf)+2*pi;
    end
    
end




end
