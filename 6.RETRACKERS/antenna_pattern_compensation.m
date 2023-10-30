function [AoA_corr] = antenna_pattern_compensation(AoA,antenna_pattern)
%figure;
%plot(-10:0.001:10,angle(antenna_pattern.el_pattern_ref_complex)*180/pi);
%plot(L1B.phase_diff_meas_ku_l1b_echo(index_L1b,:)*180/pi)

%aixo es la phase dif quan hi ha un mispointing de 0.1 deg (10 samples),
%llavors s'haura de fer l'angle d'això
%plot(unwrap(angle(antenna_pattern.el_pattern_ref_complex(1:end-mispointing/ang_res+1))*180/pi) - unwrap(angle(antenna_pattern.el_pattern_ref_complex((mispointing/ang_res+1):end))*180/pi));
%plot(-10:0.001:9.990,unwrap(angle(antenna_pattern.az_pattern_ref_complex(1:end-10))*180/pi) - unwrap(angle(antenna_pattern.az_pattern_ref_complex(11:end))*180/pi));
%plot(ant_missp_x, chd.wv_length/(2*pi*B)* (unwrap(angle(antenna_pattern.el_pattern_ref_complex(1:end-mispointing/ang_res+1))*180/pi) - unwrap(angle(antenna_pattern.el_pattern_ref_complex((mispointing/ang_res+1):end))*180/pi)) )

B=1.3300;
chd.wv_length=0.022206848740741;
mispointing=0.1; %[deg]
ang_res=mean(diff(antenna_pattern.angle_el)); %[DEG]
mispointing_direction='right';

switch mispointing_direction
    case 'left'
        ant_missp_x=[-10:0.001:9.900];
        ant_missp_y=1*(unwrap(angle(antenna_pattern.el_pattern_ref_complex(1:end-mispointing/ang_res+1))) - unwrap(angle(antenna_pattern.el_pattern_ref_complex((mispointing/ang_res+1):end))));
    case 'right'
        ant_missp_x=[-9.900:0.001:10];
        ant_missp_y=1*(unwrap(angle(antenna_pattern.el_pattern_ref_complex((mispointing/ang_res+1):end))) - unwrap(angle(antenna_pattern.el_pattern_ref_complex(1:end-mispointing/ang_res+1))));
end
%ant_missp = [ant_missp_x.' ant_missp_y];
%ant_missp_y_AoA2=chd.wv_length *(ant_missp_y./(2*pi*B))/pi*180;
ant_missp_y_AoA=asin(chd.wv_length *(ant_missp_y./(2*pi*B)).')*180/pi;
ant_missp_AoA=[ant_missp_x.' ant_missp_y_AoA.'];

for i=1:length(AoA)
    [minValue(i),closestIndex(i)]=min( abs(AoA(i)- (ant_missp_AoA(:,1))) );
    %AoA_corr(i)=AoA(i)-ant_missp_AoA(closestIndex(i),2); % correct antenna misspointing
    AoA_corr(i)=AoA(i); % don't correct 
end
