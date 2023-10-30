cnf.scenario = 'TLI'; % TPT (Point Target) 
                   % TLI (Land Ice) 
                   % TSI (Sea Ice) 
                   % TOO (Open Ocean)
                   % TIW (Inland Waters) 
cnf.id = 1; 
cnf.input_product = 'L1B HR'; % L1B HR, L1B-S HR, L1B SL, L1B LR, GR HR, 
                              % GR Swath, GR LR, GR FF, GR LR OS, GR LR RMC
cnf.mode = 'SARin OB'; % SARin OB, SARin CB, SAR OB, SAR CB
% 
% switch cnf.scenario
%     case 'TPT'
%         switch cnf.input_product
%             case 'L1B HR'
%                 cnf.azimuth_processing = 'exact';
%                 cnf.azimuth_windowing = 0; % switched off
%                 cnf.antenna_weighting = 0; % switched off
%                 cnf.range_ZP = 512;
%                 cnf.res_step = 10^-5;
%             case 'L1B 
cnf.mode_improved_approximate = 0;