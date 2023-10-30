% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% This code generates the characterization file input to isardSATs
% re-tracker as described in isardSAT_Retracker_DPM_V1.a.docx internal
% document
%
% ---------------------------------------------------------
% gen_l2_chd_file: function that generates the configuration file input to
% isardSAT's retracker
%
% Calling
%   gen_l2_chd_file( Mission )
%
% Inputs
%   configuration parameters 
%
% Output
%   configuration file:     
%   FileName = chd_file_'Mission'_L2.bin
% ----------------------------------------------------------
% 
% Author:   Eduard Makhoul / isardSAT
%
%
%
%
% This software is built with internal funding
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c_cst
global freq_ku_chd wv_length_ku prf_chd bw_rx_ku_chd fs_clock_ku_chd brf_chd
global N_total_pulses_b_chd N_samples_chd N_bursts_cycle_chd N_samples_rmc_chd
global antenna_beamwidth_act_ku_chd antenna_beamwidth_alt_ku_chd antenna_gain_ku_chd
global roll_bias_chd pitch_bias_chd
global alpha_ga_chd alpha_gr_chd A_s2Gr_chd A_s2Ga_chd
global PTR_width_chd power_tx_ant_ku_chd
 switch cnf_p.mission
    
    case 'CS2'
        % -----------------------------------------------------------------
        % OPERATION CONFIGURATION
        % -----------------------------------------------------------------
        % Total number of pulses in a burst this is equal for SAR and SARin
        % mode (uint8)
        N_total_pulses_b_chd    =   64;
        % Number of samples per each SAR or SARin pulse (uint8)
        N_samples_chd       =   128;
        %Number of burst in a radar cycle
        switch cnf_p.mode
            case 'SAR'
                % Number of burst in a SARM tracking cycle ('uint8')
                N_bursts_cycle_chd  =   4;
                
            case 'SARin'
                % Number of burst in a SARM tracking cycle ('uint8')
                N_bursts_cycle_chd  =   1;
                
            otherwise
                error('not a valid operational mode')
        end
        
        % bias in pitch and roll (BE CAREFUL WITH SUCH PARAMETERS: IF WE READ DIRECTLY L1B ESA OR L1B ISD GENERATED WITH FBR BASELINE-C SUCH BIAS ARE ALREADY CORRECTED)
        % Pitch bias from NOAA on 17th May 2013
        pitch_bias_chd                  =   0.096*pi/180;
        
        % Roll bias
        roll_bias_chd                   =   0.086*pi/180;
        
        %PTR width 
        PTR_width_chd=2.819e-9;
        
        %TX power
        power_tx_ant_ku_chd = 20;
        switch cnf_p.L1proc
            case {'ESA','ISD'}% processed by ESA or isardSAT
                % -----------------------------------------------------------------
                % INSTRUMENT PARAMETERS
                % -----------------------------------------------------------------
                % Central frequency in Ku band (float64)
                
                freq_ku_chd         =   13.575e9; %from SCOOP POCCD
                wv_length_ku        =   c_cst/freq_ku_chd;
                
                % PRF FOR SAR and SARin modes (float64)
                prf_chd             =   1/(55e-6);
                brf_chd             =   1/(11.7929625 * 1e-3);
                
                
                % Ku received Bandwidth SAR mode (float64)
                bw_rx_ku_chd        =   320e6 ;
                fs_clock_ku_chd     =   320e6; % previous 395MHz 1/1.25e-8   
                
                % -----------------------------------------------------------------
                % ANTENNA Characterization information
                % -----------------------------------------------------------------
                %antenna_beamwidth_act_ku_chd    =   (1.1992*pi/180);      % from ESA according to Cristina
                %antenna_beamwidth_alt_ku_chd    =   (1.06*pi/180);      % from ESA according to Cristina
                antenna_beamwidth_act_ku_chd    =   (1.22*pi/180);      % from SCOOP POCCD
                antenna_beamwidth_alt_ku_chd    =   (1.095*pi/180);      % from SCOOP POCCD
                antenna_gain_ku_chd             = 42.6; %dB
        
            otherwise %CNES
                % -----------------------------------------------------------------
                % INSTRUMENT PARAMETERS
                % -----------------------------------------------------------------
                % Central frequency in Ku band (float64)
                freq_ku_chd         =   13.575e9;
                
                % PRF FOR SAR and SARin modes (float64)
                prf_chd             =   18181.818;
                
                
                % Ku received Bandwidth SAR mode (float64)
                bw_rx_ku_chd        =   320e6;
                fs_clock_ku_chd     =   bw_rx_ku_chd; % previous 395MHz 1/1.25e-8
                
                % -----------------------------------------------------------------
                % ANTENNA Characterization information
                % -----------------------------------------------------------------
                antenna_beamwidth_act_ku_chd    =   (1.16*pi/180);     % CNES as in V13
                antenna_beamwidth_alt_ku_chd    =   (1.095*pi/180);    %CNES
                antenna_gain_ku_chd             = 42.6; %dB
        end
        
     case 'S3'
            %TBD
         
     case 'S6'

            % -----------------------------------------------------------------
            % INSTRUMENT PARAMETERS
            % -----------------------------------------------------------------
            % Central frequency in Ku band (float64)
            wv_length_ku            =   0.022084;
            freq_ku_chd             =   c_cst/wv_length_ku; 
            

            % PRF FOR SAR and SARin modes (float64)
           % prf_chd                 =   1/(0.108962025*1e-3); % from where??
            prf_chd                  =  9178; %mean value according to P4 model issue 4.2

            % Ku received Bandwidth SAR mode (float64)
            bw_rx_ku_chd            =   320e6; 
            fs_clock_ku_chd         =   395e6; % previous 395MHz 1/1.25e-8

            
            % Total number of pulses in a burst this is equal for SAR and SARin
            % mode (uint8)
            N_total_pulses_b_chd    =   64;
           

            % Number of samples per each SAR or SARin pulse (uint8)
            switch cnf_p.mode
                case 'RMC'
                    N_samples_chd           =   128;
                    N_samples_rmc_chd       =   128;
                    
                case {'HR','SAR'}
                    N_samples_chd           =   256;                    
            end
            
            % Number of burst in a SARM tracking cycle ('uint8')
            N_bursts_cycle_chd      =   7;
            
            % -----------------------------------------------------------------
            % ANTENNA Characterization information
            % -----------------------------------------------------------------
            antenna_beamwidth_act_ku_chd    =   (1.35*pi/180);      % NOAA
            antenna_beamwidth_alt_ku_chd    =   (1.35*pi/180);      % NOAA


            % Pitch bias from NOAA on 17th May 2013
            pitch_bias_chd                  =   0;              
            
            % Roll bias
            roll_bias_chd                   =   0;              
           
    
    otherwise
        error('Not a valid Mission. Chose between: CS2 or JCS');
 end % end switch
 
% -----------------------------------------------------------------
% PTR Characterization information (sinc to Gaussian approx.) this is
% independent of the mission, but it is a sinc to Gaussian issue
% -----------------------------------------------------------------
% EM: need to be reviewed
if any(strcmp(cnf_p.retracker_name,'ANALYTICAL') | strcmp(cnf_p.retracker_name,'SAMOSA'))
    switch cnf_p.window_type_a
        case 'Boxcar'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
%             alpha_chd     =   1;
%             [A_s2Ga_chd, alpha_ga_chd]    =   sinc2Gauss(alpha_chd, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
%             alpha_ga_chd=1/(2*0.513^2);          
            A_s2Ga_chd=1.0196;
            alpha_ga_chd=1.0/(2*(0.36012).^2);
        case 'Hanning'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
%             alpha_chd     =   0.5;
%             [A_s2Ga_chd, alpha_ga_chd]    =   sinc2Gauss(alpha_chd, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
            A_s2Ga_chd=1.0101;
            alpha_ga_chd=1.0/(2*(0.59824).^2);
        case 'Hamming' % the Gaussian approximation for a hamming window filtered waveform previous to FFT
%             alpha_chd     =   0.54;
%             [A_s2Ga_chd, alpha_ga_chd]    =   sinc2Gauss(alpha_chd, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
            A_s2Ga_chd=1.0081;
            alpha_ga_chd=1.0/(2.0*(0.54351).^2);
        case 'forced'			
			alpha_ga_chd=1.0/(2.0*(0.545333774157098).^2);            
        otherwise
            error('invalid window type. Chose Hamming, Hanning or Boxcar')
    end
    disp(strcat('Azimuth PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_ga_chd)))))
    
    switch cnf_p.window_type_r
        case 'Boxcar'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
%             alpha_chd     =   1;
%             [A_s2Gr_chd, alpha_gr_chd]    =   sinc2Gauss(alpha_chd, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
%             alpha_gr_chd =1/(2*0.513^2);            
            A_s2Gr_chd=1.0196;
            alpha_gr_chd=1.0/(2*(0.36012).^2);
        case 'Hanning'  % the Gaussian approximation for a boxcar window filtered waveform previous to FFT
            alpha_chd     =   0.5;
            [A_s2Gr_chd, alpha_gr_chd]    =   sinc2Gauss(alpha_chd, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
            A_s2Gr_chd=1.0101;
            alpha_gr_chd=1.0/(2*(0.59824).^2);
        case 'Hamming' % the Gaussian approximation for a hamming window filtered waveform previous to FFT
            alpha_chd     =   0.54;
            [A_s2Gr_chd, alpha_gr_chd]    =   sinc2Gauss(alpha_chd, [-2+1/64:4/64:2]);    % note that alpha_g varies with the number of samples not the rank of values
            A_s2Gr_chd=1.0081;
            alpha_gr_chd=1.0/(2.0*(0.54351).^2);
		case 'forced'
			alpha_gr_chd=1.0/(2.0*(0.545333774157098).^2);            
        otherwise
            error('invalid window type. Chose Hamming, Hanning or Boxcar');
    end
    disp(strcat('Range PTR approx:',{''},num2str(sqrt(1.0/(2.0*alpha_gr_chd)))))
end









