% SPECTRAL WINDOW COMPUTATION: ASSUMING A ZERO- FREQUENCY-CENTERED WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
% JasonCS 
% This code implements the computation of the spectral window for a zero-frequency 
% centered input array
%
% ---------------------------------------------------------
% Objective: Perform the SECONDARY RANGE COMPRESSION
% 
% INPUTs: 
%
%     MANDATORY
%
%     Name                    Description  
%     win_type:              string: indicating window type 
%                                    ('Hanning','Hamming','Kaiser','Blackman','Boxcar')
%     BW_to_process:         Double: indicating the bandwidth in Hz to be processed
%     f_vec:                 Double vector: indicating the vector of input
%                                           frequencies (spectral representation)
%     OPTIONAL
%
%     Name                    Description  
%     win_param:             Double: Value of the window parameter for the
%                                    Kaiser window for instance (or any future window)
%     normalize:             flag: 0: no normalization, 1: normalization of
%                                  weights of the window
%     normalize_type:        flag: type of normalziation to be used 
%                                  0: energy normalization
%                                  1: amplitude normalization
%     
% OUTPUTs:  
%     Name                    Description
%     win                     Output generated spectral weights of the
%                             window
%     
% ----------------------------------------------------------
% Author:    Eduard Makhoul  / isardSAT
%
% Reviewer:  Monica Roca   / isardSAT
% Last rev.: Monica Roca   / isardSAT (14/11/2017)

function [win]...
      = window_computation (win_type,BW_to_process,f_vec,f_centroid,f_s,...
                                  varargin)
%% Global/initial/optional variables
%optional variable inputs
p = inputParser;
p.addParamValue('normalize',0); 
p.addParamValue('normalize_type',0); 
p.addParamValue('win_param',0.5); 
p.parse(varargin{:});
normalize=p.Results.normalize;
normalize_type=p.Results.normalize_type;
win_param=p.Results.win_param;

N_samples = length(f_vec);

win = zeros(1,N_samples);

%% 1. COMPUTATION OF VALID INDEXES FOR THE PROCESSING BANDWIDTH
idx_PBW = find(f_vec>= -BW_to_process/2 & f_vec<= BW_to_process/2);
N_samples_PBW = length(idx_PBW);

%% 2. COMPUTATION OF THE WINDOW WEIGHTS 
switch lower(win_type)
    case {'boxcar'}
        win(idx_PBW) = ones(1,N_samples_PBW);
    case {'hamming'}
        win(idx_PBW) = hamming(N_samples_PBW);
    case {'hanning'}
        win(idx_PBW) = hann(N_samples_PBW);
    case {'kaiser'}  
        win(idx_PBW) = kaiser(N_samples_PBW,win_param);
    case {'blackman'}
        win(idx_PBW) = blackman(N_samples_PBW,win_param);
end

%% 3. Center to the corresponding centroid
win = circshift(win,[0,round(N_samples/2+f_centroid*N_samples/f_s)]);


%% 4. NORMALIZE THE OUTPUT WINDOW
if normalize 
    switch normalize_type
        case 0
            % ENERGY NORMALIZATION
            win = win./(sqrt(mean(win.^2)));
        case 1            
            % AMPLITUDE NORMALIZATION
            win = win./((mean(win)));
    end
end




end