function delta_h = compute_snow_depth_correction(PP, chd)

%% Correction computed following Isobel publication
% Lawrence, I., Tsamados, M., Stroeve, J., Armitage, T., and Ridout, A.: 
% Estimating snow depth over Arctic sea ice from calibrated dual-frequency radar freeboards, The Cryosphere Discuss., 
% https://doi.org/10.5194/tc-2018-54, in review, 2018.

switch chd.band
    case 'Ku'
        delta_h = 0.07.*PP-0.51;
        % paper Isobel exact
        delta_h = 0.06.*PP-0.46;
        % test values (physical retracker)
        delta_h = 0.02.*PP+0.23;
    case 'Ka'
        % paper Isobel exact
        delta_h = -0.16.*PP+0.76;
        % test values (physical retracker)
        delta_h = -0.10.*PP+0.90;
    otherwise 
        disp('Nor Band defined')
end
%tests
%PP=ones(1,4158);
%PP=(4.0+randn(1,4158)).*PP;

end