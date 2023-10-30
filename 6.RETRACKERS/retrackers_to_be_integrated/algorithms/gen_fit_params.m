function ini_p = gen_fit_params (data,cnf_L2) 
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This function initializates fitting parameters structure
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Albert Garcia / isardSAT v1.1 10/06/2019
%
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       data    =   input data structure for the L2 processor
%       cnf_L2           =   configuration parameters
% OUTPUT:
%       ini_p           =   Fitting parameters seed

% -------------------------------------------------------------------------

% v1.1 changed cnf struct

    % initialitzate fitting parameters    
    if ~cnf_L2.seed
        ini_p.Epoch     =   cnf_L2.ini_Epoch;
    else
        ini_p.Epoch     =   data.MEA.seed(1);    
    end
    
    ini_p.Hs        =   cnf_L2.ini_Hs;
    ini_p.sigmaz    =   ini_p.Hs/4;
    if (strcmp(cnf_L2.multilook_option,'Cris_old'))
        ini_p.Pu        =   cnf_L2.ini_Pu;  
    end
    ini_p.rou       =   cnf_L2.init_rou;
%     %added by EM 16.03.2016 include optimization of the noise thermal
%     %estimation window
%     if cnf_L2.fit_noise
%        ini_p.Thn_1st=cnf_L2.ini_Thn_1st;
%        ini_p.Thn_wdth=cnf_L2.ini_Thn_wdth;
%     end
    
%     % if LUT initializate LUT tables
%     if (cnf_L2.lut_flag)
%         % Look up tables initialization 
%         % Generate LUT configuration file
%         switch cnf_L2.pdf_surf
%             case 'Gaussian'
%                 order = 1; % this generates up to f_order(xi)
%                 gen_LUT(order,[cnf_L2.LUT_ximin:((cnf_L2.LUT_ximax-cnf_L2.LUT_ximin)/(cnf_L2.LUT_NLUT)):cnf_L2.LUT_ximax],cnf_L2);
%             case 'non-Gaussian'
%                 print('WARNING: Code not optimized for non-Gaussian')
%                 order = 4; % this generates up to f_order(xi)
%                 gen_LUT(order,[cnf_L2.LUT_ximin:((cnf_L2.LUT_ximax-cnf_L2.LUT_ximin)/(cnf_L2.LUT_NLUT)):cnf_L2.LUT_ximax],cnf_L2);
%             otherwise
%                 error('Not a valid surface pdf')
%         end
%     end

        
end