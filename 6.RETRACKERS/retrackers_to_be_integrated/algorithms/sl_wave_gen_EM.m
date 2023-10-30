function sl_waveform = sl_wave_gen_EM(x, l, fit_p, nf_p, cnf_p,func_f0,func_f1)%varargin)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% Routine to generate SAMOSA 3 addapted. The code is optimized for the
% fitting of epoch, sigmaz(Hs) and Pu
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Cristina Martin-Puig / isardSAT
%
% Last revision:    Cristina Martin-Puig / isardSAT V9 04/07/2014
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%       x           =   Range bin index 
%       l           =   Look index
%       fit_p       =   Parameters to be fitted 
%       nf_p        =   parameters of need for the waveform generation, but not to be
%                       fitted
%       cnf_p       =   configuration parameters of L2 processor
% OUTPUT:
%       sl_waveform =   single look waveform for a given look and all range
%                       bins
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% CALLED FUNCTIONS/ROUTINES
% - sl_wave_gen: in charge of genrating the single look waveform based on
% the model proposed by Chris in IEEE TGRS "SAR Altimeter Backscattered Waveform Model"
% DOI:10.1109/TGRS.2014.23330423
% -
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% COMMENTS/RESTRICTIONS
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Versions control:
% Based on the V9 version defined by Cristina Martin-Puig and updated
% according to revisit of the different expressions for antenna and
% radiation patterns and the f0 and f1 functions (first and second order power waveform model), 
% inclusion of the option to use external values based on LUTs for the f0 and f1: 

    %**********************************************************************
    %*************** Input optional parameters ****************************
    %**********************************************************************
%     p = inputParser;
%     defaultff0=0.0;
%     addOptional(p,'func_f0',defaultff0,@isnumeric);
%     addOptional(p,'func_f1',defaultff0,@isnumeric);
%     parse(p,varargin{:});
%     func_f0=p.Results.func_f0;
%     func_f1=p.Results.func_f1;
%     clear p
    
    if nargin < 7
        func_f1=0; % initialize altitude to 0;
    end
    if nargin <6
        func_f0=0; % initialize altitude to 0;
    end
    
    epoch=fit_p(1);
 
    if cnf_p.rou_flag
        sigmaz=0/4; rou=fit_p(2);
    else
        rou=nf_p.rou; sigmaz=fit_p(2);
    end
    
    
    %**********************************************************************
    %% ******************** Samosa3 ***************************************
    %**********************************************************************
    
    %------------dilation term---------------------------------------------
    g       =   sqrt(2*nf_p.alphag_a*nf_p.alphag_r/(nf_p.alphag_a + nf_p.alphag_r * 4 * (nf_p.Lx/nf_p.Ly)^4 * l^2 + 2 * nf_p.alphag_a * nf_p.alphag_r * (sigmaz/nf_p.Lz)^2 ));    
    
    %--------------Varibales  ---------------------------------------------
    k=x - epoch; xl=nf_p.Lx*l; alpha_sigma=1/(nf_p.h^2*rou); yk=nf_p.Ly.*abs(sqrt(k)); gk=g.*k;
    %            
    %**********************************************************************
    %************************ Antenna & Surface ***************************
    %**********************************************************************
    %Constant Term
    %Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2).*exp(-alpha_sigma * xl.^2).*exp(-nf_p.alphay * nf_p.yp^2).*exp(-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);
    Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2-alpha_sigma * xl.^2-nf_p.alphay * nf_p.yp^2-(nf_p.alphay + alpha_sigma).*(yk).^2).*cosh(2*nf_p.alphay*nf_p.yp*yk);
    %Bkl=2.0*exp(-nf_p.alphax *(xl-nf_p.xp).^2-alpha_sigma * xl.^2-(alpha_sigma).*(yk).^2);
    
    % Linear Term
%     Tkl(k~=0)=(nf_p.Ly./abs(sqrt(k(k~=0)))).*(nf_p.alphay*nf_p.yp).*tanh(2*nf_p.alphay*nf_p.yp*yk(k~=0))-...
%                 (nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
    %small angle approximation for 0
    Tkl=nf_p.Ly^2.*(nf_p.alphay*nf_p.yp)^2*2-(nf_p.alphay + alpha_sigma)*nf_p.Ly^2;
    
    %**********************************************************************
    %************************ Dilation functions **************************
    %**********************************************************************
    funcf0=0.*gk;
    if cnf_p.lut_flag  
        funcf0(gk==0)=1.077900274770464; 
        indexes=ceil((gk(gk>=cnf_p.LUT_ximin & gk<=cnf_p.LUT_ximax)-cnf_p.LUT_ximin)./cnf_p.LUT_step); 
        funcf0(gk>=cnf_p.LUT_ximin & gk<=cnf_p.LUT_ximax)=func_f0(indexes); 
        funcf0(gk>cnf_p.LUT_ximax)=sqrt(pi./(2.0*gk(gk>cnf_p.LUT_ximax))).*(1+3./(8*(gk(gk>cnf_p.LUT_ximax)).^2)+105./(16*8*(gk(gk>cnf_p.LUT_ximax)).^4));
        switch cnf_p.power_wfm_model
            case 'complete'
                funcf1=0.*gk; 
                funcf1(gk==0)=0.515224256147498; 
                funcf1(gk>=cnf_p.LUT_ximin & gk<=cnf_p.LUT_ximax)=func_f1(indexes); 
                funcf1(gk>cnf_p.LUT_ximax)=-1.0*sqrt(pi.*gk(gk>cnf_p.LUT_ximax)/8).*(1./((gk(gk>cnf_p.LUT_ximax)).^2)+15./(8.*(gk(gk>cnf_p.LUT_ximax)).^4));
        end
    else   
        funcf0(gk~=0)=pi/4.0*sqrt(abs(gk(gk~=0))).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*besseli(1/4,1/4*(gk(gk~=0)).^2,1)); funcf0(gk==0)=1.077900274770464;%2^(1/4)*gamma(5/4);
        switch cnf_p.power_wfm_model
            case 'complete'
              funcf1=zeros(1,length(x)); 
              funcf1(gk~=0)=-1.0*pi/8.0*(abs(gk(gk~=0)).^(3/2)).*(besseli(1/4,1/4*(gk(gk~=0)).^2,1)-besseli(-3/4,1/4*(gk(gk~=0)).^2,1)+sign(gk(gk~=0)).*(besseli(-1/4,1/4*(gk(gk~=0)).^2,1)-besseli(3/4,1/4*(gk(gk~=0)).^2,1))); 
              funcf1(gk==0)=0.515224256147498;%gamma(3/4)/(2.0*(2)^(1/4));
        end
        
    end
    %**********************************************************************
    %************************ Power Waveform ******************************
    %**********************************************************************
    switch cnf_p.power_wfm_model
        case 'simple'
            sl_waveform         =   sqrt(g).*Bkl.*funcf0;
        case 'complete'
            sl_waveform         =   sqrt(g).*Bkl.*(funcf0+Tkl.*g.*((sigmaz/nf_p.Lz)^2).*funcf1);
    end
       
end
