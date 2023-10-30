% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
% FITGAUS_CRY_STACK_2: function that fits the main lobe of the input signal
% (wvf_in: a sinc) with a Gaussian
%
% Calling
%   parabola = fitgaus_CRY_stack_2 (wvf_in)
%
% Inputs
%   wvf_in: input signal 
%
% Outputs
%   parabola(1:nbeams,1): sample with the maximum of the Gaussian
%   parabola(1:nbeams,2): number of beam
%
% ----------------------------------------------------------
% 
% Author:   Albert Garcia / isardSAT
%           Mercedes Reche / Pildo Labs
%           Josep Montolio / Pildo Labs
% Reviewer: Monica Roca / isardSAT
%
% Last revision: Albert Garcia  / Pildo Labs (27/06/12)
%               Minor changes to avoid errors
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parabola = fitgaus (wvf_in)

s = size (wvf_in);

for db =1:s(1) 
    
    % --- retreive only the center of the sinc +-1/4 ---
          [max_v, max_p] = max(wvf_in(db,:));
          
        if(max_v ==0)
            parabola(db,1) = s(2)/2;
            parabola(db,2)= db;
            
        else
            
              %Main lobe is find defining left_lim and right_lim limits
              if(max_p==s(2))
                pos=max_p;
                aux=max_p;
              else
                pos=max_p+1;
                aux=max_p;
              end
              while ((wvf_in(db,pos)<wvf_in(db,aux))||(aux==s(2))),
                  pos=pos+1;
                  aux=aux+1;
                  %introduce to detect rigth limit at the final samples.
                  if(pos>size(wvf_in,2))
                      if(aux>size(wvf_in,2))
                          aux=aux-1;    
                      end
                      break;
                  end
              end
                  

              right_lim=aux;


%               pos=max_p-1;
%               if(max_p==1)
%                   pos=1;
%               end
%               aux=max_p; 

%               while ((wvf_in(db,pos)<wvf_in(db,aux))||(pos==1)),
%                   pos=pos-1;
%                   aux=aux-1;
%                   if(pos<1||aux<1)
%                       aux=aux+1;
%                       break;
%                   end
%               end
% 
%               left_lim=aux;
             margin = right_lim-max_p;
             left_lim=max_p-margin;
%               ind = (left_lim+10:right_lim-10);
              %Removed +10 and -10 to avoid empty indexing
              ind = (left_lim:right_lim);
              %JPLZ: added condition to avoid 0 or negative values in the waveform index for badly behaving waveforms
              if any(ind==0 | ind<0)
                 indxs_dump=find(ind<1);
                 ind=ind(length(indxs_dump)+1:length(ind));
              end
              wvf = wvf_in(db,ind);   

            % --- Computing the offset of the wvf in the x-axis and the amplitud  (to be used as C0 and A0)    ---
            [Aq, Ind] = max(wvf);  
            C0 = Ind; 
            A0 = Aq;

            % --- Computing the initial sigma ---
            B0 = 1;    

            % --- Call fminsearch_stackcal to find the parameters of the Gaussian ---
            [coeff, FVAL,EXITFLAG] = fminsearch_stackcal('echofit_CRY',[A0 B0 C0],[], wvf);
            a = coeff(1) ;
            b = coeff(2) ;
            c = coeff(3) ;

            ampl(db) = a;
            pos(db) = c; 

            
            %---------------------------------------------------

            parabola(db,1) = c + ind(1) -1; % I substract 1 because if c was the first sample, I would obtain parabola= 1+ind(1) instead of parabola = Ind(1), and this would not be correct
            parabola(db,2) = db ;
            parabola(db,3) = a;
            parabola(db,4) = b;
        end
  
end
% for i=1:45; 
% end;


