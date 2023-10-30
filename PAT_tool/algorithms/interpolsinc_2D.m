%+======================================================================
%  Created by @Eduard Makhoul 
%=======================================================================
% NAME:
% interpolsinc_2D.pro
%
% INCLUDED FUNCTIONS/PROCEDURES/STRUCTURES:
% -
%
% CALLING SEQUENCE:
%   -
%
% INPUTS:
% Mandatory:
%   -M: Input matrix to be interpolated
%   factor_c: resampling factor in columns
%   factor_r: resampling factor in rows
%    
%   
%     
% Optional:
%  -
%
% OUTPUT:
% RETURN    :
%
% KEYWORD Parameters:
% None.
%
% CALLED PROCEDURES/FUNCTIONS/CLASSES:
% None.
%
% COMMON BLOCKS:
% None.
%
% SIDE EFFECTS:
% None.
%
% RESTRICTIONS:
%  Dimensions of the matrices must be an odd number (to ensure it a zero-padding is performed)
%
% COMMENTS:
% MATRIX NOTATION: [row, column]
% 
%
% MISSING IMPLEMENTATIONS:
%  - CONSIDER OTHER OPTIONS OF DATA TYPE (DIFFERENT FROM THE DOUBLE COMPLEX)
%  - 
% MODIFICATION HISTORY:
%
% written: Eduardo Makhoul 13.04.2018
%   modifications:  

%======================================================================

function  [M_interp]=interpolsinc_2D(M,factor_c,factor_r)


s=size(M);
num1=s(1);%number of input rows
num2=s(2);%number of input columns

M_exp=M;

%ZERO-PADDING TO ENSURE THE IMAGE IS ODD
if mod(num1,2)==1 
  num1=num1+1;
  M_exp=[M_exp; zeros(1,num2)];
end
if mod(num2,2)==1
  num2=num2+1;
  M_exp=[M_exp,zeros(num1,1)];
end

num1_int=num1*factor_r;%number of output rows
num2_int=num2*factor_c;%number of output columns

IMATGE=fft2(M_exp);
%IMATGE=ifft(fft(M_exp,[],1),[],2);
clear M_exp;

%interpolated matrix
IMATGE_INTERP=zeros(num1_int,num2_int);

% Primer quadrant
IMATGE_INTERP(1:floor(num1/2),1:floor(num2/2))=IMATGE(1:floor(num1/2),1:floor(num2/2));
IMATGE_INTERP(1:floor(num1/2),floor(num2/2)+1)=IMATGE(1:floor(num1/2),floor(num2/2)+1)/2;
IMATGE_INTERP(floor(num1/2)+1,1:floor(num2/2))=IMATGE(floor(num1/2)+1,1:floor(num2/2))/2;
IMATGE_INTERP(floor(num1/2)+1,floor(num2/2)+1)=IMATGE(floor(num1/2)+1,floor(num2/2)+1)/4;

% Segon quadrant
IMATGE_INTERP(1:floor(num1/2),(num2_int+1)-(floor(num2/2)-1):num2_int)=IMATGE(1:floor(num1/2),floor(num2/2)+2:num2);
IMATGE_INTERP(1:floor(num1/2),(num2_int+1)-floor(num2/2))=IMATGE(1:floor(num1/2),floor(num2/2)+1)/2;
IMATGE_INTERP(floor(num1/2)+1,(num2_int+1)-(floor(num2/2)-1):num2_int)=IMATGE(floor(num1/2)+1,floor(num2/2)+2:num2)/2;
IMATGE_INTERP(floor(num1/2)+1,(num2_int+1)-floor(num2/2))=IMATGE(floor(num1/2)+1,floor(num2/2)+1)/4;

%Tercer quadrant
IMATGE_INTERP((num1_int+1)-(floor(num1/2)-1):num1_int,1:floor(num2/2))=IMATGE(floor(num1/2)+2:num1,1:floor(num2/2));
IMATGE_INTERP((num1_int+1)-(floor(num1/2)-1):num1_int,floor(num2/2)+1)=IMATGE(floor(num1/2)+2:num1,floor(num2/2)+1)/2;
IMATGE_INTERP((num1_int+1)-floor(num1/2),1:floor(num2/2))=IMATGE(floor(num1/2)+1,1:floor(num2/2))/2;
IMATGE_INTERP((num1_int+1)-floor(num1/2),floor(num2/2)+1)=IMATGE(floor(num1/2)+1,floor(num2/2)+1)/4;

% Quart quadrant
IMATGE_INTERP((num1_int+1)-(floor(num1/2)-1):num1_int,(num2_int+1)-(floor(num2/2)-1):num2_int)=IMATGE(floor(num1/2)+2:num1,floor(num2/2)+2:num2);
IMATGE_INTERP((num1_int+1)-(floor(num1/2)-1):num1_int,(num2_int+1)-floor(num2/2))=IMATGE(floor(num1/2)+2:num1,floor(num2/2)+1)/2;
IMATGE_INTERP((num1_int+1)-floor(num1/2),(num2_int+1)-(floor(num2/2)-1):num2_int)=IMATGE(floor(num1/2)+1,floor(num2/2)+2:num2)/2;
IMATGE_INTERP((num1_int+1)-floor(num1/2),(num2_int+1)-floor(num2/2))=IMATGE(floor(num1/2)+1,floor(num2/2)+1)/4;

%print, IMATGE_INTERP

%print,''
%print,'***zero-padding-->FINISHED'

%imatge_interpt=dcomplexarr(num2_int,num1_int)
%imatge_interpt=FFT(IMATGE_INTERP,/INVERSE)

M_interp=ifft2(IMATGE_INTERP.*factor_c*factor_r);
%M_interp=fft(ifft(IMATGE_INTERP,[],1),[],2);

end