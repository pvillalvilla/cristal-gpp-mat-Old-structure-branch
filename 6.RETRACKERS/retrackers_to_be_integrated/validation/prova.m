function prova(a,varargin)
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<1 || nargin>(1+2*2))
    error('Wrong number of input parameters');   
end
p = inputParser;
p.addParameter('filename_L2_STL',{''},@(x)ischar(x));
p.addParameter('retrackers',{''},@(x)iscellstr(x));
p.parse(varargin{:});
filename_L2_STL=char(p.Results.filename_L2_STL);
retrackers=p.Results.retrackers;
clear p;

disp(a)
disp(filename_L2_STL)

end