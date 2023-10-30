%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
%
% ---------------------------------------------------------
% Objective: Format the paths to work in unix and windows OS.   
% 
% Calling: path = fix_paths (path)
% INPUTs: path string 
%           
% OUTPUTs: reformatted path string
%   
% ----------------------------------------------------------
% Author:    Albert Garcia  / isardSAT          
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path = fix_paths (path)

if(strfind( path, '\'))
    path(strfind( path, '\')) = '/';    % Change '\' for '/'
    if(~strcmp(path(end), '\'))
        path = [path '/'];              % Add / at the end if needed
    end
end

end