function strout=read_xml_aux_sarincb(xmlfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT 
% --------------------------------------------------------
% CRISTAL mission
%
% ---------------------------------------------------------
% Objective: read CRISTAL SARIN CB CHD/CNF/CST XML auxiliary file
%
% Calling: 
% INPUTs: xmlfile
%   xmlfile: xml aux file to be read (can be CST, CNF or CHD)
%  (note that there is a unique CST aux file: same output for all modes)
%
% OUTPUTs:
%   strout: Matlab structure with all xml variables
%           Fields are: Name / Description / Value / Units / Type      
%
% ----------------------------------------------------------
% Author:    Pablo Garcia  / isardSAT
%            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

structxml = xml2struct(xmlfile);

numparam=length(structxml.Children(4).Children);

par_count=0;
for i_par = 1:numparam
    if strcmp(structxml.Children(4).Children(i_par).Name,'Parameter')
        par_count = par_count + 1;
        for i_atr = 1:length(structxml.Children(4).Children(i_par).Children)
            if ~strcmp(structxml.Children(4).Children(i_par).Children(i_atr).Name,'#text') && ~isempty(structxml.Children(4).Children(i_par).Children(i_atr).Children) % if there is data to read
                if isempty(structxml.Children(4).Children(i_par).Children(i_atr).Attributes) % if there are not Attributes to read (same for all modes)
                    strtemp(par_count).(structxml.Children(4).Children(i_par).Children(i_atr).Name) = structxml.Children(4).Children(i_par).Children(i_atr).Children.Data;
                elseif (strcmp(structxml.Children(4).Children(i_par).Children(i_atr).Attributes.Name,'mode') && strcmp(structxml.Children(4).Children(i_par).Children(i_atr).Attributes.Value,'CB')) || ... % mode = CB
                        (strcmp(structxml.Children(4).Children(i_par).Children(i_atr).Attributes.Name,'mode') && strcmp(structxml.Children(4).Children(i_par).Children(i_atr).Attributes.Value,'SARIN CB'))    % mode = SARIN CB
                    strtemp(par_count).(structxml.Children(4).Children(i_par).Children(i_atr).Name) = structxml.Children(4).Children(i_par).Children(i_atr).Children.Data;
                end
            end
        end
    end
end

% clean from long "only space" xml strings (xml line jumps)
for i_par=1:length(strtemp)
    strtemp(i_par).Description = strrep(strtemp(i_par).Description, '                ', ' ');
end

% Build the output matlab file with 2 groups: 
% 1- data (values)
% 2- attributes (description / units / type)
for i=1:length(strtemp)
strout.data.(strtemp(i).Name)=str2double(strtemp(i).Value);
strout.attributes.(strtemp(i).Name).description=strtemp(i).Description;
strout.attributes.(strtemp(i).Name).units=strtemp(i).Units;
strout.attributes.(strtemp(i).Name).type=strtemp(i).Type;
end

end