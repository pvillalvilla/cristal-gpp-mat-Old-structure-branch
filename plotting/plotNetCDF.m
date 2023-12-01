% only working with two levels of groups /data/band 

% Specify the path to your NetCDF file
current_path=pwd;
files_netcdf=dir([current_path '/*.nc']);
for i_file=1:length(files_netcdf)
    
ncFilePath = files_netcdf(i_file).name;
load('G:\My Drive\Matlab\altimetry_tools_matlab\plotting\colormaps\colormap_blues.mat');
set_default_plot;

%Modify plot size
mida(3:4)=[1920,1080]./2;
set(0,'defaultfigurePosition',mida);
set(0,'defaultfigureVisible','off');


% Create a subfolder for saving plots
outputFolder = fullfile(pwd, ncFilePath(1:end-20));
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
else
    continue;
end

fidResults=fopen([outputFolder '/' ncFilePath(1:end-20) 'plot_log.csv'],'w');

% Get information about the NetCDF file
info = ncinfo(ncFilePath);

% Loop through and display the variables and top-level groups
for i = 1:numel(info.Variables)
    varname = info.Variables(i).Name;
    data = ncread(ncFilePath, varname);
    
    % Retrieve units from the Attributes structure
    units = '';
    for attrIndex = 1:numel(info.Variables(i).Attributes)
        if strcmp(info.Variables(i).Attributes(attrIndex).Name, 'units')
            units = info.Variables(i).Attributes(attrIndex).Value;
            units = strrep(units, '_', ' ');
            break;
        end
    end
    
    % Handle variable names with "_" symbol
    varname = strrep(varname, '_', ' ');
    disp(['Variable Name: ', varname, ' [', units, ']']);
    if(size(data,1)==3)
        figure('Name','2','visible','off'); 
        subplot(3,1,1);
        plot(data(1,:), 'o-');
        title([varname, ' X component [', units, ']']);
        xlabel('X-axis');
        ylabel(['[' units ']']);
        subplot(3,1,2);
        plot(data(2,:), 'o-');
        title([varname, ' Y component [', units, ']']);
        xlabel('X-axis');
        ylabel(['[' units ']']);
        subplot(3,1,3);
        plot(data(3,:), 'o-');
        title([varname, ' Z component [', units, ']']);
        xlabel('X-axis');
        ylabel(['[' units ']']);
        pngFilename = fullfile(outputFolder, [varname, '.png']);
        saveas(gcf, pngFilename);
        fprintf(fidResults,'%s,  [%s];  see values in plot \n', varname, units);
        
    elseif(length(data)>1)
        if isvector(data) % Check if the variable is a 1D vector
            % Create a 1D plot with '.-' line symbol
            figure('Name','1','visible','off');  
            plot(data, '.-');
            title([varname, ' [', units, ']']);
            xlabel('X-axis');
            ylabel(['[' units ']']);
            pngFilename = fullfile(outputFolder, [varname, '.png']);
            saveas(gcf, pngFilename);
            fprintf(fidResults,'%s,  [%s];  see values in plot \n', varname, units);
            
        elseif length(size(data))==2 % Check if the variable is a 2D matrix
            % Create an imagesc plot
            figure('Name','1','visible','off');  
            imagesc(data);
            colorbar;
            title([varname, ' [', units, ']']);
            xlabel('X-axis');
            ylabel('Y-axis');
            pngFilename = fullfile(outputFolder, [varname, '.png']);
            saveas(gcf, pngFilename);
            fprintf(fidResults,'%s,  [%s];  see values in plot \n', varname, units);
        else
            
        end
    else
        %write single value variables in a csv file
        fprintf(fidResults,'%s,  [%s]; %d\n', varname, units, data);
        
        
    end
end

% Loop through and display the groups
for i = 1:numel(info.Groups)
    groupname = info.Groups(i).Name;
    
    % Handle group names with "_" symbol
    groupname_txt = strrep(groupname, '_', ' ');
    
    disp(['Group Name: ', groupname]);
    fprintf(fidResults,'Group Name:%s;\n', groupname);
    % Get information about the variables within the group
    groupInfo = ncinfo(ncFilePath, groupname);
    
    
    
    % Loop through and display the variables in the group
    for j = 1:numel(groupInfo.Variables)
        varname = groupInfo.Variables(j).Name;
        data = ncread(ncFilePath, [groupname, '/', varname]);
        
        % Retrieve units from the Attributes structure
        units = '';
        for attrIndex = 1:numel(groupInfo.Variables(j).Attributes)
            if strcmp(groupInfo.Variables(j).Attributes(attrIndex).Name, 'units')
                units = groupInfo.Variables(j).Attributes(attrIndex).Value;
                units = strrep(units, '_', ' ');
                break;
            end
        end
        
        % Handle variable names with "_" symbol
        varname = strrep(varname, '_', ' ');
        
        disp(['   Variable Name: ', varname, ' [', units, ']']);
        if(size(data,1)==3)
            figure('Name','2','visible','off'); 
            subplot(3,1,1);
            plot(data(1,:), 'o-');
            title([varname, ' X component [', units, ']']);
            xlabel('X-axis');
            ylabel(['[' units ']']);
            subplot(3,1,2);
            plot(data(2,:), 'o-');
            title([varname, ' Y component [', units, ']']);
            xlabel('X-axis');
            ylabel(['[' units ']']);
            subplot(3,1,3);
            plot(data(3,:), 'o-');
            plot(data, '.-');
            title([varname, ' Z component [', units, ']']);
            xlabel('X-axis');
            ylabel(['[' units ']']);
            pngFilename = fullfile(outputFolder, [groupname '_' varname, '.png']);
            saveas(gcf, pngFilename);
            fprintf(fidResults,'    %s,  [%s];  see values in plot \n', varname, units);
            
        elseif(length(data)>1)
            if isvector(data) % Check if the variable is a 1D vector
                % Create a 1D plot with '.-' line symbol
                figure('Name','1','visible','off');  
                plot(data, 'o-');
                title([varname, ' [', units, ']']);
                xlabel('X-axis');
                ylabel(['[' units ']']);
                pngFilename = fullfile(outputFolder, [groupname '_' varname, '.png']);
                saveas(gcf, pngFilename);
                fprintf(fidResults,'    %s,  [%s];  see values in plot \n', varname, units);
            elseif length(size(data))==2  % Check if the variable is a 2D matrix
                % Create an imagesc plot
                figure('Name','1','visible','off');  
                imagesc(data);
                colormap(colormap_blues);
                hcb=colorbar;
                hcb.Label.String = ['[' units ']'];
                title([varname, ' [', units, ']']);
                xlabel('X-axis');
                ylabel('Y-axis');
                pngFilename = fullfile(outputFolder, [groupname '_' varname, '.png']);
                saveas(gcf, pngFilename);
                fprintf(fidResults,'    %s,  [%s];  see values in plot \n', varname, units);
            elseif length(size(data))==3  % Check if the variable is a ?D matrix
                fprintf(fidResults,'    %s,  [%s];  information not plotted \n', varname, units);
            else
                
            end
        else
            %write single value variables in a csv file
            fprintf(fidResults,'    %s,  [%s]; %d\n', varname, units, data);
        end
    end
    
    % Explore one more level of subgroups within the group
    for j = 1:numel(groupInfo.Groups)
        subgroupname = groupInfo.Groups(j).Name;
        
        disp(['   Subgroup Name: ', subgroupname]);
        fprintf(fidResults,'    Subgroup Name::%s;\n', subgroupname);
        % Get information about the variables within the subgroup
        subgroupInfo = ncinfo(ncFilePath, [groupname, '/', subgroupname]);
        
        % Handle subgroup names with "_" symbol
        subgroupname_txt = strrep(subgroupname, '_', ' ');
        
        % Loop through and display the variables in the subgroup
        for k = 1:numel(subgroupInfo.Variables)
            varname = subgroupInfo.Variables(k).Name;
            data = ncread(ncFilePath, [groupname, '/', subgroupname, '/', varname]);
            
            % Retrieve units from the Attributes structure
            units = '';
            for attrIndex = 1:numel(subgroupInfo.Variables(k).Attributes)
                if strcmp(subgroupInfo.Variables(k).Attributes(attrIndex).Name, 'units')
                    units = subgroupInfo.Variables(k).Attributes(attrIndex).Value;
                    units = strrep(units, '_', ' ');
                    break;
                end
            end
            
            % Handle variable names with "_" symbol
            varname = strrep(varname, '_', ' ');
            
            disp(['      Variable Name: ', varname, ' [', units, ']']);
            if(size(data,1)==3)
                figure('Name','2','visible','off'); 
                subplot(3,1,1);
                plot(data(1,:), 'o-');
                title([varname, ' X component [', units, ']']);
                xlabel('X-axis');
                ylabel(['[' units ']']);
                subplot(3,1,2);
                plot(data(2,:), 'o-');
                title([varname, ' Y component [', units, ']']);
                xlabel('X-axis');
                ylabel(['[' units ']']);
                subplot(3,1,3);
                plot(data(3,:), 'o-');
                title([varname, ' Z component [', units, ']']);
                xlabel('X-axis');
                ylabel(['[' units ']']);
                pngFilename = fullfile(outputFolder, [groupname '_' subgroupname '_' varname, '.png']);
                saveas(gcf, pngFilename);
                fprintf(fidResults,'        %s,  [%s];  see values in plot \n', varname, units);
            elseif(length(data)>1)
                if isvector(data) % Check if the variable is a 1D vector
                    % Create a 1D plot with '.-' line symbol
                    figure('Name','1','visible','off');  
                    plot(data, 'o-');
                    title([groupname_txt, '/', subgroupname_txt, '/', varname, ' [', units, ']']);
                    xlabel('X-axis');
                    ylabel(['[' units ']']);
                    % Save the plot as a PNG file
                    pngFilename = fullfile(outputFolder, [groupname '_' subgroupname '_' varname, '.png']);
                    saveas(gcf, pngFilename);
                    fprintf(fidResults,'        %s,  [%s];  see values in plot \n', varname, units);
                elseif length(size(data))==2  % Check if the variable is a 2D matrix
                    % Create an imagesc plot
                    figure('Name','1','visible','off');  
                    imagesc(data);
                    colormap(colormap_blues);
                    hcb=colorbar;
                    hcb.Label.String = ['[' units ']'];
                    title([groupname_txt, '/', subgroupname_txt, '/', varname, ' [', units, ']']);
                    xlabel('X-axis');
                    ylabel('Y-axis');
                    % Save the plot as a PNG file
                    pngFilename = fullfile(outputFolder, [groupname '_' subgroupname '_' varname, '.png']);
                    saveas(gcf, pngFilename);
                    fprintf(fidResults,'        %s,  [%s];  see values in plot \n', varname, units);
                elseif length(size(data))==3  % Check if the variable is a 2D matrix
                    
                end
            else
                
                %write single value variables in a csv file
                fprintf(fidResults,'        %s,  [%s]; %d\n', varname, units, data);
                
            end
        end
    end
end
close all
fclose (fidResults);
end