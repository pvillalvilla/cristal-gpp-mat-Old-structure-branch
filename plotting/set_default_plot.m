% close all;

% call in L2_processing.m

global defaultLineLineWidth defaultLineMarkerSize colors legend_fontsize

defaultLineLineWidth=2;
defaultLineMarkerSize=5;

text_fontsize=14;
textbox_fontsize=12;
colorbar_fontsize=10;
legend_fontsize=14;
text_interpreter = 'tex'; % 'latex' for LateX font

colors=[0, 0, 1; 1, 0, 0;  0, 1, 0];
% colors=[20/255, 115/255, 175/255; 209/255, 96/255, 47/255;  229/255, 158/255, 55/255; 25/255, 158/255, 116/255; 204/255, 121/255, 167/255; 240/255, 228/255, 66/255; 0/255, 0/255, 0/255];

% width = 8;
% height = 6;
VISIBLE=1;

% load('.\plotting\colormapDEM.mat');
% colormapATDD = colormap('hot');
% colormapATDD(1,:)=[1 1 1];
set(0,'defaultLineLineWidth',defaultLineLineWidth);   % set the default line width
set(0,'defaultLineMarkerSize',defaultLineMarkerSize);  % set the default line marker size
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',text_fontsize);
set(0,'defaultTextFontName','Arial');
set(0,'defaultTextFontSize',text_fontsize);
set(0, 'DefaultTextInterpreter', text_interpreter);
set(0, 'defaultAxesTickLabelInterpreter',text_interpreter); set(0, 'defaultLegendInterpreter',text_interpreter);

% Make ALL the plots invisible
if VISIBLE == 0
    set(0,'defaultFigureVisible','off');
else
    set(0,'defaultFigureVisible','on');
end

% Set the default Size for display
mida = get(0,'ScreenSize');
mida(3:4)=[1920,1080];
% set(0,'defaultFigurePosition', [defpos(1) defpos(2)-50 width*100, height*100]);
set(0,'defaultFigurePosition',mida);
%close

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
% set(0,'defaultFigurePaperUnits','centimeters');
% defsize = get(gcf, 'PaperSize');
% left = (defsize(1)- width)/2;
% bottom = (defsize(2)- height)/2;
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);
set(0,'defaultFigurePaperUnits','points');
set(0,'defaultFigurePaperPosition', mida);

% plot(...)
% print('figName','-r100','-dpng')