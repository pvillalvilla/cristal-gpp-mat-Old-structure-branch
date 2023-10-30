% savevideo(workingDir)

%workingDir = '/net/n8800pro/share/roger/Jason-CS/Plots/GPP_Ph3/Scenario3/Videos/4_Ambiguity_Analysis';
% workingDir = './results/stacks2/';
workingDir = ['/data/PICE/SC62_OC_OB_KA_SAR/results/PICE_SIRS__KAHR_1B_20180824T120034_20180824T12003_isd_long/stack_analysis/'];

% mkdir(workingDir);
% mkdir(workingDir,'images');

outputVideo = VideoWriter(fullfile(workingDir,'Stack_analysis_SC62_OC_OB_KA_SAR.avi'));
outputVideo.FrameRate = 2;
open(outputVideo);

imageNames = dir(fullfile(workingDir,'*.png'));
imageNames = {imageNames.name}';

for i_imag = 1:length(imageNames)
    img = imread(fullfile(workingDir,imageNames{i_imag}));

    writeVideo(outputVideo,img);
end
close(outputVideo);