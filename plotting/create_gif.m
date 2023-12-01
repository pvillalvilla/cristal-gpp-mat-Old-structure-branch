% Specify the folder containing your images
imageFolder = 'G:\My Drive\S3NGT\codes\S3NGT_TOO4\inputs\';


% Get a list of all image files in the folder
imageFiles = dir(fullfile(imageFolder, '*.png')); % Change '*.png' to the appropriate file extension of your images
[~, order] = sort( [imageFiles.datenum]);
imageFiles = imageFiles(order);

% Create a GIF file
outputGIF = [imageFolder 'S3NGT_TOO4.gif'];

load('G:\My Drive\Matlab\altimetry_tools_matlab\plotting\colormaps\colormap_blues.mat');
cmap = colormap_blues;

% Loop through each image and add it to the GIF
for i = 1:numel(imageFiles)
    % Read the image
    currentImage = imread(fullfile(imageFolder, imageFiles(i).name));
    [indexedImage, map] = rgb2ind(currentImage, size(cmap, 1));

    % Write the image to the GIF file
    if i == 1
        % If it's the first image, create the GIF file
        imwrite(indexedImage, map, outputGIF, 'gif', 'Loopcount', inf, 'DelayTime', 0.12);
    else
        % If it's not the first image, append to the existing GIF file
        imwrite(indexedImage, map, outputGIF, 'gif', 'WriteMode', 'append', 'DelayTime', 0.12);
    end
end

disp('GIF creation complete.');
