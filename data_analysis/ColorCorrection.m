%% Convert images to sRGB
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Description
%
% Reads in spectral and/or colour images, and converts them to the sRGB
% colour space. There is no output aside from TIFF image files.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 10, 2019

%% Input data and parameters

% Path and filename of a '.mat' file containing the sampling information needed
% to convert spectral bands to raw colour channels
sampling_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_dispersion/RunOnDataset_20190208_ComputarLens_rawCaptured_dispersion.mat';

% Path and filename of a '.mat' file containing the colour space conversion data
% needed to convert spectral bands to raw colour channels
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/SonyColorMapData.mat';

% Path and filename of a '.mat' file containing the conversion matrix for
% raw colour channels to XYZ
xyz_weights_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_correctedRGB/CalibrateColorCorrectionData.mat';
xyz_weights_variable = 'M_homog'; % Variable name in the above file

% Whitepoint to use for XYZ to sRGB conversion
whitepoint = [1, 1, 1];

% Wildcard for 'ls()' to find the spectral images to process.
% '.mat' or image files can be loaded
spectral_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_MATFiles/*_latent*.mat';
spectral_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the colour images to process.
% '.mat' or image files can be loaded
color_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_MATFiles/*_rgb*.mat';
color_variable_name = 'I_rgb'; % Used only when loading '.mat' files

% Output directory
output_directory = '/home/llanos/Downloads';

%% Process the images

n_spectral = 0;
if ~isempty(spectral_wildcard)
    filenames.spectral = listFiles(spectral_wildcard);
    n_spectral = length(filenames.spectral);
end

n_color = 0;
if ~isempty(color_wildcard)
    filenames.color = listFiles(color_wildcard);
    n_color = length(filenames.color);
end

n_images = n_spectral + n_color;

if n_spectral ~= 0
    variables_required = { 'sensor_map', 'bands' };
    load(color_map_filename, variables_required{:});
    if ~all(ismember(variables_required, who))
        error('One or more of the required colour space conversion variables is not loaded.')
    end
    bands_color = bands;
    clear bands

    variables_required = { 'findSamplingOptions', 'bands' };
    load(sampling_filename, variables_required{:});
    if ~all(ismember(variables_required, who))
        error('One or more of the required sampling variables is not loaded.')
    end
    color_weights = colorWeights(...
        sensor_map, bands_color, bands, findSamplingOptions...
    );
end

if n_images ~= 0
    load(xyz_weights_filename, xyz_weights_variable);
    if exist(xyz_weights_variable, 'var')
        xyz_weights = eval(xyz_weights_variable);
    end
    if ~exist(xyz_weights_variable, 'var') || isempty(xyz_weights)
        error('No raw colour to XYZ conversion matrix loaded.')
    end
end

for i = 1:n_images
    if i <= n_spectral
        [~, filename] = fileparts(filenames.spectral{i});
        I = loadImage(filenames.spectral{i}, spectral_variable_name);
        I = channelConversion(I, color_weights);
    else
        [~, filename] = fileparts(filenames.color{i - n_spectral});
        I = loadImage(filenames.color{i - n_spectral}, color_variable_name);
    end
    
    I = channelConversion(I, xyz_weights);
    I = xyz2rgb(I, 'ColorSpace', 'srgb', 'WhitePoint', whitepoint, 'OutputType', 'uint8');
    saveImages(...
        'image', output_directory, filename,...
        I, '_correctedRGB', []...
    );
end
