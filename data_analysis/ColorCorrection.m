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
% to convert spectral bands to raw colour channels. Can be empty (`[]`) if there
% are no spectral images to process.
sampling_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/RunOnDataset_20190208_ComputarLens_rawCaptured_dispersion.mat';

% Path and filename of a '.mat' file containing the colour space conversion
% data needed to convert spectral bands to raw colour channels. Can be
% empty (`[]`) if there are no spectral images to process.
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/SonyColorMapData.mat';

% Whether to use white balancing (`true`), or a more advanced colour
% conversion method (`false`)
use_chromadapt = false;

if use_chromadapt
    % Path and filename of a '.mat' file containing the illuminant colour
    % estimated from a neutral patch, as the variable 'wb_illum'
    wb_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/CalibrateColorCorrectionData_unfiltered.mat';
else
    % Path and filename of a '.mat' file containing the conversion data
    % structure for raw colour channels to XYZ
    xyz_weights_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/CalibrateColorCorrectionData_unfiltered.mat';
    xyz_weights_variable = 'M_homog'; % Variable name in the above file
    
    % Whitepoint to use for XYZ to sRGB conversion
    whitepoint = [1, 1, 1];
end

% Wildcard for 'ls()' to find the spectral images to process (can be empty (`[]`)).
% '.mat' or image files can be loaded
spectral_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/*_latent*.mat';
spectral_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the colour images to process (can be empty (`[]`)).
% '.mat' or image files can be loaded
color_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/*_rgb*.mat';
color_variable_name = 'I_rgb'; % Used only when loading '.mat' files

% Output directory
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/correctedRGB';

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
    [sensor_map, ~, bands_color] = loadColorMap(color_map_filename);
    [findSamplingOptions, bands] = loadVariables(sampling_filename, {'findSamplingOptions', 'bands'});
    color_weights = colorWeights(...
        sensor_map, bands_color, bands, findSamplingOptions...
    );
end

if n_images ~= 0
    if use_chromadapt
        [wb_illum, wb_scale] = loadVariables(wb_filename, {'wb_illum', 'wb_scale'});
        postfix = '_wb';
    else
        xyz_weights = loadVariables(xyz_weights_filename, xyz_weights_variable);
        postfix = ['_' xyz_weights_variable];
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
    
    if use_chromadapt
        I = I ./ wb_scale;
        I = chromadapt(I, wb_illum, 'ColorSpace','linear-rgb');
        I = lin2rgb(I, 'OutputType', 'uint8');
    else
        if strcmp(xyz_weights_variable, 'M_rp')
            sz = size(I);
            I = reshape(I, [], sz(3));
            I = reshape(xyz_weights.cfun(I, xyz_weights.matrix, xyz_weights.terms), sz);
        else
            I = channelConversion(I, xyz_weights);
        end
        I = xyz2rgb(I, 'ColorSpace', 'srgb', 'WhitePoint', whitepoint, 'OutputType', 'uint8');
    end
    saveImages(...
        'image', output_directory, filename,...
        I, postfix, []...
    );
end
