%% Generate input data for the modified version of Choi et al. 2017's code
%
% Output reflectance images in '.mat' format, and information on how to
% convert those images to simulated RAW images.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Notes
% - If the spectral sampling of the reference images does not match the
%   spectral sampling used by Choi et al. 2017, then this script will
%   resample in the spectral space when generating input images for Choi et
%   al. 2017. The conversion matrices mapping between the two spaces will
%   pad spectral signals with zeros outside of their domains.
%
% ## References
% - Choi, I., Jeon, D. S., Gutierrez, D., & Kim, M. H. (2017).
%   "High-Quality Hyperspectral Reconstruction Using a Spectral Prior." ACM
%   Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 218:1-13.
%   10.1145/3130800.3130810

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 26, 2018

% List of parameters to save with results
parameters_list = {
    'crop',...
    'color_map_filename',...
    'bands_filename',...
    'bands_variable',...
    'desired_bands_filename',...
    'desired_bands_variable',...
    'illuminant_filename',...
    'illuminant_temperature',...
    'normalization_channel',...
};

%% Input data and parameters

% Wildcard for 'ls()' to find the input *reflectance* images
% '.mat' or image files can be loaded
input_images_wildcard = '${FILEPATH}';
input_images_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% A crop box for the images (leave empty for no cropping)
crop = []; % [top-left x, top-left y, width, height]

% Colour space conversion data
color_map_filename = '${FILEPATH}';

% Path and filename of a '.mat' file containing the wavelengths or colour
% channel indices corresponding to the images
bands_filename = '${FILEPATH}';
bands_variable = 'bands'; % Variable name in the above file

% Path and filename of a '.mat' file containing the wavelengths or colour
% channel indices corresponding to the images that will be estimated by the
% code from Choi et al. 2017
desired_bands_filename = '${FILEPATH}';
desired_bands_variable = 'wvls2b'; % Variable name in the above file

% ## Parameters for creating radiance images

% CIE D-illuminant
illuminant_filename = '${FILEPATH}';
illuminant_temperature = 6504; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D

% Colour channel to use for radiance normalization
normalization_channel = 2;

% Output directory for all images and saved parameters
output_directory = '${DIRPATH}';

% Set the Bayer pattern code, and spectral resampling options
run('SetFixedParameters.m')

%% Load wavelengths

bands_spectral = loadVariables(bands_filename, bands_variable);
desired_bands_spectral = loadVariables(desired_bands_filename, desired_bands_variable);

%% Load sensor data

[sensor_map, ~, bands_color] = loadColorMap(color_map_filename, false);

%% Load the illuminant

illuminant_data = csvread(illuminant_filename);
bands_illuminant = illuminant_data(:, 1);
S_illuminant = illuminant_data(:, 2:end);
spd_illuminant = ciedIlluminant(...
    illuminant_temperature, bands_illuminant, S_illuminant, bands_illuminant...
);

%% Generate colour and spectral conversion matrices

resample_images = (length(desired_bands_spectral) ~= length(bands_spectral)) || any(desired_bands_spectral ~= bands_spectral);
if resample_images
    spectral_weights_output = resamplingWeights(...
        bands_spectral, desired_bands_spectral, findSamplingOptions.interpolant_ref, 0 ...
    );
    spectral_weights_input = resamplingWeights(...
        desired_bands_spectral, bands_spectral, findSamplingOptions.interpolant_ref, 0 ...
    );
else
    spectral_weights_output = eye(length(desired_bands_spectral));
    spectral_weights_input = spectral_weights_output;
end

[bands_radiance, ~, radiance_normalized_weights] = reflectanceToRadiance(...
    bands_illuminant, spd_illuminant,...
    desired_bands_spectral, eye(length(desired_bands_spectral)),...
    bands_color, sensor_map.',...
    normalization_channel, findSamplingOptions.int_method...
);
spectral_weights_output = spectral_weights_output * resampleArrays(...
    bands_radiance, radiance_normalized_weights, desired_bands_spectral,...
    'spline', 0 ...
);

color_weights_options = findSamplingOptions;
color_weights_options.interpolant = findSamplingOptions.interpolant_ref;
color_weights_radiance = colorWeights(...
    sensor_map, bands_color, bands_radiance, color_weights_options...
);

% Columns in `radiance_normalized_weights` represent radiances
% corresponding to the impulse reflectances. So
% `radiance_normalized_weights` is a conversion matrix from reflectances to
% radiances. Then `radiance_normalized_weights.'` is a matrix where the
% rows represent spectral radiances corresponding to the impulse
% reflectances. Converting the rows to colours will give a mapping from
% impulse reflectances to colours, after taking the transpose.
color_weights = channelConversion(...
    radiance_normalized_weights.', color_weights_radiance, 2 ...
);
color_weights = color_weights.';

%% Output images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

if ~isempty(crop) && (any(mod(crop(1:2), 2) ~= 1) || any(mod(crop(3:4), 2) ~= 0))
    error('The crop region''s top-left corner should have odd-integer coordinates, and its size should be even integers, to preserve the colour filter array.');
end

for i = 1:n_images
    [I, name] = loadImage(image_filenames{i}, input_images_variable_name);
    image_height = size(I, 1);
    image_width = size(I, 2);
    if any(mod([image_height, image_width], 2) ~= 0)
        error('The image dimensions must be even integers in order for the image to be a valid colour filter array.');
    end
    postfix = '_choiIn';
    if ~isempty(crop)
        if crop(2) < 1 || crop(2) > image_height ||...
           crop(1) < 1 || crop(1) > image_width
            error('The top-left corner of the crop region is outside image %d.', i);
        end
        crop_end = crop(1:2) + crop(3:4) - 1;
        if crop_end(2) < 1 || crop_end(2) > image_height ||...
           crop_end(1) < 1 || crop_end(1) > image_width
            error('The bottom-right corner of the crop region is outside image %d.', i);
        end
        I = I(crop(2):crop_end(2), crop(1):crop_end(1), :);
        postfix = [postfix, sprintf('_X%dY%dW%dH%d', crop(1), crop(2), crop(3), crop(4))];
    end
    
    if resample_images
        I = channelConversion(...
            I, spectral_weights_input, 3 ...
        );
    end

    saveImages(...
        'data', output_directory, name,...
        I, postfix, 'img_hs'...
    );
end

%% Save parameters and common data to a file
save_variables_list = [ parameters_list, {...
        'bands_spectral',...
        'desired_bands_spectral',...
        'resample_images',...
        'spectral_weights_output',...
        'spectral_weights_input',...
        'bands_radiance',...
        'radiance_normalized_weights',...
        'color_weights',...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'Choi2017Input.mat');
save(save_data_filename, save_variables_list{:});