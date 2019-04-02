%% Assemble images for spectral bands into spectral and colour images
% Load images taken under bandpass optical filters, and use them to create a
% dataset of spectral images.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
% 
% A set of RAW images taken under different bandpass optical filters. The
% image filenames are expected to contain the wavelengths of the spectral
% bands to which they correspond. Images are expected to have been taken
% for one or more scenes, where, for all scenes, images for the same
% spectral bands were captured.
%
% Images are expected to have been preprocessed, such as using
% 'PreprocessRAWImages.m', so that they do not need to be linearized after
% being loaded, and so that all image pixels are properly exposed (being
% neither saturated nor zero). As such, the individual images might have
% been formed from sets of images taken under different exposures. For
% image format files, images will simply be loaded with the Image
% Processing Toolbox 'imread()' function. For '.mat' files, the variable to
% be loaded must be provided in the script parameters.
%
% These images are assembled into spectral images.
%
% ### Calibration images
%
% A set of RAW images taken under different bandpass optical filters; possibly a
% subset of the input images. These images must satisfy the same criteria
% described for the input images above, but are used for calibrating the
% relative spectral sensitivities of the colour channels of the camera, instead
% of being combined into spectral images.
%
% ## Output
%
% ### Relative spectral sensitivities, and wavelength information
%
% A 'sensor.mat' file containing the following variables:
% - 'channel_mode': A Boolean value set to `false` to indicate that the
%   data in `sensor_map` represents spectral sensitivities, not colour
%   channel mappings
% - 'sensor_map_relative': A 2D array, where `sensor_map_relative(i, j)` is the
%   sensitivity of the i-th colour channel (Red, Green, or Blue) to light of the
%   j-th wavelength (in the j-th spectral band), relative to the sensitivity of
%   the Green channel. Consequently, `sensor_map_relative(2, :)` is a vector of
%   ones.
% - 'channel_bias': A 2D array with the same dimensions as
%   'sensor_map_relative'. 'channel_bias' is the 'bias' output argument of
%   'relativeSensitivity()', and provides corrections to conversions between
%   colour channels computed using only 'sensor_map_relative'. The bias will not
%   be used if `use_bias`, set in the script parameters below, is `false`.
% - 'sensor_map': A version of `sensor_map_relative` that has been normalized
%   so that the Green channel value resulting from a spectral image with an
%   intensity of one in all bands is one.
% - 'bands': A vector containing the wavelengths corresponding to the
%   columns of `sensor_map`.
%
% The relative spectral sensitivities of the Red and Blue channels with respect
% to the Green channel are computed, by 'relativeSensitivity()', using the set
% of calibration images.
%
% ### Output images
%
% One of each of the following types of images will be created for each
% scene (group of images differing only by spectral band). The scene name,
% is represented by '*' below.
%
% #### RAW images
%
% '*_raw.tif' and '*_raw.mat' contain a colour filter array image (stored
% in the variable 'I_raw'). This image is created by mosaicing a quasi-RGB
% image (described below). While it is a simulated image, it is effectively
% a ground truth image. The only aspect of its generation that is not
% directly measured is the assumption that the spectral data between the
% measured spectral bands in the quasi-hyperspectral images (described
% below) is well-approximated by interpolating the measured spectral bands.
% There is no error from demosaicing, because the hyperspectral images were
% not obtained by demosaicing. Nor is there error from the estimated
% spectral sensitivities of the sensor, because the process of creating the
% quasi-hyperspectral image, and then converting it to a colour image,
% divides, and then multiplies, by the relative spectral sensitivities.
%
% #### Quasi RGB images
%
% '*_q3.tif' and '*_q3.mat' contain a 3-channel colour image (stored in the
% variable 'I_3'). This image is created by converting the
% quasi-hyperspectral image (described below) to colour, according to the
% relative spectral sensitivities of the sensor (described above).
%
% #### Quasi hyperspectral images
%
% '*_qHyper.mat' contains a spectral radiance image, stored in the variable
% 'I_hyper'. This image is created from the spectral band images taken for the
% scene. Each spectral band image's colour channels are divided by their
% spectral sensitivities (in 'sensor_map_relative'), to recover an image from a
% hypothetical monochromatic camera with the spectral sensitivity of the Green
% channel. The resulting images are then concatenated to form a hyperspectral
% image.
%
% #### Demosaiced hyperspectral images
%
% '*_dHyper.mat' contains a spectral radiance image, stored in the variable
% 'I_hyper'. This image is created from the demosaiced Green channels of the
% spectral band images taken for the scene. Demosaicing is an alternative method
% to the channel scaling approach for creating images from a hypothetical
% monochromatic camera with the spectral sensitivity of the Green channel. While
% the demosaicing algorithm introduces some error, the resulting images should
% be less noisy.
%
% #### Demosaiced RGB images
%
% '*_d3.tif' and '*_d3.mat' contain a 3-channel colour image (stored in the
% variable 'I_3'). This image is created by converting the demosaiced
% hyperspectral image (described above) to colour, according to the
% relative spectral sensitivities of the sensor (described above).
%
% ### Parameters and intermediate variables
%
% The 'sensor.mat' file also contains the values of all parameters in the
% first section of the script below, for reference. (Specifically, those
% listed in `parameters_list`, which should be updated if the set of
% parameters is changed.)
%
% Additionally, there are some intermediate variables in the file:
% - 'grouped_input_filenames': A cell vector of cell vectors of input image
%   filenames retrieved based on the wildcard provided in the parameters
%   section of the script. Each of the inner cell vectors groups the images
%   of the same scene taken for the different spectral bands.
% - 'grouped_reference_filenames': A cell vector of cell vectors of calibration
%   image filenames retrieved based on the wildcard provided in the parameters
%   section of the script. Each of the inner cell vectors groups images taken
%   for the same spectral band.
% - 'reference_index': For reference, the index of the colour channel (Green)
%   with respect to which the relative spectral sensitivities were computed.
%
% ## Notes
% - The approach for computing relative spectral sensitivities
%   ('relativeSensitivity()') works well only in the event that all colour
%   channels have non-negligible sensitivities in all spectral bands. If some
%   colour channels have zero spectral sensitivity in certain spectral bands, a
%   better approach may be to compute hyperspectral images without using
%   spectral sensitivities. For example, each spectral band of a hyperspectral
%   image could be set to an interpolated version of the most well-exposed
%   colour channel for that spectral band. Unfortunately, this new approach
%   would introduce error from demosaicing.
% - Output images are clipped to zero.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 12, 2018

% List of parameters to save with results
parameters_list = {
        'input_bands_regex',...
        'reference_bands_regex',...
        'range',...
        'use_bias',...
        'bayer_pattern',...
        'quantiles'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process. All images are
% expected to be in one directory.
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/preprocessed_images/exposure_blended/*nm.mat';
input_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the images used to calibrate the relative spectral
% sensitivity. All images are expected to be in one directory.
reference_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/preprocessed_images/exposure_blended/*nm.mat';
reference_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Wavelengths will be expected within filenames, extracted using this
% regular expression
input_bands_regex = '_(\d+)nm'; % In input images
reference_bands_regex = input_bands_regex; % In calibration images

% Range of pixel values used to calibrate scaling factors between colour
% channels
range = [0, 0.95];

% Convert between colour channels using affine (`true`) or proportional
% (`false`) linear models
use_bias = false;

% ## Output directory
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling';

% ## Colour conversion parameters
run('SetFixedParameters.m')

% Colour-filter pattern (overrides the version in 'SetFixedParameters.m')
bayer_pattern = 'gbrg';

% Quantiles used for clipping to produce nice output images (for display, not
% for calculation)
quantiles = [0.01, 0.99];

% ## Debugging Flags

relativeSensitivityVerbose = true;
% If `true`, image format files will be output that contain rescaled versions of
% the spectral and colour images suitable for viewing.
outputScaledImages = true;

%% Find the images

[...
    grouped_input_filenames, ~, group_names, bands...
] = findAndGroupImages(input_images_wildcard, input_bands_regex);
n_groups = length(grouped_input_filenames);
n_bands = length(bands);

[...
    filenames, ~, ~, bands_reference...
] = findAndGroupImages(reference_images_wildcard, reference_bands_regex);
if n_bands ~= length(bands_reference) || ~all(bands == bands_reference)
    error('The input and calibration files were captured under different bands.');
end
grouped_reference_filenames = repmat({cell(length(filenames), 1)}, n_bands, 1);
for g = 1:length(filenames)
    for i = 1:n_bands
        grouped_reference_filenames{i}{g} = filenames{g}{i};
    end
end

%% Relative spectral sensitivity calibration

n_channels_rgb = 3;
channel_mode = false;
[sensor_map_relative, channel_bias] = relativeSensitivity(...
    reference_images_variable_name, grouped_reference_filenames,...
    range, bayer_pattern, relativeSensitivityVerbose...
);
% Use the Green colour channel as the reference spectral sensitivity
% Note that the choice of the Green channel is hardcoded in 'relativeSensitivity()'.
sensor_map_relative_sum = sum(sensor_map_relative, 2);
reference_index = find(sensor_map_relative_sum == n_bands);

color_weights_options = findSamplingOptions;
color_weights_options.interpolant = findSamplingOptions.interpolant_ref;
color_weights = colorWeights(...
    sensor_map_relative, bands, bands, color_weights_options...
);

color_weights_sum = sum(color_weights(reference_index, :));
color_weights = color_weights ./ color_weights_sum;
sensor_map = sensor_map_relative ./ color_weights_sum;

%% Visualization

figure;
hold on
plot(bands, sensor_map_relative(1, :), 'r');
plot(bands, sensor_map_relative(2, :), 'g');
plot(bands, sensor_map_relative(3, :), 'b');
hold off
xlabel('Wavelength [nm]')
ylabel('Spectral sensitivity relative to Green')
title('Relative spectral sensitivity')
legend('Red', 'Green', 'Blue')

figure;
hold on
plot(bands, sensor_map(1, :), 'r');
plot(bands, sensor_map(2, :), 'g');
plot(bands, sensor_map(3, :), 'b');
hold off
xlabel('Wavelength [nm]')
ylabel('Spectral sensitivity')
title('Normalized relative spectral sensitivity')
legend('Red', 'Green', 'Blue')

%% Process the images

demosaic_channels = false(n_channels_rgb, 1);
demosaic_channels(reference_index) = true;

for g = 1:n_groups
    for i = 1:n_bands
        I = loadImage(grouped_input_filenames{g}{i}, input_images_variable_name);
        if size(I, 3) ~= 1
            error('The input image "%s" is not a RAW image.', grouped_input_filenames{g}{i});
        end
        
        if i == 1
            image_height = size(I, 1);
            image_width = size(I, 2);
            channel_mask = bayerMask(image_height, image_width, bayer_pattern);
            I_hyper_q = zeros(image_height, image_width, n_bands);
            I_hyper_d = zeros(image_height, image_width, n_bands);
        end
        
        I_d = bilinearDemosaic(I, bayer_pattern, demosaic_channels);
        I_hyper_d(:, :, i) = I_d;
            
        for c = 1:n_channels_rgb
            if use_bias
                I(channel_mask(:, :, c)) =...
                    (I(channel_mask(:, :, c)) - channel_bias(c, i)) ./ ...
                    sensor_map_relative(c, i);
            else
                I(channel_mask(:, :, c)) = I(channel_mask(:, :, c)) ./ ...
                    sensor_map_relative(c, i);
            end
        end
        I_hyper_q(:, :, i) = I;
    end
    
    % Convert to colour
    [I_3_q, ~, I_raw] = imageFormation(...
        I_hyper_q, bands, sensor_map, bands,...
        imageFormationSamplingOptions, imageFormationPatchOptions,...
        [], bayer_pattern...
    );
    I_3_d = imageFormation(...
        I_hyper_d, bands, sensor_map, bands,...
        imageFormationSamplingOptions, imageFormationPatchOptions...
    );

    % Clip to zero
    I_hyper_q(I_hyper_q < 0) = 0;
    I_hyper_d(I_hyper_d < 0) = 0;
    I_3_q(I_3_q < 0) = 0;
    I_3_d(I_3_d < 0) = 0;
    I_raw(I_raw < 0) = 0;

    % Save the results
    saveImages(...
        'data', output_directory, group_names{g},...
        I_hyper_q, '_qHyper', 'I_hyper',...
        I_hyper_d, '_dHyper', 'I_hyper',...
        I_3_q, '_q3', 'I_3',...
        I_3_d, '_d3', 'I_3',...
        I_raw, '_raw', 'I_raw'...
    );

    if outputScaledImages
        type_symbols = {'q', 'd'};
        for t = 1:2
            if t == 1
                I_hyper_t = I_hyper_q;
                I_3_t = I_3_q;
            elseif t == 2
                I_hyper_t = I_hyper_d;
                I_3_t = I_3_d;
            else
                error('Unrecognized value %d of `t`.', t);
            end
            
            for i = 1:n_bands
                I_out_debug = clipAndRemap(I_hyper_t(:, :, i), 'uint8', 'quantiles', quantiles);
                saveImages(...
                    'image', output_directory, group_names{g},...
                    I_out_debug, sprintf('_%sHyper_%dnm', type_symbols{t}, bands(i)), []...
                );
            end
                        
            I_out_debug = clipAndRemap(I_3_t, 'uint8', 'quantiles', quantiles);
            saveImages(...
                'image', output_directory, group_names{g},...
                I_out_debug, sprintf('_%s3_01', type_symbols{t}), []...
            );
        end
        
        I_out_debug = clipAndRemap(I_raw, 'uint8', 'quantiles', quantiles);
        saveImages(...
            'image', output_directory, group_names{g},...
            I_out_debug, '_raw01', []...
        );
    end
end

%% Save parameters, intermediate variables, and output data
save_variables_list = [ parameters_list, {...
        'grouped_input_filenames',...
        'grouped_reference_filenames',...
        'bands',...
        'reference_index'...
        'channel_mode',...
        'sensor_map_relative',...
        'channel_bias',...
        'sensor_map'...
    } ];
save_data_filename = fullfile(output_directory, 'sensor.mat');
save(save_data_filename, save_variables_list{:});