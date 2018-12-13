%% Assemble images for spectral bands into spectral and colour images
% Load images taken under bandpass optical filters and use them to create a
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
% ## Output
%
% ### Relative spectral sensitivities, and wavelength information
%
% A 'sensor.mat' file containing the following variables:
% - 'channel_mode': A Boolean value set to `false` to indicate that the
%   data in `sensor_map` represents spectral sensitivities, not colour
%   channel mappings
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel (Red, Green, or Blue) to light of the j-th
%   wavelength (in the j-th spectral band), relative to the sensitivity of
%   the Green channel. Consequently, `sensor_map(2, :)` is a vector of
%   ones.
% - 'bands': A vector containing the wavelengths corresponding to the
%   columns of `sensor_map`.
% - 'color_weights': A matrix for converting pixels in the hyperspectral
%   images to colour, as determined by 'sensor_map', and by the type of
%   numerical intergration to perform (set in 'SetFixedParameters.m').
%
% The relative spectral sensitivities of the Red and Blue channels with
% respect to the Green channel are computed, for each spectral band, as the
% ratios of the average colour channel values in the images captured for
% that spectral band. The averaging is done over all scenes, but is
% optionally done only over a specific window within the images.
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
% divides, and then multiplies, by the spectral sensitivities.
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
% 'I_hyper'. This image is created from the spectral band images taken for
% the scene. Each spectral band image's Red and Blue pixels are divided by
% the relative spectral sensitivities of these colour channels for the
% spectral band, to recover an image from a hypothetical monochromatic
% camera with the spectral sensitivity of the Green channel. The resulting
% images are then concatenated to form a hyperspectral image.
%
% #### Demosaiced hyperspectral images
%
% '*_dHyper.mat' contains a spectral radiance image, stored in the variable
% 'I_hyper'. This image is created from the demosaiced green channels of
% the spectral band images taken for the scene. Demosaicing is an
% alternative method to the channel scaling approach for creating images
% from a hypothetical monochromatic camera with the spectral sensitivity of
% the Green channel. While the demosaicing algorithm introduces some error,
% the resulting images should be less noisy.
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
% - 'grouped_filenames': A cell vector of cell vectors of input image
%   filenames retrieved based on the wildcard provided in the parameters
%   section of the script. Each of the inner cell vectors groups the images
%   of the same scene taken for the different spectral bands.
%
% ## Notes
% - The approach for computing relative spectral sensitivities works well
%   only in the event that all colour channels have non-negligible
%   sensitivities in all spectral bands. If some colour channels have zero
%   spectral sensitivity in certain spectral bands, a better approach may
%   be to compute hyperspectral images without using spectral
%   sensitivities. For example, each spectral band of a hyperspectral image
%   could be set to an interpolated version of the most well-exposed colour
%   channel for that spectral band. Unfortunately, this new approach would
%   introduce error from demosaicing.
% - Input images are clipped to zero when computing spectral sensitivities,
%   and also when generating output images.
%
% ## References
%
% The approach taken in this script for computing relative spectral
% sensitivities is inspired by:
%
%   Darrodi, M. M., Finlayson, G., Goodman, T., & Mackiewicz, M. (2015).
%   Reference data set for camera spectral sensitivity estimation. Journal
%   of the Optical Society of America A: Optics and Image Science, and
%   Vision, 32(3), 381-391. doi:10.1364/JOSAA.32.000381

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 12, 2018

% List of parameters to save with results
parameters_list = {
        'crop',...
        'bands_regex',...
        'bayer_pattern',...
        'reference_index'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process. All images are
% expected to be in one directory.
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/preprocessed/blended_averaged/pos1_1mmDots_*.mat';
input_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% A crop box for the region to use for relative spectral sensitivity
% estimation (leave empty for no cropping)
crop = []; % [top-left x, top-left y, width, height]

% Wavelengths will be expected within filenames, extracted using this
% regular expression
bands_regex = '_(\d+)nm';

% Colour-filter pattern
bayer_pattern = 'gbrg';

% Use the Green colour channel as the reference spectral sensitivity
reference_index = 2;

% ## Output directory
output_directory = '/home/llanos/Downloads';

% ## Colour conversion parameters
run('SetFixedParameters.m')

%% Find the images

[...
    grouped_filenames, ~, group_names, bands...
] = findAndGroupImages(input_images_wildcard, bands_regex);
n_groups = length(grouped_filenames);
n_bands = length(bands);

%% Relative spectral sensitivity calibration

n_channels_rgb = 3;
px_sums = zeros(n_channels_rgb, n_bands);

image_size = [];
for g = 1:n_groups
    for i = 1:n_bands
        I = loadImage(grouped_filenames{g}{i}, input_images_variable_name);
        if isempty(image_size)
            % First time setup
            image_size = size(I);
            image_height = image_size(1);
            image_width = image_size(2);
            channel_mask = bayerMask(image_height, image_width, bayer_pattern);
            if isempty(crop)
                channel_mask_cropped = channel_mask;
            else
                if crop(2) < 1 || crop(2) > image_height ||...
                   crop(1) < 1 || crop(1) > image_width
                    error('The top-left corner of the crop region is outside the image.');
                end
                crop_end = crop(1:2) + crop(3:4) - 1;
                if crop_end(2) < 1 || crop_end(2) > image_height ||...
                   crop_end(1) < 1 || crop_end(1) > image_width
                    error('The bottom-right corner of the crop region is outside the image.');
                end
                channel_mask_cropped = channel_mask_cropped(crop(2):crop_end(2), crop(1):crop_end(1), :);
            end
            
            channel_counts = reshape(sum(sum(channel_mask_cropped, 1), 2), [], 1);
        elseif any(image_size ~= size(I))
            error('Not all images have the same dimensions.');
        end
        
        if isempty(crop)
            patch = I;
        else
            patch = I(crop(2):crop_end(2), crop(1):crop_end(1), :);
        end
        patch(patch < 0) = 0;
        for c = 1:n_channels_rgb
            px_sums(c, i) = px_sums(c, i) + sum(patch(channel_mask_cropped(:, :, c)));
        end
    end
end

channel_mode = false;
sensor_map = px_sums ./ repmat(channel_counts * n_groups * n_bands, 1, n_bands);
sensor_map = sensor_map ./ repmat(sensor_map(reference_index, :), n_channels_rgb, 1);

samplingWeightsOptions.power_threshold = 1;
samplingWeightsOptions.n_bands = 0;
samplingWeightsOptions.interpolant = @triangle;
[...
    ~, ~, ~, color_weights...
] = samplingWeights(...
  sensor_map, bands, bands, samplingWeightsOptions...
);

%% Visualization

figure;
hold on
plot(bands, sensor_map(1, :), 'r');
plot(bands, sensor_map(2, :), 'g');
plot(bands, sensor_map(3, :), 'b');
hold off
xlabel('Wavelength [nm]')
ylabel('Spectral sensitivity relative to Green')
title('Relative spectral sensitivity computed from bandpass-filtered images')
legend('Red', 'Green', 'Blue')

%% Process the images

imageFormationOptions.patch_size = patch_sizes(1, :);
imageFormationOptions.padding = paddings(1);
demosaic_channels = false(n_channels_rgb, 1);
demosaic_channels(reference_index) = true;

for g = 1:n_groups
    I_hyper_q = zeros(image_height, image_width, n_bands);
    I_hyper_d = zeros(image_height, image_width, n_bands);
    for i = 1:n_bands
        I = loadImage(grouped_filenames{g}{i}, input_images_variable_name);
        I(I < 0) = 0;
        
        I_d = bilinearDemosaic(I, bayer_pattern, demosaic_channels);
        I_hyper_d(:, :, i) = I_d;
            
        for c = 1:n_channels_rgb
            I(channel_mask(:, :, c)) = I(channel_mask(:, :, c)) ./ sensor_map(c, i);
        end
        I_hyper_q(:, :, i) = I;
    end
    
    % Convert to colour
    [I_3_q, ~, I_raw] = imageFormation(...
        I_hyper_q, color_weights, imageFormationOptions,...
        [], bands, bayer_pattern...
    );
    I_3_d = imageFormation(...
        I_hyper_d, color_weights, imageFormationOptions...
    );

    % Save the results
    saveImages(...
        'data', output_directory, group_names{g},...
        I_hyper_q, '_qHyper', 'I_hyper',...
        I_hyper_d, '_dHyper', 'I_hyper',...
        I_3_q, '_q3', 'I_3',...
        I_3_d, '_d3', 'I_3',...
        I_raw, '_raw', 'I_raw'...
    );
end

%% Save parameters, intermediate variables, and output data
save_variables_list = [ parameters_list, {...
        'grouped_filenames',...
        'bands',...
        'channel_mode',...
        'sensor_map',...
        'color_weights',...
    } ];
save_data_filename = fullfile(output_directory, 'sensor.mat');
save(save_data_filename, save_variables_list{:});