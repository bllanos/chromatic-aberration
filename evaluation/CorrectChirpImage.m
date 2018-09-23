%% Evaluate hyperspectral ADMM-based reconstruction of a chirp image
% This script investigates the following questions:
% - How does reconstruction accuracy change with the spatial and spectral
%   frequencies of the image?
% - How does reconstruction accuracy change with the signal-to-noise ratio?
% - How do the spatial and spectral frequencies of the image affect the
%   choice of regularization weights?
% - How is reconstruction accuracy affected by the patch size?
% - How much overlap should patches have?
% - How does the amount of dispersion affect reconstruction accuracy?
% - How fast does image reconstruction finish with different patch and
%   padding sizes?
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Colour space conversion data
% A '.mat' file containing several variables, which is the output of
% 'SonyColorMap.m', for example. The following variables are required:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel or spectral band in the input images to the
%   j-th colour channel or spectral band of the latent images. For example,
%   `sensor_map` is a matrix mapping discretized spectral power
%   distributions to RGB colours.
% - 'channel_mode': A Boolean value indicating whether the latent colour
%   space is a set of colour channels (true) or a set of spectral bands
%   (false). 'channel_mode' must be false.
% - 'bands': A vector containing the wavelengths or colour channel indices
%   corresponding to the second dimension of 'sensor_map'. 'bands' is
%   required to resample 'sensor_map' so that it maps from the colour space
%   of the latent image to the colour space of the input RAW image.
%
% ## Output
%
% ### Graphical output
%
% Figures are created to answer the questions given above, and are saved to
% MATLAB '.fig' files.
%
% ### Images
%
% The following types of images are created for each set of parameters used
% to create the test chirp image, and to control image estimation. True
% images are created for each noise and dispersion setting, whereas
% estimated images are created for each noise setting, dispersion setting,
% patch size, and patch padding size combination. Therefore, true images
% can be distinguished from estimated images by the parameter information
% in their filenames. The string of parameter information is represented by
% '*' below:
% - '*_latent.tif' and '*_latent.mat': The spectral image, stored in the
%   variable 'I_latent'. The '.tif' image is only output if the spectral
%   image is a greyscale or 3-channel image.
% - '*_warped.tif' and '*_warped.mat': A version of the spectral image
%   (stored in the variable 'I_warped') created by warping the spectral
%   image according to the dispersion model. The '.tif' image is only
%   output if if the spectral image is a greyscale or 3-channel image.
% - '*_rgb.tif' and '*_rgb.mat': A colour image (stored in the variable
%   'I_rgb') created by converting the spectral image to the RGB colour
%   space of the camera.
% - '*_rgb_warped.tif' and '*_rgb_warped.mat': A colour image (stored in
%   the variable 'I_full') created by warping the spectral image according
%   to the dispersion model, then converting the image to the RGB colour
%   space of the camera.
% - '*_raw_noNoise.tif' and '*_raw_noNoise.mat': A simulation of The RAW
%   image corresponding to the spectral image, stored in the variable
%   'I_raw_noNoise'.
% - '*_raw.tif' and '*_raw.mat': A simulation of The RAW image
%   corresponding to the spectral image, with added noise, stored in the
%   variable 'I_raw'.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'bands': The value of the 'bands' variable defined in
%   'SetFixedParameters.m', which determines the spectral resolution of the
%   image. This script requires 'bands' to contain equally-spaced
%   wavelengths.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference.
% - 'sensor_map_resampled': The resampled version of the 'sensor_map'
%   variable, generated for compatibility with the true latent image.
%   'sensor_map' may have been normalized so that colour images will have
%   values in the appropriate range, as set in the parameters section
%   below.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## Notes
% - This script uses the first row of `weights` defined in
%   'SetFixedParameters.m' to determine which regularization weights to
%   set. Elements of `weights(1, :)` set to zero disable the corresponding
%   regularization terms.
% - This script does not estimate downsampled images, and so it ignores
%   `downsampling_factor` in 'SetFixedParameters.m'.
% - The body of this script is largely copied from 'solvePatchesAligned()'.
%
% ## References
% This experiment is inspired by Figure 20 of:
%
%   Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%     single-shot hyperspectral imaging using a prism." ACM Transactions
%     on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%     doi:10.1145/3130800.3130896
%
% The method for automatically selecting regularization weights without
% knowledge of the true image is from:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 21, 2018

% List of parameters to save with results
parameters_list = {
    'color_map_filename',...
    'output_directory'...
};

%% Input data and parameters

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180923_TestingChirpImageGeneration/CIE1931ColorMapData.mat';

% Whether or not to normalize spectral sensitivity functions, assuming an
% illuminant which is uniform across the spectrum. The normalization code
% is based on 'reflectanceToRadiance()'
normalize_color_map = true;
% Which spectral sensitivity function to use for computing the
% normalization
normalize_color_map_reference_index = 2; % Use '2' for the CIE tristimulus functions (since 2 denotes Y) 

% Image dimensions (height, width)
image_sampling = [128, 128];

% Number of samples for antialiasing during image generation
n_samples = 1;

% Number of patch sizes to test. Patches will be square.
n_patch_sizes = 5;

% Number of patch padding sizes to test
n_padding = 4;

% Size of the largest padding size as a multiple of the maximum dispersion
padding_ratio_max = 2;

% Number of dispersion magnitudes to test
n_dispersion = 3;

% Noise fractions
noise_fractions = [0, 0.05, 0.1, 0.25, 0.5];

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% Tweaks to the parameters set by 'SetFixedParameters.m'
selectWeightsGridOptions.parallel = false;
trainWeightsOptions.parallel = false;
% Not all patches will have full borders, so do the following for
% simplicity
baek2017Algorithm2Options.l_err_border = [0, 0];
trainWeightsOptions.border = 0;

% Debugging flags
baek2017Algorithm2Verbose = false;
verbose_progress = true;

%% Validate parameters and construct intermediate parameters

if ~isempty(downsampling_factor)
    if downsampling_factor ~= 1
        warning('`downsampling_factor` is ignored.');
    end
end

diff_bands = diff(bands);
if any(diff_bands ~= diff_bands(1))
    error('Expected `bands` to contain equally-spaced wavelengths.');
end
n_bands = length(bands);
delta_lambda = (bands(end) - bands(1)) / (n_bands - 1);
lambda_range = [bands(1) - delta_lambda / 2, bands(end) + delta_lambda / 2];
image_sampling_3 = [image_sampling, n_bands];

min_spatial_size = min(image_sampling);
patch_sizes = repmat(...
    round(logspace(0, log10(min_spatial_size), n_patch_sizes)), 1, 2 ...
);

[~, dispersion_max] = chirpImage(...
  image_sampling_3, lambda_range, 1, [], 'params'...
);
paddings = [0, round(logspace(0, log10(padding_ratio_max * dispersion_max), n_padding))];

dispersion_fractions = linspace(0, 1, n_dispersion);

imageFormationOptions.patch_size = [100, 100];
imageFormationOptions.padding = 0; % No need to accommodate dispersion

%% Load calibration data

bands_script = bands;
bands = [];
model_variables_required = { 'sensor_map', 'channel_mode', 'bands' };
load(color_map_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end
if isempty(bands)
    error('No (non-empty) variable `bands` loaded from colour space conversion data.');
end
if channel_mode
    error('The input space of the colour conversion data must be a spectral space, not a space of colour channels.')
end
bands_color = bands;
bands = bands_script;

% Normalize spectral sensitivities to avoid out-of-range colour values.
if normalize_color_map
    ybar = sensor_map(normalize_color_map_reference_index, :);
    weights = integrationWeights(bands_color, int_method);
    N = dot(ybar, weights);
    sensor_map = sensor_map / N;
end

baek2017Algorithm2Options.int_method = int_method;
selectWeightsGridOptions.int_method = int_method;
trainWeightsOptions.int_method = int_method;
imageFormationOptions.int_method = int_method;

% Resample colour space conversion data if necessary
if n_bands ~= length(bands_color) || any(bands ~= bands_color)
    [sensor_map_resampled, bands] = resampleArrays(...
        bands_color, sensor_map.', bands,...
        bands_interp_method...
        );
    if n_bands ~= length(bands)
        error('The colour space conversion data does not cover a sufficiently large range of wavelengths.');
    end
    sensor_map_resampled = sensor_map_resampled.';
else
    sensor_map_resampled = sensor_map;
end

%% Collect results

baek2017Algorithm2Options.add_border = false;
baek2017Algorithm2Options.l_surface = true;
n_auxiliary_images = 4;
auxiliary_images = cell(n_auxiliary_images, 1);
n_noise_fractions = length(noise_fractions);

for d = 1:n_dispersion
    dispersion_fraction = dispersion_fractions(d);

    % Synthesize the true image
    [I_spectral_gt, I_warped_gt, dispersionfun] = chirpImage(...
        image_sampling_3, lambda_range, dispersion_fraction, n_samples...
    );
    I_rgb_gt = imageFormation(...
        I_spectral_gt, sensor_map_resampled, bands,...
        imageFormationOptions...
    );
    I_rgb_warped_gt = imageFormation(...
        I_warped_gt, sensor_map_resampled, bands,...
        imageFormationOptions...
    );
    I_raw_noNoise_gt = mosaic(I_rgb_warped_gt, bayer_pattern);
    
    name_params_gt = sprintf('dispersion%e_', dispersion_fraction);
    saveImages(...
        output_directory, name_params_gt,...
        I_raw_noNoise_gt, 'raw_noNoise', 'I_raw_noNoise',...
        I_spectral_gt, 'latent', 'I_latent',...
        I_rgb_gt, 'rgb', 'I_rgb',...
        I_rgb_warped_gt, 'rgb_warped', 'I_full',...
        I_warped_gt, 'warped', 'I_warped'...
    );
        
    for no = 1:n_noise_fractions
        noise_fraction = noise_fractions(no);
        snr = noiseFractionToSNR(noise_fraction);
        I_raw_gt = addNoise(I_raw_noNoise_gt, snr);
                
        name_params_gt = sprintf(...
            'dispersion%e_noise%e_', dispersion_fraction, noise_fraction...
        );
        saveImages(...
            output_directory, name_params_gt,...
            I_raw_gt, 'raw', 'I_raw'...
        );
                
        has_dispersion = ~isempty(dispersionfun);
        
        for ps = 1:n_patch_sizes
            patch_size = patch_sizes(ps, :);
            for pp = 1:n_padding
                padding = paddings(pp);
                
                tic
                 
                if verbose_progress
                    disp('Splitting the input image into patches...');
                end
                i_vector = 1:patch_size(1):image_sampling_3(1);
                n_i_vector = length(i_vector);
                j_vector = 1:patch_size(2):image_sampling_3(2);
                n_j_vector = length(j_vector);
                n_patches = n_i_vector * n_j_vector;
                patches_J = cell(1, n_j_vector);
                patches_I = cell(1, n_j_vector);
                patches_auxiliary = cell(n_i_vector, n_j_vector, n_auxiliary_images);
                patch_limits = zeros(n_i_vector, n_j_vector, 4);
                patch_trim = zeros(n_i_vector, n_j_vector, 4);
                corners = zeros(n_i_vector, n_j_vector, 2);

                % Divide the input image into patches to be sent to individual parallel
                % workers
                for j = 1:n_j_vector
                    for i = 1:n_i_vector
                        corners(i, j, :) = [i_vector(i), j_vector(j)];
                        [ patch_lim, trim ] = patchBoundaries(...
                            image_sampling, patch_size, options.padding, corners(i, j, :)...
                        );
                        patch_trim(i, j, :) = reshape(trim, 1, 1, 4);
                        patch_limits(i, j, :) = reshape(patch_lim, 1, 1, 4);

                        if i == 1
                            patches_J{j} = J(:, patch_lim(1, 2):patch_lim(2, 2), :);
                        end
                    end
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('Parallel processing of patches...');
                end

                % Process each patch
                parfor j = 1:n_j_vector
                    patch_limits_j = patch_limits(:, j, :);
                    patch_trim_j = patch_trim(:, j, :);
                    corners_j = corners(:, j, :);
                    patches_J_j = patches_J{j};
                    patches_I_j = zeros(...
                        size(patches_J_j, 1),...
                        patch_trim_j(1, 1, 4) - patch_trim_j(1, 1, 3) + 1,...
                        n_bands...
                    );
                    patches_auxiliary_j = cell(n_i_vector, 1, n_auxiliary_images);
                    for i = 1:n_i_vector
                        patch_lim = reshape(patch_limits_j(i, 1, :), 2, 2);
                        if isempty(bayer_pattern)
                            align_f = [];
                        else
                            align_f = offsetBayerPattern(patch_lim(1, :), bayer_pattern);
                        end
                        image_sampling_f = diff(patch_lim, 1, 1) + 1;
                        trim = reshape(patch_trim_j(i, 1, :), 2, 2);
                        
                        if has_dispersion
                            dispersion_matrix_patch = dispersionfunToMatrix(...
                                dispersionfun, bands, image_sampling_f, image_sampling_f,...
                                [0, 0, image_sampling_f(2), image_sampling_f(1)], true,...
                                [corners_j(i, 1, 2), corners_j(i, 1, 1)] - 1 ...
                            );
                        else
                            dispersion_matrix_patch = [];
                        end

                        % Process the patch
                        % Most of the options to selectWeightsGrid() are set in 'SetFixedParameters.m'
                        [ ~, ~, patches_I_ij ] = selectWeightsGrid(...
                            patches_J_j(patch_lim(1, 1):patch_lim(2, 1), :, :),...
                            align_f, dispersion_matrix_patch, sensor_map_resampled, bands,...
                            rho, baek2017Algorithm2Options, selectWeightsGridOptions,...
                            corners_j(i, 1, :), selectWeightsGridVerbose...
                        );

                        padding_filter = false(image_sampling_f);
                        padding_filter((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2))) = true;
                        patch_trimmed_size = diff(trim, 1, 1) + 1;
                        [patches_I_ij, patches_auxiliary_j(i, 1, :)] = estimateAuxiliaryImages(...
                                patches_I_ij, dispersion_matrix_patch, padding_filter,...
                                patch_trimmed_size,...
                                sensor_map_resampled, bands, int_method,...
                                align_f, n_auxiliary_images...
                        );
                        patches_I_j(...
                                ((i - 1) * patch_size(1) + 1):((i - 1) * patch_size(1) + patch_trimmed_size(1)),...
                                :, :...
                            ) = patches_I_ij;
                        if verbose_progress
                            fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i_vector, n_patches);
                        end
                    end
                    patches_I{j} = patches_I_j;
                    patches_auxiliary(:, j, :) = patches_auxiliary_j;
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('Recombining results from patches...');
                end

                % Recombine patches
                I_3D = cell2mat(patches_I);
                for im = 1:n_auxiliary_images
                    auxiliary_images{im} = cell2mat(patches_auxiliary(:, :, im));
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    toc
                end

                % Save the results
                name_params = sprintf(...
                    'noise%e_dispersion%e_patch%dx%d_pad%d_',...
                    noise_fraction, dispersion_fraction, patch_size(1),...
                    patch_size(2), padding...
                );
            
                saveImages(...
                    output_directory, name_params,...
                    I_3D, 'latent', 'I_latent',...
                    auxiliary_images{1}, 'rgb', 'I_rgb',...
                    auxiliary_images{2}, 'rgb_warped', 'I_full',...
                    auxiliary_images{3}, 'reestimated', 'I_raw',...
                    auxiliary_images{4}, 'warped', 'I_warped'...
                );
            end
        end
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'bands_color',...
        'bands',...
        'sensor_map_resampled'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectChirpImage.mat');
save(save_data_filename, save_variables_list{:});