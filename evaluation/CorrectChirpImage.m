%% Evaluate hyperspectral ADMM-based reconstruction of a chirp image
%
% This script investigates the following questions:
% - How does image reconstruction accuracy change with the spatial and
%   spectral frequencies of the image?
% - How does image reconstruction accuracy change with the signal-to-noise
%   ratio?
% - How do the spatial and spectral frequencies of the image affect the
%   choice of regularization weights?
% - How do the regularization weights chosen by 'selectWeightsGrid()' and
%   'trainWeights()' differ, and how does this difference affect image
%   reconstruction time and accuracy?
% - How is image reconstruction accuracy affected by the patch size?
% - How much overlap should patches have?
% - How does the amount of dispersion affect image reconstruction accuracy?
% - How long does image reconstruction take with different patch and
%   padding sizes?
%
% Presently, this script evaluates two image estimation methods:
% - Regularization weight selection, per image patch, using
%   'selectWeightsGrid()', followed by image estimation by an ADMM-family
%   algorithm, using the selected per-patch regularization weights.
% - Regularization weight selection, per image patch, using
%   'trainWeights()', followed by image estimation by an ADMM-family
%   algorithm, using the selected per-patch regularization weights.
%
% The same ADMM-family algorithm is used for both methods:
% 'baek2017Algorithm2()', configured according to the parameters in
% 'SetFixedParameters.m'. Note that 'baek2017Algorithm2()' is embedded in
% 'selectWeightsGrid()', but explicitly passed to 'trainWeights()' as an
% input argument. If 'selectWeightsGrid()' is updated such that it runs a
% different image estimation algorithm, this script would need to be
% modified to ensure the two methods remain comparable.
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
% Additionally, some of the per-pixel raw data used to create figures is
% saved to separate '.mat' files. The string of parameter information
% describing the context of the data is represented by '*' below:
% - '*_weight${i}Image.mat': The 'I_weights' variable is a 2D image-sized
%   array, where each element stores the i-th regularization weight for the
%   patch covering the corresponding pixel location.
% - '*_penalty${i}.mat': The 'I_penalty' variable is a 2D image-sized
%   array, where each element stores the i-th regularization penalty for
%   the patch covering the corresponding pixel location. Note that the
%   penalty is evaluated on the true image, but the patch dimensions are a
%   parameter of the image estimation process, so these images are specific
%   to the image estimation parameters.
%
% ### Images
%
% The following types of images are created for each set of parameters used
% to create the test chirp image, and to control image estimation. True
% images are created for each noise and dispersion setting, whereas
% estimated images are created for each noise setting, dispersion setting,
% patch size, and patch padding size combination. Therefore, true images
% can be distinguished from estimated images by the amount of parameter
% information in their filenames. The string of parameter information is
% represented by '*' below:
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
% #### Parameters, and execution time results
%
% A '.mat' file containing the following variables:
%
% - 'bands': The value of the 'bands' variable defined in
%   'SetFixedParameters.m', which determines the spectral resolution of the
%   image. This script requires 'bands' to contain equally-spaced
%   wavelengths.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference.
% - 'dataset_params': A structure of the form of the `dataset_params`
%   output argument of 'describeDataset()'. 'dataset_params' allows the
%   chirp images created by this script to be used as a dataset for
%   evaluating other image estimation algorithms. Presently,
%   'dataset_params' is missing dispersion information, as dispersion is
%   not the same between images, and so cannot be described with a single
%   dispersion model.
% - 'sensor_map_resampled': The resampled version of the 'sensor_map'
%   variable, generated for compatibility with the true latent image.
%   'sensor_map' may have been normalized so that colour images will have
%   values in the appropriate range, as configured in the parameters
%   section below.
% - 'selectWeightsGrid_time': Execution timing information, stored as a 4D
%   array, for image estimation and regularization weights selection by
%   'selectWeightsGrid()'. `selectWeightsGrid_time(pp, ps, no, d)` is the
%   time taken with the pp-th patch padding size, the ps-th patch size, the
%   no-th noise fraction in the input image, and the d-th level of
%   dispersion in the input image. Time values include the time taken to
%   generate output images beyond the estimated spectral image, including
%   visualizations of selected regularization weights, and regularization
%   penalty values. (In other words, all of the types of images listed
%   above which are dependent on the image estimation algorithm's
%   parameters.) Time taken to save images to disk, and to generate figures
%   from images, is excluded.
% - 'trainWeightsGrid_time': An array with the same format and
%   interpretation as 'selectWeightsGrid_time', containing execution timing
%   information for image estimation and regularization weights selection
%   by 'trainWeights()'.
% - 'n_patches': A vector containing the number of patches resulting from
%   each patch size used in the experiment.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% #### Evaluation results
%
% For each synthesized image, RGB error metrics and spectral error metrics
% are output in the form of CSV files. Each CSV file contains results for
% all image estimation methods. The RGB error metrics are saved as
% '*_evaluateRGB.csv', whereas the spectral error metrics are saved as
% '*_evaluateSpectral.csv'.
%
% Error metrics are also aggregated across images, and saved as
% 'chirp_evaluateRGB.csv' and 'chirp_evaluateSpectral.csv'. The aggregation
% may not be so meaningful, because the synthetic images differ in
% dispersion magnitude and signal-to-noise ratio.
%
% ## Notes
% - This script uses the first row of `weights` defined in
%   'SetFixedParameters.m' to determine which regularization weights to
%   set. Elements of `weights(1, :)` set to zero disable the corresponding
%   regularization terms.
% - This script does not estimate downsampled images, and so it ignores
%   `downsampling_factor` in 'SetFixedParameters.m'.
% - The image estimation portion of this script is largely copied from
%   'solvePatchesAligned()'.
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
%
% Dispersion magnitudes are set based on the following article, which
% states that 0.1 pixels is the threshold at which dispersion becomes
% perceivable, whereas dispersion of 0.3 pixels and above is noticeable:
%
%   V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%     Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%     12–22, 2014.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 21, 2018

% List of parameters to save with results
parameters_list = {
    'color_map_filename',...
    'normalize_color_map',...
    'normalize_color_map_reference_index',...
    'image_sampling',...
    'n_samples',...
    'n_patch_sizes',...
    'patch_size_max',...
    'n_padding',...
    'padding_ratio_max',...
    'dispersion_px',...
    'n_dispersion_additional',...
    'noise_fractions',...
    'output_directory'...
};

%% Input data and parameters

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180923_TestingChirpImageGeneration/CIE1931ColorMapData.mat';

% Whether or not to normalize spectral sensitivity functions, assuming an
% illuminant which has a uniform spectral power distribution. The
% normalization code is based on 'reflectanceToRadiance()'
normalize_color_map = true;
% Which spectral sensitivity function to use for computing the
% normalization
normalize_color_map_reference_index = 2; % Use '2' for the CIE tristimulus functions (since 2 denotes Y) 

% Image dimensions (height, width)
image_sampling = [32, 32];

% Number of samples for antialiasing during image generation
n_samples = 1;

% Number of patch sizes to test. Patches will be square.
% The smallest patch size will be 2 x 2 pixels.
n_patch_sizes = 5;
% Maximum patch side length. This value will be clipped to the largest
% image dimension before use.
patch_size_max = 50;

% Number of patch padding sizes to test, including no padding.
n_padding = 4;

% Size of the largest padding size as a multiple of the maximum dispersion
padding_ratio_max = 2;

% Dispersion magnitudes in pixels to test. Note that zero dispersion will
% always be tested (and so will be added to the list if it is not specified
% here). Negative dispersion values are not allowed.
dispersion_px = [0.1, 0.3, 1, 2, 3];
% Number of additional dispersion magnitudes to test, provided that the
% largest value in `dispersion_px` is below the suggested maximum
% dispersion value output by 'chirpImage()'. (Otherwise, no additional
% dispersion magnitudes will be tested.) The additional dispersion
% magnitudes are logarithmically-spaced.
n_dispersion_additional = 3;

% Noise fractions: The standard deviation of the noise added to a given
% image value is these fractions of the value
noise_fractions = [0, 0.05, 0.1, 0.25, 0.5];

% Number of patches to show spectral error plots for
n_eval_patches_x = 4;
n_eval_patches_y = 4;

% Spectral evaluation patch size (Must be an odd integer)
eval_patch_size = 5;

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% ## Debugging flags
baek2017Algorithm2Verbose = false;
verbose_progress = true;

%% Validate parameters, and construct intermediate parameters

% Tweaks to the parameters set by 'SetFixedParameters.m'
selectWeightsGridOptions.parallel = false;
trainWeightsOptions.parallel = false;
% Not all patches will have full borders, so do the following for
% simplicity
baek2017Algorithm2Options.l_err_border = [0, 0];
trainWeightsOptions.border = 0;

if add_border
    % Estimating a border area results in images which are usually not
    % registered with the ground truth.
    error('Estimating a border around images prevents evaluation of the mean square error criterion.');
end
baek2017Algorithm2Options.add_border = false;
baek2017Algorithm2Options.l_surface = true;

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

patch_sizes = repmat(...
    round(logspace(log10(2), log10(min([max(image_sampling), patch_size_max])), n_patch_sizes)).', 1, 2 ...
);

if n_padding < 1
    error('At least one padding size must be tested. The first padding size is always zero.');
end
[~, dispersion_max] = chirpImage(...
  image_sampling_3, lambda_range, 1, [], 'params'...
);
paddings = [0, round(logspace(0, log10(padding_ratio_max * dispersion_max), n_padding - 1))];

dispersion_px = sort(dispersion_px);
if dispersion_px(1) < 0
    error('Dispersion magnitudes to test must be non-negative.');
elseif dispersion_px(1) == 0
    dispersion_px_all = dispersion_px;
else
    dispersion_px_all = [0 dispersion_px];
end
if dispersion_px(end) < dispersion_max
    dispersion_px_all = [ dispersion_px_all, logspace(...
        log10(dispersion_px(end)), log10(dispersion_max), n_dispersion_additional...
    )];
end
n_dispersion = length(dispersion_px_all);

n_noise_fractions = length(noise_fractions);

imageFormationOptions.patch_size = [100, 100];
% There is no need to use a patch border, because 'imageFormation()' will
% not be asked to perform image warping.
imageFormationOptions.padding = 0;

enabled_weights = trainWeightsOptions.enabled_weights;
if any(enabled_weights ~= selectWeightsGridOptions.enabled_weights)
    error('Expected `trainWeightsOptions` and `selectWeightsGridOptions` to have the same lists of enabled weights.');
end
n_active_weights = sum(enabled_weights);
n_weights = length(enabled_weights);
to_all_weights = find(enabled_weights);

%% Generate a dataset description for future use

raw_filename_postfix = 'raw';
dataset_params.raw_images_wildcard = fullfile(output_directory, ['*_', raw_filename_postfix, '.mat']);
dataset_params.raw_images_variable = 'I_raw';
rgb_filename_postfix = 'rgb';
dataset_params.rgb_images_wildcard = fullfile(output_directory, ['*_', rgb_filename_postfix, '.mat']);
dataset_params.rgb_images_variable = 'I_rgb';
spectral_filename_postfix = 'latent';
dataset_params.spectral_images_wildcard = fullfile(output_directory, ['*_', spectral_filename_postfix, '.mat']);
dataset_params.spectral_images_variable = 'I_latent';
dataset_params.spectral_reflectances = false;
dataset_params.dispersion_rgb_forward = []; % Presently, datasets cannot have per-image dispersion
dataset_params.dispersion_rgb_reverse = []; % Presently, datasets cannot have per-image dispersion
dataset_params.dispersion_spectral_reverse = []; % Presently, datasets cannot have per-image dispersion
dataset_params.color_map = color_map_filename;
save_data_filename = fullfile(output_directory, 'CorrectChirpImage.mat'); 
dataset_params.wavelengths = save_data_filename;
dataset_params.evaluation = struct(...
    'global_rgb', struct('error_map', true),...
    'custom_rgb', struct,...
    'global_spectral', struct(...
        'error_map', true,...
        'mi_bands', [1, n_bands],...
        'bands_diff', [1, n_bands]...
    ),...
    'custom_spectral', struct...
);

[eval_patches_x, eval_patches_y] = meshgrid(...
    round(linspace(ceil(eval_patch_size / 2), floor(image_sampling(2) - eval_patch_size / 2), n_eval_patches_x)),...
    round(linspace(ceil(eval_patch_size / 2), floor(image_sampling(1) - eval_patch_size / 2), n_eval_patches_y))...
);
eval_patches_x = reshape(eval_patches_x, [], 1);
eval_patches_y = reshape(eval_patches_y, [], 1);
dataset_params.evaluation.global_spectral.radiance = [...
    eval_patches_x, eval_patches_y,...
    repmat(eval_patch_size, n_eval_patches_x * n_eval_patches_y, 2)...
];

% Use two diagonal lines: One shows increasing variation in both spatial
% and spectral roughness. The other transitions from extreme spectral
% roughness and no spatial roughness, to extreme spatial roughness and no
% spectral roughness
dataset_params.evaluation.global_spectral.scanlines = [
    1, 1, image_sampling(2), image_sampling(1);...
    1, image_sampling(1), image_sampling(2), 1
];

warped_filename_postfix = 'warped';
warped_images_variable = 'I_warped';
rgb_warped_filename_postfix = 'rgb_warped';
rgb_warped_images_variable = 'I_full';

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

% Initialize output variables
n_weights_functions = 2; % 'trainWeights()' and 'selectWeightsGrid()'
estimated_images = cell(1, n_weights_functions);
n_auxiliary_images = 4;
auxiliary_images = cell(n_auxiliary_images, n_weights_functions);
weights_images = cell(n_active_weights, n_weights_functions);
err_images = cell(n_active_weights, 1);

selectWeightsGrid_time = zeros(n_padding, n_patch_sizes, n_noise_fractions, n_dispersion);
trainWeightsGrid_time = zeros(n_padding, n_patch_sizes, n_noise_fractions, n_dispersion);
n_patches = zeros(1, n_patch_sizes);

% Variables for plotting
mdc_color = [1, 0, 0];
mse_color = [0, 1, 0];
weights_functions_names = {'Minimum distance criterion', 'Mean square error'};
weights_functions_abbrev = {'MDC', 'MSE'};
weights_functions_filenames = {'selectWeightsGrid_', 'trainWeights_'};
[patch_sizes_grid, paddings_grid] = meshgrid(patch_sizes, paddings);

% Evaluation variables
n_images = n_dispersion * n_noise_fractions;
n_spectral_evaluations = (n_patch_sizes * n_padding) * n_weights_functions + 1; % Add one for the aberrated image
evaluation_plot_colors = jet(n_spectral_evaluations);
evaluation_plot_colors_admm = evaluation_plot_colors(2:end, :);
evaluation_plot_markers_admm = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
evaluation_plot_styles_admm = {'--', ':', '-.'};
e_rgb_tables = cell(n_images, 1);
e_spectral_tables = cell(n_images, 1);

all_name_params_tables = cell(n_spectral_evaluations, 1);
all_name_params_tables{1} = 'Aberrated';
dataset_name = 'chirp';

image_number = 1;
for d = 1:n_dispersion
    dispersion_px_d = dispersion_px_all(d);

    % Synthesize the true image
    [I_spectral_gt, I_warped_gt, dispersionfun] = chirpImage(...
        image_sampling_3, lambda_range, dispersion_px_d, n_samples...
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
    
    name_params_gt = sprintf('dispersion%e_', dispersion_px_d);
    saveImages(...
        output_directory, name_params_gt,...
        I_raw_noNoise_gt, 'raw_noNoise', 'I_raw_noNoise',...
        I_spectral_gt, spectral_filename_postfix, dataset_params.spectral_images_variable,...
        I_rgb_gt, rgb_filename_postfix, dataset_params.rgb_images_variable,...
        I_rgb_warped_gt, rgb_warped_filename_postfix, rgb_warped_images_variable,...
        I_warped_gt, warped_filename_postfix, warped_images_variable...
    );

    % Compare the aberrated image to the original
    name_params_table_gt = sprintf('Dispersion %g', dispersion_px_d);
    e_rgb_table = evaluateAndSaveRGB(...
        I_rgb_warped_gt, I_rgb_gt, dataset_params, name_params_table_gt, all_name_params_tables{1},...
        fullfile(output_directory, [name_params_gt '_aberrated'])...
    );
    dataset_params.evaluation.global_spectral.plot_color = evaluation_plot_colors(1, :);
    dataset_params.evaluation.global_spectral.plot_marker = 'none';
    dataset_params.evaluation.global_spectral.plot_style = '-';
    [e_spectral_table, fg_spectral] = evaluateAndSaveSpectral(...
        I_warped_gt, I_spectral_gt, bands,...
        dataset_params, name_params_table_gt, all_name_params_tables{1},...
        fullfile(output_directory, [name_params_gt '_aberrated'])...
    );
        
    for no = 1:n_noise_fractions
        noise_fraction = noise_fractions(no);
        snr = noiseFractionToSNR(noise_fraction);
        I_raw_gt = addNoise(I_raw_noNoise_gt, snr);
                
        name_params_gt = sprintf(...
            'dispersion%e_noise%e_', dispersion_px_d, noise_fraction...
        );
        name_params_table_gt = sprintf(...
            'Dispersion %g, noise %g', dispersion_px_d, noise_fraction...
        );
        saveImages(...
            output_directory, name_params_gt,...
            I_raw_gt, raw_filename_postfix, dataset_params.raw_images_variable...
        );
                
        has_dispersion = ~isempty(dispersionfun);
        
        color_ind = 1;
        for ps = 1:n_patch_sizes
            patch_size = patch_sizes(ps, :);
            for pp = 1:n_padding
                padding = paddings(pp);
                
                whole_time_start = tic;
                 
                if verbose_progress
                    fprintf(...
                        'Running dispersion %d, noise %d, patch size %d x %d, padding %d.\n',...
                        dispersion_px_d, noise_fraction, patch_size(1),...
                        patch_size(2), padding...
                    );
                    disp('Splitting the input image into patches...');
                end
                i_vector = 1:patch_size(1):image_sampling_3(1);
                n_i_vector = length(i_vector);
                j_vector = 1:patch_size(2):image_sampling_3(2);
                n_j_vector = length(j_vector);
                n_patches_ps = n_i_vector * n_j_vector;
                n_patches(ps) = n_patches_ps;
                patches_J = cell(1, n_j_vector);
                patches_I = cell(1, n_j_vector);
                patches_I_gt = cell(1, n_j_vector);
                patches_I_rgb_gt = cell(1, n_j_vector);
                patches_auxiliary = cell(n_i_vector, n_j_vector, n_auxiliary_images, n_weights_functions);
                patches_weights = cell(1, n_j_vector);
                points_weights = cell(1, n_j_vector);
                patches_err = cell(1, n_j_vector);
                points_err = cell(1, n_j_vector);
                patch_limits = zeros(n_i_vector, n_j_vector, 4);
                patch_trim = zeros(n_i_vector, n_j_vector, 4);
                corners = zeros(n_i_vector, n_j_vector, 2);

                % Divide the input image into patches to be sent to individual parallel
                % workers
                for j = 1:n_j_vector
                    for i = 1:n_i_vector
                        corners(i, j, :) = [i_vector(i), j_vector(j)];
                        [ patch_lim, trim ] = patchBoundaries(...
                            image_sampling, patch_size, padding, corners(i, j, :)...
                        );
                        patch_trim(i, j, :) = reshape(trim, 1, 1, 4);
                        patch_limits(i, j, :) = reshape(patch_lim, 1, 1, 4);

                        if i == 1
                            patches_J{j} = I_raw_gt(:, patch_lim(1, 2):patch_lim(2, 2), :);
                            patches_I_gt{j} = I_spectral_gt(:, patch_lim(1, 2):patch_lim(2, 2), :);
                            patches_I_rgb_gt{j} = I_rgb_gt(:, patch_lim(1, 2):patch_lim(2, 2), :);
                        end
                    end
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('Parallel processing of patches...');
                end

                % Process each patch
                time_trainWeights = 0;
                time_selectWeightsGrid = 0;
                parfor j = 1:n_j_vector
                    patch_limits_j = patch_limits(:, j, :);
                    patch_trim_j = patch_trim(:, j, :);
                    corners_j = corners(:, j, :);
                    patches_J_j = patches_J{j};
                    patches_I_gt_j = patches_I_gt{j};
                    patches_I_rgb_gt_j = patches_I_rgb_gt{j};
                    rows_I_j = size(patches_J_j, 1);
                    cols_I_j = patch_trim_j(1, 1, 4) - patch_trim_j(1, 1, 3) + 1;
                    patches_I_j = zeros(...
                        rows_I_j, cols_I_j, n_bands, n_weights_functions...
                    );
                    patches_auxiliary_j = cell(n_i_vector, 1, n_auxiliary_images, n_weights_functions);
                    patches_weights_j = zeros(...
                        rows_I_j, cols_I_j, n_active_weights, n_weights_functions...
                    );
                    points_weights_j = zeros(...
                        n_i_vector, 1, n_active_weights, n_weights_functions...
                    );
                    patches_err_j = zeros(...
                        rows_I_j, cols_I_j, n_active_weights...
                    );
                    points_err_j = zeros(...
                        n_i_vector, 1, n_active_weights...
                    );
                    for i = 1:n_i_vector
                        patch_lim = reshape(patch_limits_j(i, 1, :), 2, 2);
                        if isempty(bayer_pattern)
                            align_f = [];
                        else
                            align_f = offsetBayerPattern(patch_lim(1, :), bayer_pattern);
                        end
                        image_sampling_f = diff(patch_lim, 1, 1) + 1;
                        trim = reshape(patch_trim_j(i, 1, :), 2, 2);
                        padding_filter = false(image_sampling_f);
                        padding_filter((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2))) = true;
                        patch_trimmed_size = diff(trim, 1, 1) + 1;
                        rows_I_ij = ((i - 1) * patch_size(1) + 1):((i - 1) * patch_size(1) + patch_trimmed_size(1));
                        
                        if has_dispersion
                            dispersion_matrix_patch = dispersionfunToMatrix(...
                                dispersionfun, bands, image_sampling_f, image_sampling_f,...
                                [0, 0, image_sampling_f(2), image_sampling_f(1)], true,...
                                [corners_j(i, 1, 2), corners_j(i, 1, 1)] - 1 ...
                            );
                        else
                            dispersion_matrix_patch = [];
                        end
                        
                        % Analyze the true image patch
                        image_sampling_f_3 = [image_sampling_f, n_bands];
                        patches_I_gt_ij = patches_I_gt_j(patch_lim(1, 1):patch_lim(2, 1), :, :);
                        patches_I_rgb_gt_ij = patches_I_rgb_gt_j(patch_lim(1, 1):patch_lim(2, 1), :, :);
                        for w = 1:n_active_weights
                            aw = to_all_weights(w);
                            if aw == 1 || aw == 2
                                G = spatialGradient(image_sampling_f_3);
                            end
                            if aw == 2
                                G_lambda = spectralGradient(image_sampling_f_3, baek2017Algorithm2Options.full_GLambda);
                                G_lambda_sz1 = size(G_lambda, 1);
                                G_lambda_sz2 = size(G_lambda, 2);
                                % The product `G_lambda * G_xy` must be defined, so `G_lambda` needs to be
                                % replicated to operate on both the x and y-gradients.
                                G_lambda = [
                                    G_lambda, sparse(G_lambda_sz1, G_lambda_sz2);
                                    sparse(G_lambda_sz1, G_lambda_sz2), G_lambda
                                    ];
                                G = G_lambda * G;
                            end
                            if aw == 3
                                G = antiMosaicMatrix(image_sampling_f, align_f);
                                err_vector = G * reshape(patches_I_rgb_gt_ij, [], 1);
                            else
                                err_vector = G * reshape(patches_I_gt_ij, [], 1);
                            end
                            if baek2017Algorithm2Options.norms(aw)
                                penalty = mean(abs(err_vector));
                            else
                                penalty = dot(err_vector, err_vector) / length(err_vector);
                            end
                            penalty = log10(penalty);
                            patches_err_j(rows_I_ij, :, w) = penalty;
                            points_err_j(i, 1, w) = penalty;
                        end
                        
                        patches_J_ij = patches_J_j(patch_lim(1, 1):patch_lim(2, 1), :, :);

                        % Process the patch with selectWeightsGrid()
                        weights_function_index = 1;
                        selectWeightsGrid_time_start = tic;
                        [ weights, ~, patches_I_ij ] = selectWeightsGrid(...
                            patches_J_ij, align_f, dispersion_matrix_patch,...
                            sensor_map_resampled, bands, rho,...
                            baek2017Algorithm2Options, selectWeightsGridOptions,...
                            corners_j(i, 1, :), selectWeightsGridVerbose...
                        );
                        [patches_I_ij, patches_auxiliary_j(i, 1, :, weights_function_index)] = estimateAuxiliaryImages(...
                                patches_I_ij, dispersion_matrix_patch, padding_filter,...
                                patch_trimmed_size,...
                                sensor_map_resampled, bands, int_method,...
                                align_f, n_auxiliary_images...
                        );
                        patches_I_j(rows_I_ij, :, :, weights_function_index) = patches_I_ij;
                        for w = 1:n_active_weights
                            weight = log10(weights(to_all_weights(w)));
                            patches_weights_j(rows_I_ij, :, w, weights_function_index) = weight;
                            points_weights_j(i, 1, w, weights_function_index) = weight;
                        end
                        time_selectWeightsGrid = time_selectWeightsGrid + toc(selectWeightsGrid_time_start);
                        
                        % Process the patch with trainWeights()
                        weights_function_index = 2;
                        trainWeights_time_start = tic;
                        [ weights, ~, patches_I_ij ] = trainWeights(...
                            patches_I_gt_ij, patches_J_ij, align_f,...
                            dispersion_matrix_patch, sensor_map_resampled,...
                            bands, trainWeightsOptions, @baek2017Algorithm2,...
                            {rho, baek2017Algorithm2Options, baek2017Algorithm2Verbose},...
                            corners_j(i, 1, :), trainWeightsVerbose...
                        );
                        [patches_I_ij, patches_auxiliary_j(i, 1, :, weights_function_index)] = estimateAuxiliaryImages(...
                                patches_I_ij, dispersion_matrix_patch, padding_filter,...
                                patch_trimmed_size,...
                                sensor_map_resampled, bands, int_method,...
                                align_f, n_auxiliary_images...
                        );
                        patches_I_j(rows_I_ij, :, :, weights_function_index) = patches_I_ij;
                        for w = 1:n_active_weights
                            weight = log10(weights(to_all_weights(w)));
                            patches_weights_j(rows_I_ij, :, w, weights_function_index) = weight;
                            points_weights_j(i, 1, w, weights_function_index) = weight;
                        end
                        time_trainWeights = time_trainWeights + toc(trainWeights_time_start);
                        
                        if verbose_progress
                            fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i_vector, n_patches_ps);
                        end
                    end
                    patches_I{j} = patches_I_j;
                    patches_auxiliary(:, j, :, :) = patches_auxiliary_j;
                    patches_weights{j} = patches_weights_j;
                    points_weights{j} = points_weights_j;
                    patches_err{j} = patches_err_j;
                    points_err{j} = points_err_j;
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('Recombining results from patches...');
                end

                % Recombine patches
                I_4D = cell2mat(patches_I);
                weights_images_4D = cell2mat(patches_weights);
                err_images_3D = cell2mat(patches_err);
                for f = 1:n_weights_functions
                    estimated_images{f} = I_4D(:, :, :, f);
                    for im = 1:n_auxiliary_images
                        auxiliary_images{im, f} = cell2mat(patches_auxiliary(:, :, im, f));
                    end
                    for w = 1:n_active_weights
                        weights_images{w, f} = weights_images_4D(:, :, w, f);
                    end
                end
                for w = 1:n_active_weights
                    err_images{w} = err_images_3D(:, :, w);
                end
                points_weights_4D = cell2mat(points_weights);
                points_err_3D = cell2mat(points_err);

                whole_time_elapsed = toc(whole_time_start);
                if verbose_progress
                    fprintf('\tDone in %g seconds.\n', whole_time_elapsed);
                    disp('Performing evaluation and file output...')
                end
                selectWeightsGrid_time(pp, ps, no, d) = whole_time_elapsed - time_trainWeights;
                trainWeightsGrid_time(pp, ps, no, d) = whole_time_elapsed - time_selectWeightsGrid;

                % Save the output images
                name_params = sprintf(...
                    'noise%e_dispersion%e_patch%dx%d_pad%d_',...
                    noise_fraction, dispersion_px_d, patch_size(1),...
                    patch_size(2), padding...
                );
                name_params_tables = sprintf(...
                    'patch %d x %d, padding %d',...
                    patch_size(1), patch_size(2), padding...
                );
                name_params_files = fullfile(output_directory, name_params);
            
                for f = 1:n_weights_functions
                    name_params_f = [name_params, weights_functions_filenames{f}];
                    saveImages(...
                        output_directory, name_params_f,...
                        estimated_images{f}, spectral_filename_postfix, dataset_params.spectral_images_variable,...
                        auxiliary_images{1, f}, rgb_filename_postfix, dataset_params.rgb_images_variable,...
                        auxiliary_images{2, f}, rgb_warped_filename_postfix, rgb_warped_images_variable,...
                        auxiliary_images{3, f}, raw_filename_postfix, dataset_params.raw_images_variable,...
                        auxiliary_images{4, f}, warped_filename_postfix, warped_images_variable...
                    );
                    for w = 1:n_active_weights
                        aw = to_all_weights(w);
                        saveImages(...
                            'data', output_directory, name_params_f,...
                            weights_images{w, f}, sprintf('weight%dImage', aw), 'I_weights'...
                        );
                    
                        fg = figure;
                        imagesc(weights_images{w, f});
                        c = colorbar;
                        c.Label.String = sprintf('log_{10}(weight %d)', aw);
                        xlabel('Image x-coordinate')
                        ylabel('Image y-coordinate')
                        title(sprintf('Per-patch %s weight %d', weights_functions_names{f}, aw));
                        savefig(...
                            fg,...
                            [name_params_files sprintf('weight%dImage.fig', aw)], 'compact'...
                            );
                        close(fg);
                    end
                end
                for w = 1:n_active_weights
                    aw = to_all_weights(w);
                    saveImages(...
                        'data', output_directory, name_params_f,...
                        err_images{w}, sprintf('penalty%d', aw), 'I_penalty'...
                    );

                    fg = figure;
                    imagesc(err_images{w});
                    c = colorbar;
                    c.Label.String = sprintf('log_{10}(penalty %d)', aw);
                    xlabel('Image x-coordinate')
                    ylabel('Image y-coordinate')
                    title(sprintf('Per-patch regularization penalty %d', aw));
                    savefig(...
                        fg,...
                        [name_params_files sprintf('penalty%dPatches.fig', aw)], 'compact'...
                        );
                    close(fg);
                end
                
                % Comparative visualization of weights                
                mdc_weights = reshape(points_weights_4D(:, :, :, 1), n_patches_ps, n_active_weights);
                mse_weights = reshape(points_weights_4D(:, :, :, 2), n_patches_ps, n_active_weights);
                
                fg = figure;
                hold on
                if n_active_weights == 1
                    scatter(...
                        mse_weights, mdc_weights, 'filled'...
                        );
                    line_limits = [...
                        min(min(mse_weights, mdc_weights));
                        max(max(mse_weights, mdc_weights))
                        ];
                    line(line_limits, line_limits, 'Color', 'b');
                    legend('Weights', 'y = x');
                    xlabel(sprintf('Weight selected using the %s', weights_functions_abbrev{2}));
                    ylabel(sprintf('Weight selected using the %s', weights_functions_abbrev{1}));
                elseif n_active_weights == 2
                    scatter(...
                        mdc_weights(:, 1), mdc_weights(:, 2), [], mdc_color, 'filled'...
                        );
                    scatter(...
                        mse_weights(:, 1), mse_weights(:, 2), [], mse_color, 'filled'...
                        );
                    legend(weights_functions_abbrev{:});
                    xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
                    ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
                elseif n_active_weights == 3
                    scatter3(...
                        mdc_weights(:, 1), mdc_weights(:, 2), mdc_weights(:, 3),...
                        [], mdc_color, 'filled'...
                        );
                    scatter3(...
                        mse_weights(:, 1), mse_weights(:, 2), mse_weights(:, 3),...
                        [], mse_color, 'filled'...
                        );
                    legend(weights_functions_abbrev{:});
                    xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
                    ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
                    zlabel(sprintf('log_{10}(weight %d)', to_all_weights(3)))
                else
                    error('Unexpected number of active weights.');
                end
                title(sprintf('Agreement between %s and %s weights', weights_functions_abbrev{:}));
                axis equal
                hold off
                
                savefig(...
                    fg,...
                    [name_params_files 'weightsCorrelation.fig'], 'compact'...
                    );
                close(fg);
                
                for w = 1:n_active_weights
                    aw = to_all_weights(w);
                    penalties_w = reshape(points_err_3D(:, :, aw), n_patches_ps, 1);
                    
                    fg = figure;
                    hold on
                    scatter(penalties_w, mdc_weights(:, w), [], mdc_color, 'filled');
                    scatter(penalties_w, mse_weights(:, w), [], mse_color, 'filled');
                    legend(weights_functions_abbrev{:});
                    xlabel(sprintf('log_{10}(Penalty %d)', aw))
                    ylabel(sprintf('log_{10}(Weight %d)', aw))
                    title(sprintf('%s and %s weights as a function of the corresponding penalty term', weights_functions_abbrev{:}));
                    hold off
                    
                    savefig(...
                        fg,...
                        [name_params_files sprintf('weight%dVsPenalty.fig', aw)], 'compact'...
                        );
                    close(fg);
                end
               
                for f = 1:n_weights_functions
                    % Spectral evaluation
                    name_params_f = [name_params, weights_functions_filenames{f}'];
                    name_params_tables_f = [name_params_tables, weights_functions_abbrev{f}];
                    dataset_params.evaluation.global_spectral.plot_color =...
                        evaluation_plot_colors_admm(color_ind, :);
                    dataset_params.evaluation.global_spectral.plot_marker =...
                        evaluation_plot_markers_admm(...
                            mod(color_ind - 1, length(evaluation_plot_markers_admm)) + 1 ...
                        );
                    dataset_params.evaluation.global_spectral.plot_style =...
                        evaluation_plot_styles_admm(...
                            mod(color_ind - 1, length(evaluation_plot_styles_admm)) + 1 ...
                        );
                    color_ind = color_ind + 1;
                    all_name_params_tables{color_ind} = name_params_tables_f;
                    [e_spectral_table_current, fg_spectral] = evaluateAndSaveSpectral(...
                        estimated_images{f}, I_spectral_gt, bands, dataset_params,...
                        name_params_table_gt, name_params_tables_f,...
                        fullfile(output_directory, name_params_f(1:(end-1))),...
                        fg_spectral...
                        );
                    e_spectral_table = union(e_spectral_table_current, e_spectral_table);
                    
                    % RGB evaluation
                    e_rgb_table_current = evaluateAndSaveRGB(...
                        auxiliary_images{1, f}, I_rgb_gt, dataset_params,...
                        name_params_table_gt, name_params_tables_f,...
                        fullfile(output_directory, name_params_f(1:(end-1)))...
                        );
                    e_rgb_table = union(e_rgb_table_current, e_rgb_table);
                end
            end
        end
        
        % Write evaluations to a file
        name_params_figures_gt = fullfile(output_directory, name_params_gt);
        writetable(...
            e_rgb_table, [name_params_figures_gt, '_evaluateRGB.csv']...
        );
        e_rgb_tables{image_number} = e_rgb_table;
        writetable(...
            e_spectral_table, [name_params_figures_gt, '_evaluateSpectral.csv']...
        );
        % Also save completed figures
        evaluateAndSaveSpectral(...
            output_directory, dataset_params, name_params_table_gt, all_alg_names, fg_spectral...
        );
        e_spectral_tables{image_number} = e_spectral_table;
    
        % Show and save time results
        
        fg = figure;
        hold on
        surf(...
            patch_sizes_grid, paddings_grid, selectWeightsGrid_time(:, :, no, d),...
            'FaceAlpha', 0.5, 'FaceColor', mdc_color...
        );
        surf(...
            patch_sizes_grid, paddings_grid, trainWeights_time(:, :, no, d),...
            'FaceAlpha', 0.5, 'FaceColor', mse_color...
        );
        colorbar
        xlabel('Patch size');
        ylabel('Padding size');
        zlabel('Time to process image [s]');
        legend(weights_functions_abbrev{:});
        title(sprintf(...
            'Time for a %d x %d x %d image under dispersion %f and noise %f',...
            image_sampling_3(1), image_sampling_3(2), image_sampling_3(3),...
            dispersion_px_d, noise_fraction...
        ));
        hold off
        savefig(...
            fg,...
            [name_params_figures_gt 'time.fig'], 'compact'...
        );
        close(fg); 
        
        fg = figure;
        hold on
        surf(...
            patch_sizes_grid, paddings_grid,...
            selectWeightsGrid_time(:, :, no, d) ./ ((patch_sizes_grid + paddings_grid) .^ 2),...
            'FaceAlpha', 0.5, 'FaceColor', mdc_color...
        );
        surf(...
            patch_sizes_grid, paddings_grid,...
            trainWeights_time(:, :, no, d) ./ ((patch_sizes_grid + paddings_grid) .^ 2),...
            'FaceAlpha', 0.5, 'FaceColor', mse_color...
        );
        colorbar
        xlabel('Patch size');
        ylabel('Padding size');
        zlabel('Time per pixel in a patch [s]');
        legend(weights_functions_abbrev{:});
        title(sprintf(...
            'Time for a %d x %d x %d image under dispersion %f and noise %f, per patch pixel',...
            image_sampling_3(1), image_sampling_3(2), image_sampling_3(3),...
            dispersion_px_d, noise_fraction...
        ));
        hold off
        savefig(...
            fg,...
            [name_params_figures_gt 'timePerPixel.fig'], 'compact'...
        );
        close(fg);
        
        fg = figure;
        hold on
        n_patches_rep = repmat(n_patches, n_padding, 1);
        surf(...
            patch_sizes_grid, paddings_grid,...
            selectWeightsGrid_time(:, :, no, d) ./ n_patches_rep,...
            'FaceAlpha', 0.5, 'FaceColor', mdc_color...
        );
        surf(...
            patch_sizes_grid, paddings_grid,...
            trainWeights_time(:, :, no, d) ./ n_patches_rep,...
            'FaceAlpha', 0.5, 'FaceColor', mse_color...
        );
        colorbar
        xlabel('Patch size');
        ylabel('Padding size');
        zlabel('Time per patch [s]');
        legend(weights_functions_abbrev{:});
        title(sprintf(...
            'Time for a %d x %d x %d image under dispersion %f and noise %f, per patch',...
            image_sampling_3(1), image_sampling_3(2), image_sampling_3(3),...
            dispersion_px_d, noise_fraction...
        ));
        hold off
        savefig(...
            fg,...
            [name_params_figures_gt 'timePerPatch.fig'], 'compact'...
        );
        close(fg);
        
        image_number = image_number + 1;
        
        if verbose_progress
            disp('Finished evaluation and file output.')
        end
    end
end

%% Save evaluation results for all images

e_rgb_summary_table = mergeRGBTables(e_rgb_tables);
writetable(...
    e_rgb_summary_table,...
    fullfile(output_directory, [dataset_name, '_evaluateRGB.csv'])...
);
e_spectral_summary_table = mergeSpectralTables(e_spectral_tables);
writetable(...
    e_spectral_summary_table,...
    fullfile(output_directory, [dataset_name, '_evaluateSpectral.csv'])...
);

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'dataset_params',...
        'bands_color',...
        'bands',...
        'sensor_map_resampled',...
        'selectWeightsGrid_time',...
        'trainWeights_time',...
        'n_patches'...
    } ];
save(save_data_filename, save_variables_list{:});