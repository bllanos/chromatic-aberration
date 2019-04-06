%% Evaluate hyperspectral ADMM-based reconstruction of a chirp image
%
% This script investigates the following questions:
% - How does image reconstruction accuracy change with the spatial and
%   spectral frequencies of the image?
% - How does image reconstruction accuracy change with the signal-to-noise
%   ratio?
% - How do the spatial and spectral frequencies of the image affect the
%   choice of regularization weights?
% - How do the regularization weights chosen by the minimum distance
%   criterion differ, and by image similarity, differ, and how does this
%   difference affect image reconstruction time and accuracy?
% - How is image reconstruction accuracy affected by the patch size?
% - How much overlap should patches have?
% - How does the amount of dispersion affect image reconstruction accuracy?
% - How long does image reconstruction take with different patch and
%   padding sizes?
%
% Presently, this script evaluates three image estimation methods:
% - Regularization weight selection, per image patch, using
%   the minimum distance criterion, followed by image estimation by an
%   ADMM-family algorithm, using the selected per-patch regularization
%   weights.
% - Regularization weight selection, per image patch, using
%   image similarity, specifically mean squared error with respect to the
%   true image, followed by image estimation by an ADMM-family algorithm,
%   using the selected per-patch regularization weights.
% - Same as above, but selecting regularization weights using the mean
%   squared error with respect to a demosaicing result.
%
% The same ADMM-family algorithm is used for both methods:
% 'baek2017Algorithm2LowMemory()' (via 'solvePatchesSpectral()'),
% configured according to the parameters in 'SetFixedParameters.m'.
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
%   corresponding to the second dimension of 'sensor_map'.
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
%   output if the spectral image is a greyscale or 3-channel image.
% - '*_rgb.tif' and '*_rgb.mat': A colour image (stored in the variable
%   'I_rgb') created by converting the spectral image to the RGB colour
%   space of the camera. For estimated images, only the '.tif' file will be
%   saved.
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
% Of the above types of images, the following *estimated* images will only
% be saved if `save_all_images` set in 'SetFixedParameters.m' is `true`:
% - '*_warped.tif' and '*_warped.mat'
% - '*_rgb_warped.tif' and '*_rgb_warped.mat'
% - '*_raw.tif' and '*_raw.mat'
%
% ### Data file output
%
% #### Parameters, and execution time results
%
% In the following, 'f' indexes regularization weight selection methods:
% 1 - Minimum distance criterion
% 2 - Error with respect to the true spectral image
% 3 - Error with respect to a demosaicking result
%
% The output data file is a '.mat' file containing the following variables:
%
% - 'bands': A vector containing the wavelengths at which the estimated
%   images are sampled.
% - 'bands_spectral': A vector containing the wavelengths at which the true
%   images are sampled.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference.
% - 'dataset_params': A structure of the form of the `dataset_params`
%   output argument of 'describeDataset()'. 'dataset_params' allows the
%   chirp images created by this script to be used as a dataset for
%   evaluating other image estimation algorithms. Presently,
%   'dataset_params' is missing dispersion information, as dispersion is
%   not the same between images, and so cannot be described with a single
%   dispersion model.
% - 'time': Execution timing information, stored as a 5D array.
%   `time(pp, ps, no, d, f)` is the time taken (in seconds) with the pp-th
%   patch padding size, the ps-th patch size, the no-th noise fraction in
%   the input image, and the d-th level of dispersion in the input image.
%   Time values include the time taken to generate output images beyond the
%   estimated spectral image, including visualizations of selected
%   regularization weights.
% - 'n_patches': A vector containing the number of patches resulting from
%   each patch size used in the experiment.
% - 'num_workers': The number of workers in the parallel pool.
%   'num_workers' is one if there is no parallel pool, or if the parallel
%   pool no longer exists, for some reason, when this variable is set just
%   after the parallel for-loop.
% - 'mrae': Mean relative absolute error over all spectral bands for the
%   images output the image estimation algorithms, stored as a 5D array.
%   `mrae(pp, ps, no, d, f)` is the error for the pp-th patch padding size,
%   the ps-th patch size, the no-th noise fraction in the input image, and
%   the d-th level of dispersion in the input image.
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
% - This script uses `solvePatchesSpectralOptions.reg_options.enabled`
%   defined in 'SetFixedParameters.m' to determine which regularization
%   weights to set.
% - If `solvePatchesSpectralOptions.sampling_options.show_steps` is
%   `true`, output images will be saved only for the highest spectral
%   resolution, even though the images returned by
%   'solvePatchesSpectral()' will contain multiple spectral resolutions.
%   However, there will be additional spectral evaluation figures and CSV
%   files comparing the results between spectral resolutions (for the same
%   image estimation algorithm and image synthesis condition). Refer to the
%   documentation of 'solvePatchesSpectral.m' for more information.
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
% perceptible, whereas dispersion of 0.3 pixels and above is noticeable:
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
    'max_spatial_freq',...
    'power_threshold',...
    'n_samples',...
    'chirp_image_type',...
    'n_patch_sizes',...
    'patch_size_min',...
    'patch_size_max',...
    'n_padding',...
    'padding_min',...
    'padding_ratio_max',...
    'dispersion_px',...
    'n_dispersion_additional',...
    'noise_fractions',...
    'output_directory'...
};

%% Input data and parameters

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181127_TestingChoiEtAl2017/NikonD5100ColorMapData.mat';

% Whether or not to normalize spectral sensitivity functions, assuming an
% illuminant which has a uniform spectral power distribution. The
% normalization code is based on 'reflectanceToRadiance()'
normalize_color_map = true;
% Which spectral sensitivity function to use for computing the
% normalization
normalize_color_map_reference_index = 2; % Use '2' for the CIE tristimulus functions (since 2 denotes Y) 

% Image dimensions (height, width)
image_sampling = [128, 128];

% Maximum spatial frequency at the right border of the image
max_spatial_freq = 0.5; % Cycles per pixel displacement

% Power threshold input to 'bandlimit()' to determine the maximum spectral
% frequency at the bottom border of the image
power_threshold = 0.95;

% Number of samples for antialiasing during image generation
n_samples = 100;

% Type of chirp image to generate
chirp_image_type = 'phase';

% Number of patch sizes to test. Patches will be square.
% Padding sizes will be rounded to the nearest even integers.
n_patch_sizes = 1;

% Minimum patch side length
patch_size_min = 24;

% Maximum patch side length. This value will be clipped to the largest
% image dimension before use.
patch_size_max = 24;

% Number of patch padding sizes to test
% Padding sizes will be rounded to the nearest even integers.
n_padding = 1;

% Minimum amount of padding
padding_min = 16;

% Size of the largest padding size as a multiple of the maximum dispersion
padding_ratio_max = 1;

% Dispersion magnitudes in pixels to test. Note that zero dispersion will
% always be tested (and so will be added to the list if it is not specified
% here). Negative dispersion values are not allowed.
dispersion_px = []; %[0.1, 0.3, 1, 2, 3];
% Number of additional dispersion magnitudes to test, provided that the
% largest value in `dispersion_px` is below the suggested maximum
% dispersion value output by 'chirpImage()'. (Otherwise, no additional
% dispersion magnitudes will be tested.) The additional dispersion
% magnitudes are logarithmically-spaced.
n_dispersion_additional = 0;

% Noise fractions: The standard deviation of the noise added to a given
% image value is these fractions of the value
noise_fractions = 0; %[0, 0.05, 0.1, 0.25, 0.5];

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
verbose_progress = true;

%% Load calibration data

model_variables_required = { 'sensor_map', 'channel_mode', 'bands' };
load(color_map_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end
if channel_mode
    error('The input space of the colour conversion data must be a spectral space, not a space of colour channels.')
end
bands_color = bands;
n_channels = size(sensor_map, 1);

% Normalize spectral sensitivities to avoid out-of-range colour values.
if normalize_color_map
    ybar = sensor_map(normalize_color_map_reference_index, :);
    weights = integrationWeights(bands_color, findSamplingOptions.int_method);
    N = dot(ybar, weights);
    sensor_map = sensor_map ./ N;
end

%% Determine spectral sampling for image estimation and evaluation

max_freq = [max_spatial_freq, 0];
% Find bandlimit of spectral sensitivities
max_freq(2) = bandlimit(sensor_map, power_threshold); % Units of cycles per index
max_freq(2) = max_freq(2) / mean(diff(bands_color)); % Units of cycles per unit change in wavelength

sampling_options_gt = solvePatchesSpectralOptions.sampling_options;
sampling_options_gt.power_threshold = 1;
sampling_options_gt.n_bands = 0;
sampling_options_gt.interpolant = @triangle;
[...
    ~, ~, bands_spectral...
] = findSampling(...
  sensor_map, bands_color, {}, sampling_options_gt, findSamplingVerbose...
);

[...
    ~, spectral_weights, bands...
] = findSampling(...
  sensor_map, bands_color, bands_spectral, solvePatchesSpectralOptions.sampling_options, findSamplingVerbose...
);
n_bands = length(bands);
n_bands_spectral = length(bands_spectral);
image_sampling_3 = [image_sampling, n_bands];

%% Validate parameters, and construct intermediate parameters
if use_fixed_weights
    warning('`use_fixed_weights` is set, so weights will not be automatically selected.');
end

patch_sizes = repmat(...
    round(logspace(log10(patch_size_min), log10(min([max(image_sampling), patch_size_max])), n_patch_sizes)).', 1, 2 ...
);
patch_sizes(mod(patch_sizes, 2) ~= 0) = patch_sizes(mod(patch_sizes, 2) ~= 0) + 1;

dispersion_max = chirpImage(...
  image_sampling, max_freq, bands_spectral, 'params'...
);
paddings = round(logspace(log10(padding_min), log10(padding_ratio_max * dispersion_max), n_padding));
paddings(mod(paddings, 2) ~= 0) = paddings(mod(paddings, 2) ~= 0) + 1;

if isempty(dispersion_px)
    dispersion_px = 0;
else
    dispersion_px = sort(dispersion_px);
end
if dispersion_px(1) < 0
    error('Dispersion magnitudes to test must be non-negative.');
end
if dispersion_px(1) == 0
    dispersion_px_all = dispersion_px;
else
    dispersion_px_all = [0 dispersion_px];
end
if dispersion_px_all(end) < dispersion_max
    dispersion_additional = logspace(...
        log10(dispersion_px(end)), log10(dispersion_max), n_dispersion_additional + 1 ...
    );
    dispersion_px_all = [ dispersion_px_all, dispersion_additional(2:end)];
end
n_dispersion = length(dispersion_px_all);

n_noise_fractions = length(noise_fractions);

enabled_weights = solvePatchesSpectralOptions.reg_options.enabled;
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
    'global_rgb', struct('error_map', false),...
    'custom_rgb', struct,...
    'global_spectral', struct(...
        'metric', 'gof',...
        'error_map', true,...
        'mi_bands', [1, n_bands_spectral],...
        'bands_diff', [1, n_bands_spectral]...
    ),...
    'custom_spectral', struct...
);

if n_eval_patches_x > 0 && n_eval_patches_y > 0
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
end

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

%% Collect results

% Initialize output variables
n_criteria = length(criteria);
estimated_images = cell(1, n_criteria);
if save_all_images
    n_auxiliary_images = 4;
else
    n_auxiliary_images = 1;
end
auxiliary_images = cell(n_auxiliary_images, n_criteria);
weights_images = cell(n_criteria, 1);    
err_images = cell(n_active_weights, 1);

mrae = zeros(n_padding, n_patch_sizes, n_noise_fractions, n_dispersion, n_criteria);
time = zeros(n_padding, n_patch_sizes, n_noise_fractions, n_dispersion, n_criteria);
n_patches = zeros(1, n_patch_sizes);

% Variables for plotting
[patch_sizes_grid, paddings_grid] = meshgrid(patch_sizes(:, 1), paddings);
error_name = 'Mean relative absolute error';
error_filename = 'mrae';

% Evaluation variables
n_images = n_dispersion * n_noise_fractions;
n_spectral_evaluations = (n_patch_sizes * n_padding) * n_criteria + 1; % Add one for the aberrated image
evaluation_plot_colors = jet(n_spectral_evaluations);
evaluation_plot_colors_admm = evaluation_plot_colors(2:end, :);
evaluation_plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
evaluation_plot_styles = {'--', ':', '-.'};
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
        image_sampling, max_freq, bands_spectral, dispersion_px_d, n_samples, chirp_image_type...
    );
    I_rgb_gt = imageFormation(...
        I_spectral_gt, bands_spectral, sensor_map, bands_color,...
        imageFormationSamplingOptions, imageFormationPatchOptions...
    );
    I_rgb_warped_gt = imageFormation(...
        I_warped_gt, bands_spectral, sensor_map, bands_color,...
        imageFormationSamplingOptions, imageFormationPatchOptions...
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
        
    for no = 1:n_noise_fractions
        % Compare the aberrated image to the original
        name_params_table_gt = sprintf('Dispersion %g', dispersion_px_d);
        e_rgb_table = evaluateAndSaveRGB(...
            I_rgb_warped_gt, I_rgb_gt, dataset_params, name_params_table_gt, all_name_params_tables{1},...
            fullfile(output_directory, [name_params_gt 'aberrated'])...
        );
        dataset_params.evaluation.global_spectral.plot_color = evaluation_plot_colors(1, :);
        dataset_params.evaluation.global_spectral.plot_marker = 'none';
        dataset_params.evaluation.global_spectral.plot_style = '-';
        [e_spectral_table, fg_spectral] = evaluateAndSaveSpectral(...
            I_warped_gt, I_spectral_gt, bands_spectral, eye(n_bands_spectral),...
            dataset_params, name_params_table_gt, all_name_params_tables{1},...
            fullfile(output_directory, [name_params_gt 'aberrated'])...
        );

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
            solvePatchesSpectralOptions.patch_options.patch_size = patch_size;
            for pp = 1:n_padding
                padding = paddings(pp);
                solvePatchesSpectralOptions.patch_options.padding = padding;
                                 
                if verbose_progress
                    fprintf(...
                        'Running dispersion %d, noise %d, patch size %d x %d, padding %d.\n',...
                        dispersion_px_d, noise_fraction, patch_size(1),...
                        patch_size(2), padding...
                    );
                end
                
                % Minimum distance criterion
                criterion_index = mdc_index;
                time_start = tic;
                solvePatchesSpectralOptions.reg_options.demosaic = false;
                [...
                    ~, estimated_images{criterion_index},...
                    auxiliary_images{1, criterion_index},...
                    weights_images{criterion_index},...
                    auxiliary_images{2:end, criterion_index}...
                ] = solvePatchesSpectral(...
                  [], I_raw_gt, bayer_pattern, dispersionfun,...
                  sensor_map, bands_color,...
                  solvePatchesSpectralOptions.sampling_options,...
                  solvePatchesSpectralOptions.admm_options,...
                  solvePatchesSpectralOptions.reg_options,...
                  solvePatchesSpectralOptions.patch_options,...
                  solvePatchesSpectralVerbose...
                );
                time(pp, ps, no, d, criterion_index) = toc(time_start);
            
                % Image similarity criterion
                criterion_index = mse_index;
                time_start = tic;
                I_in.I = I_spectral_gt;
                I_in.spectral_bands = bands_spectral;
                [...
                    bands_all, estimated_images{criterion_index},...
                    auxiliary_images{1, criterion_index},...
                    weights_images{criterion_index},...
                    auxiliary_images{2:end, criterion_index}...
                ] = solvePatchesSpectral(...
                  I_in, I_raw_gt, bayer_pattern, dispersionfun,...
                  sensor_map, bands_color,...
                  solvePatchesSpectralOptions.sampling_options,...
                  solvePatchesSpectralOptions.admm_options,...
                  solvePatchesSpectralOptions.reg_options,...
                  solvePatchesSpectralOptions.patch_options,...
                  solvePatchesSpectralVerbose...
                );
                time(pp, ps, no, d, criterion_index) = toc(time_start);
                
                % Demosaicing image similarity criterion
                criterion_index = dm_index;
                time_start = tic;
                solvePatchesSpectralOptions.reg_options.demosaic = true;
                [...
                    ~, estimated_images{criterion_index},...
                    auxiliary_images{1, criterion_index},...
                    weights_images{criterion_index},...
                    auxiliary_images{2:end, criterion_index}...
                ] = solvePatchesSpectral(...
                  [], I_raw_gt, bayer_pattern, dispersionfun,...
                  sensor_map, bands_color,...
                  solvePatchesSpectralOptions.sampling_options,...
                  solvePatchesSpectralOptions.admm_options,...
                  solvePatchesSpectralOptions.reg_options,...
                  solvePatchesSpectralOptions.patch_options,...
                  solvePatchesSpectralVerbose...
                );
                time(pp, ps, no, d, criterion_index) = toc(time_start);
                
                name_params = sprintf(...
                    'noise%e_dispersion%e_patch%dx%d_pad%d_',...
                    noise_fraction, dispersion_px_d, patch_size(1),...
                    patch_size(2), padding...
                );
                name_params_tables = sprintf(...
                    'patch %d x %d, padding %d',...
                    patch_size(1), patch_size(2), padding...
                );
                
                % Spectral evaluation of intermediate images
                if solvePatchesSpectralOptions.sampling_options.show_steps
                    n_steps = length(bands_all);
                    step_plot_colors = jet(n_steps);
                    step_name_params_tables = cell(n_steps, 1);
                    for f = 1:n_criteria
                        spectral_inc = 0;
                        fg_spectral_step = struct;
                        for t = 1:n_steps
                            name_params_f = [name_params, sprintf('step%d_', t), criteria_filenames{f}];
                            name_params_tables_f = sprintf('%s, step %d, %s', name_params_tables, t, criteria_abbrev{f});
                            dataset_params.evaluation.global_spectral.plot_color =...
                                step_plot_colors(t, :);
                            dataset_params.evaluation.global_spectral.plot_marker =...
                                evaluation_plot_markers{...
                                mod(t - 1, length(evaluation_plot_markers)) + 1 ...
                                };
                            dataset_params.evaluation.global_spectral.plot_style =...
                                evaluation_plot_styles{...
                                mod(t - 1, length(evaluation_plot_styles)) + 1 ...
                                };
                            step_name_params_tables{t} = name_params_tables_f;
                            
                            n_bands_t = length(bands_all{t});
                            spectral_weights_step = resamplingWeights(...
                                bands_spectral, bands_all{t},...
                                solvePatchesSpectralOptions.sampling_options.interpolant,...
                                solvePatchesSpectralOptions.sampling_options.bands_padding...
                                );
                            [e_spectral_table_step_current, fg_spectral_step] = evaluateAndSaveSpectral(...
                                estimated_images{f}(:, :, (spectral_inc + 1):(spectral_inc + n_bands_t)),...
                                I_spectral_gt, bands_spectral, spectral_weights_step,...
                                dataset_params, name_params_table_gt, name_params_tables_f,...
                                fullfile(output_directory, name_params_f(1:(end-1))),...
                                fg_spectral_step...
                                );
                            if t == 1
                                e_spectral_table_step = e_spectral_table_step_current;
                            else
                                e_spectral_table_step = union(e_spectral_table_step_current, e_spectral_table_step);
                            end
                            spectral_inc = spectral_inc + n_bands_t;
                        end
                        writetable(...
                            e_spectral_table_step,...
                            fullfile(output_directory, [name_params, criteria_filenames{f}, 'multiStep_evaluateSpectral.csv'])...
                        );
                        evaluateAndSaveSpectral(...
                            output_directory, dataset_params, [name_params_table_gt, '_', criteria_filenames{f}, 'multiStep'], step_name_params_tables, fg_spectral_step...
                        );
                    
                        % Retain only the highest spectral resolution data
                        % for further study
                        estimated_images{f} = estimated_images{f}(:, :, (end - n_bands + 1):end);
                        auxiliary_images{1, f} = auxiliary_images{1, f}(:, :, (end - n_channels + 1):end);
                        if save_all_images
                            auxiliary_images{2, f} = auxiliary_images{2, f}(:, :, (end - n_channels + 1):end);
                            auxiliary_images{3, f} = auxiliary_images{3, f}(:, :, (end - size(I_raw_gt, 3) + 1):end);
                            auxiliary_images{4, f} = auxiliary_images{4, f}(:, :, (end - n_bands + 1):end);
                        end
                        weights_images{f} = weights_images{f}(:, :, (end - n_active_weights + 1):end);
                    end
                end
            
                for f = 1:n_criteria
                    weights_images{f} = log10(weights_images{f});
                end
                
                if verbose_progress
                    disp('[Image characterization] Splitting the input image into patches...');
                end
                i_vector = 1:patch_size(1):image_sampling_3(1);
                n_i_vector = length(i_vector);
                j_vector = 1:patch_size(2):image_sampling_3(2);
                n_j_vector = length(j_vector);
                n_patches_ps = n_i_vector * n_j_vector;
                n_patches(ps) = n_patches_ps;
                patches_I_gt = cell(1, n_j_vector);
                patches_I_rgb_gt = cell(1, n_j_vector);
                if ~use_fixed_weights
                    points_weights = zeros(n_i_vector, n_j_vector, n_active_weights, n_criteria);
                end
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
                        
                        if ~use_fixed_weights
                            for f = 1:n_criteria
                                points_weights(i, j, :, f) = weights_images{f}(...
                                    corners(i, j, 1), corners(i, j, 2), :...
                                );
                            end
                        end

                        if i == 1
                            patches_I_gt{j} = I_spectral_gt(:, patch_lim(1, 2):patch_lim(2, 2), :);
                            patches_I_rgb_gt{j} = I_rgb_gt(:, patch_lim(1, 2):patch_lim(2, 2), :);
                        end
                    end
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('[Image characterization] Parallel processing of patches...');
                end

                % Process each patch
                parfor j = 1:n_j_vector
                    patch_limits_j = patch_limits(:, j, :);
                    patch_trim_j = patch_trim(:, j, :);
                    corners_j = corners(:, j, :);
                    patches_I_gt_j = patches_I_gt{j};
                    patches_I_rgb_gt_j = patches_I_rgb_gt{j};
                    rows_I_j = size(patches_I_gt_j, 1);
                    cols_I_j = patch_trim_j(1, 1, 4) - patch_trim_j(1, 1, 3) + 1;
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
                        
                        % Analyze the true image patch
                        image_sampling_f_3 = [image_sampling_f, n_bands_spectral];
                        patches_I_gt_ij = patches_I_gt_j(patch_lim(1, 1):patch_lim(2, 1), :, :);
                        patches_I_rgb_gt_ij = patches_I_rgb_gt_j(patch_lim(1, 1):patch_lim(2, 1), :, :);
                        for w = 1:n_active_weights
                            aw = to_all_weights(w);
                            if aw == 1 || aw == 2
                                G = spatialGradient(image_sampling_f_3);
                            end
                            if aw == 2
                                G_lambda = spectralGradient(image_sampling_f_3, solvePatchesSpectralOptions.admm_options.full_GLambda);
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
                                G = spatialLaplacian(image_sampling_f_3);
                            end
                            err_vector = G * reshape(patches_I_gt_ij, [], 1);
                            if solvePatchesSpectralOptions.admm_options.norms(aw)
                                penalty = mean(abs(err_vector));
                            else
                                penalty = dot(err_vector, err_vector) / length(err_vector);
                            end
                            penalty = log10(penalty);
                            patches_err_j(rows_I_ij, :, w) = penalty;
                            points_err_j(i, 1, w) = penalty;
                        end
                        
                        if verbose_progress
                            fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i_vector, n_patches_ps);
                        end
                    end
                    patches_err{j} = patches_err_j;
                    points_err{j} = points_err_j;
                end
                
                pool = gcp('nocreate');
                if isempty(pool)
                    num_workers = 1;
                else
                    num_workers = pool.NumWorkers;
                end

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('[Image characterization] Recombining results from patches...');
                end

                % Recombine patches
                err_images_3D = cell2mat(patches_err);
                for w = 1:n_active_weights
                    err_images{w} = err_images_3D(:, :, w);
                end
                points_err_3D = cell2mat(points_err);

                if verbose_progress
                    fprintf('\tDone.\n');
                    disp('Performing evaluation and file output...')
                end

                % Save the output images
                name_params_files = fullfile(output_directory, name_params);
            
                for f = 1:n_criteria
                    name_params_f = [name_params, criteria_filenames{f}];
                    saveImages(...
                        output_directory, name_params_f,...
                        estimated_images{f}, spectral_filename_postfix, dataset_params.spectral_images_variable...
                    );
                    saveImages(...
                        'image', output_directory, name_params_f,...
                        auxiliary_images{1, f}, rgb_filename_postfix, dataset_params.rgb_images_variable...
                    );
                    if save_all_images
                        saveImages(...
                            output_directory, name_params_f,...
                            auxiliary_images{2, f}, rgb_warped_filename_postfix, rgb_warped_images_variable,...
                            auxiliary_images{3, f}, raw_filename_postfix, dataset_params.raw_images_variable,...
                            auxiliary_images{4, f}, warped_filename_postfix, warped_images_variable...
                        );
                    end

                    if ~use_fixed_weights
                        for w = 1:n_active_weights
                            aw = to_all_weights(w);
                            saveImages(...
                                'data', output_directory, name_params_f,...
                                weights_images{f}(:, :, w), sprintf('weight%dImage', aw), 'I_weights'...
                            );

                            fg = figure;
                            imagesc(weights_images{f}(:, :, w));
                            c = colorbar;
                            c.Label.String = sprintf('log_{10}(weight %d)', aw);
                            xlabel('Image x-coordinate')
                            ylabel('Image y-coordinate')
                            title(sprintf('Per-patch %s weight %d', criteria_names{f}, aw));
                            savefig(...
                                fg,...
                                [name_params_files  criteria_filenames{f} sprintf('weight%dImage.fig', aw)], 'compact'...
                                );
                            close(fg);
                        end
                    end
                end
                for w = 1:n_active_weights
                    aw = to_all_weights(w);
                    saveImages(...
                        'data', output_directory, name_params,...
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
                if ~use_fixed_weights
                    points_weights = reshape(points_weights, n_patches_ps, n_active_weights, n_criteria);

                    fg = figure;
                    hold on
                    if n_active_weights == 1
                        scatter(...
                            points_weights(:, :, mse_index), points_weights(:, :, mdc_index), 'filled'...
                            );
                        scatter(...
                            points_weights(:, :, mse_index), points_weights(:, :, dm_index), 'filled'...
                            );
                        line_limits = [...
                            min(points_weights(:));
                            max(points_weights(:))
                            ];
                        line(line_limits, line_limits, 'Color', 'k');
                        legend(...
                            sprintf('%s vs. %s', criteria_abbrev{mdc_index}, criteria_abbrev{mse_index}),...
                            sprintf('%s vs. %s', criteria_abbrev{dm_index}, criteria_abbrev{mse_index}),...
                            'y = x'...
                        );
                        xlabel(sprintf('Weight selected using the %s criterion', criteria_abbrev{mse_index}));
                        ylabel('Weight selected using the other criterion');
                    elseif n_active_weights == 2
                        for f = 1:n_criteria
                            scatter(...
                                points_weights(:, 1, f), points_weights(:, 2, f), [], criteria_colors(f, :), 'filled'...
                            );
                        end
                        legend(criteria_abbrev{:});
                        xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
                        ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
                    elseif n_active_weights == 3
                        for f = 1:n_criteria
                            scatter3(...
                                points_weights(:, 1, f), points_weights(:, 2, f), points_weights(:, 3, f),...
                                [], criteria_colors(f, :), 'filled'...
                            );
                        end
                        legend(criteria_abbrev{:});
                        xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
                        ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
                        zlabel(sprintf('log_{10}(weight %d)', to_all_weights(3)))
                    else
                        error('Unexpected number of active weights.');
                    end
                    title('Agreement between weights');
                    axis equal
                    hold off

                    savefig(...
                        fg,...
                        [name_params_files 'weightsCorrelation.fig'], 'compact'...
                        );
                    close(fg);

                    for w = 1:n_active_weights
                        aw = to_all_weights(w);
                        penalties_w = reshape(points_err_3D(:, :, w), n_patches_ps, 1);

                        fg = figure;
                        hold on
                        for f = 1:n_criteria
                            scatter(penalties_w, points_weights(:, w, f), [], criteria_colors(f, :), 'filled');
                        end
                        legend(criteria_abbrev{:});
                        xlabel(sprintf('log_{10}(Penalty %d)', aw))
                        ylabel(sprintf('log_{10}(Weight %d)', aw))
                        title('Weights as a function of the corresponding penalty term');
                        hold off

                        savefig(...
                            fg,...
                            [name_params_files sprintf('weight%dVsPenalty.fig', aw)], 'compact'...
                            );
                        close(fg);
                    end
                end
               
                for f = 1:n_criteria
                    % Spectral evaluation
                    name_params_f = [name_params, criteria_filenames{f}];
                    name_params_tables_f = sprintf('%s, %s', name_params_tables, criteria_abbrev{f});
                    dataset_params.evaluation.global_spectral.plot_color =...
                        evaluation_plot_colors_admm(color_ind, :);
                    dataset_params.evaluation.global_spectral.plot_marker =...
                        evaluation_plot_markers{...
                            mod(color_ind - 1, length(evaluation_plot_markers)) + 1 ...
                        };
                    dataset_params.evaluation.global_spectral.plot_style =...
                        evaluation_plot_styles{...
                            mod(color_ind - 1, length(evaluation_plot_styles)) + 1 ...
                        };
                    color_ind = color_ind + 1;
                    all_name_params_tables{color_ind} = name_params_tables_f;
                    [e_spectral_table_current, fg_spectral] = evaluateAndSaveSpectral(...
                        estimated_images{f}, I_spectral_gt, bands_spectral, spectral_weights,...
                        dataset_params, name_params_table_gt, name_params_tables_f,...
                        fullfile(output_directory, name_params_f(1:(end-1))),...
                        fg_spectral...
                        );
                    e_spectral_table = union(e_spectral_table_current, e_spectral_table);
                    
                    mrae(pp, ps, no, d, f) = e_spectral_table_current.MRAE_mean;
                    
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
            e_rgb_table, [name_params_figures_gt, 'evaluateRGB.csv']...
        );
        e_rgb_tables{image_number} = e_rgb_table;
        writetable(...
            e_spectral_table, [name_params_figures_gt, 'evaluateSpectral.csv']...
        );
        % Also save completed figures
        evaluateAndSaveSpectral(...
            output_directory, dataset_params, name_params_table_gt, all_name_params_tables, fg_spectral...
        );
        e_spectral_tables{image_number} = e_spectral_table;
    
        % Show and save time results
        fg = figure;
        hold on
        if n_patch_sizes > 1 && n_padding > 1
            for f = 1:n_criteria
                surf(...
                    patch_sizes_grid, paddings_grid, time(:, :, no, d, f),...
                    'FaceAlpha', 0.5, 'FaceColor', criteria_colors(f, :)...
                );
            end
            xlabel('Patch size');
            ylabel('Padding size');
            zlabel('Time to process image [s]');
        elseif n_padding > 1
            for f = 1:n_criteria
                plot(...
                    paddings_grid, time(:, :, no, d, f),...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Padding size');
            ylabel('Time to process image [s]');
        else
            for f = 1:n_criteria
                plot(...
                    patch_sizes_grid, time(:, :, no, d, f),...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Patch size');
            ylabel('Time to process image [s]');
        end
        legend(criteria_abbrev);
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
        if n_patch_sizes > 1 && n_padding > 1
            for f = 1:n_criteria
                surf(...
                    patch_sizes_grid, paddings_grid,...
                    time(:, :, no, d, f) ./ ((patch_sizes_grid + 2 * paddings_grid) .^ 2),...
                    'FaceAlpha', 0.5, 'FaceColor', criteria_colors(f, :)...
                );
            end
            xlabel('Patch size');
            ylabel('Padding size');
            zlabel('Time per pixel in a patch [s]');
        elseif n_padding > 1
            for f = 1:n_criteria
                plot(...
                    paddings_grid, time(:, :, no, d, f) ./ ((patch_sizes_grid + 2 * paddings_grid) .^ 2),...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Padding size');
            ylabel('Time per pixel in a patch [s]');
        else
            for f = 1:n_criteria
                plot(...
                    patch_sizes_grid, time(:, :, no, d, f) ./ ((patch_sizes_grid + 2 * paddings_grid) .^ 2),...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Patch size');
            ylabel('Time per pixel in a patch [s]');
        end
        legend(criteria_abbrev);
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
        if n_patch_sizes > 1 && n_padding > 1
            for f = 1:n_criteria
                surf(...
                    patch_sizes_grid, paddings_grid,...
                    time(:, :, no, d, f) ./ n_patches_rep,...
                    'FaceAlpha', 0.5, 'FaceColor', criteria_colors(f, :)...
                );
            end
            xlabel('Patch size');
            ylabel('Padding size');
            zlabel('Time per patch [s]');
        elseif n_padding > 1
            for f = 1:n_criteria
                plot(...
                    paddings_grid, time(:, :, no, d, f) ./ n_patches_rep,...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Padding size');
            ylabel('Time per patch [s]');
        else
            for f = 1:n_criteria
                plot(...
                    patch_sizes_grid, time(:, :, no, d, f) ./ n_patches_rep,...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Patch size');
            ylabel('Time per patch [s]');
        end
        legend(criteria_abbrev);
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
        
        % Show and save error results
        
        fg = figure;
        hold on
        if n_patch_sizes > 1 && n_padding > 1
            for f = 1:n_criteria
                surf(...
                    patch_sizes_grid, paddings_grid, mrae(:, :, no, d, f),...
                    'FaceAlpha', 0.5, 'FaceColor', criteria_colors(f, :)...
                );
            end
            xlabel('Patch size');
            ylabel('Padding size');
            zlabel(error_name);
        elseif n_padding > 1
            for f = 1:n_criteria
                plot(...
                    paddings_grid, mrae(:, :, no, d, f),...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Padding size');
            ylabel(error_name);
        else
            for f = 1:n_criteria
                plot(...
                    patch_sizes_grid, mrae(:, :, no, d, f),...
                    'Color', criteria_colors(f, :), 'LineWidth', 2,...
                    'Marker', 'o', 'LineStyle', '--'...
                );
            end
            xlabel('Patch size');
            ylabel(error_name);
        end
        legend(criteria_abbrev);
        title(sprintf(...
            '%s for a %d x %d x %d image under dispersion %f and noise %f',...
            error_name,...
            image_sampling_3(1), image_sampling_3(2), image_sampling_3(3),...
            dispersion_px_d, noise_fraction...
        ));
        hold off
        savefig(...
            fg,...
            [name_params_figures_gt error_filename '.fig'], 'compact'...
        );
        close(fg);
        
        image_number = image_number + 1;
        
        if verbose_progress
            disp('Finished evaluation and file output.')
        end
    end
end

%% Show and save minimum (over patch and padding sizes) error results

mrae_min = shiftdim(min(min(mrae, [], 1), [], 2), 2);
[dispersion_grid, noise_fractions_grid] = meshgrid(dispersion_px_all, noise_fractions);
fg = figure;
hold on
if n_dispersion > 1 && n_noise_fractions > 1
    for f = 1:n_criteria
        surf(...
            dispersion_grid, noise_fractions_grid, mrae_min(:, :, f),...
            'FaceAlpha', 0.5, 'FaceColor', criteria_colors(f, :)...
        );
    end
    xlabel('Dispersion magnitude [pixels]');
    ylabel('Noise fraction');
    zlabel(error_name);
elseif n_noise_fractions > 1
    for f = 1:n_criteria
        plot(...
            noise_fractions_grid, mrae_min(:, :, f),...
            'Color', criteria_colors(f, :), 'LineWidth', 2,...
            'Marker', 'o', 'LineStyle', '--'...
        );
    end
    xlabel('Noise fraction');
    ylabel(error_name);
else
    for f = 1:n_criteria
        plot(...
            dispersion_grid, mrae_min(:, :, f),...
            'Color', criteria_colors(f, :), 'LineWidth', 2,...
            'Marker', 'o', 'LineStyle', '--'...
        );
    end
    xlabel('Dispersion magnitude [pixels]');
    ylabel(error_name);
end
legend(criteria_abbrev);
title(sprintf(...
    '%s for a %d x %d x %d image (minimum over padding and patch sizes)',...
    error_name,...
    image_sampling_3(1), image_sampling_3(2), image_sampling_3(3)...
));
hold off
savefig(...
    fg,...
    fullfile(output_directory, [dataset_name, '_', error_filename, '_min.fig']), 'compact'...
);
close(fg);

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
        'bands_spectral',...
        'time',...
        'num_workers',...
        'n_patches',...
        'mrae'...
    } ];
save(save_data_filename, save_variables_list{:});