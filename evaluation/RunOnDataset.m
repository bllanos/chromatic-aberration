%% Evaluate demosaicking, spectral reconstruction, and/or chromatic aberration correction
% Run algorithms on a dataset to evaluate demosaicking, spectral
% reconstruction, and/or chromatic aberration correction
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% The dataset determines the data to be loaded, the algorithms to be
% tested, and the types of evaluations to perform, as encapsulated by the
% 'describeDataset()' function.
%
% The documentation in the scripts 'CorrectByHyperspectralADMM.m' and
% 'CorrectByWarping.m' contains more information on the formats of the
% various types of data associated with the datasets.
%
% This script also runs 'SetFixedParameters.m' to set the values of
% seldomly-changed parameters. These parameters are briefly documented in
% 'SetFixedParameters.m'.
%
% In contrast with 'CorrectByHyperspectralADMM.m', the wavelengths at which
% hyperspectral images are to be sampled are either determined from ground
% truth hyperspectral data, or are otherwise set by 'SetFixedParameters.m',
% but are not loaded from colour space conversion data, or dispersion model
% data.
%
% ## Output
%
% ### Estimated images
%
% The following types of images are created for each input image, depending
% on the image estimation algorithms. The filename of the input image,
% concatenated with a string of parameter information, is represented by
% '*' below:
% - '*_roi.tif' and '*_roi.mat': A cropped version of the input image
%   (stored in the variable 'I_raw'), containing the portion used as input.
%   This region of interest was determined using the domain of the model of
%   dispersion associated with the dataset. If no model of dispersion is
%   associated with the dataset, the cropped region is the entire input
%   image. All of the other output images listed below are limited to the
%   region shown in this output image.
% - '*_latent.mat': The estimated latent spectral image (stored in the
%   variable 'I_latent') corresponding to the input image.
% - '*_rgb.tif' and '*_rgb.mat': A colour image (stored in the variable
%   'I_rgb'). If it was not estimated directly, it was created by
%   converting the latent image to the RGB colour space of the input image.
%
% ### Data file output
%
% #### Intermediate data and parameters
% A '.mat' file containing the following variables, as appropriate:
% - 'bands': A vector containing the wavelengths of the spectral
%   bands used in hyperspectral image estimation.
% - 'bands_spectral': A vector containing the wavelengths of the spectral
%   bands associated with ground truth hyperspectral images.
% - 'sensor_map_resampled': A resampled version of the 'sensor_map'
%   variable loaded from colour space conversion data, used for
%   hyperspectral image estimation. 'sensor_map_resampled' is the spectral
%   response functions of the camera (or, more generally, of the output
%   3-channel colour space) approximated at the wavelengths in `bands`.
% - 'sensor_map_spectral': A resampled version of the 'sensor_map'
%   variable loaded from colour space conversion data, used to convert
%   ground truth spectral images to color. 'sensor_map_spectral' is the
%   spectral response functions of the camera (or, more generally, of the
%   output 3-channel colour space) approximated at the wavelengths in
%   `bands_spectral`.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`, which is initialized in this file, and then augmented
% by 'SetFixedParameters.m'.
%
% The file is saved as 'RunOnDataset_${dataset_name}_evaluateSpectral.mat'.
%
% #### Evaluation results
%
% For each image, RGB error metrics and (if applicable) spectral error
% metrics are output in the form of CSV files. Each CSV file contains
% results for all algorithms tested. The RGB error metrics are saved as
% '*_evaluateRGB.csv', whereas the spectral error metrics are saved as
% '*_evaluateSpectral.csv'.
%
% Error metrics are also aggregated across images, and saved as
% '${dataset_name}_evaluateRGB.csv' and
% '${dataset_name}_evaluateSpectral.csv'.
%
% ## Notes
% - This script only uses the first row of `patch_sizes`, and the first
%   element of `paddings`, defined in 'SetFixedParameters.m'.
% - This script ignores the `downsampling_factor` parameter defined in
%   'SetFixedParameters.m'.
%
% ## References
%
% The adaptive residual interpolation demosaicking algorithm
% ('third_party/Sensors_ARI/') was developed by Yusuke Monno and Daisuke
% Kiku, and was retrieved from
% http://www.ok.sc.e.titech.ac.jp/res/DM/RI.html
%
% It is described in:
%
%   Yusuke Monno, Daisuke Kiku, Masayuki Tanaka, and Masatoshi Okutomi,
%   "Adaptive Residual Interpolation for Color and Multispectral Image
%   Demosaicking," Sensors, vol.17, no.12, pp.2787-1-21, 2017.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 27, 2018

% List of parameters to save with results
parameters_list = {
        'dataset_name',...
        'output_directory'...
    };

%% Input data and parameters

dataset_name = 'kodak';

% ADMM family algorithms to run
admm_algorithms.spectralL1L1 = struct(...
    'name', 'Spectral L1L1',... % Pretty name for tables and graphs
    'file', 'L1L1',... % Short name for filenames
    'enabled', false,... % Whether or not to run the algorithm
    'spectral', true,... % Estimate spectral or colour images
    'priors', [true, true, false],... % Priors to enable (`false` means the corresponding weight will be zero)
    'options', struct('norms', [true, true, false], 'nonneg', false)... % Custom options for baek2017Algorithm2()
);
admm_algorithms.spectralL1L1NonNeg = struct(...
    'name', 'Spectral L1L1NonNeg',...
    'file', 'L1L1NonNeg',...
    'enabled', false,...
    'spectral', true,...
    'priors', [true, true, false],...
    'options', struct('norms', [true, true, false], 'nonneg', true)...
);
admm_algorithms.spectralL1L2NonNeg = struct(...
    'name', 'Spectral L1L2NonNeg',...
    'file', 'L1L2NonNeg',...
    'enabled', false,...
    'spectral', true,...
    'priors', [true, true, false],...
    'options', struct('norms', [true, false, false], 'nonneg', true)...
);
admm_algorithms.spectralL2L1NonNeg = struct(...
    'name', 'Spectral L2L1NonNeg',...
    'file', 'L2L1NonNeg',...
    'enabled', false,...
    'spectral', true,...
    'priors', [true, true, false],...
    'options', struct('norms', [false, true, false], 'nonneg', true)...
);
admm_algorithms.spectralL2L2NonNeg = struct(...
    'name', 'Spectral L2L2NonNeg',...
    'file', 'L2L2NonNeg',...
    'enabled', false,...
    'spectral', true,...
    'priors', [true, true, false],...
    'options', struct('norms', [false, false, false], 'nonneg', true)...
);
admm_algorithms.spectralMinNorm = struct(...
    'name', 'Spectral MinNorm',...
    'file', 'MinNorm',...
    'enabled', false,...
    'spectral', true,...
    'priors', [false, false, false],...
    'options', struct('norms', [false, false, false], 'nonneg', false)...
);
admm_algorithms.spectralL1L1L1NonNeg = struct(...
    'name', 'Spectral L1L1L1NonNeg',...
    'file', 'L1L1L1NonNeg',...
    'enabled', false,...
    'spectral', true,...
    'priors', [true, true, true],...
    'options', struct('norms', [true, true, true], 'nonneg', true)...
);

admm_algorithms.colorL1L1 = struct(...
    'name', 'Color L1L1',...
    'file', 'L1L1',...
    'enabled', false,...
    'spectral', false,...
    'priors', [true, true, false],...
    'options', struct('norms', [true, true, false], 'nonneg', false)...
);
admm_algorithms.colorL1L1NonNeg = struct(...
    'name', 'Color L1L1NonNeg',...
    'file', 'L1L1NonNeg',...
    'enabled', true,...
    'spectral', false,...
    'priors', [true, true, false],...
    'options', struct('norms', [true, true, false], 'nonneg', true)...
);
admm_algorithms.colorL1L2NonNeg = struct(...
    'name', 'Color L1L2NonNeg',...
    'file', 'L1L2NonNeg',...
    'enabled', false,...
    'spectral', false,...
    'priors', [true, true, false],...
    'options', struct('norms', [true, false, false], 'nonneg', true)...
);
admm_algorithms.colorL2L1NonNeg = struct(...
    'name', 'Color L2L1NonNeg',...
    'file', 'L2L1NonNeg',...
    'enabled', false,...
    'spectral', false,...
    'priors', [true, true, false],...
    'options', struct('norms', [false, true, false], 'nonneg', true)...
);
admm_algorithms.colorL2L2NonNeg = struct(...
    'name', 'Color L2L2NonNeg',...
    'file', 'L2L2NonNeg',...
    'enabled', false,...
    'spectral', false,...
    'priors', [true, true, false],...
    'options', struct('norms', [false, false, false], 'nonneg', true)...
);
admm_algorithms.colorMinNorm = struct(...
    'name', 'Color MinNorm',...
    'file', 'MinNorm',...
    'enabled', false,...
    'spectral', false,...
    'priors', [false, false, false],...
    'options', struct('norms', [false, false, false], 'nonneg', false)...
);
admm_algorithms.colorL1L1L1NonNeg = struct(...
    'name', 'Color L1L1L1NonNeg',...
    'file', 'L1L1L1NonNeg',...
    'enabled', false,...
    'spectral', false,...
    'priors', [true, true, true],...
    'options', struct('norms', [true, true, true], 'nonneg', true)...
);

% Demosaicking algorithms
demosaic_algorithms.bilinear = struct(...
    'name', 'Bilinear demosaicking',... % Bilinear interpolation for demosaicking
    'file', 'bilinear',...
    'enabled', false,...
    'fn', @bilinearDemosaic... % Function to call
);
demosaic_algorithms.matlab = struct(...
    'name', 'MATLAB demosaicking',... % MATLAB's built-in demosaic() function
    'file', 'MATLABdemosaic',...
    'enabled', false,...
    'fn', 'matlab'... % Special handling is required for this algorithm
);
demosaic_algorithms.ari = struct(...
    'name', 'ARI demosaicking',... % Adaptive residual interpolation
    'file', 'ARI',...
    'enabled', false,...
    'fn', 'ARI'... % Special handling is required for this algorithm
);

% ## Parameters for creating radiance images

% CIE D-illuminant
illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
illuminant_temperature = 6504; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
illuminant_name = 'd65';

% Colour channel to use for radiance normalization
normalization_channel = 2;

% ## Operational parameters

% Patch size and padding to use for converting spectral images to colour
% and RAW images
imageFormationOptions.patch_size = [100, 100];
imageFormationOptions.padding = 10;

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% Produce console output to describe the processing in this script
verbose = true;

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% Check for problematic parameters
if add_border
    % Estimating a border area results in images which are usually not
    % registered with the ground truth.
    error('Estimating a border around images prevents quantitative evaluation');
end

%% Load the illuminant

illuminant_data = csvread(illuminant_filename);
bands_illuminant = illuminant_data(:, 1);
S_illuminant = illuminant_data(:, 2:end);
spd_illuminant = ciedIlluminant(...
    illuminant_temperature, bands_illuminant, S_illuminant, bands_illuminant...
);

%% Look up the dataset

dp = describeDataset(dataset_name);

%% Find and/or prepare to generate the dataset images

has_raw = ~isempty(dp.raw_images_wildcard);
has_rgb = ~isempty(dp.rgb_images_wildcard);
has_spectral = ~isempty(dp.spectral_images_wildcard);
has_color_map = ~isempty(dp.color_map);

if ~has_spectral && ~has_rgb
    error('Ground truth (colour or spectral images) must be provided for evaluation.')
end
if ~has_raw && ~has_rgb && (~has_spectral || (has_spectral && ~has_color_map))
    error('RAW images are not available, and cannot be generated.')
end

bands_color = [];
if has_color_map
    % Load colour conversion data
    bands_script = bands;
    bands = [];
    
    optional_variable = 'bands';
    model_variables_required = { 'sensor_map', 'channel_mode' };
    load(dp.color_map, model_variables_required{:}, optional_variable);
    if ~all(ismember(model_variables_required, who))
        error('One or more of the required colour space conversion variables is not loaded.')
    end
    bands_color = bands;
    bands = bands_script;
    
    % Compare with colour space conversion data
    n_bands = length(bands);
    n_bands_sensor_map = size(sensor_map, 2);
    resample_bands = false;
    if ~isempty(bands_color)
        bands_for_interp = bands_color;
        if n_bands ~= length(bands_color) || any(bands ~= bands_color)
            % Resampling is needed
            resample_bands = true;
        end
    elseif n_bands_sensor_map ~= n_bands
        % Resampling is needed, but will be "blind"
        resample_bands = true;
        bands_for_interp = linspace(bands(1), bands(end), n_bands_sensor_map);
    end
    % Resample colour space conversion data if necessary
    if resample_bands
        [sensor_map_resampled, bands] = resampleArrays(...
            bands_for_interp, sensor_map.', bands,...
            bands_interp_method...
            );
        n_bands = length(bands);
        sensor_map_resampled = sensor_map_resampled.';
    else
        sensor_map_resampled = sensor_map;
    end
end

bands_spectral = [];
sensor_map_spectral = [];
if has_spectral
    spectral_filenames = listFiles(dp.spectral_images_wildcard);
    n_images_spectral = length(spectral_filenames);
    n_images = n_images_spectral;
    names = trimCommon(spectral_filenames);
    
    bands = [];
    load(dp.wavelengths, optional_variable);
    if isempty(bands)
        error('No wavelength band information is associated with the spectral images.')
    end
    bands_spectral = bands;
    bands = bands_script;
    
    if has_color_map
        % Check if quantitative evaluation of spectral images is possible
        can_evaluate_spectral = (length(bands_spectral) == n_bands) && all(bands_spectral == bands); 
        if can_evaluate_spectral
            sensor_map_spectral = sensor_map_resampled;
        else
            warning(['Quantitative evaluation of latent images is not possibl'...
                'e because they will be produced at different wavelength'...
                's from the true latent images.']...
            );
        
            % Allow for conversion to colour images
            [sensor_map_spectral, bands_spectral] = resampleArrays(...
                bands_for_interp, sensor_map.', bands_spectral,...
                bands_interp_method...
                );
            sensor_map_spectral = sensor_map_spectral.';
        end
    else
        can_evaluate_spectral = false;
    end
end
if has_rgb
    rgb_filenames = listFiles(dp.rgb_images_wildcard);
    n_images_rgb = length(rgb_filenames);
    rgb_names = trimCommon(rgb_filenames);
    if has_spectral
        if (n_images_spectral ~= n_images_rgb)
            error('Mismatched number of spectral and colour images.');
        end
        for i = 1:n_images_rgb
            if ~strcmp(rgb_names{i}, names{i})
                error('Not all spectral image filenames and colour image filenames match.');
            end
        end
    end
    n_images = n_images_rgb;
    names = rgb_names;
end
if has_raw
    raw_filenames = listFiles(dp.raw_images_wildcard);
    n_images_raw = length(raw_filenames);
    raw_names = trimCommon(raw_filenames);
    if (n_images ~= n_images_raw)
        error('Mismatched number of spectral/colour and RAW images.');
    end
    for i = 1:n_images
        if ~strcmp(raw_names{i}, names{i})
            error('Not all spectral/colour image filenames and RAW image filenames match.');
        end
    end
end

%% Load dispersion models

has_dispersion_rgb = ~isempty(dp.dispersion_rgb_forward) && ~isempty(dp.dispersion_rgb_reverse);
has_dispersion_spectral = ~isempty(dp.dispersion_spectral_reverse);

if has_dispersion_rgb
    [...
        dd_rgb_forward, ~, td_rgb_forward...
    ] = loadDispersionModel(dp.dispersion_rgb_forward, true, false);
    [...
        dd_rgb_reverse, ~, td_rgb_reverse...
    ] = loadDispersionModel(dp.dispersion_rgb_reverse, false, false);
end
if has_dispersion_spectral
    [...
        dd_spectral_reverse, ~, td_spectral_reverse...
    ] = loadDispersionModel(dp.dispersion_spectral_reverse, false, false);
end

%% Process the images

n_channels_rgb = 3;
bands_rgb = 1:n_channels_rgb;
sensor_map_rgb = eye(n_channels_rgb);
n_weights = size(weights, 1);
patch_size = patch_sizes(1, :);
padding = paddings(1);
imageFormationOptions.int_method = int_method;

e_rgb_tables = cell(n_images, 1);
e_spectral_tables = cell(n_images, 1);

% Fixed options for ADMM
solvePatchesOptions.add_border = add_border;
baek2017Algorithm2Options.add_border = false;
solvePatchesOptions.patch_size = patch_size;
solvePatchesOptions.padding = padding;

for i = 1:n_images
    if verbose
        fprintf('[RunOnDataset, image %d] Starting\n', i);
    end

    % Generate or load input images, and instantiate dispersion information
    if has_spectral
        I_spectral_gt = loadImage(spectral_filenames{i}, dp.spectral_images_variable);
        
        % Convert to radiance images, if required
        if dp.spectral_reflectances
            original_size = size(I_spectral_gt);
            reflectances = reshape(I_spectral_gt, [], length(bands_spectral)).';
            
            % Resample the illuminant, to avoid resampling the image (which
            % is computationally more expensive)
            [spd_illuminant_i, bands_illuminant_i] = resampleArrays(...
                bands_illuminant, spd_illuminant, bands_spectral,...
                'spline', 'none'...
            );
            if length(bands_illuminant_i) ~= length(bands_spectral) ||...
                    any(bands_illuminant_i ~= bands_spectral)
                error(['Cannot convert spectral image %d to a radiance image',...
                    ', because the illuminant is not sampled at a sufficientl',...
                    'y wide range of wavelengths.'], i);
            end
        
            [~, ~, rad_normalized] = reflectanceToRadiance(...
                bands_illuminant_i, spd_illuminant_i,...
                bands_spectral, reflectances,...
                bands_color, sensor_map.',...
                normalization_channel, int_method...
            );
            
            I_spectral_gt = reshape(rad_normalized.', original_size);
        end
    
        if has_dispersion_spectral
            [df_spectral_reverse, I_spectral_gt] = makeDispersionForImage(...
                dd_spectral_reverse, I_spectral_gt, td_spectral_reverse...
            );
        else
            df_spectral_reverse = [];
        end
        image_sampling = [size(I_spectral_gt, 1), size(I_spectral_gt, 2)];
        I_spectral_lin = reshape(I_spectral_gt, [], 1);
    elseif has_color_map && has_dispersion_spectral
        df_spectral_reverse = makeDispersionForImage(dd_spectral_reverse);
    elseif has_color_map
        df_spectral_reverse = [];
    end
    
    if has_spectral && has_color_map
        [...
            I_rgb_gt_simulated, I_rgb_gt_warped, I_raw_gt_simulated,...
            I_spectral_gt_warped...
        ] = imageFormation(...
            I_spectral_gt, sensor_map_spectral, bands_spectral,...
            imageFormationOptions, df_spectral_reverse, bayer_pattern...
        );
    else
        I_rgb_gt_warped = [];
        I_spectral_gt_warped = [];
    end

    if has_rgb
        I_rgb_gt = loadImage(rgb_filenames{i}, dp.rgb_images_variable);
        if has_dispersion_rgb
            [df_rgb_reverse, I_rgb_gt] = makeDispersionForImage(...
                    dd_rgb_reverse, I_rgb_gt, td_rgb_reverse...
                );
            df_rgb_forward = makeDispersionForImage(...
                    dd_rgb_forward, I_rgb_gt, td_rgb_forward...
                );
        else
            df_rgb_reverse = [];
            df_rgb_forward = [];
        end
        if has_spectral
            if any([size(I_rgb_gt, 1), size(I_rgb_gt, 2)] ~= image_sampling)
                error('The colour version of %s has different spatial dimensions from its spectral version.',...
                    names{i}...
                );
            end
        else
            image_sampling = [size(I_rgb_gt, 1), size(I_rgb_gt, 2)];
        end
        I_rgb_lin = reshape(I_rgb_gt, [], 1);
    elseif has_spectral && has_color_map
        I_rgb_gt = I_rgb_gt_simulated;
        I_raw_gt = I_raw_gt_simulated;
    end
    if ~has_rgb && has_dispersion_rgb
        df_rgb_reverse = makeDispersionForImage(dd_rgb_reverse);
        df_rgb_forward = makeDispersionForImage(dd_rgb_forward);
    elseif ~has_dispersion_rgb
        df_rgb_reverse = [];
        df_rgb_forward = [];
    end
    
    if has_raw
        I_raw_gt = loadImage(raw_filenames{i}, dp.raw_images_variable);
        if ~ismatrix(I_raw_gt)
            error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
        end
        
        % Crop to the region of valid dispersion
        roi = [];
        if has_spectral && has_dispersion_spectral
            roi = modelSpaceTransform(...
                [size(I_raw_gt, 1), size(I_raw_gt, 2)], td_spectral_reverse.model_space, td_spectral_reverse.fill...
                );
        elseif has_rgb && has_dispersion_rgb
            roi = modelSpaceTransform(...
                [size(I_raw_gt, 1), size(I_raw_gt, 2)], td_rgb_reverse.model_space, td_rgb_reverse.fill...
                );
        end
        if ~isempty(roi)
            I_raw_gt = I_raw_gt(roi(1):roi(2), roi(3):roi(4), :);
        end
        
        if any([size(I_raw_gt, 1), size(I_raw_gt, 2)] ~= image_sampling)
            error('The RAW version of %s has different spatial dimensions from its other versions.',...
                names{i}...
                );
        end
    else
        if has_rgb
            if has_dispersion_rgb
                if verbose
                    fprintf('[RunOnDataset, image %d] Calculating the reverse colour dispersion matrix...\n', i);
                end
                W_reverse = dispersionfunToMatrix(...
                    df_rgb_reverse, bands_rgb, image_sampling, image_sampling,...
                    [0, 0, image_sampling(2),  image_sampling(1)], true...
                    );
                if verbose
                    fprintf('\t...done\n');
                end
                I_rgb_warped = warpImage(I_rgb_gt, W_reverse, image_sampling);
            else
                I_rgb_warped = I_rgb_gt;
            end
            I_raw_gt = mosaic(I_rgb_warped, bayer_pattern);
        end
    end        
        
    saveImages(...
        output_directory, names{i},...
        I_raw_gt, '_roi', 'I_raw'...
    );
    
    % Compare the aberrated image to the original
    
    if isempty(I_rgb_gt_warped)
        e_rgb_table = [];
    else
        % Evaluate the aberrated image as a baseline
        e_rgb_table = evaluateAndSaveRGB(...
            I_rgb_gt_warped, I_rgb_gt, dp, names{i}, 'Aberrated',...
            fullfile(output_directory, [names{i} '_aberrated'])...
        );
    end
    
    admm_algorithm_fields = fieldnames(admm_algorithms);
    n_admm_algorithms = length(admm_algorithm_fields);
    n_spectral_evaluations = 0;
    if has_spectral && can_evaluate_spectral
        for f = 1:n_admm_algorithms
            algorithm = admm_algorithms.(admm_algorithm_fields{f});
            if algorithm.enabled && algorithm.spectral
                n_spectral_evaluations = n_spectral_evaluations + 1;
            end
        end
        n_spectral_evaluations = n_spectral_evaluations * n_weights;
    end
    if ~isempty(I_spectral_gt_warped)
        n_spectral_evaluations = n_spectral_evaluations + 1;
    end
    if n_spectral_evaluations > 0
        evaluation_plot_colors = jet(n_spectral_evaluations);
        if isempty(I_spectral_gt_warped)
            evaluation_plot_colors_admm = evaluation_plot_colors;
        else
            evaluation_plot_colors_admm = evaluation_plot_colors(2:end, :);
        end
    end
    if isempty(I_spectral_gt_warped)
        e_spectral_table = [];
        fg_spectral = struct;
        all_alg_names = {};
    else
        dp.evaluation.global_spectral.plot_color = evaluation_plot_colors(1, :);
        all_alg_names = {'Aberrated'};
        [e_spectral_table, fg_spectral] = evaluateAndSaveSpectral(...
            I_spectral_gt_warped, I_spectral_gt, bands_spectral,...
            dp, names{i}, all_alg_names{1},...
            fullfile(output_directory, [names{i} '_aberrated'])...
        );
    end
    
    % Run the algorithms
    
    % ADMM
    color_ind = 1;
    for w = 1:n_weights
        for f = 1:n_admm_algorithms
            algorithm = admm_algorithms.(admm_algorithm_fields{f});
            if ~algorithm.enabled
                continue;
            end
            
            if algorithm.spectral
                if has_color_map
                    if channel_mode
                        baek2017Algorithm2Options.int_method = 'none';
                        solvePatchesOptions.int_method = 'none';
                    else
                        baek2017Algorithm2Options.int_method = int_method;
                        solvePatchesOptions.int_method = int_method;
                    end
                else
                    continue;
                end
            else
                baek2017Algorithm2Options.int_method = 'none';
                solvePatchesOptions.int_method = 'none';
            end
            
            weights_f = weights(w, :);
            weights_f(~algorithm.priors) = 0;
            
            baek2017Algorithm2Options = mergeStructs(...
                baek2017Algorithm2Options, algorithm.options, false, true...
            );
    
            name_params = sprintf(...
                '%s_patch%dx%d_pad%d_weights%ew%ew%e_',...
                algorithm.file, patch_size(1), patch_size(2), padding,...
                weights_f(1), weights_f(2), weights_f(3)...
                );
            alg_name_params = sprintf(...
                '%s, patch %d x %d, padding %d, weights (%g, %g, %g)',...
                algorithm.file, patch_size(1), patch_size(2), padding,...
                weights_f(1), weights_f(2), weights_f(3)...
                );
            if algorithm.spectral
                name_params = [...
                    names{i}, sprintf('_bands%d_', n_bands), name_params...
                    ];
                alg_name_params = [...
                    alg_name_params, sprintf(', %d bands', n_bands)...
                    ];
                [...
                    I_latent, ~, I_rgb...
                ] = solvePatchesAligned(...
                    I_raw_gt, bayer_pattern, df_spectral_reverse,...
                    sensor_map_resampled,...
                    bands, solvePatchesOptions, @baek2017Algorithm2,...
                    {...
                        weights_f, rho,...
                        baek2017Algorithm2Options, baek2017Algorithm2Verbose...
                    }...
                );
                saveImages(...
                    output_directory, name_params,...
                    I_latent, 'latent', 'I_latent',...
                    I_rgb, 'latent_rgb', 'I_rgb'...
                );
            
                % Spectral evaluation
                if can_evaluate_spectral
                    dp.evaluation.global_spectral.plot_color =...
                        evaluation_plot_colors_admm(color_ind, :);
                    color_ind = color_ind + 1;
                    all_alg_names{end + 1} = algorithm.name;
                    [e_spectral_table_current, fg_spectral] = evaluateAndSaveSpectral(...
                        I_latent, I_spectral_gt, bands, dp, names{i},...
                        alg_name_params,...
                        fullfile(output_directory, name_params(1:(end-1))),...
                        fg_spectral...
                    );
                    if ~isempty(e_spectral_table)
                        e_spectral_table = union(e_spectral_table_current, e_spectral_table);
                    else
                        e_spectral_table = e_spectral_table_current;
                    end
                end
            else
                name_params = [...
                    names{i}, '_RGB_', name_params...
                ];
                alg_name_params = [...
                    alg_name_params, ', RGB'...
                ];
                I_rgb = solvePatchesAligned(...
                    I_raw_gt, bayer_pattern, df_rgb_reverse,...
                    sensor_map_rgb,...
                    bands_rgb, solvePatchesOptions, @baek2017Algorithm2,...
                    {...
                        weights_f, rho,...
                        baek2017Algorithm2Options, baek2017Algorithm2Verbose...
                    }...
                );
                saveImages(...
                    output_directory, name_params,...
                    I_rgb, 'rgb', 'I_rgb'...
                );
            end
            
            % RGB evaluation
            e_rgb_table_current = evaluateAndSaveRGB(...
                I_rgb, I_rgb_gt, dp, names{i}, alg_name_params,...
                fullfile(output_directory, name_params(1:(end-1)))...
            );
            if ~isempty(e_rgb_table)
                e_rgb_table = union(e_rgb_table_current, e_rgb_table);
            else
                e_rgb_table = e_rgb_table_current;
            end
        end
    end
    
    % Demosaicking and colour channel warping
    demosaic_algorithm_fields = fieldnames(demosaic_algorithms);
    W_forward = [];
    for f = 1:length(demosaic_algorithm_fields)
        algorithm = demosaic_algorithms.(demosaic_algorithm_fields{f});
        if ~algorithm.enabled
            continue;
        end
    
        if ischar(algorithm.fn)
            if strcmp(algorithm.fn, 'matlab')
                I_raw_int = im2uint16(I_raw_gt);
                I_rgb_warped = im2double(demosaic(I_raw_int, bayer_pattern));
            elseif strcmp(algorithm.fn, 'ARI')
                I_rgb_warped = demosaic_ARI(...
                    repmat(I_raw_gt, 1, 1, n_channels_rgb), bayer_pattern...
                );
            else
                error('Unrecognized demosaicking algorithm name.');
            end
        else
            I_rgb_warped = algorithm.fn(I_raw_gt, bayer_pattern);
        end
        saveImages(...
            output_directory, names{i},...
            I_rgb_warped, sprintf('_%s', algorithm.file), 'I_rgb'...
        );
    
        % RGB evaluation
        e_rgb_table_current = evaluateAndSaveRGB(...
            I_rgb_warped, I_rgb_gt, dp, names{i}, algorithm.name,...
            fullfile(output_directory, [names{i} '_' algorithm.file])...
        );
        if ~isempty(e_rgb_table)
            e_rgb_table = union(e_rgb_table_current, e_rgb_table);
        else
            e_rgb_table = e_rgb_table_current;
        end
            
        if has_dispersion_rgb
            
            if isempty(W_forward)
                if verbose
                    fprintf('[RunOnDataset, image %d] Calculating the forward colour dispersion matrix...\n', i);
                end
                W_forward = dispersionfunToMatrix(...
                    df_rgb_forward, bands_rgb, image_sampling, image_sampling,...
                    [0, 0, image_sampling(2),  image_sampling(1)], false...
                    );
                if verbose
                    fprintf('\t...done\n');
                end
            end
    
            I_rgb = warpImage(I_rgb_warped, W_forward, image_sampling);
            saveImages(...
                output_directory, names{i},...
                I_rgb, sprintf('_%s_channelWarp', algorithm.file), 'I_rgb'...
            );
        
            % RGB evaluation
            e_rgb_table_current = evaluateAndSaveRGB(...
                I_rgb, I_rgb_gt, dp, names{i},...
                sprintf('%s, warp-corrected', algorithm.name),...
                fullfile(output_directory, [names{i} '_' algorithm.file '_channelWarp'])...
            );
            if ~isempty(e_rgb_table)
                e_rgb_table = union(e_rgb_table_current, e_rgb_table);
            else
                e_rgb_table = e_rgb_table_current;
            end
        end
    end

    % Write evaluations to a file
    if ~isempty(e_rgb_table)
        writetable(...
            e_rgb_table,...
            fullfile(output_directory, [names{i}, '_evaluateRGB.csv'])...
        );
        e_rgb_tables{i} = e_rgb_table;
    end
    if ~isempty(e_spectral_table)
        writetable(...
            e_spectral_table,...
            fullfile(output_directory, [names{i}, '_evaluateSpectral.csv'])...
        );
        % Also save completed figures
        evaluateAndSaveSpectral(output_directory, dp, names{i}, all_alg_names, fg_spectral);
        e_spectral_tables{i} = e_spectral_table;
    end

    if verbose
        fprintf('[RunOnDataset, image %d] Finished\n', i);
    end
end

%% Save results for all images

e_rgb_tables = e_rgb_tables(~cellfun(@isempty, e_rgb_tables, 'UniformOutput', true));
if ~isempty(e_rgb_tables)
    e_rgb_summary_table = mergeRGBTables(e_rgb_tables);
    writetable(...
        e_rgb_summary_table,...
        fullfile(output_directory, [dataset_name, '_evaluateRGB.csv'])...
    );
end

e_spectral_tables = e_spectral_tables(~cellfun(@isempty, e_spectral_tables, 'UniformOutput', true));
if ~isempty(e_spectral_tables)
    e_spectral_summary_table = mergeSpectralTables(e_spectral_tables);
    writetable(...
        e_spectral_summary_table,...
        fullfile(output_directory, [dataset_name, '_evaluateSpectral.csv'])...
    );
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {'bands'} ];
if has_spectral
    save_variables_list = [save_variables_list, {'bands_spectral'}];
end
if has_color_map
    save_variables_list = [save_variables_list, {'sensor_map_resampled'}];
    if has_spectral
        save_variables_list = [save_variables_list, {'sensor_map_spectral'}];
    end
end
save_data_filename = fullfile(output_directory, ['RunOnDataset_' dataset_name '.mat']);
save(save_data_filename, save_variables_list{:});