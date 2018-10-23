%% Select regularization weights for images in a dataset
% Choose regularization weights for random image patches to determine the
% best weights for each algorithm, and how much the optimal weights vary
% within and between images.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% The dataset determines the data to be loaded, and algorithms to be
% tested, as encapsulated by the 'describeDataset()' function.
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
% ### Data and parameters
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
% - 'admm_algorithms': A structure describing the algorithms for which
%   regularization weights were selected. 'admm_algorithms' is created by
%   'SetAlgorithms.m', and then each algorithm is given two additional
%   fields by this script:
%   - 'mdc_weights': Regularization weights selected using the minimum
%     distance criterion of Song et al. 2016. 'mdc_weights' is a 2D array
%     with one row per image in the dataset.
%   - 'mse_weights': Regularization weights selected to minimize the mean
%     square error with respect to the true images. 'mse_weights' is a 2D
%     array with one row per image in the dataset.
% - 'corners': A three-dimensional array containing the top-left corners of
%   the image patches used to select regularization weights. `corners(pc,
%   :, i)` is a two-element vector containing the row and column indices,
%   respectively, of the i-th image's pc-th patch's top-left corner.
% - 'patch_penalties_spectral_L1': A three-dimensional array, containing
%   the L1 norms of the regularization penalties evaluated on the true
%   spectral image patches. `patch_penalties_spectral_L1(pc, w, i)`
%   corresponds to the w-th enabled regularization term for the pc-th patch
%   of the i-th image.
% - 'patch_penalties_spectral_L2': An array similar to
%   'patch_penalties_spectral_L1', but which contains L2 norms of
%   regularization penalties.
% - 'patch_penalties_rgb_L1': An array similar to
%   'patch_penalties_spectral_L1', but which corresponds to regularization
%   penalties evaluated on colour versions of the images.
% - 'patch_penalties_rgb_L2': An array similar to 'patch_penalties_rgb_L1',
%   but which contains L2 norms of regularization penalties.
% - 'all_mdc_weights': Regularization weights selected using the method of
%   Song et al. 2016. `all_mdc_weights(pc, :, i, f)` is the vector of
%   regularization weights selected for the f-th algorithm on the pc-th
%   patch of the i-th image
% - 'all_mse_weights': An array with the same layout as 'all_mdc_weights',
%   but which contains regularization weights selected by minimizing the
%   true mean square error.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`, which is initialized in this file, and then augmented
% by 'SetFixedParameters.m'.
%
% The file is saved as 'SelectWeightsForDataset_${dataset_name}.mat'.
%
% ### Graphical output
%
% Figures are generated for each ADMM algorithm showing the distributions
% of weights selected by the minimum distance and mean square error
% criteria, and are saved to '.fig' files.
%
% ## Notes
% - This script only uses the first row of `patch_sizes`, and the first
%   element of `paddings`, defined in 'SetFixedParameters.m', by using
%   `solvePatchesADMMOptions.patch_options`.
% - Regularization weights will not be selected for algorithms which are
%   disabled (using the 'enabled' fields in the structures describing the
%   algorithms). Instead, values of zero will be output for disabled
%   algorithms in 'all_mdc_weights' and 'all_mse_weights', and the
%   'mdc_weights' and 'mse_weights' fields will not be added to their
%   structures in 'admm_algorithms'.
% - Only image patches that are not clipped by the image edges will be
%   chosen for selecting regularization weights.
% - Image patch spectral derivatives, used for graphical output, not for
%   regularization weights selection, are computed according to the
%   'full_GLambda' field of the 'solvePatchesADMMOptions.admm_options'
%   structure defined in 'SetFixedParameters.m', rather than according to
%   per-algorithm options.
%
% ## References
%
% The following article discusses the grid-search method for minimizing the
% minimum distance function:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 10, 2018

% List of parameters to save with results
parameters_list = {
        'dataset_name',...
        'n_patches',...
        'output_directory'...
    };

%% Input data and parameters

dataset_name = 'kaist-crop';

% Describe algorithms to run
run('SetAlgorithms.m')

% Number of patches to select for each image
n_patches = 10;

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% Produce console output to describe the processing in this script
verbose = true;

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Check for problematic parameters

if use_fixed_weights
    error('Weights should be fixed by running ''SelectWeightsForDataset.m'', not using the `use_fixed_weights` parameter in ''SetFixedParameters.m''');
end

%% Preprocess the dataset

dp = describeDataset(dataset_name);

run('PreprocessDataset.m')

if has_spectral && ~can_evaluate_spectral
    error('The estimated spectral images must have the same channels as the true spectral images.');
end

%% Prepare for patch processing

n_weights = length(solvePatchesADMMOptions.reg_options.enabled);

patch_size = solvePatchesADMMOptions.patch_options.patch_size;
padding = solvePatchesADMMOptions.patch_options.padding;
full_patch_size = patch_size + padding * 2;

if has_spectral
    image_sampling_patch_spectral = [patch_size, n_bands];
    n_spectral_weights = n_weights - 1;
    patch_operators_spectral = cell(n_spectral_weights, 1);
    for w = 1:n_weights
        if w > n_spectral_weights
            % The anti-mosaicking prior is applied in RGB space
            continue;
        end
        if w == 1 || w == 2
            G = spatialGradient(image_sampling_patch_spectral);
        end
        if w == 2
            G_lambda = spectralGradient(image_sampling_patch_spectral, solvePatchesADMMOptions.admm_options.full_GLambda);
            G_lambda_sz1 = size(G_lambda, 1);
            G_lambda_sz2 = size(G_lambda, 2);
            % The product `G_lambda * G_xy` must be defined, so `G_lambda` needs to be
            % replicated to operate on both the x and y-gradients.
            G_lambda = [
                G_lambda, sparse(G_lambda_sz1, G_lambda_sz2);
                sparse(G_lambda_sz1, G_lambda_sz2), G_lambda
                ]; %#ok<AGROW>
            G = G_lambda * G;
        end
        patch_operators_spectral{w} = G;
    end
end

patch_operators_rgb = cell(n_weights, 1);
image_sampling_patch_rgb = [patch_size, n_channels_rgb];
for w = 1:n_weights
    if w == 1 || w == 2
        G = spatialGradient(image_sampling_patch_rgb);
    end
    if w == 2
        G_lambda = spectralGradient(image_sampling_patch_rgb, solvePatchesADMMOptions.admm_options.full_GLambda);
        G_lambda_sz1 = size(G_lambda, 1);
        G_lambda_sz2 = size(G_lambda, 2);
        % The product `G_lambda * G_xy` must be defined, so `G_lambda` needs to be
        % replicated to operate on both the x and y-gradients.
        G_lambda = [
            G_lambda, sparse(G_lambda_sz1, G_lambda_sz2);
            sparse(G_lambda_sz1, G_lambda_sz2), G_lambda
            ]; %#ok<AGROW>
        G = G_lambda * G;
    end
    if w == 3
        % The Bayer pattern code changes depending on the location of the
        % patch, but only if the patch has an odd integer vertical or
        % horizontal offset relative to the image origin
        G = antiMosaicMatrix(patch_size, bayer_pattern);
    end
    patch_operators_rgb{w} = G;
end

admm_algorithm_fields = fieldnames(admm_algorithms);
n_admm_algorithms = length(admm_algorithm_fields);
for f = 1:n_admm_algorithms
    algorithm = admm_algorithms.(admm_algorithm_fields{f});
    if ~algorithm.enabled
        continue;
    end
    admm_algorithms.(admm_algorithm_fields{f}).mdc_weights = zeros(n_images, n_weights); 
    admm_algorithms.(admm_algorithm_fields{f}).mse_weights = zeros(n_images, n_weights);
end

%% Process the images

if has_spectral
    patch_penalties_spectral_L1 = zeros(n_patches, n_weights, n_images);
    patch_penalties_spectral_L2 = zeros(n_patches, n_weights, n_images);
end
patch_penalties_rgb_L1 = zeros(n_patches, n_weights, n_images);
patch_penalties_rgb_L2 = zeros(n_patches, n_weights, n_images);

all_mdc_weights = zeros(n_patches, n_weights, n_images, n_admm_algorithms);
all_mse_weights = zeros(n_patches, n_weights, n_images, n_admm_algorithms);

for i = 1:n_images
    if verbose
        fprintf('[SelectWeightsForDataset, image %d] Starting\n', i);
    end

    % Generate or load input images, and instantiate dispersion information
    run('LoadAndConvertImage.m');
    
    % Select image patches
    if any(full_patch_size > image_sampling)
        error('Image %d is smaller than the patch size', i);
    end
    
    corners = [
        randi(image_sampling(1) - full_patch_size(1), n_patches, 1),...
        randi(image_sampling(2) - full_patch_size(2), n_patches, 1)
    ];
    % Avoid changing the Bayer pattern
    corners(mod(corners, 2) == 0) = corners(mod(corners, 2) == 0) - 1;

    % Characterize the patches
    if has_spectral
        for pc = 1:n_patches
            patch = I_spectral_gt(...
                corners(pc, 1):(corners(pc, 1) + (full_patch_size(1) - 1)),...
                corners(pc, 2):(corners(pc, 2) + (full_patch_size(2) - 1)), : ...
            );
            for w = 1:n_weights
                if w > n_spectral_weights
                    continue;
                end
                err_vector = patch_operators_spectral{w} * reshape(patch, [], 1);
                patch_penalties_spectral_L1(pc, w, i) = mean(abs(err_vector));
                patch_penalties_spectral_L2(pc, w, i) = dot(err_vector, err_vector) / length(err_vector);
            end
        end
    end
    for pc = 1:n_patches
        patch = I_rgb_gt(...
            corners(pc, 1):(corners(pc, 1) + (full_patch_size(1) - 1)),...
            corners(pc, 2):(corners(pc, 2) + (full_patch_size(2) - 1)), : ...
        );
        for w = 1:n_weights
            err_vector = patch_operators_rgb{w} * reshape(patch, [], 1);
            patch_penalties_rgb_L1(pc, w, i) = mean(abs(err_vector));
            patch_penalties_rgb_L2(pc, w, i) = dot(err_vector, err_vector) / length(err_vector);
            if has_spectral && w == n_weights
                patch_penalties_spectral_L1(pc, w, i) = patch_penalties_rgb_L1(pc, w, i);
                patch_penalties_spectral_L2(pc, w, i) = patch_penalties_rgb_L2(pc, w, i);
            end
        end
    end

    % Run the algorithms   
    for f = 1:n_admm_algorithms
        algorithm = admm_algorithms.(admm_algorithm_fields{f});
        if ~algorithm.enabled
            continue;
        end

        solvePatchesADMMOptions.admm_options.int_method = 'none';
        if algorithm.spectral && has_color_map
            if ~channel_mode
                solvePatchesADMMOptions.admm_options.int_method = int_method;
            end
        elseif algorithm.spectral
            continue;
        end
        
        reg_options_f = mergeStructs(...
            solvePatchesADMMOptions.reg_options, solvePatchesADMMOptions.reg_options, true, false...
        );
        reg_options_f.enabled = algorithm.priors;

        admm_options_f = mergeStructs(...
            solvePatchesADMMOptions.admm_options, algorithm.options, false, true...
        );
      
        mdc_weights_patches = zeros(n_patches, n_weights);
        mse_weights_patches = zeros(n_patches, n_weights);
        if algorithm.spectral
            for pc = 1:n_patches
                solvePatchesADMMOptions.patch_options.target_patch = corners(pc, :);
                [...
                    ~, ~, weights_images...
                ] = solvePatchesADMM(...
                  [], I_raw_gt, bayer_pattern, df_spectral_reverse,...
                  sensor_map_resampled, bands,...
                  admm_options_f, reg_options_f,...
                  solvePatchesADMMOptions.patch_options,...
                  solvePatchesADMMVerbose...
                );
                mdc_weights_patches(pc, reg_options_f.enabled) = weights_images(1, 1, :);
                
                [...
                    ~, ~, weights_images...
                ] = solvePatchesADMM(...
                  I_spectral_gt, I_raw_gt, bayer_pattern, df_spectral_reverse,...
                  sensor_map_resampled, bands,...
                  admm_options_f, reg_options_f,...
                  solvePatchesADMMOptions.patch_options,...
                  solvePatchesADMMVerbose...
                );
                mse_weights_patches(pc, reg_options_f.enabled) = weights_images(1, 1, :);                   
            end
        else
            for pc = 1:n_patches
                solvePatchesADMMOptions.patch_options.target_patch = corners(pc, :);
                [...
                    ~, ~, weights_images...
                ] = solvePatchesADMM(...
                  [], I_raw_gt, bayer_pattern, df_rgb_reverse,...
                  sensor_map_rgb, bands_rgb,...
                  admm_options_f, reg_options_f,...
                  solvePatchesADMMOptions.patch_options,...
                  solvePatchesADMMVerbose...
                );
                mdc_weights_patches(pc, reg_options_f.enabled) = weights_images(1, 1, :);
                
                [...
                    ~, ~, weights_images...
                ] = solvePatchesADMM(...
                  I_rgb_gt, I_raw_gt, bayer_pattern, df_rgb_reverse,...
                  sensor_map_rgb, bands_rgb,...
                  admm_options_f, reg_options_f,...
                  solvePatchesADMMOptions.patch_options,...
                  solvePatchesADMMVerbose...
                );
                mse_weights_patches(pc, reg_options_f.enabled) = weights_images(1, 1, :);
            end
        end
        mdc_weights = geomean(mdc_weights_patches, 1);
        mse_weights = geomean(mse_weights_patches, 1);
        all_mdc_weights(:, :, i, f) = mdc_weights_patches;
        all_mse_weights(:, :, i, f) = mse_weights_patches;
        admm_algorithms.(admm_algorithm_fields{f}).mdc_weights(i, :) = mdc_weights;
        admm_algorithms.(admm_algorithm_fields{f}).mse_weights(i, :) = mse_weights;
    end

    if verbose
        fprintf('[SelectWeightsForDataset, image %d] Finished\n', i);
    end
end

%% Visualization of weights selected for the dataset

mdc_color = [1, 0, 0];
mse_color = [0, 1, 0];

min_nz_weight = min(min(all_mdc_weights(all_mdc_weights ~= 0)), min(all_mse_weights(all_mse_weights ~= 0)));
max_nz_weight = max(max(all_mdc_weights(all_mdc_weights ~= 0)), max(all_mse_weights(all_mse_weights ~= 0)));
log_min_nz_weight = log10(min_nz_weight);
log_max_nz_weight = log10(max_nz_weight);
plot_limits = [log_min_nz_weight - 1, log_max_nz_weight + 1];

if has_spectral
    log_patch_penalties_spectral_L1 = log10(patch_penalties_spectral_L1);
    log_patch_penalties_spectral_L2 = log10(patch_penalties_spectral_L2);
end
log_patch_penalties_rgb_L1 = log10(patch_penalties_rgb_L1);
log_patch_penalties_rgb_L2 = log10(patch_penalties_rgb_L2);

for f = 1:n_admm_algorithms
    algorithm = admm_algorithms.(admm_algorithm_fields{f});
    if ~algorithm.enabled || (algorithm.spectral && ~has_color_map)
        continue;
    end
    
    name_params = sprintf(...
        '%s_patch%dx%d_pad%d_',...
        algorithm.file, patch_size(1), patch_size(2), padding...
    );
    if algorithm.spectral
        name_params = [...
            sprintf('bands%d_', n_bands), name_params...
        ];
    else
        name_params = [...
            'RGB_', name_params...
        ];
    end
    name_params = fullfile(output_directory, name_params);
    
    admm_options_f = mergeStructs(...
        solvePatchesADMMOptions.admm_options, algorithm.options, false, true...
    );

    enabled_weights = algorithm.priors;
    n_active_weights = sum(enabled_weights);
    to_all_weights = find(enabled_weights);

    mdc_weights = reshape(permute(all_mdc_weights(:, :, :, f), [1, 3, 2, 4]), [], n_weights);
    mdc_weights = mdc_weights(:, enabled_weights);
    log_mdc_weights = log10(mdc_weights);
    mse_weights = reshape(permute(all_mse_weights(:, :, :, f), [1, 3, 2, 4]), [], n_weights);
    mse_weights = mse_weights(:, enabled_weights);
    log_mse_weights = log10(mse_weights);

    fg = figure;
    hold on
    if n_active_weights == 1
        scatter(...
            log_mse_weights, log_mdc_weights, 'filled'...
        );
        line_limits = [...
            min(min(log_mse_weights, log_mdc_weights));
            max(max(log_mse_weights, log_mdc_weights))
        ];
        line(line_limits, line_limits, 'Color', 'b');
        legend('Weights', 'y = x');
        xlabel('Weight selected using the mean square error');
        ylabel('Weight selected using the minimum distance criterion');
        xlim(plot_limits)
        ylim(plot_limits)
    elseif n_active_weights == 2
        scatter(...
            log_mdc_weights(:, 1), log_mdc_weights(:, 2), [], mdc_color, 'filled'...
        );
        scatter(...
            log_mse_weights(:, 1), log_mse_weights(:, 2), [], mse_color, 'filled'...
        );
        legend('MDC', 'MSE');
        xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
        ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
        xlim(plot_limits)
        ylim(plot_limits)
    elseif n_active_weights == 3
        scatter3(...
            log_mdc_weights(:, 1), log_mdc_weights(:, 2), log_mdc_weights(:, 3),...
            [], mdc_color, 'filled'...
        );
        scatter3(...
            log_mse_weights(:, 1), log_mse_weights(:, 2), log_mse_weights(:, 3),...
            [], mse_color, 'filled'...
        );
        legend('MDC', 'MSE');
        xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
        ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
        zlabel(sprintf('log_{10}(weight %d)', to_all_weights(3)))
        xlim(plot_limits)
        ylim(plot_limits)
        zlim(plot_limits)
    else
        error('Unexpected number of active weights.');
    end
    title(sprintf(...
        'Agreement between MDC and MSE weights selected for %s',...
        algorithm.name...
    ));
    hold off
    
    savefig(...
        fg,...
        [name_params 'weightsCorrelation.fig'], 'compact'...
    );
    close(fg);

    for w = 1:n_active_weights
        aw = to_all_weights(w);
        if algorithm.spectral && has_color_map
            if admm_options_f.norms(aw)
                patch_penalties = log_patch_penalties_spectral_L1(:, w, :);
            else
                patch_penalties = log_patch_penalties_spectral_L2(:, w, :);
            end
        else
            if admm_options_f.norms(aw)
                patch_penalties = log_patch_penalties_rgb_L1(:, w, :);
            else
                patch_penalties = log_patch_penalties_rgb_L2(:, w, :);
            end
        end
        patch_penalties = reshape(patch_penalties, [], 1);
        
        fg = figure;
        hold on
        scatter(patch_penalties, log_mdc_weights(:, w), [], mdc_color, 'filled');
        scatter(patch_penalties, log_mse_weights(:, w), [], mse_color, 'filled');
        legend('MDC', 'MSE');
        xlabel(sprintf('log_{10}(Penalty %d)', aw))
        ylabel(sprintf('log_{10}(Weight %d)', aw))
        title(sprintf(...
            'Agreement between MDC and MSE weights selected for %s',...
            algorithm.name...
        ));
        ylim(plot_limits)
        hold off
        
        savefig(...
            fg,...
            [name_params sprintf('weight%d.fig', aw)], 'compact'...
        );
        close(fg);
    end
end

%% Save parameters and data to a file
save_variables_list = [ parameters_list, {...
    'bands', 'admm_algorithms', 'corners',...
    'patch_penalties_rgb_L1', 'patch_penalties_rgb_L2',...
    'all_mdc_weights', 'all_mse_weights'...
} ];
if has_spectral
    save_variables_list = [save_variables_list, {...
        'bands_spectral', 'patch_penalties_spectral_L1', 'patch_penalties_spectral_L2'...
    }];
end
if has_color_map
    save_variables_list = [save_variables_list, {'sensor_map_resampled'}];
    if has_spectral
        save_variables_list = [save_variables_list, {'sensor_map_spectral', }];
    end
end
save_data_filename = fullfile(output_directory, ['SelectWeightsForDataset_' dataset_name '.mat']);
save(save_data_filename, save_variables_list{:});