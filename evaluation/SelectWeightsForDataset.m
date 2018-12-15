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
% The documentation in the script 'CorrectByHyperspectralADMM.m' contains
% more information on the formats of the various types of data associated
% with the datasets.
%
% This script also runs 'SetFixedParameters.m' to set the values of
% seldomly-changed parameters. These parameters are briefly documented in
% 'SetFixedParameters.m'.
%
% ## Output
%
% ### Data and parameters
% A '.mat' file containing the following variables, as appropriate:
% - 'bands': A vector containing the wavelengths of the spectral
%   bands used in hyperspectral image estimation.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference.
% - 'bands_spectral': A vector containing the wavelengths of the spectral
%   bands associated with ground truth hyperspectral images.
% - 'color_weights': A matrix for converting pixels in the estimated
%   hyperspectral images to colour, as determined by the 'sensor_map'
%   variable loaded from the colour space conversion data file, and by the
%   type of numerical intergration to perform.
% - 'color_weights_reference': A matrix for converting pixels in the ground
%   truth hyperspectral images to colour, as determined by the 'sensor_map'
%   variable loaded from the colour space conversion data file, and by the
%   type of numerical intergration to perform.
% - 'spectral_weights': A matrix for converting pixels in the spectral
%   space of the estimated hyperspectral images to the spectral space of
%   the true hyperspectral images.
% - 'admm_algorithms': A structure describing the algorithms for which
%   regularization weights were selected. 'admm_algorithms' is created by
%   'SetAlgorithms.m', and then each algorithm is given additional fields
%   by this script. Of the following fields, only those corresponding to
%   `true` values in the `criteria` vector set in 'SetFixedParameters.m'
%   will be added.
%   - 'mdc_weights': Regularization weights selected using the minimum
%     distance criterion of Song et al. 2016. 'mdc_weights' is a 3D array
%     with one frame per image in the dataset. The first dimension indexes
%     spectral resolutions (applicable only for 'solvePatchesMultiADMM()'),
%     whereas the second dimension indexes regularization weights.
%   - 'mse_weights': Regularization weights selected to minimize the mean
%     square error with respect to the true images. 'mse_weights' is an
%     array with the same format as 'mdc_weights'.
%   - 'dm_weights': Regularization weights selected to minimize the mean
%     square error with respect to a demosaicing result. 'dm_weights' is an
%     array with the same format as 'mdc_weights'.
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
% - 'all_weights': Regularization weights selected using the minimum
%   distance criterion, the mean squared error with respect to the true
%   image, and using the mean squared error with respect to a demosaicking
%   result. `all_weights{f}(pc, :, t, i, cr)` is the vector of
%   regularization weights selected for the f-th algorithm on the pc-th
%   patch of the i-th image. `t` is the spectral resolution, which is
%   relevant only for 'solvePatchesMultiADMM()', but is one for single
%   spectral resolution approaches (i.e. 'solvePatchesADMM()'). `cr` is the
%   criterion (from the ordered list given immediately above) used to
%   select regularization weights. If `cr` corresponds to a disabled
%   criterion, the corresponding regularization weights are all zeros.
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
%   algorithms in 'all_weights', and the 'mdc_weights', etc., fields will
%   not be added to their structures in 'admm_algorithms'.
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

dataset_name = '20181212_RealData_spectralAsRAW';

% Describe algorithms to run
run('SetAlgorithms.m')

% Number of patches to select for each image
n_patches = 10;

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads/results';

% Produce console output to describe the processing in this script
verbose = true;

% ## Parameters which do not usually need to be changed
% Note that this sets the value of `criteria`.
run('SetFixedParameters.m')

%% Check for problematic parameters

if use_fixed_weights
    error('Weights should be fixed by running ''SelectWeightsForDataset.m'', not using the `use_fixed_weights` parameter in ''SetFixedParameters.m''');
end

if sum(criteria) == 0
    error('All regularization weight selection criteria are disabled.');
end

%% Preprocess the dataset

dp = describeDataset(dataset_name);

run('PreprocessDataset.m')

%% Prepare for patch processing

n_weights = length(solvePatchesADMMOptions.reg_options.enabled);

patch_size = solvePatchesADMMOptions.patch_options.patch_size;
padding = solvePatchesADMMOptions.patch_options.padding;
full_patch_size = patch_size + padding * 2;

solvePatchesMultiADMMOptions.sampling_options.show_steps = true;

if has_spectral
    image_sampling_patch_spectral = [full_patch_size, length(bands_spectral)];
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
image_sampling_patch_rgb = [full_patch_size, n_channels_rgb];
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
        G = antiMosaicMatrix(full_patch_size, bayer_pattern);
    end
    patch_operators_rgb{w} = G;
end

admm_algorithm_fields = fieldnames(admm_algorithms);
n_admm_algorithms = length(admm_algorithm_fields);

%% Process the images

if has_spectral
    patch_penalties_spectral_L1 = zeros(n_patches, n_weights, n_images);
    patch_penalties_spectral_L2 = zeros(n_patches, n_weights, n_images);
end
patch_penalties_rgb_L1 = zeros(n_patches, n_weights, n_images);
patch_penalties_rgb_L2 = zeros(n_patches, n_weights, n_images);

all_weights = cell(n_admm_algorithms, 1);
n_bands_all = cell(n_admm_algorithms, 1);
n_criteria = length(criteria);

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
        if ~algorithm.enabled || (algorithm.spectral && ~has_color_map)
            continue;
        end
        for cr = 1:n_criteria
            if ~criteria(cr)
                continue;
            end
        
            true_spectral = algorithm.spectral && ~channel_mode;

            if true_spectral
                reg_options_f = solvePatchesMultiADMMOptions.reg_options;
            else
                reg_options_f = solvePatchesADMMOptions.reg_options;
            end
            reg_options_f.enabled = algorithm.priors;
            enabled_weights = reg_options_f.enabled;
            n_active_weights = sum(enabled_weights);

            if true_spectral
                admm_options_f = mergeStructs(...
                    solvePatchesMultiADMMOptions.admm_options, algorithm.options, false, true...
                );
            else
                admm_options_f = mergeStructs(...
                    solvePatchesADMMOptions.admm_options, algorithm.options, false, true...
                );
            end
            
            if cr == dm_index
                reg_options_f.demosaic = true;
            else
                reg_options_f.demosaic = false;
            end
            
            weights_patches = zeros(n_patches, n_weights);
            n_steps = 1;
            have_steps = false;
            if algorithm.spectral
                if true_spectral
                    for pc = 1:n_patches
                        solvePatchesMultiADMMOptions.patch_options.target_patch = corners(pc, :);
                        if cr == mse_index
                            I_in.I = I_spectral_gt;
                            I_in.spectral_bands = bands_spectral;
                            [...
                                bands_all, ~, ~, weights_images...
                            ] = solvePatchesMultiADMM(...
                              I_in, I_raw_gt, bayer_pattern, df_spectral_reverse,...
                              sensor_map, bands_color,...
                              solvePatchesMultiADMMOptions.sampling_options,...
                              admm_options_f, reg_options_f,...
                              solvePatchesMultiADMMOptions.patch_options,...
                              solvePatchesMultiADMMVerbose...
                            );
                        else
                            [...
                                bands_all, ~, ~, weights_images...
                            ] = solvePatchesMultiADMM(...
                              [], I_raw_gt, bayer_pattern, df_spectral_reverse,...
                              sensor_map, bands_color,...
                              solvePatchesMultiADMMOptions.sampling_options,...
                              admm_options_f, reg_options_f,...
                              solvePatchesMultiADMMOptions.patch_options,...
                              solvePatchesMultiADMMVerbose...
                            );
                        end
                        
                        if ~have_steps
                            n_steps = size(weights_images, 3) / n_active_weights;
                            weights_patches = repmat(weights_patches, 1, 1, n_steps);
                            if n_steps > 1
                                n_bands_all{f} = zeros(length(bands_all), 1);
                                for b = 1:length(bands_all)
                                    n_bands_all{f}(b) = length(bands_all{b});
                                end
                            end
                            have_steps = true;
                        end
                        weights_patches(pc, reg_options_f.enabled, :) = weights_images(1, 1, :);
                    end
                else
                    for pc = 1:n_patches
                        solvePatchesADMMOptions.patch_options.target_patch = corners(pc, :);
                        if cr == mse_index
                            I_in.I = I_spectral_gt;
                            I_in.spectral_weights = spectral_weights;
                            [...
                                ~, ~, weights_images...
                            ] = solvePatchesADMM(...
                              I_in, I_raw_gt, bayer_pattern, df_spectral_reverse,...
                              color_weights, bands,...
                              admm_options_f, reg_options_f,...
                              solvePatchesADMMOptions.patch_options,...
                              solvePatchesADMMVerbose...
                            );
                        else
                            [...
                                ~, ~, weights_images...
                            ] = solvePatchesADMM(...
                              [], I_raw_gt, bayer_pattern, df_spectral_reverse,...
                              color_weights, bands,...
                              admm_options_f, reg_options_f,...
                              solvePatchesADMMOptions.patch_options,...
                              solvePatchesADMMVerbose...
                            );
                        end
                        weights_patches(pc, reg_options_f.enabled) = weights_images(1, 1, :);
                    end
                end
            else
                for pc = 1:n_patches
                    solvePatchesADMMOptions.patch_options.target_patch = corners(pc, :);
                    if cr == mse_index
                        I_in.I = I_rgb_gt;
                        I_in.spectral_weights = sensor_map_rgb;
                        [...
                            ~, ~, weights_images...
                        ] = solvePatchesADMM(...
                          I_in, I_raw_gt, bayer_pattern, df_rgb_reverse,...
                          sensor_map_rgb, bands_rgb,...
                          admm_options_f, reg_options_f,...
                          solvePatchesADMMOptions.patch_options,...
                          solvePatchesADMMVerbose...
                        );
                    else
                        [...
                            ~, ~, weights_images...
                        ] = solvePatchesADMM(...
                          [], I_raw_gt, bayer_pattern, df_rgb_reverse,...
                          sensor_map_rgb, bands_rgb,...
                          admm_options_f, reg_options_f,...
                          solvePatchesADMMOptions.patch_options,...
                          solvePatchesADMMVerbose...
                        );
                    end
                    weights_patches(pc, reg_options_f.enabled) = weights_images(1, 1, :);
                end
            end
            field_weights = geomean(weights_patches, 1);
            if true_spectral && n_steps > 1
                field_weights = squeeze(field_weights).';
            end
            if i == 1 && cr == 1
                all_weights{f} = zeros(n_patches, n_weights, n_steps, n_images, n_criteria);
                admm_algorithms.(admm_algorithm_fields{f}).(criteria_fields{cr}) = zeros(n_steps, n_weights, n_images); 
            end
            all_weights{f}(:, :, :, i, cr) = weights_patches;
            admm_algorithms.(admm_algorithm_fields{f}).(criteria_fields{cr})(:, :, i) = field_weights;
        end
    end

    if verbose
        fprintf('[SelectWeightsForDataset, image %d] Finished\n', i);
    end
end

%% Visualization of weights selected for the dataset

if criteria(2)
    reference_criteria = 2;
elseif criteria(3)
    reference_criteria = 3;
elseif criteria(1)
    reference_criteria = 1;
end

min_nz_weight = Inf;
max_nz_weight = -Inf;
for f = 1:n_admm_algorithms
    if admm_algorithms.(admm_algorithm_fields{f}).enabled
        for cr = 1:n_criteria
            if criteria(cr)
                all_weights_fcr = all_weights{f}(:, :, :, :, cr);
                min_nz_weight = min(...
                    min_nz_weight, min(all_weights_fcr(all_weights_fcr ~= 0))...
                );
                max_nz_weight = max(...
                    max_nz_weight, max(all_weights_fcr(all_weights_fcr ~= 0))...
                );
            end
        end
    end
end
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
    
    true_spectral = algorithm.spectral && ~channel_mode;
    
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
    
    if true_spectral
        admm_options_f = mergeStructs(...
            solvePatchesADMMOptions.admm_options, algorithm.options, false, true...
            );
    else
        admm_options_f = mergeStructs(...
            solvePatchesMultiADMMOptions.admm_options, algorithm.options, false, true...
            );
    end
    
    enabled_weights = algorithm.priors;
    n_active_weights = sum(enabled_weights);
    to_all_weights = find(enabled_weights);
    
    n_steps = size(all_weights{f}, 3);            
    for t = 1:n_steps
        if true_spectral
            name_params_t = [...
                name_params, sprintf('step%d_', t), ...
                ];
        else
            name_params_t = name_params;
        end
        name_params_t = fullfile(output_directory, name_params_t);

        ref_weights = reshape(permute(all_weights{f}(:, :, t, :, reference_criteria), [1, 4, 2, 3, 5]), [], n_weights);
        ref_weights = ref_weights(:, enabled_weights);
        log_ref_weights = log10(ref_weights);
        
        fg = figure;
        hold on
        if n_active_weights == 1
            line(plot_limits, plot_limits, 'Color', 'b');
        end
        for cr = 1:n_criteria
            if ~criteria(cr)
                continue;
            end
            field_weights = reshape(permute(all_weights{f}(:, :, t, :, cr), [1, 4, 2, 3, 5]), [], n_weights);
            field_weights = field_weights(:, enabled_weights);
            log_weights = log10(field_weights);

            if n_active_weights == 1
                scatter(...
                    log_ref_weights, log_weights, [], criteria_colors(cr, :), 'filled'...
                );
                xlabel('Weight selected using the reference method');
                ylabel('Weight selected using the other method');
                xlim(plot_limits)
                ylim(plot_limits)
            elseif n_active_weights == 2
                scatter(...
                    log_weights(:, 1), log_weights(:, 2), [], criteria_colors(cr, :), 'filled'...
                );
                xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
                ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
                xlim(plot_limits)
                ylim(plot_limits)
            elseif n_active_weights == 3
                scatter3(...
                    log_weights(:, 1), log_weights(:, 2), log_weights(:, 3),...
                    [], criteria_colors(cr, :), 'filled'...
                );
                xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
                ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
                zlabel(sprintf('log_{10}(weight %d)', to_all_weights(3)))
                xlim(plot_limits)
                ylim(plot_limits)
                zlim(plot_limits)
            else
                error('Unexpected number of active weights.');
            end
        end
        hold off
        if true_spectral
            title(sprintf(...
                'Agreement between weights selected for %s, step %d',...
                algorithm.name, t...
            ));
        else
            title(sprintf(...
                'Agreement between weights selected for %s',...
                algorithm.name...
            ));
        end
        savefig(...
            fg,...
            [name_params_t 'weightsCorrelation.fig'], 'compact'...
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
            for cr = 1:n_criteria
                if ~criteria(cr)
                    continue;
                end
                field_weights = reshape(permute(all_weights{f}(:, :, t, :, cr), [1, 4, 2, 3, 5]), [], n_weights);
                field_weights = field_weights(:, enabled_weights);
                log_weights = log10(field_weights);
                scatter(patch_penalties, log_weights(:, w), [], criteria_colors(cr, :), 'filled');
            end
            hold off
            if true_spectral
                title(sprintf(...
                    'Agreement between weights selected for %s, step %d',...
                    algorithm.name, t...
                ));
            else
                title(sprintf(...
                    'Agreement between weights selected for %s',...
                    algorithm.name...
                ));
            end
            ylim(plot_limits)
            xlabel(sprintf('log_{10}(Penalty %d)', aw))
            ylabel(sprintf('log_{10}(Weight %d)', aw))
            
            savefig(...
                fg,...
                [name_params_t sprintf('weight%d.fig', aw)], 'compact'...
                );
            close(fg);
        end
    end
    
    % Plot weights vs. image estimation step
    if true_spectral && n_steps > 1
        n_bands_plot = repelem(n_bands_all{f}, n_patches * n_images);
        for w = 1:n_active_weights
            aw = to_all_weights(w);

            fg = figure;
            hold on
            for cr = 1:n_criteria
                if ~criteria(cr)
                    continue;
                end
                field_weights = reshape(permute(all_weights{f}(:, :, :, :, cr), [1, 4, 3, 2, 5]), [], n_weights);
                field_weights = field_weights(:, enabled_weights);
                log_weights = log10(field_weights);
            
                scatter(n_bands_plot, log_weights(:, w), [], criteria_colors(cr, :), 'filled');
            end
            hold off
            xlabel('Number of bands in image estimation step')
            ylabel(sprintf('log_{10}(Weight %d)', aw))
            title(sprintf(...
                'Agreement between weight %d selected for %s',...
                aw, algorithm.name...
            ));
            ylim(plot_limits)
            legend(criteria_abbrev{criteria});
            
            savefig(...
                fg,...
                fullfile(output_directory, [name_params, sprintf('weight%d.fig', aw)]), 'compact'...
                );
            close(fg);
        end
    end
end

%% Save parameters and data to a file
save_variables_list = [ parameters_list, {...
    'admm_algorithms', 'corners',...
    'patch_penalties_rgb_L1', 'patch_penalties_rgb_L2',...
    'all_weights'...
} ];
if has_spectral
    save_variables_list = [save_variables_list, {...
        'bands_spectral', 'spectral_weights', 'color_weights_reference',...
        'patch_penalties_spectral_L1', 'patch_penalties_spectral_L2'...
    }];
end
if has_color_map
    save_variables_list = [save_variables_list, {'bands', 'bands_color', 'color_weights'}];
end
save_data_filename = fullfile(output_directory, ['SelectWeightsForDataset_' dataset_name '.mat']);
save(save_data_filename, save_variables_list{:});