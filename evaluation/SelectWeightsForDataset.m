%% Select regularization weights for images in a dataset
% Choose regularization weights for image patches to determine the best weights
% for each algorithm, and how much the optimal weights vary within and between
% images.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% The dataset determines the data to be loaded, and algorithms to be tested, as
% encapsulated by the 'describeDataset()' function. It also optionally provides
% a list of image patches for each image to use for selecting regularization
% weights. If no patches are specified for an image, this script will select
% `n_patches` random patches.
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
%     spectral resolutions (applicable only for 'solvePatchesSpectral()'),
%     whereas the second dimension indexes regularization weights.
%   - 'mse_weights': Regularization weights selected to minimize the mean
%     square error with respect to the true images. 'mse_weights' is an
%     array with the same format as 'mdc_weights'.
%   - 'dm_weights': Regularization weights selected to minimize the mean
%     square error with respect to a demosaicing result. 'dm_weights' is an
%     array with the same format as 'mdc_weights'.
% - 'corners': A cell vector containing the top-left corners of the image
%   patches used to select regularization weights. `corners{i}` is a two-column
%   matrix, where the columns contain the row and column indices, respectively,
%   of the top left corners of patches used for the i-th image. Note that
%   'corners' does not account for the padding around each patch that is used to
%   reduce patch boundary effects.
% - 'patch_penalties_spectral_L1': A cell vector containing the L1 norms of the
%   regularization penalties evaluated on the true spectral image patches.
%   `patch_penalties_spectral_L1{i}` is a matrix where the columns correspond to
%   enabled regularization terms, and the rows correspond to image patches, for
%   the i-th image.
% - 'patch_penalties_spectral_L2': A cell vector similar to
%   'patch_penalties_spectral_L1', but which contains L2 norms of
%   regularization penalties.
% - 'patch_penalties_rgb_L1': A cell vector similar to
%   'patch_penalties_spectral_L1', but which corresponds to regularization
%   penalties evaluated on colour versions of the images.
% - 'patch_penalties_rgb_L2': A cell vector similar to 'patch_penalties_rgb_L1',
%   but which contains L2 norms of regularization penalties.
% - 'all_weights': Regularization weights selected using the minimum
%   distance criterion, the mean squared error with respect to the true image,
%   and using the mean squared error with respect to a demosaicking result.
%   `all_weights{f}{i, cr}` is a 3D array of regularization weights selected for
%   the f-th algorithm on the i-th image. `cr` is the criterion (from the
%   ordered list given immediately above) used to select regularization weights.
%   If `cr` corresponds to a disabled criterion, the corresponding matrix is
%   empty. Rows of `all_weights{f}{i, cr}` correspond to patches, columns
%   correspond to regularization terms, and the third dimension indexes the
%   spectral resolution, which is relevant only for 'solvePatchesSpectral()',
%   and has size one otherwise.
% - 'time_admm': Execution timing information, stored as a 3D array.
%   `time_admm(f, i, cr)` is the average time taken (in seconds) to process a
%   patch of the i-th image with the f-th ADMM-family algorithm defined in
%   'SetAlgorithms.m', according to the cr-th regularization weight selection
%   criterion. Entries corresponding to disabled algorithms or disabled weight
%   selection criterion will be set to `NaN`.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`, which is initialized in this file, and then augmented
% by 'SetFixedParameters.m'.
%
% The file is saved as 'SelectWeightsForDataset_${dataset_name}.mat'.
%
% ### Graphical output
%
% Figures are generated for each ADMM algorithm showing the distributions of
% weights selected by the different criteria, and are saved to '.fig' files.
%
% A figure is also generated for each image showing the locations of the patches
% used to select weights.
%
% ## Notes
% - This script uses 'patch_size' and 'padding' defined in the dataset
%   description, not those set in 'SetFixedParameters.m'.
% - Regularization weights will not be selected for algorithms which are
%   disabled (using the 'enabled' fields in the structures describing the
%   algorithms). Instead, values of zero will be output for disabled
%   algorithms in 'all_weights', and the 'mdc_weights', etc., fields will
%   not be added to their structures in 'admm_algorithms'.
% - Only image patches that are not clipped by the image edges will be
%   chosen for selecting regularization weights. An warning will be thrown if
%   there are patches defined by describeDataset() which are clipped by the
%   image edges. Note that the clipping test comes after cropping the image to
%   the domain of the model of dispersion to be used for image estimation.
% - Image patch spectral derivatives, used for graphical output, not for
%   regularization weights selection, are computed according to the
%   'full_GLambda' field of the 'solvePatchesColorOptions.admm_options'
%   structure defined in 'SetFixedParameters.m', rather than according to
%   per-algorithm options.
% - If the true images are affected by dispersion, but image estimation
%   will involve dispersion correction, then the 'MSE' regularization
%   weight selection criterion defined in 'SetFixedParameters.m' is
%   evaluated against versions of the true images which have been
%   approximately corrected by dispersion using image warping.
%   Consequently, the images used to select regularization weights are no
%   longer ground truth images, but it is still better than evaluating the
%   criterion on images which are incomparable because of differences in
%   dispersion.
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

dataset_name = '20190421_ComputarLens_dHyper_dispersion';

% Describe algorithms to run
run('SetAlgorithms.m')

% Default number of patches to select for each image, when none are provided
n_patches = 10;

% Output directory for all images and saved parameters
output_directory = 'C:\Users\GraphicsLab\Documents\llanos\Results\weights_selection_dHyper_spectral_MSEBugFixed';

% Produce console output to describe the processing in this script
verbose = true;

% ## Parameters which do not usually need to be changed
% Note that this sets the value of `criteria`.
run('SetFixedParameters.m')

%% Check for problematic parameters

if use_fixed_weights
    error('The `use_fixed_weights` parameter in ''SetFixedParameters.m'' should be `false` when running this script.');
end

if sum(criteria) == 0
    error('All regularization weight selection criteria are disabled.');
end

%% Preprocess the dataset

dp = describeDataset(dataset_name);

run('PreprocessDataset.m')

%% Prepare for patch processing

n_weights = length(solvePatchesColorOptions.reg_options.enabled);

patch_size = dp.patch_size;
padding = dp.padding;
full_patch_size = patch_size + padding * 2;

solvePatchesSpectralOptions.sampling_options.show_steps = true;

if has_spectral
    image_sampling_patch_spectral = [full_patch_size, length(bands_spectral)];
    patch_operators_spectral = cell(n_weights, 1);
    for w = 1:n_weights
        if w == 1 || w == 2
            G = spatialGradient(image_sampling_patch_spectral);
        end
        if w == 2
            G_lambda = spectralGradient(image_sampling_patch_spectral, solvePatchesColorOptions.admm_options.full_GLambda);
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
        G_lambda = spectralGradient(image_sampling_patch_rgb, solvePatchesColorOptions.admm_options.full_GLambda);
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
        G = spatialLaplacian(image_sampling_patch_rgb);
    end
    patch_operators_rgb{w} = G;
end

admm_algorithm_fields = fieldnames(admm_algorithms);
n_admm_algorithms = length(admm_algorithm_fields);

%% Process the images

if has_spectral
    patch_penalties_spectral_L1 = cell(n_images, 1);
    patch_penalties_spectral_L2 = cell(n_images, 1);
end
patch_penalties_rgb_L1 = cell(n_images, 1);
patch_penalties_rgb_L2 = cell(n_images, 1);

all_weights = cell(n_admm_algorithms, 1);
n_bands_all = cell(n_admm_algorithms, 1);
n_criteria = length(criteria);

time_admm = nan(n_admm_algorithms, n_images, n_criteria);

corners = cell(n_images, 1);
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
    
    corners_i = [];
    if isfield(dp, 'params_patches') && isfield(dp.params_patches, names{i})
        corners_i = fliplr(dp.params_patches.(names{i})) -...
            repmat(ceil(patch_size / 2), size(dp.params_patches.(names{i}), 1), 1);
        filter = ((corners_i(:, 1) - padding) >= 1) & ((corners_i(:, 2) - padding) >= 1) &...
            ((corners_i(:, 1) + (patch_size(1) + padding - 1)) <= image_sampling(1)) & ((corners_i(:, 2) + (patch_size(2) + padding - 1)) <= image_sampling(2));
        if ~any(filter)
            warning([
                'No patches for image "%s" are fully within the image borders (accounting for cropping to any dispersion model).\n',...
                'Using random patches instead'], names{i});
        elseif ~all(filter)
            warning('Some patches for image "%s" are not fully within the image borders (accounting for cropping to any dispersion model).', names{i});
        end
        corners_i = corners_i(filter, :);
    end
    if isempty(corners_i)
        corners_i = [
            randi(image_sampling(1) - full_patch_size(1), n_patches, 1),...
            randi(image_sampling(2) - full_patch_size(2), n_patches, 1)
        ] + padding;
    end
    n_patches_i = size(corners_i, 1);
    if n_patches_i == 0
        error('No patches for image "%s".', names{i});
    end
    % Avoid changing the Bayer pattern
    corners_i(mod(corners_i, 2) == 0) = corners_i(mod(corners_i, 2) == 0) - 1;
    corners{i} = corners_i;
    
    % Generate a figure showing all patches
    fg = figure;
    labels = cell(n_patches_i, 1);
    font_size = max(12, floor(0.02 * max(image_sampling)));
    for pc = 1:n_patches_i
        labels{pc} = sprintf('%d', pc);
    end
    figure_image = insertObjectAnnotation(...
        I_raw_gt, 'rectangle', [fliplr(corners_i), fliplr(repmat(patch_size, n_patches_i, 1))],...
        labels,...
        'TextBoxOpacity', 0.9, 'FontSize', font_size, 'LineWidth', 2,...
        'Color', jet(n_patches_i)...
        );
    imshow(figure_image);
    title(sprintf('Weight selection patches for image "%s"', names{i}));
    figure_save_name = sprintf('%s_patches.fig', names{i});
    savefig(...
        fg,...
        fullfile(output_directory, figure_save_name), 'compact'...
        );
    close(fg);

    % Characterize the patches
    if has_spectral
        patch_penalties_spectral_L1{i} = zeros(n_patches_i, n_weights);
        patch_penalties_spectral_L2{i} = zeros(n_patches_i, n_weights);
        for pc = 1:n_patches_i
            patch = I_spectral_gt(...
                (corners_i(pc, 1) - padding):(corners_i(pc, 1) + (patch_size(1) + padding - 1)),...
                (corners_i(pc, 2) - padding):(corners_i(pc, 2) + (patch_size(2) + padding - 1)), : ...
            );
            for w = 1:n_weights
                err_vector = patch_operators_spectral{w} * reshape(patch, [], 1);
                patch_penalties_spectral_L1{i}(pc, w) = mean(abs(err_vector));
                patch_penalties_spectral_L2{i}(pc, w) = dot(err_vector, err_vector) / length(err_vector);
            end
        end
    end
    patch_penalties_rgb_L1{i} = zeros(n_patches_i, n_weights);
    patch_penalties_rgb_L2{i} = zeros(n_patches_i, n_weights);
    for pc = 1:n_patches_i
        patch = I_rgb_gt(...
            (corners_i(pc, 1) - padding):(corners_i(pc, 1) + (patch_size(1) + padding - 1)),...
            (corners_i(pc, 2) - padding):(corners_i(pc, 2) + (patch_size(2) + padding - 1)), : ...
        );
        for w = 1:n_weights
            err_vector = patch_operators_rgb{w} * reshape(patch, [], 1);
            patch_penalties_rgb_L1{i}(pc, w) = mean(abs(err_vector));
            patch_penalties_rgb_L2{i}(pc, w) = dot(err_vector, err_vector) / length(err_vector);
            if has_spectral && w == n_weights
                patch_penalties_spectral_L1{i}(pc, w) = patch_penalties_rgb_L1{i}(pc, w);
                patch_penalties_spectral_L2{i}(pc, w) = patch_penalties_rgb_L2{i}(pc, w);
            end
        end
    end

    % Run the algorithms
    if use_warped_spectral
        dispersion_options = struct('bands_in', bands_spectral);
        I_spectral_gt_unwarped = dispersionfunToMatrix(...
            df_spectral_forward, dispersion_options, I_spectral_gt, false...
        );
    elseif has_spectral
        I_spectral_gt_unwarped = I_spectral_gt;
    end
    if use_warped_rgb
        dispersion_options = struct('bands_in', bands_rgb);
        I_rgb_gt_unwarped = dispersionfunToMatrix(...
            df_rgb_forward, dispersion_options, I_rgb_gt, false...
        );
    elseif use_warped_spectral
        I_rgb_gt_unwarped = imageFormation(...
            I_spectral_gt_unwarped, bands_spectral, sensor_map, bands_color,...
            imageFormationSamplingOptions, imageFormationPatchOptions...
        );
    elseif has_rgb
        I_rgb_gt_unwarped = I_rgb_gt;
    end
    
    for f = 1:n_admm_algorithms
        algorithm = admm_algorithms.(admm_algorithm_fields{f});
        if ~algorithm.enabled || (algorithm.spectral && ~has_color_map) ||...
                (algorithm.spectral && channel_mode)
            continue;
        end
        for cr = 1:n_criteria
            if ~criteria(cr)
                continue;
            end
        
            if algorithm.spectral
                reg_options_f = solvePatchesSpectralOptions.reg_options;
            else
                reg_options_f = solvePatchesColorOptions.reg_options;
            end
            reg_options_f.enabled = algorithm.priors;
            enabled_weights = reg_options_f.enabled;
            n_active_weights = sum(enabled_weights);

            if algorithm.spectral
                admm_options_f = mergeStructs(...
                    solvePatchesSpectralOptions.admm_options, algorithm.options, false, true...
                );
            else
                admm_options_f = mergeStructs(...
                    solvePatchesColorOptions.admm_options, algorithm.options, false, true...
                );
            end
            
            if cr == dm_index
                reg_options_f.demosaic = true;
            else
                reg_options_f.demosaic = false;
            end
            
            weights_patches = zeros(n_patches_i, n_weights);
            n_steps = 1;
            have_steps = false;
            time_start = tic;
            if algorithm.spectral
                for pc = 1:n_patches_i
                    solvePatchesSpectralOptions.patch_options.target_patch = corners_i(pc, :);
                    if cr == mse_index
                        I_in.I = I_spectral_gt_unwarped;
                        I_in.spectral_bands = bands_spectral;
                        [...
                            bands_all, ~, ~, weights_images...
                        ] = solvePatchesSpectral(...
                          I_in, I_raw_gt, bayer_pattern, df_spectral_reverse,...
                          sensor_map, bands_color,...
                          solvePatchesSpectralOptions.sampling_options,...
                          admm_options_f, reg_options_f,...
                          solvePatchesSpectralOptions.patch_options,...
                          solvePatchesSpectralVerbose...
                        );
                    else
                        [...
                            bands_all, ~, ~, weights_images...
                        ] = solvePatchesSpectral(...
                          [], I_raw_gt, bayer_pattern, df_spectral_reverse,...
                          sensor_map, bands_color,...
                          solvePatchesSpectralOptions.sampling_options,...
                          admm_options_f, reg_options_f,...
                          solvePatchesSpectralOptions.patch_options,...
                          solvePatchesSpectralVerbose...
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
                    weights_patches(pc, reg_options_f.enabled, :) = reshape(weights_images(1, 1, :), 1, n_active_weights, []);
                end
            else
                for pc = 1:n_patches_i
                    solvePatchesColorOptions.patch_options.target_patch = corners_i(pc, :);
                    if cr == mse_index
                        I_in.I = I_rgb_gt_unwarped;
                        [...
                            ~, weights_images...
                        ] = solvePatchesColor(...
                          I_in, I_raw_gt, bayer_pattern, df_rgb_reverse,...
                          admm_options_f, reg_options_f,...
                          solvePatchesColorOptions.patch_options,...
                          solvePatchesColorVerbose...
                        );
                    else
                        [...
                            ~, weights_images...
                        ] = solvePatchesColor(...
                          [], I_raw_gt, bayer_pattern, df_rgb_reverse,...
                          admm_options_f, reg_options_f,...
                          solvePatchesColorOptions.patch_options,...
                          solvePatchesColorVerbose...
                        );
                    end
                    weights_patches(pc, reg_options_f.enabled) = reshape(weights_images(1, 1, :), 1, []);
                end
            end
            time_admm(f, i, cr) = toc(time_start) / n_patches_i;
            field_weights = geomean(weights_patches, 1);
            if algorithm.spectral && n_steps > 1
                field_weights = squeeze(field_weights).';
            end
            if i == 1
                if cr == 1
                    all_weights{f} = cell(n_images, n_criteria);
                end
                admm_algorithms.(admm_algorithm_fields{f}).(criteria_fields{cr}) = zeros(n_steps, n_weights, n_images); 
            end
            all_weights{f}{i, cr} = weights_patches;
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
all_weights_concat = cell(n_admm_algorithms, n_criteria);
for f = 1:n_admm_algorithms
    if admm_algorithms.(admm_algorithm_fields{f}).enabled &&...
            ~(admm_algorithms.(admm_algorithm_fields{f}).spectral && ~has_color_map)
        for cr = 1:n_criteria
            if criteria(cr)
                all_weights_concat{f, cr} = vertcat(all_weights{f}{:, cr});
                min_nz_weight = min(...
                    min_nz_weight, min(all_weights_concat{f, cr}(all_weights_concat{f, cr} ~= 0))...
                );
                max_nz_weight = max(...
                    max_nz_weight, max(all_weights_concat{f, cr}(all_weights_concat{f, cr} ~= 0))...
                );
            end
        end
    end
end
log_min_nz_weight = log10(min_nz_weight);
log_max_nz_weight = log10(max_nz_weight);
plot_limits = [log_min_nz_weight - 1, log_max_nz_weight + 1];

if has_spectral
    log_patch_penalties_spectral_L1 = log10(vertcat(patch_penalties_spectral_L1{:}));
    log_patch_penalties_spectral_L2 = log10(vertcat(patch_penalties_spectral_L2{:}));
end
log_patch_penalties_rgb_L1 = log10(vertcat(patch_penalties_rgb_L1{:}));
log_patch_penalties_rgb_L2 = log10(vertcat(patch_penalties_rgb_L2{:}));

for f = 1:n_admm_algorithms
    algorithm = admm_algorithms.(admm_algorithm_fields{f});
    if ~algorithm.enabled || (algorithm.spectral && ~has_color_map) ||...
            (algorithm.spectral && channel_mode)
        continue;
    end
        
    name_params = sprintf('%s_', algorithm.file);
    if algorithm.spectral
        name_params = [...
            sprintf('bands%d_', n_bands), name_params...
        ];
    else
        name_params = [...
            'RGB_', name_params...
        ];
    end
    
    if algorithm.spectral
        admm_options_f = mergeStructs(...
            solvePatchesSpectralOptions.admm_options, algorithm.options, false, true...
            );
    else
        admm_options_f = mergeStructs(...
            solvePatchesColorOptions.admm_options, algorithm.options, false, true...
            );
    end
    
    enabled_weights = algorithm.priors;
    n_active_weights = sum(enabled_weights);
    to_all_weights = find(enabled_weights);
    
    n_steps = max(cellfun(@(x) size(x, 3), all_weights{f, :}(:), 'UniformOutput', true));            
    for t = 1:n_steps
        if algorithm.spectral
            name_params_t = [...
                name_params, sprintf('step%d_', t), ...
                ];
        else
            name_params_t = name_params;
        end
        name_params_t = fullfile(output_directory, name_params_t);

        ref_weights = all_weights_concat{f, reference_criteria}(:, :, t);
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
            field_weights = all_weights_concat{f, cr}(:, :, t);
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
        if algorithm.spectral
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
        if n_active_weights == 1
            legend({'y = x', criteria_abbrev{criteria}});
        else
            legend(criteria_abbrev{criteria});
        end
        
        savefig(...
            fg,...
            [name_params_t 'weightsCorrelation.fig'], 'compact'...
            );
        close(fg);
        
        for w = 1:n_active_weights
            aw = to_all_weights(w);
            if algorithm.spectral
                if admm_options_f.norms(aw)
                    patch_penalties = log_patch_penalties_spectral_L1(:, w);
                else
                    patch_penalties = log_patch_penalties_spectral_L2(:, w);
                end
            else
                if admm_options_f.norms(aw)
                    patch_penalties = log_patch_penalties_rgb_L1(:, w);
                else
                    patch_penalties = log_patch_penalties_rgb_L2(:, w);
                end
            end

            fg = figure;
            hold on
            for cr = 1:n_criteria
                if ~criteria(cr)
                    continue;
                end
                field_weights = all_weights_concat{f, cr}(:, :, t);
                field_weights = field_weights(:, enabled_weights);
                log_weights = log10(field_weights);
                scatter(patch_penalties, log_weights(:, w), [], criteria_colors(cr, :), 'filled');
            end
            hold off
            if algorithm.spectral
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
            legend(criteria_abbrev{criteria});
            savefig(...
                fg,...
                [name_params_t sprintf('weight%d.fig', aw)], 'compact'...
                );
            close(fg);
        end
    end
    
    % Plot weights vs. image estimation step
    if algorithm.spectral && n_steps > 1
        n_patches_total = max(cellfun(@(x) size(x, 1), all_weights_concat(f, :), 'UniformOutput', true));  
        n_bands_plot = repelem(n_bands_all{f}, n_patches_total);
        for w = 1:n_active_weights
            aw = to_all_weights(w);

            fg = figure;
            hold on
            for cr = 1:n_criteria
                if ~criteria(cr)
                    continue;
                end
                field_weights = reshape(permute(all_weights_concat{f, cr}, [1, 3, 2]), [], n_weights);
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
    'all_weights', 'time_admm'...
} ];
if has_spectral
    save_variables_list = [save_variables_list, {...
        'bands_spectral', 'spectral_weights',...
        'patch_penalties_spectral_L1', 'patch_penalties_spectral_L2'...
    }];
end
if has_color_map
    save_variables_list = [save_variables_list, {'bands', 'bands_color'}];
end
save_data_filename = fullfile(output_directory, ['SelectWeightsForDataset_' dataset_name '.mat']);
save(save_data_filename, save_variables_list{:});
