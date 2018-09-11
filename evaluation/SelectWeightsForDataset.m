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
%   - 'mdf_weights': Regularization weights selected using the minimum
%     distance criterion of Song et al. 2016. 'mdf_weights' is a 2D array
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
% - 'all_mdf_weights': Regularization weights selected using the method of
%   Song et al. 2016. `all_mdf_weights(pc, :, i, f)` is the vector of
%   regularization weights selected for the f-th algorithm on the pc-th
%   patch of the i-th image
% - 'all_mse_weights': An array with the same layout as 'all_mdf_weights',
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
% - This script uses the first row of `patch_sizes`, and the first element
%   of `paddings`, defined in 'SetFixedParameters.m' to set patch sizes for
%   selecting regularization weights.
% - This script ignores the `downsampling_factor` parameter defined in
%   'SetFixedParameters.m'.
% - Regularization weights will not be selected for algorithms which are
%   disabled (using the 'enabled' fields in the structures describing the
%   algorithms). Instead, values of zero will be output for disabled
%   algorithms in 'all_mdf_weights' and 'all_mse_weights', and the
%   'mdf_weights' and 'mse_weights' fields will not be added to their
%   structures in 'admm_algorithms'.
% - Only image patches that are not clipped by the image edges will be
%   chosen for selecting regularization weights.
% - Image patch spectral derivatives, used for graphical output, not for
%   regularization weights selection, are computed according to the
%   'full_GLambda' field of the 'baek2017Algorithm2Options' structure
%   defined in 'SetFixedParameters.m', rather than according to
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

dataset_name = '20180817_TestSpectralDataset';

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

% Check for problematic parameters
if add_border
    % Estimating a border area results in images which are usually not
    % registered with the ground truth.
    error('Estimating a border around images prevents evaluation of the mean square error criterion.');
end

%% Preprocess the dataset

dp = describeDataset(dataset_name);

run('PreprocessDataset.m')

if has_spectral && ~can_evaluate_spectral
    error('The estimated spectral images must have the same channels as the true spectral images.');
end

%% Prepare for patch processing

n_weights = length(trainWeightsOptions.enabled_weights);

patch_size = patch_sizes(1, :);
padding = paddings(1);
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
            G_lambda = spectralGradient(image_sampling_patch_spectral, baek2017Algorithm2Options.full_GLambda);
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
        G_lambda = spectralGradient(image_sampling_patch_rgb, baek2017Algorithm2Options.full_GLambda);
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
    admm_algorithms.(admm_algorithm_fields{f}).mdf_weights = zeros(n_images, n_weights); 
    admm_algorithms.(admm_algorithm_fields{f}).mse_weights = zeros(n_images, n_weights);
end

%% Process the images

% Fixed options for ADMM
baek2017Algorithm2Options.add_border = false;
baek2017Algorithm2Options.l_surface = true;

if has_spectral
    patch_penalties_spectral_L1 = zeros(n_patches, n_weights, n_images);
    patch_penalties_spectral_L2 = zeros(n_patches, n_weights, n_images);
end
patch_penalties_rgb_L1 = zeros(n_patches, n_weights, n_images);
patch_penalties_rgb_L2 = zeros(n_patches, n_weights, n_images);

corners = zeros(n_patches, 2, n_images);

all_mdf_weights = zeros(n_patches, n_weights, n_images, n_admm_algorithms);
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
    
    corners(:, :, i) = [
        randi(image_sampling(1) - full_patch_size(1), n_patches, 1),...
        randi(image_sampling(2) - full_patch_size(2), n_patches, 1)
    ];

    % Characterize the patches
    if has_spectral
        for pc = 1:n_patches
            patch = I_spectral_gt(...
                corners(pc, 1, i):(corners(pc, 1, i) + (full_patch_size(1) - 1)),...
                corners(pc, 2, i):(corners(pc, 2, i) + (full_patch_size(2) - 1)), : ...
            );
            patch_cropped = patch(...
                (padding + 1):(end - padding), (padding + 1):(end - padding), : ...
            );
            for w = 1:n_weights
                if w > n_spectral_weights
                    continue;
                end
                err_vector = patch_operators_spectral{w} * reshape(patch_cropped, [], 1);
                patch_penalties_spectral_L1(pc, w, i) = mean(abs(err_vector));
                patch_penalties_spectral_L2(pc, w, i) = dot(err_vector, err_vector) / length(err_vector);
            end
        end
    end
    for pc = 1:n_patches
        patch = I_rgb_gt(...
            corners(pc, 1, i):(corners(pc, 1, i) + (full_patch_size(1) - 1)),...
            corners(pc, 2, i):(corners(pc, 2, i) + (full_patch_size(2) - 1)), : ...
        );
        patch_cropped = patch(...
            (padding + 1):(end - padding), (padding + 1):(end - padding), : ...
        );
        for w = 1:n_weights
            err_vector = patch_operators_rgb{w} * reshape(patch_cropped, [], 1);
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

        if algorithm.spectral
            if has_color_map
                if channel_mode
                    baek2017Algorithm2Options.int_method = 'none';
                    selectWeightsGridOptions.int_method = 'none';
                    trainWeightsOptions.int_method = 'none';
                else
                    baek2017Algorithm2Options.int_method = int_method;
                    selectWeightsGridOptions.int_method = int_method;
                    trainWeightsOptions.int_method = int_method;
                end
            else
                continue;
            end
        else
            baek2017Algorithm2Options.int_method = 'none';
            selectWeightsGridOptions.int_method = 'none';
            trainWeightsOptions.int_method = 'none';
        end

        baek2017Algorithm2Options_f = mergeStructs(...
            baek2017Algorithm2Options, algorithm.options, false, true...
        );
    
        selectWeightsGridOptions_f = selectWeightsGridOptions;
        selectWeightsGridOptions_f.enabled_weights = algorithm.priors;
        trainWeightsOptions_f = trainWeightsOptions;
        trainWeightsOptions_f.enabled_weights = algorithm.priors;
        
        mdf_weights_patches = zeros(n_patches, n_weights);
        mse_weights_patches = zeros(n_patches, n_weights);
        if algorithm.spectral
            for pc = 1:n_patches
                mdf_weights_patches(pc, :) = selectWeightsGrid(...
                    I_raw_gt, bayer_pattern, df_spectral_reverse,...
                    sensor_map_resampled, bands,...
                    rho, baek2017Algorithm2Options_f, selectWeightsGridOptions_f,...
                    corners(pc, :, i), selectWeightsGridVerbose...
                );
                mse_weights_patches(pc, :) = trainWeights(...
                    I_spectral_gt, I_raw_gt, bayer_pattern, df_spectral_reverse,...
                    sensor_map_resampled, bands, trainWeightsOptions_f,...
                    @baek2017Algorithm2, {...
                        rho, baek2017Algorithm2Options_f, false...
                    }, corners(pc, :, i), trainWeightsVerbose...
                );
            end 
        else
            for pc = 1:n_patches
                mdf_weights_patches(pc, :) = selectWeightsGrid(...
                    I_raw_gt, bayer_pattern, df_rgb_reverse,...
                    sensor_map_rgb, bands_rgb,...
                    rho, baek2017Algorithm2Options_f, selectWeightsGridOptions_f,...
                    corners(pc, :, i), selectWeightsGridVerbose...
                );
                mse_weights_patches(pc, :) = trainWeights(...
                    I_rgb_gt, I_raw_gt, bayer_pattern, df_rgb_reverse,...
                    sensor_map_rgb, bands_rgb, trainWeightsOptions_f,...
                    @baek2017Algorithm2, {...
                        rho, baek2017Algorithm2Options_f, false...
                    }, corners(pc, :, i), trainWeightsVerbose...
                );
            end
        end
        mdf_weights = geomean(mdf_weights_patches, 1);
        mse_weights = geomean(mse_weights_patches, 1);
        all_mdf_weights(:, :, i, f) = mdf_weights_patches;
        all_mse_weights(:, :, i, f) = mse_weights_patches;
        admm_algorithms.(admm_algorithm_fields{f}).mdf_weights(i, :) = mdf_weights;
        admm_algorithms.(admm_algorithm_fields{f}).mse_weights(i, :) = mse_weights;
    end

    if verbose
        fprintf('[SelectWeightsForDataset, image %d] Finished\n', i);
    end
end

%% Visualization of weights selected for the dataset

mdf_color = [1, 0, 0];
mse_color = [0, 1, 0];

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
    
    baek2017Algorithm2Options_f = mergeStructs(...
        baek2017Algorithm2Options, algorithm.options, false, true...
    );

    enabled_weights = algorithm.priors;
    n_active_weights = sum(enabled_weights);
    to_all_weights = find(enabled_weights);

    mdf_weights = reshape(permute(all_mdf_weights(:, :, :, f), [1, 3, 2, 4]), [], n_weights);
    mdf_weights = mdf_weights(:, enabled_weights);
    log_mdf_weights = log10(mdf_weights);
    mse_weights = reshape(permute(all_mse_weights(:, :, :, f), [1, 3, 2, 4]), [], n_weights);
    mse_weights = mse_weights(:, enabled_weights);
    log_mse_weights = log10(mse_weights);

    fg = figure;
    hold on
    if n_active_weights == 1
        scatter(...
            log_mse_weights, log_mdf_weights, 'filled'...
        );
        line_limits = [...
            min(min(log_mse_weights, log_mdf_weights));
            max(max(log_mse_weights, log_mdf_weights))
        ];
        line(line_limits, line_limits, 'Color', 'b');
        legend('Weights', 'y = x');
        xlabel('Weight selected using the mean square error');
        ylabel('Weight selected using the minimum distance criterion');
    elseif n_active_weights == 2
        scatter(...
            log_mdf_weights(:, 1), log_mdf_weights(:, 2), [], mdf_color, 'filled'...
        );
        scatter(...
            log_mse_weights(:, 1), log_mse_weights(:, 2), [], mse_color, 'filled'...
        );
        legend('MDF', 'MSE');
        xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
        ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
    elseif n_active_weights == 3
        scatter3(...
            log_mdf_weights(:, 1), log_mdf_weights(:, 2), log_mdf_weights(:, 3),...
            [], mdf_color, 'filled'...
        );
        scatter3(...
            log_mse_weights(:, 1), log_mse_weights(:, 2), log_mse_weights(:, 3),...
            [], mse_color, 'filled'...
        );
        legend('MDF', 'MSE');
        xlabel(sprintf('log_{10}(weight %d)', to_all_weights(1)))
        ylabel(sprintf('log_{10}(weight %d)', to_all_weights(2)))
        zlabel(sprintf('log_{10}(weight %d)', to_all_weights(3)))
    else
        error('Unexpected number of active weights.');
    end
    title(sprintf(...
        'Agreement between MDF and MSE weights selected for %s',...
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
            if baek2017Algorithm2Options_f.norms(aw)
                patch_penalties = log_patch_penalties_spectral_L1(:, w, :);
            else
                patch_penalties = log_patch_penalties_spectral_L2(:, w, :);
            end
        else
            if baek2017Algorithm2Options_f.norms(aw)
                patch_penalties = log_patch_penalties_rgb_L1(:, w, :);
            else
                patch_penalties = log_patch_penalties_rgb_L2(:, w, :);
            end
        end
        patch_penalties = reshape(patch_penalties, [], 1);
        
        fg = figure;
        hold on
        scatter(patch_penalties, log_mdf_weights(:, w), [], mdf_color, 'filled');
        scatter(patch_penalties, log_mse_weights(:, w), [], mse_color, 'filled');
        legend('MDF', 'MSE');
        xlabel(sprintf('log_{10}(Penalty %d)', aw))
        ylabel(sprintf('log_{10}(Weight %d)', aw))
        title(sprintf(...
            'Agreement between MDF and MSE weights selected for %s',...
            algorithm.name...
        ));
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
    'all_mdf_weights', 'all_mse_weights'...
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