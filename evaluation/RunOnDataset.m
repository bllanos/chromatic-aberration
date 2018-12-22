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
% The documentation in the script 'CorrectByHyperspectralADMM.m' contains
% more information on the formats of the various types of data associated
% with the datasets.
%
% This script also runs 'SetFixedParameters.m' to set the values of
% seldomly-changed and/or common parameters. These parameters are briefly
% documented in 'SetFixedParameters.m'. Regularization weights set in
% 'SetFixedParameters.m' are ignored; Regularization weights are
% automatically selected if there are no regularization weights associated
% with the individual ADMM-family algorithms, loaded from a file generated
% by 'SelectWeightsForDataset.m'. An output file from
% 'SelectWeightsForDataset.m' will also override the set of ADMM-family
% algorithms to run, which is otherwise determined by running
% 'SetAlgorithms.m'.
%
% ## Output
%
% ### Estimated images
%
% The following types of images are created for each input image, depending
% on the image estimation algorithms. The filename of the input image,
% concatenated with a string of parameter information, is represented by
% '*' below:
% - '*_roi.tif': A cropped version of the input image, containing the
%   portion used as input. This region of interest was determined using the
%   domain of the model of dispersion associated with the dataset. If no
%   model of dispersion is associated with the dataset, the cropped region
%   is the entire input image. All of the other output images listed below
%   are limited to the region shown in this output image. Note that this
%   image is normalized by its maximum value - It is intended for viewing,
%   not for quantitative analysis.
% - '*_latent.mat': The estimated latent spectral image (stored in the
%   variable 'I_latent') corresponding to the input image.
% - '*_rgb.tif': A colour image. If it was not estimated directly, it was
%   created by converting the latent image to the RGB colour space of the
%   input image. Images which are estimated directly are also saved as
%   '.mat' files (under the variable name 'I_rgb').
%
% If the dataset contains images affected by dispersion, and if there are
% models of dispersion for use during image estimation, additional images
% are saved:
% - '*_latent_ab.mat': A version of '*_latent.mat' subject to dispersion.
% - '*_rgb_ab.tif': A version of '*_rgb.tif' subject to dispersion.
%
% For demosaicking algorithms, the colour images are saved under the names
% of the demosaicking algorithms, with no '_rgb' suffix. Furthermore, if
% there is a model of colour dispersion associated with the dataset,
% additional images, '_channelWarp.tif' and '_channelWarp.mat' (using the
% variable name 'I_rgb') are output. These images are corrections of the
% demosaicking results for chromatic aberration using the colour dispersion
% model.
%
% ### Regularization weights images
%
% If automatic regularization weight selection is enabled (see
% `admm_algorithms_filename` in the parameters below), then the image
% estimation algorithm will automatically choose weights on the
% regularization terms in the ADMM optimization problem. For the i-th
% enabled regularization term in the ADMM optimization problem, an image
% will be output, as the variable 'I_weights', in the file
% '*_weight${i}Image.mat', where '*' represents the filename of the input
% image concatenated with a string of parameter information. A pixel in the
% image will contain the weight on the i-th regularization term used when
% estimating the pixel.
%
% ### Data file output
%
% #### Intermediate data and parameters
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
% - 'radiance_normalized_weights': This variable is only output for
%   datasets of reflectance images. It is a matrix for converting pixels in
%   the ground truth reflectance images to normalized radiances. As such,
%   it accounts for both the illuminant used to generate radiances from
%   reflectances, and the image sensor for which the radiances are
%   normalized.
% - 'spectral_weights': A matrix for converting pixels in the spectral
%   space of the estimated hyperspectral images to the spectral space of
%   the true hyperspectral images.
% - 'admm_algorithms': A structure describing the ADMM algorithms being
%   evaluated, created by 'SetAlgorithms.m' and possibly updated by
%   'SelectWeightsForDataset.m'.
% - 'demosaic_algorithms': A structure describing the demosaicking
%   algorithms being evaluated, created by 'SetAlgorithms.m'.
% - 'time': A structure containing execution timing information, measured in
%   seconds. 'time' has the following fields:
%   - 'admm': Execution timing information, stored as a 3D array, for
%     ADMM-family algorithms. `time.admm(f, i, cr)` is the time taken to process
%     the i-th image with the f-th ADMM-family algorithm defined in
%     'SetAlgorithms.m', according to weights selected using the cr-th
%     regularization weight selection criterion. Entries corresponding to
%     disabled algorithms or disabled weight selection criterion will be set to
%     `NaN`.
%   - 'demosaic': Execution timing information, stored as a 2D array, for
%     demosaicing algorithms. `time.demosaic(f, i)` is the time taken to process
%     the i-th image with the f-th demosaicing algorithm defined in
%     'SetAlgorithms.m'. Entries corresponding to disabled algorithms will be
%     set to `NaN`.
%   - 'warp': This field exists for datasets which have models of colour
%     channel-space dispersion. 'time.warp' is a vector where the elements are
%     the times taken to calculate warping matrices for correcting dispersion in
%     the images.
%   - 'warp_apply': This field exists for datasets which have models of
%     colour channel-space dispersion. 'time.warp_apply' is a vector where the
%     elements are the times taken to apply warping matrices for correcting
%     dispersion to the images.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`, which is initialized in this file, and then augmented
% by 'SetFixedParameters.m'.
%
% The file is saved as 'RunOnDataset_${dataset_name}.mat'.
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
% If the dataset contains images affected by dispersion, and if there are
% models of dispersion for use during image estimation, versions of the
% above files are saved for evaluations performed on versions of the
% estimated images subject to dispersion. The resulting files are given
% names ending in '_ab'.
%
% ## Notes
% - This script only uses the first row of `patch_sizes`, and the first
%   element of `paddings`, defined in 'SetFixedParameters.m', by using
%   `solvePatchesADMMOptions.patch_options`.
% - If `solvePatchesMultiADMMOptions.sampling_options.show_steps` is
%   `true`, output images will be saved only for the highest spectral
%   resolution, even though the images returned by
%   'solvePatchesMultiADMM()' will contain multiple spectral resolutions.
%   However, there will be additional spectral evaluation figures and CSV
%   files comparing the results between spectral resolutions (for the same
%   image estimation algorithm). Refer to the documentation of
%   'solvePatchesMultiADMM.m' for more information.
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
%
% The RGB-to-spectral image estimation comparison methods are:
%
%   Choi, I., Jeon, D. S., Gutierrez, D., & Kim, M. H. (2017).
%   "High-Quality Hyperspectral Reconstruction Using a Spectral Prior." ACM
%   Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 218:1-13.
%   10.1145/3130800.3130810

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 27, 2018

% List of parameters to save with results
parameters_list = {
        'dataset_name',...
        'admm_algorithms_filename',...
        'output_directory'...
    };

%% Input data and parameters

dataset_name = '20181212_RealData_RGBAsRAW';

% Describe algorithms to run
run('SetAlgorithms.m')

% Optionally override the list of ADMM-family algorithms to run, and the
% regularization weights to run them with, from the output file of
% 'SelectWeightsForDataset.m'. (Leave empty otherwise)
admm_algorithms_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/results/rgb_raw_images/weights_selection/SelectWeightsForDataset_20181212_RealData_RGBAsRAW.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads/results/RunOnDataset';

% Produce console output to describe the processing in this script
verbose = true;

% ## Parameters which do not usually need to be changed
% Note that this sets the value of `criteria`.
run('SetFixedParameters.m')

%% Check for problematic parameters

if use_fixed_weights
    error('Weights should be fixed by running ''SelectWeightsForDataset.m'', not using the `use_fixed_weights` parameter in ''SetFixedParameters.m''');
end

use_automatic_weights = isempty(admm_algorithms_filename);
if use_automatic_weights && sum(criteria) == 0
    error('All regularization weight selection criteria are disabled.');
end

%% Preprocess the dataset

dp = describeDataset(dataset_name);

run('PreprocessDataset.m')

%% Finalize the set of algorithms to run

if ~use_automatic_weights
    load(admm_algorithms_filename, 'admm_algorithms');
end
admm_algorithm_fields = fieldnames(admm_algorithms);
n_admm_algorithms = length(admm_algorithm_fields);
n_criteria = length(criteria_fields);
if ~use_automatic_weights
    for f = 1:n_admm_algorithms
        algorithm = admm_algorithms.(admm_algorithm_fields{f});
        if algorithm.enabled 
            for cr = 1:n_criteria
                if f == 1
                    criteria(cr) = isfield(algorithm, criteria_fields{cr});
                elseif criteria(cr) ~= isfield(algorithm, criteria_fields{cr})
                    error('Different algorithms have regularization weights selected using different methods.');
                end
                if criteria(cr) && (size(algorithm.(criteria_fields{cr}), 3) ~= n_images)
                    error('Algorithm "%s" does not have a number of "%s" weights matching the number of images.',...
                        algorithm.name, criteria_names{cr})
                end
            end
        end
    end
end
n_active_criteria = sum(criteria);

%% Process the images

patch_size = solvePatchesMultiADMMOptions.patch_options.patch_size;
padding = solvePatchesMultiADMMOptions.patch_options.padding;

e_rgb_tables = cell(n_images, 1);
if evaluate_aberrated_rgb
    e_rgb_tables_ab = cell(n_images, 1);
end
e_spectral_tables = cell(n_images, 1);
if evaluate_aberrated_spectral
    e_spectral_tables_ab = cell(n_images, 1);
end

time.admm = nan(n_admm_algorithms, n_images, n_criteria);
demosaic_algorithm_fields = fieldnames(demosaic_algorithms);
n_demosaic_algorithms = length(demosaic_algorithm_fields);
time.demosaic = nan(n_demosaic_algorithms, n_images);
if has_dispersion_rgb
    time.warp = zeros(1, n_images);
    time.warp_apply = zeros(1, n_images);
end

for i = 1:n_images
    if verbose
        fprintf('[RunOnDataset, image %d] Starting\n', i);
    end

    % Generate or load input images, and instantiate dispersion information
    run('LoadAndConvertImage.m');
        
    saveImages(...
        'image', output_directory, names{i},...
        I_raw_gt ./ max(I_raw_gt(:)), '_roi', 'I_raw'...
    );
    
    % Compare the aberrated image to the original
    
    if isempty(I_rgb_gt_warped) || dp.is_aberrated
        e_rgb_table = [];
    else
        e_rgb_table = evaluateAndSaveRGB(...
            I_rgb_gt_warped, I_rgb_gt, dp, names{i}, 'Aberrated',...
            fullfile(output_directory, [names{i} '_aberrated'])...
        );
    end
    if evaluate_aberrated_rgb
        e_rgb_table_ab = [];
    end
    
    n_spectral_evaluations = 0;
    if has_spectral
        for f = 1:n_admm_algorithms
            algorithm = admm_algorithms.(admm_algorithm_fields{f});
            if algorithm.enabled && algorithm.spectral
                n_spectral_evaluations = n_spectral_evaluations + 1;
            end
        end
        n_spectral_evaluations = n_spectral_evaluations * n_active_criteria;
        n_spectral_evaluations_admm = n_spectral_evaluations;
    end
    evaluation_ind = 0;
    if ~isempty(I_spectral_gt_warped) && ~dp.is_aberrated
        n_spectral_evaluations = n_spectral_evaluations + 1;
        evaluation_ind = evaluation_ind + 1;
        aberrated_evaluation_ind = evaluation_ind;
    end
    if has_choi_spectral
        n_spectral_evaluations = n_spectral_evaluations + 1;
        evaluation_ind = evaluation_ind + 1;
        choi_evaluation_ind = evaluation_ind;
    end
    if n_spectral_evaluations > 0
        evaluation_plot_colors = jet(n_spectral_evaluations);
        evaluation_plot_colors_admm = evaluation_plot_colors((end - n_spectral_evaluations_admm + 1):end, :);
        evaluation_plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
        evaluation_plot_styles = {'--', ':', '-.'};
    end
    if isempty(I_spectral_gt_warped) || dp.is_aberrated
        e_spectral_table = [];
        fg_spectral = struct;
        all_alg_names = {};
    else
        dp.evaluation.global_spectral.plot_color = evaluation_plot_colors(aberrated_evaluation_ind, :);
        dp.evaluation.global_spectral.plot_marker = 'none';
        dp.evaluation.global_spectral.plot_style = '-';
        all_alg_names = {'Aberrated'};
        [e_spectral_table, fg_spectral] = evaluateAndSaveSpectral(...
            I_spectral_gt_warped, I_spectral_gt, bands_spectral,...
            eye(length(bands_spectral)), dp, names{i}, all_alg_names{aberrated_evaluation_ind},...
            fullfile(output_directory, [names{i} '_aberrated'])...
        );
    end
    
    if evaluate_aberrated_spectral
        e_spectral_table_ab = [];
        fg_spectral_ab = struct;
    end
    
    % Evaluate comparison methods
    if has_choi_spectral
        dp.evaluation.global_spectral.plot_color = evaluation_plot_colors(choi_evaluation_ind, :);
        dp.evaluation.global_spectral.plot_marker = 'none';
        dp.evaluation.global_spectral.plot_style = '--';
        all_alg_names{choi_evaluation_ind} = 'Choi et al. 2017';
        choi_img = loadImage(choi_spectral_filenames{i}, 'I_latent');
        [e_spectral_table_current, fg_spectral] = evaluateAndSaveSpectral(...
            choi_img, I_spectral_gt, bands_spectral,...
            eye(length(bands_spectral)), dp, names{i}, all_alg_names{choi_evaluation_ind},...
            fullfile(output_directory, [names{i} '_choi']), fg_spectral...
        );
        if ~isempty(e_spectral_table)
            e_spectral_table = union(e_spectral_table_current, e_spectral_table);
        else
            e_spectral_table = e_spectral_table_current;
        end
        
        if evaluate_aberrated_spectral
            [e_spectral_table_current, fg_spectral_ab] = evaluateAndSaveSpectral(...
                choi_img, I_spectral_gt_warped, bands_spectral,...
                eye(length(bands_spectral)), dp, names{i}, all_alg_names{choi_evaluation_ind},...
                fullfile(output_directory, [names{i} '_choi_ab']), fg_spectral_ab...
            );
            if ~isempty(e_spectral_table_ab)
                e_spectral_table_ab = union(e_spectral_table_current, e_spectral_table_ab);
            else
                e_spectral_table_ab = e_spectral_table_current;
            end
        end
    end
    
    if has_choi_rgb
        choi_img = loadImage(choi_rgb_filenames{i}, 'I_rgb');
        e_rgb_table_current = evaluateAndSaveRGB(...
            choi_img, I_rgb_gt, dp, names{i}, all_alg_names{choi_evaluation_ind},...
            fullfile(output_directory, [names{i} '_choi'])...
        );
        if ~isempty(e_rgb_table)
            e_rgb_table = union(e_rgb_table_current, e_rgb_table);
        else
            e_rgb_table = e_rgb_table_current;
        end
        
        if evaluate_aberrated_rgb
            e_rgb_table_current = evaluateAndSaveRGB(...
                choi_img, I_rgb_gt_warped, dp, names{i}, all_alg_names{choi_evaluation_ind},...
                fullfile(output_directory, [names{i} '_choi_ab'])...
            );
            if ~isempty(e_rgb_table_ab)
                e_rgb_table_ab = union(e_rgb_table_current, e_rgb_table_ab);
            else
                e_rgb_table_ab = e_rgb_table_current;
            end
        end
    end
    
    % Run own algorithms
    
    % ADMM
    color_ind = 1;
    if evaluate_aberrated_spectral || evaluate_aberrated_rgb
        extra_images = cell(3, 1);
    else
        extra_images = cell(0, 1);
    end
    for cr = 1:n_criteria
        if ~criteria(cr)
            continue;
        end
        for f = 1:n_admm_algorithms
            algorithm = admm_algorithms.(admm_algorithm_fields{f});
            if ~algorithm.enabled || (algorithm.spectral && ~has_color_map)
                continue;
            end
            
            true_spectral = algorithm.spectral && ~channel_mode;
            if true_spectral
                reg_options_f = solvePatchesMultiADMMOptions.reg_options;
            else
                reg_options_f = solvePatchesADMMOptions.reg_options;
            end
            
            if ~use_automatic_weights
                weights_f = algorithm.(criteria_fields{cr})(:, :, i);
                if true_spectral
                    reg_options_f.multi_weights = weights_f;
                else
                    reg_options_f.minimum_weights = weights_f;
                    reg_options_f.maximum_weights = weights_f;
                end
            end
            reg_options_f.enabled = algorithm.priors;
            
            if true_spectral
                admm_options_f = mergeStructs(...
                    solvePatchesMultiADMMOptions.admm_options, algorithm.options, false, true...
                );
            else
                admm_options_f = mergeStructs(...
                    solvePatchesADMMOptions.admm_options, algorithm.options, false, true...
                );
            end

            name_params = sprintf(...
                '%s_patch%dx%d_pad%d',...
                algorithm.file, patch_size(1), patch_size(2), padding...
                );
            alg_name_params = sprintf(...
                '%s, patch %d x %d, padding %d',...
                algorithm.name, patch_size(1), patch_size(2), padding...
                );
            if use_automatic_weights
                weights_filepart = ['_', criteria_filenames{cr}];
                name_params = [name_params, weights_filepart];
                alg_name_params = [alg_name_params, sprintf(', %s', criteria_abbrev{cr})];
                enabled_weights = reg_options_f.enabled;
                n_active_weights = sum(enabled_weights);
                to_all_weights = find(enabled_weights);
            else
                weights_filepart = ['_', criteria_abbrev{cr} 'fw_'];
                name_params = [name_params, weights_filepart];
                alg_name_params = [alg_name_params, sprintf(', %s fixed weights', criteria_abbrev{cr})];
            end
            if use_automatic_weights && cr == dm_index
                reg_options_f.demosaic = true;
            else
                reg_options_f.demosaic = false;
            end
            
            time_start = tic;
            if algorithm.spectral                
                if true_spectral
                    if cr == mse_index && use_automatic_weights
                        I_in.I = I_spectral_gt;
                        I_in.spectral_bands = bands_spectral;
                        [...
                            bands_all,...
                            I_latent,...
                            I_rgb,...
                            weights_images,...
                            extra_images{:}...
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
                            bands_all,...
                            I_latent,...
                            I_rgb,...
                            weights_images,...
                            extra_images{:}...
                        ] = solvePatchesMultiADMM(...
                            [], I_raw_gt, bayer_pattern, df_spectral_reverse,...
                            sensor_map, bands_color,...
                            solvePatchesMultiADMMOptions.sampling_options,...
                            admm_options_f, reg_options_f,...
                            solvePatchesMultiADMMOptions.patch_options,...
                            solvePatchesMultiADMMVerbose...
                        );
                    end
                    
                    time.admm(f, i, cr) = toc(time_start);
                    
                    if solvePatchesMultiADMMOptions.sampling_options.show_steps
                        n_bands_t = length(bands_all{end});
                    else
                        n_bands_t = length(bands_all);
                    end
                    name_params = [...
                        names{i}, sprintf('_bands%d_', n_bands_t), name_params...
                    ];
                    alg_name_params = [...
                        alg_name_params, sprintf(', %d bands', n_bands_t)...
                    ];
                
                    % Spectral evaluation of intermediate images
                    if solvePatchesMultiADMMOptions.sampling_options.show_steps
                        n_steps = length(bands_all);
                        step_plot_colors = jet(n_steps);
                        step_name_params_tables = cell(n_steps, 1);
                        
                        spectral_inc = 0;
                        fg_spectral_step = struct;
                        if evaluate_aberrated_spectral
                            fg_spectral_step_ab = struct;
                        end
                        for t = 1:n_steps
                            name_params_t = [name_params, sprintf('step%d_', t)];
                            name_params_tables_t = sprintf('%s, step %d', alg_name_params, t);
                            dp.evaluation.global_spectral.plot_color =...
                                step_plot_colors(t, :);
                            dp.evaluation.global_spectral.plot_marker =...
                                evaluation_plot_markers{...
                                mod(t - 1, length(evaluation_plot_markers)) + 1 ...
                                };
                            dp.evaluation.global_spectral.plot_style =...
                                evaluation_plot_styles{...
                                mod(t - 1, length(evaluation_plot_styles)) + 1 ...
                                };
                            step_name_params_tables{t} = name_params_tables_t;

                            n_bands_t = length(bands_all{t});
                            spectral_weights_step = resamplingWeights(...
                                bands_spectral, bands_all{t},...
                                solvePatchesMultiADMMOptions.sampling_options.interpolant,...
                                solvePatchesMultiADMMOptions.sampling_options.bands_padding...
                            );
                            [e_spectral_table_step_current, fg_spectral_step] = evaluateAndSaveSpectral(...
                                I_latent(:, :, (spectral_inc + 1):(spectral_inc + n_bands_t)),...
                                I_spectral_gt, bands_spectral, spectral_weights_step,...
                                dp, names{i}, name_params_tables_t,...
                                fullfile(output_directory, name_params_t(1:(end-1))),...
                                fg_spectral_step...
                                );
                            if t == 1
                                e_spectral_table_step = e_spectral_table_step_current;
                            else
                                e_spectral_table_step = union(e_spectral_table_step_current, e_spectral_table_step);
                            end
                            
                            if evaluate_aberrated_spectral
                                [e_spectral_table_step_current, fg_spectral_step_ab] = evaluateAndSaveSpectral(...
                                    extra_images{3}(:, :, (spectral_inc + 1):(spectral_inc + n_bands_t)),...
                                    I_spectral_gt_warped, bands_spectral, spectral_weights_step,...
                                    dp, names{i}, name_params_tables_t,...
                                    fullfile(output_directory, [name_params_t, 'ab']),...
                                    fg_spectral_step_ab...
                                    );
                                if t == 1
                                    e_spectral_table_step_ab = e_spectral_table_step_current;
                                else
                                    e_spectral_table_step_ab = union(e_spectral_table_step_current, e_spectral_table_step_ab);
                                end
                            end
                            spectral_inc = spectral_inc + n_bands_t;
                        end
                        writetable(...
                            e_spectral_table_step,...
                            fullfile(output_directory, [name_params, 'multiStep_evaluateSpectral.csv'])...
                        );
                        if evaluate_aberrated_spectral
                            writetable(...
                                e_spectral_table_step_ab,...
                                fullfile(output_directory, [name_params, 'multiStep_evaluateSpectral_ab.csv'])...
                            );
                        end
                        name_step = [name_params, 'MS'];
                        dp.evaluation.custom_spectral.(name_step) = dp.evaluation.custom_spectral.(names{i});
                        evaluateAndSaveSpectral(...
                            output_directory, dp, name_step, step_name_params_tables, fg_spectral_step...
                        );
                        if evaluate_aberrated_spectral
                            name_step_ab = [name_params, 'MSab'];
                            dp.evaluation.custom_spectral.(name_step_ab) = dp.evaluation.custom_spectral.(names{i});
                            evaluateAndSaveSpectral(...
                                output_directory, dp, name_step_ab, step_name_params_tables, fg_spectral_step_ab...
                            );
                        end

                        % Retain only the highest spectral resolution data
                        % for further study
                        I_latent = I_latent(:, :, (end - n_bands_t + 1):end);
                        I_rgb = I_rgb(:, :, (end - n_channels_rgb + 1):end);
                        if ~isempty(extra_images)
                            extra_images{3} = extra_images{3}(:, :, (end - n_bands_t + 1):end);
                            extra_images{1} = extra_images{1}(:, :, (end - n_channels_rgb + 1):end);
                        end
                        if use_automatic_weights
                            weights_images = weights_images(:, :, (end - n_active_weights + 1):end);
                        end
                    end
                else
                    if cr == mse_index && use_automatic_weights
                        I_in.I = I_spectral_gt;
                        I_in.spectral_weights = spectral_weights;
                        [...
                            I_latent,...
                            I_rgb,...
                            weights_images,...
                            extra_images{:}...
                        ] = solvePatchesADMM(...
                            I_in, I_raw_gt, bayer_pattern, df_spectral_reverse,...
                            color_weights, bands,...
                            admm_options_f, reg_options_f,...
                            solvePatchesADMMOptions.patch_options,...
                            solvePatchesADMMVerbose...
                        );
                    else
                        [...
                            I_latent,...
                            I_rgb,...
                            weights_images,...
                            extra_images{:}...
                        ] = solvePatchesADMM(...
                            [], I_raw_gt, bayer_pattern, df_spectral_reverse,...
                            color_weights, bands,...
                            admm_options_f, reg_options_f,...
                            solvePatchesADMMOptions.patch_options,...
                            solvePatchesADMMVerbose...
                        );
                    end
                    
                    time.admm(f, i, cr) = toc(time_start);
                
                    name_params = [...
                        names{i}, sprintf('_bands%d_', n_bands), name_params...
                    ];
                    alg_name_params = [...
                        alg_name_params, sprintf(', %d bands', n_bands)...
                    ];
                end
                
                saveImages(...
                    'data', output_directory, name_params,...
                    I_latent, 'latent', 'I_latent'...
                );
                saveImages(...
                    'image', output_directory, name_params,...
                    I_rgb, 'rgb', 'I_rgb'...
                );
                if ~isempty(extra_images)
                    saveImages(...
                        'data', output_directory, name_params,...
                        extra_images{3}, 'latent_ab', 'I_latent'...
                    );
                    saveImages(...
                        'image', output_directory, name_params,...
                        extra_images{1}, 'rgb_ab', 'I_rgb'...
                    );
                end
            
                % Spectral evaluation
                if has_spectral
                    dp.evaluation.global_spectral.plot_color =...
                        evaluation_plot_colors_admm(color_ind, :);
                    dp.evaluation.global_spectral.plot_marker =...
                        evaluation_plot_markers{...
                            mod(color_ind - 1, length(evaluation_plot_markers)) + 1 ...
                        };
                    dp.evaluation.global_spectral.plot_style =...
                        evaluation_plot_styles{...
                            mod(color_ind - 1, length(evaluation_plot_styles)) + 1 ...
                        };
                    color_ind = color_ind + 1;
                    all_alg_names{end + 1} = alg_name_params;
                    [e_spectral_table_current, fg_spectral] = evaluateAndSaveSpectral(...
                        I_latent, I_spectral_gt, bands_spectral, spectral_weights,...
                        dp, names{i}, alg_name_params,...
                        fullfile(output_directory, name_params(1:(end-1))),...
                        fg_spectral...
                    );
                    if ~isempty(e_spectral_table)
                        e_spectral_table = union(e_spectral_table_current, e_spectral_table);
                    else
                        e_spectral_table = e_spectral_table_current;
                    end
                    
                    if evaluate_aberrated_spectral
                        [e_spectral_table_current, fg_spectral_ab] = evaluateAndSaveSpectral(...
                            extra_images{3}, I_spectral_gt_warped, bands_spectral, spectral_weights,...
                            dp, names{i}, alg_name_params,...
                            fullfile(output_directory, [name_params, 'ab']),...
                            fg_spectral_ab...
                        );
                        if ~isempty(e_spectral_table_ab)
                            e_spectral_table_ab = union(e_spectral_table_current, e_spectral_table_ab);
                        else
                            e_spectral_table_ab = e_spectral_table_current;
                        end
                    end
                end
            else
                name_params = [...
                    names{i}, '_RGB_', name_params...
                ];
                alg_name_params = [...
                    alg_name_params, ', RGB'...
                ];
            
                if cr == mse_index && use_automatic_weights
                    I_in.I = I_rgb_gt;
                    I_in.spectral_weights = sensor_map_rgb;
                    [...
                        I_rgb,...
                        ~,...
                        weights_images,...
                        extra_images{:}...
                    ] = solvePatchesADMM(...
                      I_in, I_raw_gt, bayer_pattern, df_rgb_reverse,...
                      sensor_map_rgb, bands_rgb,...
                      admm_options_f, reg_options_f,...
                      solvePatchesADMMOptions.patch_options,...
                      solvePatchesADMMVerbose...
                    );
                else
                    [...
                        I_rgb,...
                        ~,...
                        weights_images,...
                        extra_images{:}...
                    ] = solvePatchesADMM(...
                      [], I_raw_gt, bayer_pattern, df_rgb_reverse,...
                      sensor_map_rgb, bands_rgb,...
                      admm_options_f, reg_options_f,...
                      solvePatchesADMMOptions.patch_options,...
                      solvePatchesADMMVerbose...
                    );
                end
                
                time.admm(f, i, cr) = toc(time_start);

                saveImages(...
                    output_directory, name_params,...
                    I_rgb, 'rgb', 'I_rgb'...
                );
                if ~isempty(extra_images)
                    saveImages(...
                        output_directory, name_params,...
                        extra_images{3}, 'rgb_ab', 'I_rgb'...
                    );
                end
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
            
            if evaluate_aberrated_rgb
                if algorithm.spectral
                    I_rgb_ab = extra_images{1};
                else
                    I_rgb_ab = extra_images{3};
                end
                e_rgb_table_current = evaluateAndSaveRGB(...
                    I_rgb_ab, I_rgb_gt_warped, dp, names{i}, alg_name_params,...
                    fullfile(output_directory, [name_params, 'ab'])...
                );
                if ~isempty(e_rgb_table_ab)
                    e_rgb_table_ab = union(e_rgb_table_current, e_rgb_table_ab);
                else
                    e_rgb_table_ab = e_rgb_table_current;
                end
            end
            
            % Save the selected regularization weights
            if use_automatic_weights
                for w = 1:n_active_weights
                    aw = to_all_weights(w);
                    saveImages(...
                        'data', output_directory, name_params,...
                        weights_images(:, :, w), sprintf('weight%dImage', aw), 'I_weights'...
                    );

                    fg = figure;
                    imagesc(log10(weights_images(:, :, w)));
                    c = colorbar;
                    c.Label.String = sprintf('log_{10}(weight %d)', aw);
                    xlabel('Image x-coordinate')
                    ylabel('Image y-coordinate')
                    title(sprintf('Per-patch regularization weight %d', aw));
                    savefig(...
                        fg,...
                        fullfile(output_directory, [name_params  sprintf('weight%dImage.fig', aw)]),...
                        'compact'...
                        );
                    close(fg);
                end
            end
        end
    end
    
    % Demosaicking and colour channel warping
    W_forward = [];
    for f = 1:n_demosaic_algorithms
        algorithm = demosaic_algorithms.(demosaic_algorithm_fields{f});
        if ~algorithm.enabled
            continue;
        end
    
        time_start = tic;
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
        time.demosaic(f, i) = toc(time_start);
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
                time_start = tic;
                W_forward = dispersionfunToMatrix(...
                    df_rgb_forward, bands_rgb, image_sampling, image_sampling,...
                    [0, 0, image_sampling(2),  image_sampling(1)], false...
                    );
                time.warp(i) = toc(time_start);
                if verbose
                    fprintf('\t...done\n');
                end
            end
    
            time_start = tic;
            I_rgb = warpImage(I_rgb_warped, W_forward, image_sampling);
            time.warp_apply(i) = time.warp_apply(i) + toc(time_start);
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
    if has_dispersion_rgb
        time.warp_apply = time.warp_apply ./ sum(all(isfinite(time.demosaic), 2));
    end

    % Write evaluations to a file
    if ~isempty(e_rgb_table)
        writetable(...
            e_rgb_table,...
            fullfile(output_directory, [names{i}, '_evaluateRGB.csv'])...
        );
        e_rgb_tables{i} = e_rgb_table;
    end
    if evaluate_aberrated_rgb && ~isempty(e_rgb_table_ab)
        writetable(...
            e_rgb_table_ab,...
            fullfile(output_directory, [names{i}, '_evaluateRGB_ab.csv'])...
        );
        e_rgb_tables_ab{i} = e_rgb_table_ab;
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
    if evaluate_aberrated_spectral && ~isempty(e_spectral_table_ab)
        writetable(...
            e_spectral_table_ab,...
            fullfile(output_directory, [names{i}, '_evaluateSpectral_ab.csv'])...
        );
        % Also save completed figures
        name_ab = [names{i}, '_ab'];
        dp.evaluation.custom_spectral.(name_ab) = dp.evaluation.custom_spectral.(names{i});
        evaluateAndSaveSpectral(output_directory, dp, name_ab, all_alg_names, fg_spectral_ab);
        e_spectral_tables_ab{i} = e_spectral_table_ab;
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
if evaluate_aberrated_rgb
    e_rgb_tables_ab = e_rgb_tables_ab(~cellfun(@isempty, e_rgb_tables_ab, 'UniformOutput', true));
    if ~isempty(e_rgb_tables_ab)
        e_rgb_summary_table = mergeRGBTables(e_rgb_tables_ab);
        writetable(...
            e_rgb_summary_table,...
            fullfile(output_directory, [dataset_name, '_evaluateRGB_ab.csv'])...
        );
    end
end

e_spectral_tables = e_spectral_tables(~cellfun(@isempty, e_spectral_tables, 'UniformOutput', true));
if ~isempty(e_spectral_tables)
    e_spectral_summary_table = mergeSpectralTables(e_spectral_tables);
    writetable(...
        e_spectral_summary_table,...
        fullfile(output_directory, [dataset_name, '_evaluateSpectral.csv'])...
    );
end
if evaluate_aberrated_spectral
    e_spectral_tables_ab = e_spectral_tables_ab(~cellfun(@isempty, e_spectral_tables_ab, 'UniformOutput', true));
    if ~isempty(e_spectral_tables_ab)
        e_spectral_summary_table = mergeSpectralTables(e_spectral_tables_ab);
        writetable(...
            e_spectral_summary_table,...
            fullfile(output_directory, [dataset_name, '_evaluateSpectral_ab.csv'])...
        );
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
    'admm_algorithms', 'demosaic_algorithms', 'time'...
} ];
if has_spectral
    save_variables_list = [save_variables_list, {'bands_spectral', 'spectral_weights', 'color_weights_reference'}];
    if dp.spectral_reflectances
        save_variables_list = [save_variables_list, {'radiance_normalized_weights'}];
    end
end
if has_color_map
    save_variables_list = [save_variables_list, {'bands', 'bands_color', 'color_weights'}];
end
save_data_filename = fullfile(output_directory, ['RunOnDataset_' dataset_name '.mat']);
save(save_data_filename, save_variables_list{:});