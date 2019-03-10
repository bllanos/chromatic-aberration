%% Evaluate demosaicking, spectral reconstruction, and chromatic aberration correction
% For an image of an X-Rite ColorChecker CLASSIC colour chart.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### ColorChecker image
% 
% A raw colour filter array image of a ColorChecker chart must be provided. The
% image is expected to have been preprocessed, such as using
% 'PreprocessRAWImages.m', so that it does not need to be linearized after being
% loaded.
%
% Along with the captured image, two single-channel label images are required.
% The first is used for colour calibration, for vignetting correction, and for
% evaluation of colour-based image estimation algorithms. For the 24 patches of
% the chart, the corresponding pixels in the label image must have values equal
% to the patch indices. The uniformly-coloured portions of the frame surrounding
% the patches, used to calibrate vignetting, must be labelled 25. Lastly, the
% rest of the image should have a different label (e.g. zero).
%
% The second label image is used for evaluation of spectral-based image
% estimation algorithms. It should have all 24 patches labelled with their
% indices, but need not have other pixels labelled. The boundaries of the
% patches in the first label image should be accurate in the reference channel
% of colour-based models of chromatic aberration (e.g. the Green channel). In
% contrast, the boundaries of the patches in the second label image should be
% accurate in the spectral band used as the reference band for spectral models
% of chromatic aberration.
%
% ### ColorChecker spectral image
%
% An image of a ColorChecker chart captured under bandpass-filtered
% illumination, used for plotting a comparison line in figures of error around
% image edges. The image should be a "quasi hyperspectral image" output by
% 'RAWBandImagesToDataset.m'.
%
% ### Spectral reflectances and conversion to colour
%
% #### Spectral reflectances
% A '.csv' file containing a header row, a first column for wavelength
% values, and remaining columns for relative spectral reflectances of each
% patch in the ColorChecker.
%
% #### CIE tristimulus functions
% A '.mat' file containing a variable 'xyzbar', which can be used as the
% `C` input argument of 'cieSpectralToColor()'.
%
% ### Estimated images
%
% Two filepath wildcards in the script parameters below are used to locate the
% reconstructed RGB and spectral images to be evaluated. These images should
% have been generated using the raw image of the ColorChecker chart as input. In
% the output files describing the evaluation results, the names of the
% algorithms being evaluated will be the unique portions of the image filenames.
%
% ## Output
%
% ### Data file output
%
% #### Intermediate data and parameters
% A '.mat' file containing the following variables, as appropriate:
% - 'bands_estimated': A vector containing the wavelengths of the spectral
%   bands used in the estimated images.
% - 'bands_measured': A vector containing the wavelengths of the spectral
%   bands used by the reference reflectance data.
% - 'bands_qhyper': A vector containing the wavelengths of the spectral
%   bands in the reference bandpass-filtered illumination image.
% - 'vignetting_data': The model of vignetting describing the ColorChecker
%   image. 'vignetting_data' is the output argument of the same name of
%   'vignettingPolyfit()'. Refer to the documentation of 'vignettingPolyfit.m'
%   for details.
% - 'spectral_weights': A matrix for converting pixels in the spectral
%   space of the estimated spectral images to the spectral space of the
%   reference reflectance data.
% - 'spectral_weights_qhyper': A matrix for converting pixels in the spectral
%   space of the bandpass-filtered image to the spectral space of the reference
%   reflectance data.
% - 'algorithms': A structure containing the names of the algorithms which were
%   evaluated, extracted from the names of the input estimated images.
%   'algorithms' has the following fields:
%   - 'spectral': A cell vector of character vectors containing the names of the
%     algorithms which produced the estimated spectral images.
%   - 'color': A cell vector of character vectors containing the names of the
%     algorithms which produced the estimated colour images.
% - 'algorithm_filenames': A structure of the same form as 'algorithms'
%   containing the names and paths of the input estimated images.
% - 'color_correction': A 3 x 3 matrix for mapping from the raw RGB colours of
%   the ColorChecker image to the XYZ colours of the corresponding measured
%   reflectances.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`.
%
% The file is saved as 'EvaluateColorCheckerData.mat'.
%
% #### Evaluation results
%
% For each estimated ColorChecker patch, RGB error metrics and spectral error
% metrics, are output in the form of CSV files. Each CSV file contains results
% for all (RGB or spectral) estimated images. The estimated images will be
% listed in the files under the unique portions of their filenames. RGB error
% metrics are saved as '*_evaluateRGB.csv', whereas the spectral error metrics
% are saved as '*_evaluateSpectral.csv', where '*' contains the number of the
% patch.
%
% Error metrics are also aggregated across images, and saved as
% '*_evaluateRGB.csv' and '*_evaluateSpectral.csv', where '*' now does not
% contain the number of the patch.
%
% There are separate evaluation files for the central and edge regions of each
% patch.
%
% ## Detailed description of processing
%
% Ideal XYZ colours corresponding to the ColorChecker patches are created by
% integrating over the products of the measured spectral reflectances and the
% CIE tristimulus functions. These colours could be converted to RGB, if
% desired, using a whitepoint of [1, 1, 1] (for an equal energy radiator).
%
% A model of vignetting in the image of the ColorChecker is created from the
% variation in intensity of the ColorChecker's frame. A 3 x 3 colour correction
% matrix is then fit between the XYZ colours of the measured spectral
% reflectances and the vignetting-corrected raw image of the ColorChecker, and
% output for later use in converting the estimated images to the sRGB colour
% space, for example. The inverse of the vignetting model is applied to the
% measured spectral reflectances to obtain a "ground truth" image for spectral
% image evaluation.
%
% For RGB image evaluation, a ground truth image is created from the
% vignetting-corrected raw image of the ColorChecker as follows:
% - For each ColorChecker patch, the averages of the Red, Green, and Blue pixels
%   in a square around the centroid of the patch are computed.
% - These average values are then replicated to every pixel in the patch, and
%   the inverse of the vignetting model is applied to them to simulate an image
%   without spectral dispersion at the edges of the ColorChecker patches.
%
% For each ColorChecker patch, an independent evaluation is performed for all
% reconstructed images as follows:
% - Spectral images:
%   - For a spectral image, as the illuminant is unknown, the average spectral
%     signal is calculated within a reference patch. The image is then
%     multiplied by the ratio of the corresponding average signal in the true
%     image to this average signal.
%   - The centroid of the patch is located, and a square is cropped from both
%     the true and reconstructed spectral images. 'evaluateAndSaveSpectral()' is
%     called on these patches to evaluate spectral reconstruction of the
%     homogenous area.
%   - The border region of the patch is located, and its pixels are reshaped
%     into 1D images from the true and reconstructed spectral images.
%     'evaluateAndSaveSpectral()' is called on these patches to evaluate
%     spectral reconstruction of the border region.
% - RGB images
%   - The above two evaluations are performed for an RGB image, using
%     'evaluateAndSaveRGB()'. In this case, the reference image is created from
%     the raw image of the ColorChecker as described above, and so there is no
%     need to rescale one image's channels to match the other's.
%
% The evaluation results for all patches are merged into summaries.
%
% Note that the extraction of patch centre and border regions uses the
% appropriate label image for the reconstructed image being evaluated (spectral
% vs. RGB).
%
% ## References
%
% The colour correction matrix is calculated based on the description in:
%
%   Karaimer, Hakki C., & Brown, Michael S. (2018). "Improving Color
%   Reproduction Accuracy on Cameras," In IEEE Conference on Computer Vision and
%   Pattern Recognition (CVPR) (pp. 6440–6449).

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 7, 2019

% List of parameters to save with results
parameters_list = {
    'color_label_filename',...
    'spectral_label_filename',...
    'reference_filename',...
    'reference_variable_name',...
    'qhyper_filename',...
    'qhyper_variable_name',...
    'qhyper_bands_filename',...
    'qhyper_bands_variable',...
    'xyzbar_filename',...
    'reflectances_filename',...
    'reflectance_data_patch_columns',...
    'wavelength_range',...
    'max_degree_vignetting',...
    'n_patches',...
    'reference_patch_index',...
    'centroid_patch_width',...
    'edge_width',...
    'true_image_name',...
    'true_image_filename',...
    'bands_filename',...
    'bands_variable'...
    };

%% Input data and parameters

% Label image for colour-based dispersion correction and for vignetting
% calibration (an image file, not a '.mat' file)
color_label_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/colorChecker_preprocessed/unfiltered/d2_colorChecker30cm_unfiltered_background0_patches1to24_frame25.png';

% Label image for spectral dispersion correction
spectral_label_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/colorChecker_preprocessed/600nm/d2_colorChecker30cm_600nm_background0_patches1to24_frame25.png';

% Raw image of the ColorChecker taken under unfiltered light ('.mat' or image
% files can be loaded)
reference_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/exposure_blending/d2_colorChecker30cm_unfiltered.mat';
reference_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Image of the ColorChecker taken under bandpass-filtered light
qhyper_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/channel_scaling/d2_colorChecker30cm_qHyper.mat';
qhyper_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% Path and filename of a '.mat' file containing the wavelengths corresponding to
% the bandpass-filtered image
qhyper_bands_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/channel_scaling/sensor.mat';
qhyper_bands_variable = 'bands'; % Variable name in the above file

% CIE tristimulus functions
xyzbar_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180614_ASTM_E308/Table1_CIE1931_2DegStandardObserver.csv';

% Sample spectral reflectances
reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180626_SpectralCharacterizationOfSetup/spectra_averaged.csv';

% Categorization of the samples in the file of reflectances
reflectance_data_patch_columns = 13:36;

% Wavelength range to truncate the measured reflectances to
wavelength_range = [350, 750];

% Maximum degree of the polynomial model of vignetting
max_degree_vignetting = 5;

% Number of ColorChecker patches
n_patches = 24;

% The reference patch for spectral evaluation
reference_patch_index = 20;

% Odd-integer size of the regions to extract around the centroids of the
% ColorChecker patches
centroid_patch_width = 15;

% Width of the edge region to extract inside, and outside, the patch boundaries
edge_width = 10;

true_image_name = 'GT'; % Name used in figures. Must not contain spaces.
true_image_filename = 'colorChecker'; % Name used in filenames. Must not contain spaces.

% Wildcard for 'ls()' to find the estimated spectral images to process.
% '.mat' or image files can be loaded
spectral_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_MATFiles/d2_colorChecker30cm_*_latent.mat';
spectral_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Path and filename of a '.mat' file containing the wavelengths corresponding to
% the estimated spectral images
bands_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_dispersion/RunOnDataset_20190208_ComputarLens_rawCaptured_dispersion.mat';
bands_variable = 'bands'; % Variable name in the above file

% Wildcard for 'ls()' to find the estimated colour images to process.
% '.mat' or image files can be loaded
color_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_MATFiles/d2_colorChecker30cm_*_rgb.mat';
color_variable_name = 'I_rgb'; % Used only when loading '.mat' files

% Output directory
output_directory = '/home/llanos/Downloads/colorChecker_evaluation';

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% ## Debugging Flags
vignettingPolyfitVerbose = true;

%% Load wavelengths first, to avoid accidentally overwriting a variable

load(bands_filename, bands_variable);
if exist(bands_variable, 'var')
    bands_estimated = eval(bands_variable);
end
if ~exist(bands_variable, 'var') || isempty(bands_estimated)
    error('No wavelengths loaded.')
end
n_bands_estimated = length(bands_estimated);

load(qhyper_bands_filename, qhyper_bands_variable);
if exist(qhyper_bands_variable, 'var')
    bands_qhyper = eval(qhyper_bands_variable);
end
if ~exist(qhyper_bands_variable, 'var') || isempty(bands_qhyper)
    error('No wavelengths loaded.')
end
n_bands_qhyper = length(bands_qhyper);

%% Calibrate vignetting

I_reference = loadImage(reference_filename, reference_variable_name);
if size(I_reference, 3) ~= 1
    error('Expected the ColorChecker image, "%s", to have only one channel.', reference_filename);
end

I_color_label = imread(color_label_filename);
if size(I_color_label, 3) ~= 1
    error('Expected the colour-based dispersion label image, "%s", to have only one channel.', color_label_filename);
end
mask_vignetting = I_color_label == (n_patches + 1);

[vignettingfun, vignetting_data] = vignettingPolyfit(...
    I_reference, mask_vignetting, max_degree_vignetting, bayer_pattern, vignettingPolyfitVerbose...
);

% Make the vignetting relative to the center of the image, which is valid if we
% assume that the center of the image is within the domain of samples used to
% fit the vignetting model
image_sampling = size(I_reference);
vignetting_correction_factor = vignettingfun(fliplr(image_sampling) / 2);

%% Load reflectance data

sample_table = readtable(reflectances_filename);
variable_names = sample_table.Properties.VariableNames;
bands_measured = sample_table.(variable_names{1});
bands_filter = (bands_measured >= wavelength_range(1)) & (bands_measured <= wavelength_range(2));
bands_measured = bands_measured(bands_filter);
n_bands_measured = length(bands_measured);
reflectances = sample_table{bands_filter, reflectance_data_patch_columns};

% Find sample colours
xyzbar_table = readtable(xyzbar_filename);
lambda_xyzbar = xyzbar_table{:, 1};
xyzbar = xyzbar_table{:, 2:end};

[reflectances_rgb, reflectances_xyz] = reflectanceToColor(...
    bands_measured, ones(size(bands_measured)),... % Equal energy radiator
    bands_measured, reflectances,...
    lambda_xyzbar, xyzbar...
    );
reflectances_rgb_integer = floor(256 * reflectances_rgb);

figure;
hold on
names_legend = cell(n_patches, 1);
disp('Measured reflectance sRGB colours under an equal energy illuminant:');
for pc = 1:n_patches
    plot(...
        bands_measured, reflectances(:, pc),...
        'Color', reflectances_rgb(pc, :), 'LineWidth', 2, 'Marker', 'none'...
    );
    % Recover original variable names, which contained spaces
    names_legend{pc} = strsplit(sample_table.Properties.VariableDescriptions{reflectance_data_patch_columns(pc)}, ':');
    names_legend{pc} = names_legend{pc}{end};
    if isempty(names_legend{pc})
        names_legend{pc} = sample_table.Properties.VariableNames{reflectance_data_patch_columns(pc)};
    end

    fprintf(...
        '\t%s: %d, %d, %d\n', names_legend{pc},...
        reflectances_rgb_integer(pc, 1),...
        reflectances_rgb_integer(pc, 2),...
        reflectances_rgb_integer(pc, 3)...
    );
end
hold off
title('Measured reflectance spectral signals')
xlabel('\lambda [nm]')
ylabel('Relative spectral signal')
legend(names_legend);
ax = gca;
ax.Color = [0.5 0.5 0.5];

spectral_weights = resamplingWeights(...
    bands_measured, bands_estimated, findSamplingOptions.interpolant, findSamplingOptions.bands_padding...
);
spectral_weights_evaluation = eye(n_bands_measured);

spectral_weights_qhyper = resamplingWeights(...
    bands_measured, bands_qhyper, findSamplingOptions.interpolant, findSamplingOptions.bands_padding...
);

%% Prepare for image evaluation

% Images to be evaluated, and evaluation parameters
algorithm_filenames.spectral = listFiles(spectral_wildcard);
n_spectral_algorithms = length(algorithm_filenames.spectral);

algorithm_filenames.color = listFiles(color_wildcard);
algorithm_names = trimCommon([algorithm_filenames.spectral; algorithm_filenames.color]);
algorithms.spectral = algorithm_names(1:n_spectral_algorithms);
algorithms.color = algorithm_names((n_spectral_algorithms + 1):end);
algorithms_escaped.spectral = strrep(algorithms.spectral, '_', '\_');
algorithms_escaped.color = strrep(algorithms.color, '_', '\_');
n_color_algorithms = length(algorithms.color);

I_spectral_label = imread(spectral_label_filename);
if size(I_spectral_label, 3) ~= 1
    error('Expected the spectral dispersion label image, "%s", to have only one channel.', spectral_label_filename);
end

e_centroid_tables_rgb = cell(n_patches, 1);
e_edge_tables_rgb = cell(n_patches, 1);
e_centroid_tables_spectral = cell(n_patches, 1);
e_edge_tables_spectral = cell(n_patches, 1);

evaluation_plot_colors_spectral = jet(n_spectral_algorithms);
evaluation_plot_colors_color = jet(n_color_algorithms);
evaluation_plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
evaluation_plot_styles = {'--', ':', '-.'};

% Common variables for patch processing
n_channels_rgb = 3;
if mod(centroid_patch_width, 2) == 0
    error('`centroid_patch_width` should be an odd integer.');
end
half_width = floor(centroid_patch_width / 2);
se_clean = strel('square', 3); % Structuring element for cleaning up the label image
se_edge = strel('disk', edge_width);
channel_mask = bayerMask(image_sampling(1), image_sampling(2), bayer_pattern);

% Enumerate the positions of all pixels
[X, Y] = meshgrid(1:image_sampling(2), 1:image_sampling(1));
X = X - 0.5; % Place coordinates at pixel centres
Y = Y - 0.5;

% For calculating a colour correction matrix
patch_means_rgb = zeros(n_patches, n_channels_rgb);

%% Evaluate each patch of the ColorChecker

% For spectral evaluation, first use the reference patch to calibrate the
% intensity alignment between the reference and estimated images
is_reference_patch = true;
for pc = [reference_patch_index, 1:n_patches]
    for color_flag = [true, false]
        if color_flag
            if is_reference_patch
                continue
            end
            I_label = I_color_label;
        else
            I_label = I_spectral_label;
        end
        
        patch_mask = imopen(I_label == pc, se_clean);

        % Find the centroid of the patch
        cc_stats = regionprops(patch_mask, 'Centroid');
        if numel(cc_stats) ~= 1
            if color_flag
                error('Expected only one binary mask region for patch %d in the colour label image.', pc);
            else
                error('Expected only one binary mask region for patch %d in the spectral label image.', pc);
            end
        end
        patch_centroid = floor(cc_stats.Centroid + 0.5);

        % Define a region around the centroid
        roi = [
            patch_centroid(2) - half_width, patch_centroid(2) + half_width,...
            patch_centroid(1) - half_width, patch_centroid(1) + half_width
        ];
        if is_reference_patch
            roi_reference = roi;
        end
        patch_xy = [
            reshape(X(roi(1):roi(2), roi(3):roi(4)), [], 1),...
            reshape(Y(roi(1):roi(2), roi(3):roi(4)), [], 1)
        ];
        factors = reshape(...
            vignettingfun(patch_xy) / vignetting_correction_factor,...
            centroid_patch_width, centroid_patch_width...
        );
    
        bayer_pattern_centroid = offsetBayerPattern(roi([1, 3]), bayer_pattern);
        channel_mask_centroid = bayerMask(centroid_patch_width, centroid_patch_width, bayer_pattern_centroid);
    
        if color_flag
            I_patch = I_reference(roi(1):roi(2), roi(3):roi(4));
            I_patch = I_patch ./ factors;

            % Find the mean colour in this region, corrected for vignetting
            patch_mean = zeros(1, 1, n_channels_rgb);
            for c = 1:n_channels_rgb
                patch_mean(c) = mean(I_patch(channel_mask_centroid(:, :, c)));
            end
            patch_means_rgb(pc, :) = patch_mean;
        else
            patch_mean = reshape(reflectances(:, pc), 1, 1, []);
        end
        
        I_patch = repmat(patch_mean, centroid_patch_width, centroid_patch_width, 1);
        I_patch = I_patch .* repmat(factors, 1, 1, numel(patch_mean));
        if is_reference_patch
            patch_mean_reference = mean(mean(I_patch, 1), 2);
            is_reference_patch = false;
            continue
        end
        
        % Find the edge of the patch
        
        % Outside the edge
        I_distance = -bwdist(patch_mask);
        edge_mask = (I_distance < 0 & I_distance >= -edge_width);
        edge_ind_outer = find(edge_mask);
        edge_distance = I_distance(edge_ind_outer);
        
        % Inside the edge
        I_distance = bwdist(~patch_mask);
        edge_mask = (I_distance > 0 & I_distance <= edge_width);
        edge_ind_inner = find(edge_mask);
        edge_ind = [edge_ind_outer; edge_ind_inner];
        edge_distance = [edge_distance; I_distance(edge_ind_inner)]; %#ok<AGROW>
        
        % Synthesize the ideal edge region
        patch_xy = [X(edge_ind_inner), Y(edge_ind_inner)];
        factors = vignettingfun(patch_xy) / vignetting_correction_factor;
        n_edge_px_inner = length(edge_ind_inner);
        I_edge_inner = repmat(patch_mean, n_edge_px_inner, 1, 1);
        I_edge_inner = I_edge_inner .* repmat(factors, 1, 1, numel(patch_mean));
        
        patch_xy = [X(edge_ind), Y(edge_ind)];
        factors = vignettingfun(patch_xy) / vignetting_correction_factor;
        n_edge_px = length(edge_ind);
        I_edge = repmat(patch_mean, n_edge_px, 1, 1);
        I_edge = I_edge .* repmat(factors, 1, 1, numel(patch_mean));
        
        % To smooth plots, points will be averaged for each distance from the
        % edge
        [edge_distance_unique, ~, edge_distance_ind] = unique(edge_distance);
        n_edge_distances = length(edge_distance_unique);
        edge_distance_filters = false(n_edge_px, n_edge_distances);
        for k = 1:n_edge_distances
            edge_distance_filters(:, k) = (edge_distance_ind == k);
        end
        
        % Evaluate the ColorChecker image with respect to edge error
        fg_edge = [];
        err_averaged = zeros(n_edge_distances, 1);
        if color_flag
            err = abs(I_edge - repmat(I_reference(edge_ind), 1, 1, n_channels_rgb)) ./ abs(I_edge);
            for c = 1:n_channels_rgb
                mask_c = channel_mask(:, :, c);
                filter_c = mask_c(edge_ind);
                for k = 1:n_edge_distances
                    filter_c_edge = filter_c & edge_distance_filters(:, k);
                    err_averaged(k) = mean(squeeze(err(filter_c_edge, :, c)));
                end
                err_filter = isfinite(err_averaged);
                fg_edge(c) = figure;
                plot(edge_distance_unique(err_filter), err_averaged(err_filter), 'Color', [0, 0, 0], 'LineWidth', 2);
            end
        else
            I = loadImage(qhyper_filename, qhyper_variable_name);
            I_roi_reference = reshape(I(roi_reference(1):roi_reference(2), roi_reference(3):roi_reference(4), :), [], n_bands_qhyper);
            I_edge_qhyper = zeros(n_edge_px, 1, n_bands_qhyper);
            for c = 1:n_bands_qhyper
                I_c = I(:, :, c);
                I_edge_qhyper(:, :, c) = I_c(edge_ind);
            end
            I_edge_qhyper = channelConversion(I_edge_qhyper, spectral_weights_qhyper);
            
            % This image has a Bayer pattern
            err_averaged_c = cell(n_channels_rgb, n_edge_distances);
            for c = 1:n_channels_rgb
                I_roi_reference_c = I_roi_reference(reshape(channel_mask_centroid(:, :, c), [], 1), :);
                I_mean = reshape(mean(I_roi_reference_c, 1), 1, 1, n_bands_qhyper);
                I_mean = channelConversion(I_mean, spectral_weights_qhyper);
                radiance_ratio = patch_mean_reference ./ I_mean;
                % Prevent division by zero, and indefinite values
                radiance_ratio(I_mean == patch_mean_reference) = 1;
                radiance_ratio(~isfinite(radiance_ratio)) = 1;
                
                mask_c = channel_mask(:, :, c);
                filter_c = mask_c(edge_ind);
                I_edge_qhyper_c = I_edge_qhyper(filter_c, :, :);
                n_edge_px_c = size(I_edge_qhyper_c, 1);
                I_edge_qhyper_c = I_edge_qhyper_c .* repmat(radiance_ratio, n_edge_px_c, 1, 1);
                err = mean(abs(I_edge(filter_c, :, :) -  I_edge_qhyper_c) ./ I_edge(filter_c, :, :) , 3);
                for k = 1:n_edge_distances
                    err_averaged_c{c, k} = squeeze(err(edge_distance_filters(filter_c, k)));
                end
            end
            for k = 1:n_edge_distances
                err_averaged(k) = mean(cell2mat(err_averaged_c(:, k)));
            end
            err_filter = isfinite(err_averaged);
            
            fg_edge = figure;
            plot(...
                edge_distance_unique(err_filter), err_averaged(err_filter),...
                'Color', [0, 0, 0], 'LineWidth', 2 ...
            );     
        end

        % Evaluation
        for evaluating_center = [true, false]
            if evaluating_center
                I_reference_eval = I_patch;
                name_params = fullfile(output_directory, [true_image_filename, sprintf('_patch%dCenter', pc)]);
            else
                I_reference_eval = I_edge_inner;
                name_params = fullfile(output_directory, [true_image_filename, sprintf('_patch%dEdge', pc)]);
            end
            fg_spectral = struct;
            e_table = [];
            if color_flag
                dp.evaluation.global_rgb = struct;
                dp.evaluation.custom_rgb = struct;
            else
                evaluation_options = struct(...
                    'metric', 'mrae',...
                    'mi_bands', [1, n_bands_measured]...
                );
                if evaluating_center
                    evaluation_options_radiance = [...
                        half_width + 1, half_width + 1, centroid_patch_width, centroid_patch_width
                    ];
                else
                    center_edge_ind = floor(n_edge_px_inner / 2);
                    evaluation_options_radiance = [
                        1, center_edge_ind,...
                        1, 2 * center_edge_ind - 1
                    ];
                end
                true_image_filename_eval = [true_image_filename, sprintf('_patch%d', pc)];
                dp.evaluation.global_spectral = evaluation_options;
                dp.evaluation.custom_spectral.(true_image_name).radiance = evaluation_options_radiance;
                dp.evaluation.custom_spectral.(true_image_filename_eval).radiance = evaluation_options_radiance;
            end

            if color_flag
                for i = 1:n_color_algorithms
                    I = loadImage(algorithm_filenames.color{i}, color_variable_name);
                    if evaluating_center
                        I = I(roi(1):roi(2), roi(3):roi(4), :);
                    else
                        I_edge_estimated_inner = zeros(n_edge_px_inner, 1, n_channels_rgb);
                        for c = 1:n_channels_rgb
                            I_c = I(:, :, c);
                            I_edge_estimated_inner(:, 1, c) = I_c(edge_ind_inner);
                        end
                        for c = 1:n_channels_rgb
                            I_c = I(:, :, c);
                            I_edge_estimated_c = I_c(edge_ind);
                            figure(fg_edge(c))
                            hold on
                            err = abs(I_edge(:, :, c) -  I_edge_estimated_c) ./ abs(I_edge(:, :, c));
                            for k = 1:n_edge_distances
                                err_averaged(k) = mean(squeeze(err(edge_distance_filters(:, k))));
                            end
                            err_filter = isfinite(err_averaged);
                            plot(...
                                edge_distance_unique(err_filter), err_averaged(err_filter),...
                                'Color', evaluation_plot_colors_color(i, :),...
                                'LineWidth', 2,...
                                'Marker', evaluation_plot_markers{mod(i - 1, length(evaluation_plot_markers)) + 1},...
                                'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
                                );
                            hold off
                        end
                        I = I_edge_estimated_inner;
                    end
                    name_params_i = [name_params, '_', algorithms.color{i}];
                    e_table_current = evaluateAndSaveRGB(...
                        I, I_reference_eval, dp, true_image_name, algorithms_escaped.color{i},...
                        name_params_i...
                    );
                    if ~isempty(e_table)
                        e_table = union(e_table_current, e_table);
                    else
                        e_table = e_table_current;
                    end
                end
                writetable(e_table, [name_params, '_evaluateRGB.csv']);
                
                if ~evaluating_center
                    for c = 1:n_channels_rgb
                        figure(fg_edge(c))
                        legend(['Raw input'; algorithms_escaped.color]);
                        title(sprintf('Absolute error around patch %d edge in channel %d', pc, c));
                        xlabel('Position relative to patch borders [pixels]')
                        ylabel('Relative error')
                        savefig(...
                            fg_edge(c),...
                            [ name_params, sprintf('_errChannel%d.fig', c)], 'compact'...
                        );
                        close(fg_edge(c));
                    end
                end
                
            else
                for i = 1:n_spectral_algorithms
                    I = loadImage(algorithm_filenames.spectral{i}, spectral_variable_name);
                    name_params_i = [name_params, '_', algorithms.spectral{i}];
                    
                    I_mean = mean(mean(...
                        I(roi_reference(1):roi_reference(2), roi_reference(3):roi_reference(4), :),...
                        1), 2 ...
                    );
                    I_mean = channelConversion(I_mean, spectral_weights);
                    radiance_ratio = patch_mean_reference ./ I_mean;
                    % Prevent division by zero, and indefinite values
                    radiance_ratio(I_mean == patch_mean_reference) = 1;
                    radiance_ratio(~isfinite(radiance_ratio)) = 1;
                    
                    if evaluating_center
                        I = channelConversion(I(roi(1):roi(2), roi(3):roi(4), :), spectral_weights)...
                            .* repmat(radiance_ratio, centroid_patch_width, centroid_patch_width, 1);
                    else
                        I_edge_estimated_inner = zeros(n_edge_px_inner, 1, n_bands_estimated);
                        for c = 1:n_bands_estimated
                            I_c = I(:, :, c);
                            I_edge_estimated_inner(:, 1, c) = I_c(edge_ind_inner);
                        end
                        
                        I_edge_estimated = zeros(n_edge_px, 1, n_bands_estimated);
                        for c = 1:n_bands_estimated
                            I_c = I(:, :, c);
                            I_edge_estimated(:, :, c) = I_c(edge_ind);
                        end
                        I_edge_estimated = channelConversion(I_edge_estimated, spectral_weights)...
                            .* repmat(radiance_ratio, n_edge_px, 1, 1);
                        I = channelConversion(I_edge_estimated_inner, spectral_weights)...
                            .* repmat(radiance_ratio, n_edge_px_inner, 1, 1);
                        figure(fg_edge)
                        hold on
                        err = mean(abs(I_edge -  I_edge_estimated) ./ abs(I_edge), 3);
                        for k = 1:n_edge_distances
                            err_averaged(k) = mean(squeeze(err(edge_distance_filters(:, k))));
                        end
                        err_filter = isfinite(err_averaged);
                        plot(...
                            edge_distance_unique(err_filter), err_averaged(err_filter),...
                            'Color', evaluation_plot_colors_spectral(i, :),...
                            'LineWidth', 2,...
                            'Marker', evaluation_plot_markers{mod(i - 1, length(evaluation_plot_markers)) + 1},...
                            'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
                            );
                        hold off
                    end

                    dp.evaluation.global_spectral.plot_color = evaluation_plot_colors_spectral(i, :);
                    dp.evaluation.global_spectral.plot_marker = 'none';
                    dp.evaluation.global_spectral.plot_style =...
                        evaluation_plot_styles{...
                        mod(i - 1, length(evaluation_plot_styles)) + 1 ...
                        };  
                    [e_table_current, fg_spectral] = evaluateAndSaveSpectral(...
                        I, I_reference_eval, bands_measured, spectral_weights_evaluation,...
                        dp, true_image_name, algorithms_escaped.spectral{i},...
                        name_params_i, fg_spectral...
                        );
                    if ~isempty(e_table)
                        e_table = union(e_table_current, e_table);
                    else
                        e_table = e_table_current;
                    end
                end
                writetable(e_table, [name_params, '_evaluateSpectral.csv']);
                % Also save completed figures
                evaluateAndSaveSpectral(output_directory, dp, true_image_filename_eval, algorithms_escaped.spectral, fg_spectral);
                
                if ~evaluating_center
                    figure(fg_edge)
                    legend(['Bandpass-filtered image'; algorithms_escaped.spectral]);
                    title(sprintf('Relative error around patch %d edge', pc));
                    xlabel('Position relative to patch borders [pixels]')
                    ylabel('Relative error (MRAE)')
                    savefig(...
                        fg_edge,...
                        [ name_params, '_errSpectral.fig'], 'compact'...
                    );
                    close(fg_edge);
                end
            end
            
            if color_flag
                if evaluating_center
                    e_centroid_tables_rgb{pc} = e_table;
                else
                    e_edge_tables_rgb{pc} = e_table;
                end
            else
                if evaluating_center
                    e_centroid_tables_spectral{pc} = e_table;
                else
                    e_edge_tables_spectral{pc} = e_table;
                end
            end
        end
    end
end

%% Create a colour correction matrix

color_correction = (reflectances_xyz.') / (patch_means_rgb.');

%% Save results for all patches

e_rgb_summary_table = mergeRGBTables(e_centroid_tables_rgb);
writetable(...
    e_rgb_summary_table,...
    fullfile(output_directory, [true_image_filename, '_centers_evaluateRGB.csv'])...
);
e_rgb_summary_table = mergeRGBTables(e_edge_tables_rgb);
writetable(...
    e_rgb_summary_table,...
    fullfile(output_directory, [true_image_filename, '_edges_evaluateRGB.csv'])...
);

e_spectral_summary_table = mergeSpectralTables(e_centroid_tables_spectral);
writetable(...
    e_spectral_summary_table,...
    fullfile(output_directory, [true_image_filename, '_centers_evaluateSpectral.csv'])...
);
e_spectral_summary_table = mergeSpectralTables(e_edge_tables_spectral);
writetable(...
    e_spectral_summary_table,...
    fullfile(output_directory, [true_image_filename, '_edges_evaluateSpectral.csv'])...
);

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {
    'bands_estimated', 'bands_measured', 'bands_qhyper',...
    'vignetting_data', 'spectral_weights', 'spectral_weights_qhyper',...
    'algorithms', 'algorithm_filenames', 'color_correction'...
} ];
save_data_filename = fullfile(output_directory, 'EvaluateColorCheckerData.mat');
save(save_data_filename, save_variables_list{:});