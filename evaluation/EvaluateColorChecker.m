%% Evaluate demosaicking, spectral reconstruction, and chromatic aberration correction
% For an image of an X-Rite ColorChecker CLASSIC colour chart.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### ColorChecker image and labels
% 
% A raw colour filter array image of a ColorChecker chart must be provided. The
% image is expected to have been preprocessed, such as using
% 'PreprocessRAWImages.m', so that it does not need to be linearized after being
% loaded.
%
% Along with the captured image, two single-channel label images are required.
% The first is used for vignetting correction, and for evaluation of colour
% images. For the 24 patches of the chart, the corresponding pixels in the label
% image must have values equal to the patch indices. The uniformly-coloured
% portions of the frame surrounding the patches, used to calibrate vignetting,
% must be labelled 25. Lastly, the rest of the image should have a different
% label (e.g. zero).
%
% The second label image is used for evaluation of spectral images. It should
% have all 24 patches labelled with their indices, but need not have other
% pixels labelled. The boundaries of the patches in the first label image should
% be accurate in the reference channel of colour-based models of chromatic
% aberration (e.g. the Green channel). In contrast, the boundaries of the
% patches in the second label image should be accurate in the spectral band used
% as the reference band for spectral models of chromatic aberration.
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
% #### Colour space conversion data
% A '.mat' file containing several variables, which is the output of
% 'SonyColorMap.m', for example. The following variables are required:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel of the output sensor response images to the
%   j-th element of 'bands' (below).
% - 'channel_mode': A Boolean value indicating whether the input colour
%   space is a set of colour channels (true) or a set of spectral bands
%   (false). A value of `false` is required.
% - 'bands': A vector containing the wavelengths corresponding to the
%   second dimension of 'sensor_map'.
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
% - 'vignetting_correction_factor': A scalar used to map vignetting-corrected
%   intensities back to the scale of the original image intensities.
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
% for all estimated images. The estimated images will be listed in the files
% under the unique portions of their filenames. RGB error metrics are saved as
% '*_evaluateRGB.csv', whereas the spectral error metrics are saved as
% '*_evaluateSpectral.csv', where '*' contains the number of the patch.
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
% A model of vignetting in the image of the ColorChecker is created from the
% variation in intensity of the ColorChecker's frame. The inverse of the
% vignetting model is applied to the measured spectral reflectances to obtain a
% "ground truth" image for spectral image evaluation.
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
%   - For a spectral image, as the illuminant is unknown, the spectral pixels in
%     the patches of the true and test images are globally registered by finding
%     best-fit scaling factors for each spectral band. This alignment process is
%     the same as the 'global' alignment process used in
%     PointSpectralEvaluation.m. Note that patches are weighted inversely by
%     their pixel areas in the registration process.
%   - The centroid of the patch is located, and a square is cropped from both
%     the true and reconstructed spectral images. 'evaluateAndSaveSpectral()' is
%     called on these patches to evaluate spectral reconstruction of the
%     homogenous area.
%   - The border region of the patch is located, and its pixels are reshaped
%     into 1D images from the true and reconstructed spectral images.
%     'evaluateAndSaveSpectral()' is called on these patches to evaluate
%     spectral reconstruction of the border region.
%   - The spectral image (without spectral registration to the true spectral
%     image) is converted to colour, and also evaluated as an RGB image (see
%     below).
% - RGB images
%   - The above two evaluations are performed for an RGB image, using
%     'evaluateAndSaveRGB()'. In this case, the reference image is created from
%     the raw image of the ColorChecker as described above, and so there is no
%     need to rescale one image's channels to match the other's.
%
% The evaluation results for all patches are merged into summaries.
%
% Note that the extraction of patch centre and border regions uses the
% appropriate label image for the reconstructed image being evaluated, because
% the zero-dispersion locations differ between the spectral and colour models of
% chromatic aberration.

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
    'reflectances_filename',...
    'reflectance_data_patch_columns',...
    'wavelength_range',...
    'max_degree_vignetting',...
    'n_patches',...
    'centroid_patch_width',...
    'edge_width',...
    'smoothing_method',...
    'smoothing_span',...
    'true_image_name',...
    'true_image_filename',...
    'bands_filename',...
    'bands_variable',...
    'color_map_filename'...
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
qhyper_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/d2_colorChecker30cm_qHyper.mat';
qhyper_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% Path and filename of a '.mat' file containing the wavelengths corresponding to
% the bandpass-filtered image
qhyper_bands_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/sensor.mat';
qhyper_bands_variable = 'bands'; % Variable name in the above file

% Sample spectral reflectances
reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180626_SpectralCharacterizationOfSetup/spectra_averaged.csv';

% Categorization of the samples in the file of reflectances
reflectance_data_patch_columns = 13:36;

% Wavelength range to truncate the measured reflectances to
wavelength_range = [375, 725];

% Maximum degree of the polynomial model of vignetting
max_degree_vignetting = 5;

% Number of ColorChecker patches
n_patches = 24;

% Odd-integer size of the regions to extract around the centroids of the
% ColorChecker patches
centroid_patch_width = 15;

% Width of the edge region to extract inside, and outside, the patch boundaries
edge_width = 10;

% Method to use for smoothing curves of error around edges (Refer to the
% documentation of the Curve Fitting Toolbox's 'smooth' function)
smoothing_method = 'moving';

% Span to use for smoothing curves of error around edges (Refer to the
% documentation of the Curve Fitting Toolbox's 'smooth' function)
smoothing_span = 5;

true_image_name = 'GT'; % Name used in figures. Must not contain spaces.
true_image_filename = 'colorChecker'; % Name used in filenames. Must not contain spaces.

% Wildcard for 'ls()' to find the estimated spectral images to process (can be
% empty). '.mat' or image files can be loaded.
spectral_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/d2_colorChecker30cm_*_latent.mat';
spectral_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Path and filename of a '.mat' file containing the wavelengths corresponding to
% the estimated spectral images (`bands`), as well as the spectral resampling
% parameters governing their conversion to other spectral sampling spaces
% (`findSamplingOptions`).
sampling_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/RunOnDataset_20190208_ComputarLens_rawCaptured_dispersion.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/SonyColorMapData.mat';

% Wildcard for 'ls()' to find the estimated colour images to process (can be
% empty). '.mat' or image files can be loaded.
color_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dispersion_ignoreDispersionWeights/d2_colorChecker30cm_*_rgb.mat';
color_variable_name = 'I_rgb'; % Used only when loading '.mat' files

% Output directory
output_directory = '/home/llanos/Downloads/colorCheckerEvaluation';

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% ## Debugging Flags
vignettingPolyfitVerbose = true;

%% Load wavelengths and spectral to colour conversion information

has_spectral = ~isempty(spectral_wildcard);
has_color = ~isempty(color_wildcard);

bands_estimated = [];
bands_qhyper = [];
if has_spectral
    [findSamplingOptions, bands_estimated] = loadVariables(sampling_filename, {'findSamplingOptions', 'bands'});
    n_bands_estimated = length(bands_estimated);

    bands_qhyper = loadVariables(qhyper_bands_filename, qhyper_bands_variable);
    n_bands_qhyper = length(bands_qhyper);

    [sensor_map, ~, bands_color] = loadColorMap(color_map_filename, false);
    color_weights = colorWeights(...
        sensor_map, bands_color, bands_estimated, findSamplingOptions...
    );
end

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

spectral_weights = [];
spectral_weights_evaluation = [];
spectral_weights_qhyper = [];
if has_spectral
    spectral_weights = resamplingWeights(...
        bands_measured, bands_estimated, findSamplingOptions.interpolant, findSamplingOptions.bands_padding...
    );
    spectral_weights_evaluation = eye(n_bands_measured);

    spectral_weights_qhyper = resamplingWeights(...
        bands_measured, bands_qhyper, findSamplingOptions.interpolant_ref, findSamplingOptions.bands_padding...
    );
end

%% Prepare for image evaluation

% Images to be evaluated, and evaluation parameters
algorithm_filenames.spectral = {};
n_spectral_algorithms = 0;
if has_spectral
    algorithm_filenames.spectral = listFiles(spectral_wildcard);
    n_spectral_algorithms = length(algorithm_filenames.spectral);
end

algorithm_filenames.color = {};
if has_color
    algorithm_filenames.color = listFiles(color_wildcard);
end
if has_color || has_spectral
    algorithm_names = trimCommon([algorithm_filenames.spectral; algorithm_filenames.color]);
else
    algorithm_names = {};
end
algorithms.spectral = algorithm_names(1:n_spectral_algorithms);
algorithms.color = algorithm_names((n_spectral_algorithms + 1):end);
algorithms_escaped.spectral = strrep(algorithms.spectral, '_', '\_');
algorithms_escaped.color = strrep(algorithms.color, '_', '\_');
n_color_algorithms = length(algorithms.color);

if has_spectral
    I_spectral_label = imread(spectral_label_filename);
    if size(I_spectral_label, 3) ~= 1
        error('Expected the spectral dispersion label image, "%s", to have only one channel.', spectral_label_filename);
    end
end

e_centroid_tables_rgb = cell(n_patches, 1);
e_edge_tables_rgb = cell(n_patches, 1);
if has_spectral
    e_centroid_tables_spectral = cell(n_patches, 1);
    e_edge_tables_spectral = cell(n_patches, 1);
end

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

patch_means_rgb = zeros(n_patches, n_channels_rgb);

%% Spectral registration calibration of the reference and estimated images

if has_spectral
    % Load all patches
    patch_weights = zeros(n_patches, 1);
    patch_weights_perChannel = zeros(n_patches, n_channels_rgb);
    pixels_reference = cell(n_patches, 1);
    pixels_reference_perChannel = cell(n_patches, n_channels_rgb);
    pixels_qhyper = cell(n_patches, n_channels_rgb);
    pixels_estimated = cell(n_patches, n_spectral_algorithms);
    for pc = 1:n_patches
        patch_mask = imopen(I_spectral_label == pc, se_clean);
        patch_xy = [
            reshape(X(patch_mask), [], 1),...
            reshape(Y(patch_mask), [], 1)
        ];
        n_px_patch = size(patch_xy, 1);
        patch_weights(pc) = n_px_patch;
        factors = vignettingfun(patch_xy) / vignetting_correction_factor;

        % Reference measurements
        patch_mean = reflectances(:, pc).';
        I_patch = repmat(patch_mean, n_px_patch, 1);
        I_patch = I_patch .* repmat(factors, 1, n_bands_measured);
        pixels_reference{pc} = I_patch;

        % Captured spectral image - Account for the Bayer pattern
        I = loadImage(qhyper_filename, qhyper_variable_name);
        for c = 1:n_channels_rgb
            mask_c = channel_mask(:, :, c);
            filter_c = mask_c(patch_mask);
            n_px_patch_c = sum(filter_c);
            patch_weights_perChannel(pc, c) = n_px_patch_c;
            I_patch = zeros(n_px_patch_c, n_bands_qhyper);
            for b = 1:n_bands_qhyper
                I_b = I(:, :, b);
                I_b = I_b(patch_mask);
                I_patch(:, b) = I_b(filter_c);
            end
            I_patch = channelConversion(I_patch, spectral_weights_qhyper, 2);
            pixels_qhyper{pc, c} = I_patch;

            I_patch = repmat(patch_mean, n_px_patch_c, 1);
            I_patch = I_patch .* repmat(factors(filter_c), 1, n_bands_measured);
            pixels_reference_perChannel{pc, c} = I_patch;
        end

        % Estimated spectral images
        for i = 1:n_spectral_algorithms
            I = loadImage(algorithm_filenames.spectral{i}, spectral_variable_name);
            I_patch = zeros(n_px_patch, n_bands_estimated);
            for c = 1:n_bands_estimated
                I_c = I(:, :, c);
                I_patch(:, c) = I_c(patch_mask);
            end
            I_patch = channelConversion(I_patch, spectral_weights, 2);
            pixels_estimated{pc, i} = I_patch;
        end
    end

    % Compute scaling factors
    radiance_ratios_qhyper = zeros(n_channels_rgb, 1, n_bands_measured);
    for c = 1:n_channels_rgb
        weights = repelem(patch_weights_perChannel(:, c), patch_weights_perChannel(:, c), 1);
        px = cell2mat(pixels_qhyper(:, c));
        weighted_px = repmat(weights, 1, n_bands_measured) .* px;
        radiance_ratios_qhyper(c, 1, :) = reshape(...
            dot(weighted_px, cell2mat(pixels_reference_perChannel(:, c)), 1) ./ ...
                dot(weighted_px, px, 1),...
            1, 1, []...
        );
    end

    radiance_ratios = zeros(n_spectral_algorithms, 1, n_bands_measured);
    weights = repelem(patch_weights, patch_weights, 1);
    px_reference = cell2mat(pixels_reference);
    for i = 1:n_spectral_algorithms
        px = cell2mat(pixels_estimated(:, i));
        weighted_px = repmat(weights, 1, n_bands_measured) .* px;
        radiance_ratios(i, 1, :) = reshape(...
            dot(weighted_px, px_reference, 1) ./ dot(weighted_px, px, 1),...
            1, 1, []...
        );
    end
end

%% Evaluate each patch of the ColorChecker

color_flags = true;
if has_spectral
    color_flags(2) = false;
end
for pc = 1:n_patches
    for color_flag = color_flags
        if color_flag
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
        if ~color_flag
            I_patch_rgb = repmat(...
                reshape(patch_means_rgb(pc, :), 1, 1, n_channels_rgb),...
                centroid_patch_width, centroid_patch_width, 1 ...
            );
            I_patch_rgb = I_patch_rgb .* repmat(factors, 1, 1, n_channels_rgb);
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
        if ~color_flag
            I_edge_inner_rgb = repmat(...
                reshape(patch_means_rgb(pc, :), 1, 1, n_channels_rgb),...
                n_edge_px_inner, 1, 1 ...
            );
            I_edge_inner_rgb = I_edge_inner_rgb .* repmat(factors, 1, 1, n_channels_rgb);
        end
        
        patch_xy = [X(edge_ind), Y(edge_ind)];
        factors = vignettingfun(patch_xy) / vignetting_correction_factor;
        n_edge_px = length(edge_ind);
        I_edge = repmat(patch_mean, n_edge_px, 1, 1);
        I_edge = I_edge .* repmat(factors, 1, 1, numel(patch_mean));
        if ~color_flag
            I_edge_rgb = repmat(reshape(patch_means_rgb(pc, :), 1, 1, n_channels_rgb), n_edge_px, 1, 1);
            I_edge_rgb = I_edge_rgb .* repmat(factors, 1, 1, n_channels_rgb);
        end
        
        % To smooth plots, points will be averaged for each distance from the
        % edge
        [edge_distance_unique, ~, edge_distance_ind] = unique(edge_distance);
        n_edge_distances = length(edge_distance_unique);
        edge_distance_filters = false(n_edge_px, n_edge_distances);
        for k = 1:n_edge_distances
            edge_distance_filters(:, k) = (edge_distance_ind == k);
        end
        edge_distance_unique = double(edge_distance_unique);
        
        % Evaluate the ColorChecker image with respect to edge error
        err_averaged = zeros(n_edge_distances, 1);
        if color_flag
            fg_edge_rgb = [];
            err = abs(I_edge - repmat(I_reference(edge_ind), 1, 1, n_channels_rgb)) ./ abs(I_edge);
            for c = 1:n_channels_rgb
                mask_c = channel_mask(:, :, c);
                filter_c = mask_c(edge_ind);
                for k = 1:n_edge_distances
                    filter_c_edge = filter_c & edge_distance_filters(:, k);
                    err_averaged(k) = mean(squeeze(err(filter_c_edge, :, c)));
                end
                err_filter = isfinite(err_averaged);
                fg_edge_rgb(c) = figure; %#ok<SAGROW>
                curve = smooth(edge_distance_unique(err_filter), err_averaged(err_filter), smoothing_span, smoothing_method);
                plot(edge_distance_unique(err_filter), curve, 'Color', [0, 0, 0], 'LineWidth', 2);
            end
        else
            I = loadImage(qhyper_filename, qhyper_variable_name);
            I_edge_qhyper = zeros(n_edge_px, 1, n_bands_qhyper);
            for c = 1:n_bands_qhyper
                I_c = I(:, :, c);
                I_edge_qhyper(:, :, c) = I_c(edge_ind);
            end
            I_edge_qhyper = channelConversion(I_edge_qhyper, spectral_weights_qhyper);
            
            % This image has a Bayer pattern
            err_averaged_c = cell(n_channels_rgb, n_edge_distances);
            for c = 1:n_channels_rgb                
                mask_c = channel_mask(:, :, c);
                filter_c = mask_c(edge_ind);
                I_edge_qhyper_c = I_edge_qhyper(filter_c, :, :);
                n_edge_px_c = size(I_edge_qhyper_c, 1);
                I_edge_qhyper_c = I_edge_qhyper_c .* repmat(radiance_ratios_qhyper(c, :, :), n_edge_px_c, 1, 1);
                err = mean(abs(I_edge(filter_c, :, :) -  I_edge_qhyper_c) ./ I_edge(filter_c, :, :) , 3);
                for k = 1:n_edge_distances
                    err_averaged_c{c, k} = squeeze(err(edge_distance_filters(filter_c, k)));
                end
            end
            for k = 1:n_edge_distances
                err_averaged(k) = mean(cell2mat(err_averaged_c(:, k)));
            end
            err_filter = isfinite(err_averaged);
            
            fg_edge_spectral = figure;
            curve = smooth(edge_distance_unique(err_filter), err_averaged(err_filter), smoothing_span, smoothing_method);
            plot(...
                edge_distance_unique(err_filter), curve,...
                'Color', [0, 0, 0], 'LineWidth', 2 ...
            );     
        end

        % Evaluation
        for evaluating_center = [true, false]
            if evaluating_center
                I_reference_eval = I_patch;
                if ~color_flag
                    I_reference_eval_rgb = I_patch_rgb;
                end
                name_params = fullfile(output_directory, [true_image_filename, sprintf('_patch%dCenter', pc)]);
            else
                I_reference_eval = I_edge_inner;
                if ~color_flag
                    I_reference_eval_rgb = I_edge_inner_rgb;
                end
                name_params = fullfile(output_directory, [true_image_filename, sprintf('_patch%dEdge', pc)]);
            end
            fg_spectral = struct;
            e_table = [];
            if color_flag
                e_table_rgb = [];
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
                for i = 0:n_color_algorithms
                    if i
                        I = loadImage(algorithm_filenames.color{i}, color_variable_name);
                    else
                        I = repmat(I_reference, 1, 1, n_channels_rgb);
                        I(~channel_mask) = NaN;
                    end
                    if evaluating_center
                        I = I(roi(1):roi(2), roi(3):roi(4), :);
                    else
                        I_edge_estimated_inner = zeros(n_edge_px_inner, 1, n_channels_rgb);
                        for c = 1:n_channels_rgb
                            I_c = I(:, :, c);
                            I_edge_estimated_inner(:, 1, c) = I_c(edge_ind_inner);
                        end
                        if i
                            for c = 1:n_channels_rgb
                                I_c = I(:, :, c);
                                I_edge_estimated_c = I_c(edge_ind);
                                figure(fg_edge_rgb(c))
                                hold on
                                err = abs(I_edge(:, :, c) -  I_edge_estimated_c) ./ abs(I_edge(:, :, c));
                                for k = 1:n_edge_distances
                                    err_averaged(k) = mean(squeeze(err(edge_distance_filters(:, k))));
                                end
                                err_filter = isfinite(err_averaged);
                                curve = smooth(edge_distance_unique(err_filter), err_averaged(err_filter), smoothing_span, smoothing_method);
                                plot(...
                                    edge_distance_unique(err_filter), curve,...
                                    'Color', evaluation_plot_colors_color(i, :),...
                                    'LineWidth', 2,...
                                    'Marker', evaluation_plot_markers{mod(i - 1, length(evaluation_plot_markers)) + 1},...
                                    'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
                                    );
                                hold off
                            end
                        end
                        I = I_edge_estimated_inner;
                    end
                    if i
                        name_params_i = [name_params, '_', algorithms.color{i}];
                        escaped_name = algorithms_escaped.color{i};
                    else
                        name_params_i = [name_params, '_input'];
                        escaped_name = 'Raw input';
                    end
                    e_table_current = evaluateAndSaveRGB(...
                        I, I_reference_eval, dp, true_image_name, escaped_name,...
                        name_params_i...
                    );
                    if ~isempty(e_table)
                        e_table = union(e_table_current, e_table);
                    else
                        e_table = e_table_current;
                    end
                end
                e_table_rgb = e_table;
                
            else
                for i = 1:n_spectral_algorithms
                    I = loadImage(algorithm_filenames.spectral{i}, spectral_variable_name);
                    name_params_i = [name_params, '_', algorithms.spectral{i}];
                    
                    if evaluating_center
                        I_rgb = channelConversion(I(roi(1):roi(2), roi(3):roi(4), :), color_weights);
                        I = channelConversion(I(roi(1):roi(2), roi(3):roi(4), :), spectral_weights)...
                            .* repmat(radiance_ratios(i, :, :), centroid_patch_width, centroid_patch_width, 1);
                    else
                        I_edge_estimated_inner = zeros(n_edge_px_inner, 1, n_bands_estimated);
                        for c = 1:n_bands_estimated
                            I_c = I(:, :, c);
                            I_edge_estimated_inner(:, 1, c) = I_c(edge_ind_inner);
                        end
                        I_rgb = channelConversion(I_edge_estimated_inner, color_weights);
                        
                        I_edge_estimated = zeros(n_edge_px, 1, n_bands_estimated);
                        for c = 1:n_bands_estimated
                            I_c = I(:, :, c);
                            I_edge_estimated(:, :, c) = I_c(edge_ind);
                        end
                        I_edge_estimated_rgb = channelConversion(I_edge_estimated, color_weights);
                        
                        for c = 1:n_channels_rgb
                            figure(fg_edge_rgb(c))
                            hold on
                            err = abs(I_edge_rgb(:, :, c) -  I_edge_estimated_rgb(:, :, c)) ./ abs(I_edge_rgb(:, :, c));
                            for k = 1:n_edge_distances
                                err_averaged(k) = mean(squeeze(err(edge_distance_filters(:, k))));
                            end
                            err_filter = isfinite(err_averaged);
                            curve = smooth(edge_distance_unique(err_filter), err_averaged(err_filter), smoothing_span, smoothing_method);
                            plot(...
                                edge_distance_unique(err_filter), curve,...
                                'Color', evaluation_plot_colors_spectral(i, :),...
                                'LineWidth', 2,...
                                'Marker', evaluation_plot_markers{mod(i - 1, length(evaluation_plot_markers)) + 1},...
                                'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
                                );
                            hold off
                        end
                        
                        I_edge_estimated = channelConversion(I_edge_estimated, spectral_weights)...
                            .* repmat(radiance_ratios(i, :, :), n_edge_px, 1, 1);
                        I = channelConversion(I_edge_estimated_inner, spectral_weights)...
                            .* repmat(radiance_ratios(i, :, :), n_edge_px_inner, 1, 1);
                        figure(fg_edge_spectral)
                        hold on
                        err = mean(abs(I_edge -  I_edge_estimated) ./ abs(I_edge), 3);
                        for k = 1:n_edge_distances
                            err_averaged(k) = mean(squeeze(err(edge_distance_filters(:, k))));
                        end
                        err_filter = isfinite(err_averaged);
                        curve = smooth(edge_distance_unique(err_filter), err_averaged(err_filter), smoothing_span, smoothing_method);
                        plot(...
                            edge_distance_unique(err_filter), curve,...
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
                    
                    e_table_current = evaluateAndSaveRGB(...
                        I_rgb, I_reference_eval_rgb, dp, true_image_name, algorithms_escaped.spectral{i},...
                        name_params_i...
                        );
                    if ~isempty(e_table_rgb)
                        e_table_rgb = union(e_table_current, e_table_rgb);
                    else
                        e_table_rgb = e_table_current;
                    end
                end
                writetable(e_table, [name_params, '_evaluateSpectral.csv']);
                writetable(e_table_rgb, [name_params, '_evaluateRGB.csv']);
                % Also save completed figures
                evaluateAndSaveSpectral(output_directory, dp, true_image_filename_eval, algorithms_escaped.spectral, fg_spectral);
                
                if ~evaluating_center
                    figure(fg_edge_spectral)
                    legend(['Bandpass-filtered image'; algorithms_escaped.spectral]);
                    title(sprintf('Relative error around patch %d edge', pc));
                    xlabel('Position relative to patch borders [pixels]')
                    ylabel('Relative error (MRAE)')
                    savefig(...
                        fg_edge_spectral,...
                        [ name_params, '_errSpectral.fig'], 'compact'...
                    );
                    close(fg_edge_spectral);
                    
                    for c = 1:n_channels_rgb
                        figure(fg_edge_rgb(c))
                        legend(['Raw input'; algorithms_escaped.color; algorithms_escaped.spectral]);
                        title(sprintf('Absolute error around patch %d edge in channel %d', pc, c));
                        xlabel('Position relative to patch borders [pixels]')
                        ylabel('Relative error')
                        savefig(...
                            fg_edge_rgb(c),...
                            [ name_params, sprintf('_errChannel%d.fig', c)], 'compact'...
                        );
                        close(fg_edge_rgb(c));
                    end
                end
            end
            
            if ~color_flag || ~has_spectral
                if evaluating_center
                    if has_spectral
                        e_centroid_tables_spectral{pc} = e_table;
                    end
                    e_centroid_tables_rgb{pc} = e_table_rgb;
                else
                    if has_spectral
                        e_edge_tables_spectral{pc} = e_table;
                    end
                    e_edge_tables_rgb{pc} = e_table_rgb;
                end
                for i = 1:n_spectral_algorithms
                    algorithm_index = e_table_rgb.Algorithm == algorithms_escaped.spectral{i};
                    e_table_rgb(algorithm_index, :) = [];
                end
            end
        end
    end
end

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

if has_spectral
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
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {
    'bands_estimated', 'bands_measured', 'bands_qhyper',...
    'vignetting_data', 'vignetting_correction_factor',...
    'spectral_weights', 'spectral_weights_qhyper',...
    'algorithms', 'algorithm_filenames'...
} ];
save_data_filename = fullfile(output_directory, 'EvaluateColorCheckerData.mat');
save(save_data_filename, save_variables_list{:});