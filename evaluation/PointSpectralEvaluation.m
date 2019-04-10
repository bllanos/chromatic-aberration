%% Evaluate spectral reconstruction with respect to measured spectra
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Spectral measurements
%
% '.csv' files containing a header row, a first column for wavelength values,
% and remaining columns for relative or absolute spectral reflectances. The
% first file will be designated as the reference file for evaluation. In the
% output files describing the evaluation results, the names of the sets of
% measurements will be the unique portions of their filenames.
%
% ### Spectral images
%
% Spectral images to be evaluated, along with the following supporting
% information:
% - A set of pixel locations corresponding to the measured spectra, assumed to
%   be valid for all images. The locations should be provided as a two-column
%   matrix of patch centre integer x and y-coordinates, in a '.mat' file.
% - Either a single '.mat' file containing the wavelengths at which the spectral
%   images are defined, or one file of wavelengths per image, if the images have
%   different spectral resolutions.
%
% In the output files describing the evaluation results, the names of the images
% being evaluated will be the unique portions of their filenames.
%
% ## Output
%
% ### Data file output
%
% #### Intermediate data and parameters
% A '.mat' file containing the following variables, as appropriate:
% - 'bands': A cell vector of length 'm', where each cell contains a vector of
%   wavelengths at which the data with the same index is defined. The cells of
%   'bands' are ordered by measured data followed by spectral images.
% - 'spectra': An m x n cell array, where 'm' is the number of sets of spectral
%   data, and 'n' is the number of spectra in the sets. The rows of 'spectra'
%   index measured data first, then spectral images. Each cell contains a vector
%   of spectral intensities.
% - 'spectra_aligned_eval': An (m - 1) x n x 2 cell array, where 'm' is the
%   number of sets of spectral data to evaluate, and 'n' is the number of
%   spectra in the sets. The rows of 'spectra' index measured data first, then
%   spectral images. Each cell contains a vector of spectral intensities that
%   have been rescaled to align with the reference spectral intensities.
%   `spectra_aligned_eval(:, :, 1)` contains spectral intensities sampled
%   according to `bands{1}`, whereas `spectra_aligned_eval(i, :, 1)` contains
%   spectral intensities sampled according to `bands{i + 1}`.
% - 'spectra_aligned_reference': An (m - 1) x n cell array, where 'm' is the
%   number of sets of spectral data to evaluate, and 'n' is the number of
%   spectra in the sets. The rows of 'spectra' index measured data first, then
%   spectral images. Each cell in the i-th row contains a vector of the
%   reference spectral intensities that have been resampled to the sampling
%   space of `bands{i + 1}`.
% - 'conditions': A cell vector containing the names of the datasets which were
%   evaluated, extracted from the names of the input files. The first cell
%   contains the name of the reference data. Subsequent cells contain the names
%   of the sets of measured spectra, followed by the names of the spectral images.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`.
%
% The file is saved as 'PointSpectralEvaluationData.mat'.
%
% #### Evaluation results
%
% Spectral error metrics are output in the form of two CSV files. Each CSV file
% contains results for all measured spectra (other than the reference
% measurement) and spectral images. The spectral data being evaluated will be
% listed in the files under the unique portions of the input filenames. There
% are two entries for each test measurement, corresponding to the cases where
% either the reference, or the test measurement, was resampled when the two were
% converted to the same sampling space.
%
% Spectral error metrics computed per measurement location are saved as
% '*_perPatchEvaluation.csv'.
%
% Error metrics are also aggregated across measurement locations, and saved as
% '*_summaryEvaluation.csv'.
%
% ### Graphical output
%
% One figure is generated for each measurement location, showing the reference
% spectra and all aligned test spectra. All spectra are expressed in the
% sampling space of the reference spectra. The figures are saved to files whose
% names contain the indices of the measurement locations and the measurement
% location coordinates.
%
% ## Detailed description of processing
%
% For each non-reference set of measured spectra, and spectral image, the script
% provides three options for how it should be aligned with the reference set of
% measured spectra:
% - 'none': The spectra are treated as absolute measurements, and are compared
%   as-is.
% - 'reference': Spectra are matched to the reference set by rescaling them so
%   that they are identical to the reference set for a designated reference
%   point (e.g. a white patch).
% - 'global': Spectra are matched to the reference set using the scaling
%   transformation that provides the best least-squares fit between the two sets
%   of spectra.
%
% ## Notes
% - This script can correct spectral images for vignetting, if desired.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 8, 2019

% List of parameters to save with results
parameters_list = {
    'images_filenames',...
    'images_variable',...
    'bands_filenames',...
    'bands_variable',...
    'location_filename',...
    'location_variable',...
    'vignetting_mask_filename',...
    'vignetting_mask_label',...
    'max_degree_vignetting',...
    'vignetting_erosion_radius',...
    'vignetting_calibration_bands',...
    'spectra_filenames',...
    'spectra_data_columns',...
    'wavelength_range',...
    'reference_patch_index',...
    'alignment_method',...
    'patch_side_length',...
    'output_filename'...
};

%% Input data and parameters

% A list of filenames of spectral images to evaluate.
% '.mat' or image files can be loaded
images_filenames = { [] };
images_variable = 'img_hs'; % Used only when loading '.mat' files

% Filenames of '.mat' files containing the wavelengths corresponding to the
% bands in the spectral images. Either one file can be given, or one file per
% image.
bands_filenames = { [] };
bands_variable = 'bands'; % Variable name in the above files

% Interpolation functions to use when resampling the spectral images to
% different spectral resolutions.
images_interpolants = {
    @gaussian
};

% Filename of a '.mat' file containing the image locations corresponding to the
% measured spectra
location_filename = [];
location_variable = 'patch_centers'; % Variable name in the file

% Label image for vignetting calibration (an image file, not a '.mat' file).
% Can be empty (`[]`), in which case vignetting correction will not be used.
vignetting_mask_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/colorChecker_preprocessed/unfiltered/d2_colorChecker30cm_unfiltered_background0_patches1to24_frame25.png';
vignetting_mask_label = 25; % The value of pixels which are to be used to calibrate vignetting

% Maximum degree of the polynomial model of vignetting
max_degree_vignetting = 5;

% Radius of a disk structuring element to use for morphological erosion of the
% vignetting calibration mask
vignetting_erosion_radius = 3;

% A vector indicating, for each image, which band should be used to
% calibrate vignetting. Zero elements are interpreted to mean that vignetting
% correction should be calibrated separately for each band.
vignetting_calibration_bands = [
];

% Measured spectra CSV files, with the reference set's filename given first
spectra_filenames = {
    '/home/llanos/GoogleDrive/ThesisResearch/Data/20180626_SpectralCharacterizationOfSetup/spectra_averaged.csv';
};

% Indices of the columns of interest in the CSV files
spectra_data_columns = {
    13:36;
};

% Interpolation functions to use when resampling the measured spectra to
% different spectral resolutions.
spectra_interpolants = {
    @triangle
};

% Wavelength range to truncate spectral measurements to
wavelength_range = [375, 725];

% The reference patch for aligning spectral measurements, if such an alignment
% method ('reference') is to be used.
reference_patch_index = 20;

% Method used to align spectral measurements from different sets
alignment_method = 'global';

% Odd-integer size of the image regions to extract around the pixel location of
% each measurement
patch_side_length = 15;

% Prefix for the filenames of the output files. Must not contain spaces.
output_filename = 'colorChecker';

% Output directory
output_directory = '/home/llanos/Downloads';

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% ## Debugging Flags
vignettingPolyfitVerbose = true;

%% Load information describing the spectral images

n_images = length(images_filenames);
n_measurement_sets = length(spectra_filenames);
if n_measurement_sets == 0
    error('At least one set of spectral measurements is expected, which will be the reference set.');
elseif n_measurement_sets == 1 && n_images == 0
    error('There are no sets of spectral measurements, or images, to evaluate.');
end

n_all = n_measurement_sets + n_images;

has_images = n_images > 0;
bands_filenames_rep = bands_filenames;
if has_images && length(bands_filenames) ~= n_images
    if length(bands_filenames) == 1
        bands_filenames_rep = repelem(bands_filenames_rep, n_images, 1);
    else
        error('%d images were provided, but %d filenames defining spectral sampling information.',...
            n_images, length(bands_filenames));
    end
end
    
bands = cell(n_all, 1);
for i = 1:n_images
    load(bands_filenames_rep{i}, bands_variable);
    if exist(bands_variable, 'var')
        bands{i + n_measurement_sets} = eval(bands_variable);
    end
    if ~exist(bands_variable, 'var') || isempty(bands{i + n_measurement_sets})
        error('No wavelengths loaded from file "%s".', bands_filenames_rep{i})
    end
end

if has_images
    load(location_filename, location_variable);
    if exist(location_variable, 'var')
        patch_centers = eval(location_variable);
    end
    if ~exist(location_variable, 'var') || isempty(patch_centers)
        error('No image patch center locations loaded from file "%s".', location_filename)
    end
end

%% Load all spectral measurements and spectral images

if length(spectra_data_columns) ~= n_measurement_sets
    error('%d CSV files of spectra were provided, but %d sets of columns of interest in the files.',...
            n_measurement_sets, length(spectra_data_columns));
end

for i = 1:n_measurement_sets    
    sample_table = readtable(spectra_filenames{i});
    variable_names = sample_table.Properties.VariableNames;
    bands{i} = sample_table{:, 1};
    bands_filter = (bands{i} >= wavelength_range(1)) & (bands{i} <= wavelength_range(2));
    bands{i} = bands{i}(bands_filter);
    spectra_i = sample_table{bands_filter, spectra_data_columns{i}};
    if i == 1
        n_spectra = size(spectra_i, 2);
        if has_images && (n_spectra ~= size(patch_centers, 1))
            error('The number of reference measurements is not consistent with the number of image locations.');
        end
        spectra = cell(n_all, n_spectra);
        spectra_reference = spectra_i;
    else
        if n_spectra ~= size(spectra_i, 2)
            error('%d measurements were extracted from "%s", but there are %d in the reference set.',...
                size(spectra_i, 2), spectra_filenames{i}, n_spectra);
        end
    end
    spectra(i, :) = num2cell(spectra_i, 2);
end

correct_vignetting = ~isempty(vignetting_mask_filename);
if correct_vignetting
    I_vignetting_label = imread(vignetting_mask_filename);
    if size(I_vignetting_label, 3) ~= 1
        error('Expected the vignetting mask image, "%s", to have only one channel.', vignetting_mask_filename);
    end
    mask_vignetting = (I_vignetting_label == vignetting_mask_label);
    se_clean = strel('disk', vignetting_erosion_radius);
    mask_vignetting = imerode(mask_vignetting, se_clean);
end

if mod(patch_side_length, 2) == 0
    error('`patch_side_length` should be an odd integer.');
end
half_width = floor(patch_side_length / 2);
for i = 1:n_images
    I = loadImage(images_filenames{i}, images_variable);
    
    % Correct vignetting
    if correct_vignetting
        if vignetting_calibration_bands(i) ~= 0
            [vignettingfun, vignetting_data] = vignettingPolyfit(...
                I(:, :, vignetting_calibration_bands(i)), mask_vignetting,...
                max_degree_vignetting, bayer_pattern, vignettingPolyfitVerbose...
            );
            % Make the vignetting relative to the center of the image, which is valid if we
            % assume that the center of the image is within the domain of samples used to
            % fit the vignetting model
            vignetting_correction_factor = vignettingfun([size(I, 2), size(I, 1)] / 2);
            I = correctVignetting(I, vignettingfun) * vignetting_correction_factor;
        else
            for c = 1:size(I, 3)
                [vignettingfun, vignetting_data] = vignettingPolyfit(...
                    I(:, :, c), mask_vignetting,...
                    max_degree_vignetting, bayer_pattern, vignettingPolyfitVerbose...
                );
                vignetting_correction_factor = vignettingfun([size(I, 2), size(I, 1)] / 2);
                I(:, :, c) = correctVignetting(I(:, :, c), vignettingfun) * vignetting_correction_factor;
            end
        end
    end

    for j = 1:n_spectra
        roi = [
            patch_centers(2) - half_width, patch_centers(2) + half_width,...
            patch_centers(1) - half_width, patch_centers(1) + half_width
        ];
        spectra(i + n_measurement_sets, j) = mean(mean(I(roi(1):roi(2), roi(3):roi(4), :), 1), 2);
    end
end

%% Register all spectra with the reference spectra

n_eval = n_all - 1;
n_alignments = 2;
spectra_aligned_eval = cell(n_eval, n_spectra, n_alignments);
spectra_aligned_reference = cell(n_eval, n_spectra);
for i = 2:n_all
    spectra_i = cell2mat(spectra(i, :));
    for direction = 1:n_alignments
        align_to_reference = (direction == 1);
        
        % Match spectral sampling space
        if align_to_reference
            src_bands = bands{i};
            dst_bands = bands{1};
        else
            src_bands = bands{1};
            dst_bands = bands{i};
        end
        if align_to_reference
            if i <= n_measurement_sets
                spectral_weights = resamplingWeights(...
                    dst_bands, src_bands, spectra_interpolants{i}, findSamplingOptions.bands_padding...
                );
            else
                spectral_weights = resamplingWeights(...
                    dst_bands, src_bands, images_interpolants{i - n_measurement_sets}, findSamplingOptions.bands_padding...
                );
            end
        else
            spectral_weights = resamplingWeights(...
                dst_bands, src_bands, spectra_interpolants{1}, findSamplingOptions.bands_padding...
            );
        end
        if align_to_reference
            spectra_aligned_i = channelConversion(spectra_i, spectral_weights, 2);
            spectra_aligned_reference_i = spectra_reference;
        else
            spectra_aligned_i = spectra_i;
            spectra_aligned_reference_i = channelConversion(spectra_reference, spectral_weights, 2);
        end
        
        % Match spectral intensities
        if strcmp(alignment_method, 'global')
            scaling_i = dot(spectra_aligned_i, spectra_aligned_reference_i, 2) ./ ...
                dot(spectra_aligned_i, spectra_aligned_i, 2);
        elseif strcmp(alignment_method, 'reference')
            scaling_i = spectra_aligned_reference_i(:, reference_patch_index) ./ ...
                spectra_aligned_i(:, reference_patch_index);
        elseif strcmp(alignment_method, 'none')
            scaling_i = ones(size(spectra_aligned_i, 1));
        else
            error('Unrecognized value of `alignment_method`.');
        end
        
        spectra_aligned_eval(i - 1, :, direction) = num2cell(spectra_aligned_i * repmat(scaling_i, n_spectra, 1));
        
        if ~align_to_reference
            spectra_aligned_reference(i - 1, :) = num2cell(spectra_aligned_reference_i, 2);
        end               
    end
end

%% Prepare for evaluation

conditions = cell(n_all, 1);
for i = 1:n_all
    if i <= n_measurement_sets
        [~, conditions{i}] = fileparts(spectra_filenames{i});
    else
        [~, conditions{i}] = fileparts(image_filenames{i - n_measurement_sets});
    end
end

conditions = trimCommon(conditions);
conditions_escaped = strrep(algorithms.spectral, '_', '\_');

e_tables_reference_bands = cell(n_patches, 1);
e_tables_other_bands = cell(n_patches, 1);

plot_colors = jet(n_eval);
plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
plot_styles = {'--', ':', '-.'};

%% Evaluation and graphical output

rmse = zeros(n_eval, n_spectra, n_alignments);
mrae = zeros(n_eval, n_spectra, n_alignments);
gof = zeros(n_eval, n_spectra, n_alignments);
for s = 1:n_spectra
    fg = figure;
    hold on
    plot(...
        bands{1}, spectra_reference(:, s),...
        'Color', 'k', 'LineWidth', 2, 'Marker', 'none'...
    );
    xlabel('Wavelength [nm]');
    if strcmp(alignment_method, 'none')
        ylabel('Radiance');
    else
        ylabel('Relative radiance');
    end
    for i = 1:n_eval
        plot_color = plot_colors(i, :);
        plot_marker = 'none';
        plot_style = plot_styles{mod(i - 1, length(plot_styles)) + 1};  
        plot(...
            bands{1}, spectra_aligned_eval{i, s, 1},...
            'Color', plot_color, 'LineWidth', 2,...
            'Marker', plot_marker, 'LineStyle', plot_style...
        );
    
        for direction = 1:n_alignments
            align_to_reference = (direction == 1);
            if align_to_reference 
                [...
                    rmse(i, s, direction),...
                    mrae(i, s, direction),...
                    gof(i, s, direction)...
                ] = metrics(spectra_aligned_eval{i, s, direction}, spectra_reference(:, s), 1, 0, false);
            else
                [...
                    rmse(i, s, direction),...
                    mrae(i, s, direction),...
                    gof(i, s, direction)...
                ] = metrics(spectra_aligned_eval{i, s, direction}, spectra_aligned_reference(i, s), 1, 0, false);
            end
        end
    end
        
    title(sprintf('Spectral radiance for patch %d', s));
    hold off
    legend(conditions_escaped{:});
    savefig(...
        fg, [...
            output_filename,...
            sprintf(...
                '_evalPatch%d_X%dY%dW%dH%d.fig',...
                s,...
                patch_centers(1),...
                patch_centers(2),...
                patch_side_length,...
                patch_side_length...
            )...
        ], 'compact'...
    );
    close(fg);
end

% Output a CSV file of per-patch statistics
for s = 1:n_spectra
    e_spectral.(sprintf('Patch%d_MRAE_toReference', s)) = mrae(:, s, 1);
    e_spectral.(sprintf('Patch%d_RMSE_toReference', s)) = rmse(:, s, 1);
    e_spectral.(sprintf('Patch%d_GOF_toReference', s)) = gof(:, s, 1);
    e_spectral.(sprintf('Patch%d_MRAE_toOther', s)) = mrae(:, s, 2);
    e_spectral.(sprintf('Patch%d_RMSE_toOther', s)) = rmse(:, s, 2);
    e_spectral.(sprintf('Patch%d_GOF_toOther', s)) = gof(:, s, 2);    
end
e_spectral.Source = conditions(2:end);
e_spectral_table = struct2table(e_spectral);
writetable(e_spectral_table, [output_filename, '_perPatchEvaluation.csv']);

% Output a global CSV file
e_spectral = struct('Source', conditions(2:end));
e_spectral.MRAE_max_toReference = max(mrae(:, :, 1), [], 2);
e_spectral.MRAE_mean_toReference = mean(mrae(:, :, 1), 2);
e_spectral.MRAE_median_toReference = median(mrae(:, :, 1), 2);
e_spectral.RMSE_max_toReference = max(rmse(:, :, 1), [], 2);
e_spectral.RMSE_mean_toReference = mean(rmse(:, :, 1), 2);
e_spectral.RMSE_median_toReference = median(rmse(:, :, 1), 2);
e_spectral.GOF_min_toReference = min(gof(:, :, 1), [], 2);
e_spectral.GOF_mean_toReference = mean(gof(:, :, 1), 2);
e_spectral.GOF_median_toReference = median(gof(:, :, 1), 2);
e_spectral.MRAE_max_toOther = max(mrae(:, :, 2), [], 2);
e_spectral.MRAE_mean_toOther = mean(mrae(:, :, 2), 2);
e_spectral.MRAE_median_toOther = median(mrae(:, :, 2), 2);
e_spectral.RMSE_max_toOther = max(rmse(:, :, 2), [], 2);
e_spectral.RMSE_mean_toOther = mean(rmse(:, :, 2), 2);
e_spectral.RMSE_median_toOther = median(rmse(:, :, 2), 2);
e_spectral.GOF_min_toOther = min(gof(:, :, 2), [], 2);
e_spectral.GOF_mean_toOther = mean(gof(:, :, 2), 2);
e_spectral.GOF_median_toOther = median(gof(:, :, 2), 2);

e_spectral_table = struct2table(e_spectral);
writetable(e_spectral_table, [output_filename, '_summaryEvaluation.csv']);

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {
    'bands', 'spectra', 'spectra_aligned_eval', 'spectra_aligned_reference',...
    'conditions'...
} ];
save_data_filename = fullfile(output_directory, 'PointSpectralEvaluationData.mat');
save(save_data_filename, save_variables_list{:});