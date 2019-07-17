%% Extract colours from a ColorChecker and calibrate colour correction
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
% Along with the captured image, a single-channel label images is required. It
% is used for colour calibration and vignetting correction. For the 24 patches
% of the chart, the corresponding pixels in the label image must have values
% equal to the patch indices. The uniformly-coloured portions of the frame
% surrounding the patches, used to calibrate vignetting, must be labelled 25.
% Lastly, the rest of the image should have a different label (e.g. zero).
%
% ### Spectral reflectances and conversion to colour
%
% #### Spectral reflectances
% A '.csv' file containing a header row, a first column for wavelength values,
% and remaining columns for relative spectral reflectances of each patch in the
% ColorChecker, in the same order that they are labelled in the input label
% image.
%
% #### CIE tristimulus functions
% A '.mat' file containing a variable 'xyzbar', which can be used as the
% `C` input argument of 'cieSpectralToColor()'.
%
% ## Output
%
% ### Data file output
%
% #### Intermediate data and parameters
% A '.mat' file containing the following variables:
% - 'bands_measured': A vector containing the wavelengths of the spectral
%   bands used by the reference reflectance data.
% - 'vignetting_data': The model of vignetting describing the ColorChecker
%   image. 'vignetting_data' is the output argument of the same name of
%   'vignettingPolyfit()'. Refer to the documentation of 'vignettingPolyfit.m'
%   for details.
% - 'vignetting_correction_factor': A scalar used to map vignetting-corrected
%   intensities back to the scale of the original image intensities.
% - 'M_ls': A 3 x 3 matrix for mapping from the raw RGB colours of the
%   ColorChecker image to the XYZ colours of the corresponding measured
%   reflectances. 'C_ls' is computed using basic least squares fitting.
% - 'M_lsnonneg': Similar to 'C_ls', but computed using non-negative least
%   squares fitting (MATLAB's 'lsqnonneg' function).
% - 'M_homog': Similar to 'C_ls', but computed using the colour homography
%   technique of Finlayson, Gong, and Fisher 2019.
% - 'M_rp': A root polynomial regression colour correction structure,
%   computed according to the method of Finlayson, Mackiewicz, and Hurlbert
%   2015. 'M_rp' is the output of the third-party 'rpcal()' function.
% - 'wb_scale': A scalar intensity by which the image will be divided
%   before white balancing using MATLAB's built-in chromatic adaptation
%   functionality, in order to produce an output image with suitable
%   brightness. 'wb_scale' is the Green channel intensity of a designated
%   white patch (whether or not the white patch is saturated).
% - 'wb_illum': The colour of a neutral patch in the scene, for use with
%   MATLAB's built-in chromatic adaptation functionality, 'chromadapt()'. A
%   3-element vector.
% - 'reflectances_xyz': CIE XYZ colours of the measured spectral reflectances
%   for the patches in the ColorChecker.
% - 'patch_means_rgb': Vignetting-corrected camera RGB responses for the patches
%   in the ColorChecker.
% - 'patch_means_rgb_nonuniform': Raw camera RGB responses for the patches in
%   the ColorChecker. (A version 'patch_means_rgb', without vignetting
%   correction.)
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`.
%
% The file is saved as 'CalibrateColorCorrectionData.mat'.
%
% ## Notes
% - The whitepoint to use for XYZ to RGB conversion should be `[1, 1, 1]`,
%   because an equal energy radiator is the illuminant used to convert the
%   spectral reflectance data to XYZ.
% - To apply colour corrections to images, see 'ColorCorrection.m'.
%
% ## References
% - Finlayson, G., Gong, H., and Fisher, R.B. "Color Homography: Theory and
%   Applications." IEEE Transactions on Pattern Analysis and Machine
%   Intelligence, vol. 41, no. 1, pp. 20-33, 2019.
%   doi:10.1109/TPAMI.2017.2760833
% - Finlayson, G.D., Mackiewicz, M. and Hurlbert, A. "Color correction
%   using root-polynomial regression." IEEE Transactions on Image
%   Processing, vol. 24, no. 5, pp.1460-1470, 2015.
%   - I used the implementation by Han Gong provided in the code
%     corresponding to Finlayson, Gong, and Fisher 2019.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 15, 2019

% List of parameters to save with results
parameters_list = {
    'color_label_filename',...
    'reference_filename',...
    'reference_variable_name',...
    'xyzbar_filename',...
    'reflectances_filename',...
    'reflectance_data_patch_columns',...
    'patch_filter',...
    'max_degree_vignetting',...
    'n_patches',...
    'white_index',...
    'centroid_patch_width'...
};

%% Input data and parameters

% Label image (an image file, not a '.mat' file)
color_label_filename = fullfile('.', 'demo_data', 'colorChecker_labels', 'd2_colorChecker30cm_d3_green_labels.png');

% Raw image of the ColorChecker ('.mat' or image files can be loaded)
reference_filename = fullfile('.', 'demo_data', 'multispectral_images', 'd2_colorChecker30cm_raw.mat');
reference_variable_name = 'I_raw'; % Used only when loading '.mat' files

% CIE tristimulus functions
xyzbar_filename = '${FILEPATH}';

% Sample spectral reflectances
reflectances_filename = fullfile('.', 'demo_data', 'spectral_data', 'spectra_averaged.csv');

% Categorization of the samples in the file of reflectances
reflectance_data_patch_columns = 13:36;

% A filter indicating (elements that are `false`) which patches should be
% excluded from color correction calibration
patch_filter = [true(18, 1); false; true(5, 1)]; % The white patch is saturated

% Maximum degree of the polynomial model of vignetting
max_degree_vignetting = 3;

% Number of ColorChecker patches
n_patches = 24;

% The index of the white patch on the ColorChecker
white_index = 19;

% The index of the neutral patch used for chromatic adaptation (white balancing)
neutral_index = 20; % Use the white patch, unless it is saturated

% Odd-integer size of the regions to extract around the centroids of the
% ColorChecker patches
centroid_patch_width = 15;

% Output directory
output_directory = fullfile('.', 'demo_data', 'color_correction');

% ## Parameters which do not usually need to be changed
run('SetFixedParameters.m')

% ## Debugging Flags
vignettingPolyfitVerbose = false;

%% Parameter checking

if ~patch_filter(neutral_index)
    error('`neutral_index` should correspond to a usable patch.');
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
n_bands_measured = length(bands_measured);
reflectances = sample_table{:, reflectance_data_patch_columns};

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

%% Extract colours from the image of the ColorChecker

n_channels_rgb = 3;
if mod(centroid_patch_width, 2) == 0
    error('`centroid_patch_width` should be an odd integer.');
end
half_width = floor(centroid_patch_width / 2);
se_clean = strel('square', 3); % Structuring element for cleaning up the label image
channel_mask = bayerMask(image_sampling(1), image_sampling(2), bayer_pattern);

% Enumerate the positions of all pixels
[X, Y] = meshgrid(1:image_sampling(2), 1:image_sampling(1));
X = X - 0.5; % Place coordinates at pixel centres
Y = Y - 0.5;

patch_means_rgb = zeros(n_patches, n_channels_rgb);
patch_means_rgb_nonuniform = zeros(n_patches, n_channels_rgb);

for pc = 1:n_patches
    patch_mask = imopen(I_color_label == pc, se_clean);

    % Find the centroid of the patch
    cc_stats = regionprops(patch_mask, 'Centroid');
    if numel(cc_stats) ~= 1
        error('Expected only one binary mask region for patch %d in the label image.', pc);
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
    I_patch = I_reference(roi(1):roi(2), roi(3):roi(4));
    I_patch_vc = I_patch ./ factors;

    % Find the mean colour in this region
    for c = 1:n_channels_rgb
        patch_means_rgb(pc, c) = mean(I_patch_vc(channel_mask_centroid(:, :, c)));
        patch_means_rgb_nonuniform(pc, c) = mean(I_patch(channel_mask_centroid(:, :, c)));
    end
end

%% Create a colour correction matrix

% Least squares
rgb_filtered = patch_means_rgb(patch_filter, :);
xyz_filtered = reflectances_xyz(patch_filter, :);
n_patches_filtered = size(rgb_filtered, 1);
M_ls = (rgb_filtered \ xyz_filtered).'

% Non-negative least squares
rgb_filtered_block = zeros(n_patches_filtered * n_channels_rgb, n_channels_rgb * n_channels_rgb);
for c = 1:n_channels_rgb
    rgb_filtered_block(...
        ((c - 1) * n_patches_filtered + 1):(c * n_patches_filtered),...
        ((c - 1) * n_channels_rgb + 1):(c * n_channels_rgb)...
    ) = rgb_filtered;
end
xyz_filtered_block = reshape(xyz_filtered, [], 1);

options = optimset('Display','notify');
M_lsnonneg = reshape(...
    lsqnonneg(rgb_filtered_block, xyz_filtered_block, options),...
    n_channels_rgb, n_channels_rgb...
).'

% Colour homography
M_homog = ransachomocal_luv(...
    patch_means_rgb_nonuniform(patch_filter, :),...
    xyz_filtered,...
    reflectances_xyz(white_index, :),...
    rgb_filtered...
);
xyz_est = uea_homocvt(patch_means_rgb, M_homog);
if patch_filter(white_index)
    % Normalize by the white patch's green intensity
    M_homog = M_homog ./ xyz_est(white_index, 2)
else
    % Normalize by finding the best fit scaling
    M_homog_scale = reshape(xyz_est, [], 1) \ reshape(reflectances_xyz, [], 1);
    M_homog = M_homog * M_homog_scale
end

% Root-polynomial regression
M_rp = rpcal(...
    patch_means_rgb_nonuniform(patch_filter, :),...
    xyz_filtered...
);

%% Extract information for white balancing

wb_illum = patch_means_rgb_nonuniform(neutral_index,:)
wb_scale = patch_means_rgb_nonuniform(white_index, 2)

%% Save parameters and results to a file
save_variables_list = [ parameters_list, {
    'bands_measured', 'vignetting_data', 'vignetting_correction_factor',...
    'M_ls', 'M_lsnonneg', 'M_homog', 'M_rp', 'wb_scale', 'wb_illum',...
    'reflectances_xyz', 'patch_means_rgb', 'patch_means_rgb_nonuniform'...
} ];
save_data_filename = fullfile(output_directory, 'CalibrateColorCorrectionData.mat');
save(save_data_filename, save_variables_list{:});