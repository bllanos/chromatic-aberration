%% Chromatic aberration calibration from captured images
% Obtain a dispersion model from RAW images of disk calibration patterns.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
%
% Graphical output from 'plotXYLambdaPolyfit()'.
%
% ## References
% - Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%   single-shot hyperspectral imaging using a prism." ACM Transactions
%   on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%   doi:10.1145/3130800.3130896
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 23, 2018

%% Input data and parameters

% ## Input images
%
% Images are expected to have been preprocessed, such as using
% 'AverageRAWImages.m', so that they do not need to be linearized after
% being loaded. Images will simply be loaded with the Image Processing
% Toolbox 'imread()' function.
%
% RAW (non-demosaicked) images are expected.

% Partial filepaths containing the paths of the input images, and the
% portions of the filenames excluding the extension and the wavelengths.
% This script will append wavelength numbers and the extension to the
% partial filepaths. Each partial filepath is expected to correspond to all
% of the wavelengths in `wavelengths` below.
%
% If colour channels are to be treated as separate wavelengths, only the
% extension will be appended to the partial filepaths.
%
% This script will also search for image masks, using filenames which are
% constructed by appending '_mask' and `mask_ext` (below) to the partial
% filepaths. Masks are used to avoid processing irrelevant portions of
% images.
partial_filepaths = {
    '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20170808_OpticalTableMatrix/averaged/d44_a22_far_disksWhite';
    '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20170808_OpticalTableMatrix/averaged/d28_a22_far_disksWhite'
    };

% Filename extension (excluding the leading '.')
ext = 'tif';
% Mask filename extension
mask_ext = 'png';

% Threshold used to binarize mask images
mask_threshold = 0.5;

% ## Spectral information

% Find dispersion between colour channels, as opposed to between images
% taken under different spectral bands
rgb_mode = true;

if rgb_mode
    bands = 1:3;
    reference_index = 2; % Green colour channel
    bands_to_rgb = eye(3);
else
    % Wavelengths will be expected at the end of filenames
    bands = [];
    reference_index = 0;
    bands_to_rgb = sonyQuantumEfficiency(bands);
    % Normalize, for improved colour saturation
    bands_to_rgb = bands_to_rgb ./ max(max(bands_to_rgb));
end

% ## Disk fitting
bayer_pattern = 'gbrg'; % Colour-filter pattern
cleanup_radius = 2; % Morphological operations radius for 'findAndFitDisks()'
k0 = 0.1; % `k0` argument of 'findAndFitDisks()'
findAndFitDisks_options.bright_disks = true;
findAndFitDisks_options.mask_as_threshold = false;
findAndFitDisks_options.group_channels = ~rgb_mode;
findAndFitDisks_options.area_outlier_threshold = 2;

% ## Dispersion model generation
dispersion_fieldname = 'center';
max_degree_xy = 12;
max_degree_lambda = 12;

% ## Debugging Flags
findAndFitDisksVerbose.verbose_disk_search = false;
findAndFitDisksVerbose.verbose_disk_refinement = false;
findAndFitDisksVerbose.display_final_centers = true;

statsToDisparityVerbose.display_raw_values = true;
statsToDisparityVerbose.display_raw_disparity = true;
statsToDisparityVerbose.filter = struct(...
    dispersion_fieldname, true...
);

xylambdaPolyfitVerbose = true;
plot_polynomial_model = true;
if plot_polynomial_model && ~rgb_mode
    n_lambda_plot = min(20, length(lens_params.wavelengths));
end

%% Process the images

n_images = length(partial_filepaths);
if ~rgb_mode
    n_wavelengths = length(bands);
    wavelength_postfixes = cell(n_wavelengths, 1);
    for j = 1:n_wavelengths
        wavelength_postfixes{j} = num2str(bands(j), '%d');
    end
end

centers_cell = cell(n_images, 1);
for i = 1:n_images
    % Find any mask
    mask_filename = strcat(partial_filepaths(i), {'_mask.'}, {mask_ext});
    mask_filename = mask_filename{1};
    mask_listing = dir(mask_filename);
    if isempty(mask_listing)
        mask = [];
    else
        mask = imbinarize(imread(mask_filename), mask_threshold);
    end
    
    if rgb_mode
        filenames = strcat(partial_filepaths(i), {'.'}, {ext});
    else
        filenames = strcat(partial_filepaths(i), wavelength_postfixes, {'.'}, {ext});
    end
    
    I = imread(filenames{1});
    centers_i = findAndFitDisks(...
        I, mask, bayer_pattern, [], cleanup_radius, k0,...
        findAndFitDisks_options, findAndFitDisksVerbose...
    );
    if ~rgb_mode
        if n_wavelengths > 1
            centers_i = [
                centers_i,...
                struct(dispersion_fieldname, cell(size(centers_i, 1), n_wavelengths - 1))...
            ]; %#ok<AGROW>
            for j = 2:n_wavelengths
                I = imread(filenames{j});
                centers_i(:, j) = findAndFitDisks(...
                    I, mask, bayer_pattern, [], cleanup_radius, k0,...
                    findAndFitDisks_options, findAndFitDisksVerbose...
                );
            end
        end
    end
    centers_cell{i} = centers_i;
end
        
%% Fit a dispersion model to the results

if rgb_mode
    % Centers are already associated between colour channels
    centers = cell2mat(centers_cell);
else
    % TODO
end

x_fields = struct(...
    dispersion_fieldname, dispersion_fieldname...
);

disparity_raw_disks = statsToDisparity(...
    centers, reference_index,...
    1, 0, x_fields, bands, bands_to_rgb, statsToDisparityVerbose...
);

if rgb_mode
    polyfun_disks = xylambdaPolyfit(...
        centers, dispersion_fieldname, max_degree_xy, disparity_raw_disks, dispersion_fieldname,...
        xylambdaPolyfitVerbose...
    );
else
    polyfun_disks = xylambdaPolyfit(...
        centers, dispersion_fieldname, max_degree_xy, disparity_raw_disks, dispersion_fieldname,...
        bands, max_degree_lambda, xylambdaPolyfitVerbose...
    );
end

%% Visualization

if plot_polynomial_model
    if rgb_mode
        plotXYLambdaPolyfit(...
            centers, dispersion_fieldname, disparity_raw_disks, dispersion_fieldname,...
            reference_index, polyfun_disks...
        );
    else
        plotXYLambdaPolyfit(...
            centers, dispersion_fieldname, disparity_raw_disks, dispersion_fieldname,...
            bands, bands(reference_index), n_lambda_plot, polyfun_disks...
        );
    end
end