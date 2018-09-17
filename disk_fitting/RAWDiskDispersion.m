%% Chromatic aberration calibration from captured images
% Obtain a dispersion model from RAW images of disk calibration patterns.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% This script can either model dispersion between colour channels, or
% dispersion between images taken under different wavelength bandpass
% filters.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
%
% ### Graphical output from 'plotXYLambdaModel()'
% - Displayed if `plot_model` is `true`.
%
% ### Model fitting results
%
% A '.mat' file containing the following variables:
%
% - 'centers': The independent variables data used for fitting the
%   model of dispersion. `centers` is a structure array, with one field
%   containing the image positions of the centres of disks fitted to image
%   blobs. `centers(i, k)` is the centre of the ellipse fitted to the i-th
%   image blob, with k representing eitherthe k-th bandpass filter, or the
%   k-th colour channel (depending on `rgb_mode`).
% - 'disparity': The dependent variables data used for fitting the
%   model of dispersion. `disparity` is the first output argument of
%   'statsToDisparity()', produced when 'statsToDisparity()' was called
%   with `centers` as one of its input arguments. The format of `disparity`
%   is described in the documentation of 'statsToDisparity()'. `disparity`
%   contains the dispersion vectors between the centres of disks fit to
%   image blobs for different wavelength bands or colour channels.
% - 'dispersion_data': The model of dispersion, modeling the mapping from
%   `centers` to `disparity`. `dispersion_data` can be converted to a
%   function form using `dispersionfun = makeDispersionfun(dispersion_data)`
% - 'model_space': A structure describing the range of image coordinates
%   over which the model of dispersion is valid, having the following
%   fields:
%   - 'corners': The first and second rows contain the (x,y) image
%     coordinates of the top left and bottom right corners of the region,
%     respectively. Remember that image coordinates are 0.5 units offset
%     from pixel indices.
%   - 'image_size': A two-element vector containing the image height and
%     width in pixels.
%   - 'system': A character vector, 'image', indicating that the
%     dispersion model was constructed under image coordinate conventions,
%     wherein the y-axis is positive downards on the image plane, and the
%     origin is the top left corner of the image.
%
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.) Note that the set of parameters contains the `bands` variable
% also output by 'DoubleConvexThickLensDispersion.m' and
% 'DoubleConvexThickLensDiskDispersion.m'.
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

% List of parameters to save with results
parameters_list = {
        'partial_filepaths',...
        'ext',...
        'mask_ext',...
        'mask_threshold',...
        'rgb_mode',...
        'bands',...
        'reference_index',...
        'bands_to_rgb',...
        'bayer_pattern',...
        'cleanup_radius',...
        'k0',...
        'findAndFitDisks_options',...
        'dispersion_fieldname',...
        'max_degree_xy',...
        'max_degree_lambda',...
        'spline_smoothing_weight',...
        'model_from_reference',...
        'model_type'...
    };

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
    '/home/llanos/GoogleDrive/ThesisResearch/Results/20170808_OpticalTableMatrix/averaged/d44_a22_far_disksWhite'
    };

% Filename extension (excluding the leading '.')
ext = 'tif';
% Mask filename extension
mask_ext = 'png';

% Threshold used to binarize mask images
mask_threshold = 0.5; % In a range of intensities from 0 to 1

% ## Spectral information

% Find dispersion between colour channels, as opposed to between images
% taken under different spectral bands
rgb_mode = true;

% Note that `bands_to_rgb` is used only for visualization purposes, not for
% calculations. It does not need to be accurate.
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

% If `true`, model dispersion between bands (colour channels or spectral
% bands) as a function of positions in the reference band. If `false`,
% model dispersion as a function of positions in the non-reference bands.
% The first case is useful for warping the other bands to align with the
% reference band, such as when correcting chromatic aberration by image
% warping. The second case is useful for warping an "ideal" image to
% compare it with an observed aberrated image. In both cases, the
% dispersion vectors point from the reference band to the other bands.
model_from_reference = false;

% Spline, or global polynomial models can be fitted
model_type = 'spline';

% Parameters for polynomial model fitting
max_degree_xy = 5;
max_degree_lambda = 5;

% Parameters for spline model fitting
spline_smoothing_options = struct(...
    'n_iter', [20, 50],...
    'grid_size', [10, 4],...
    'minimum', eps,...
    'maximum', 1e10,...
    'tol', 1e-6 ...
);

% ## Debugging Flags
findAndFitDisksVerbose.verbose_disk_search = false;
findAndFitDisksVerbose.verbose_disk_refinement = false;
findAndFitDisksVerbose.display_final_centers = true;

statsToDisparityVerbose.display_raw_values = false;
statsToDisparityVerbose.display_raw_disparity = true;
statsToDisparityVerbose.filter = struct(...
    dispersion_fieldname, true...
);

xylambdaFitVerbose = true;
plot_model = true;
if plot_model && ~rgb_mode
    n_lambda_plot = min(20, length(bands));
end

%% Process the images

n_images = length(partial_filepaths);
n_bands = length(bands);
if ~rgb_mode
    wavelength_postfixes = cell(n_bands, 1);
    for j = 1:n_bands
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
        mask = imread(mask_filename);
        if size(mask, 3) ~= 1
            error('Expected the mask, %s, to have only one channel.')
        end
        mask = imbinarize(imread(mask_filename), mask_threshold);
    end
    
    if rgb_mode
        filenames = strcat(partial_filepaths(i), {'.'}, {ext});
    else
        filenames = strcat(partial_filepaths(i), wavelength_postfixes, {'.'}, {ext});
    end
    
    if rgb_mode
        I = imread(filenames{1});
        if i == 1
            image_size = size(I);
        elseif any(image_size ~= size(I))
            error('Not all images have the same dimensions.');
        end
        centers_i = findAndFitDisks(...
            I, mask, bayer_pattern, [], cleanup_radius, k0,...
            findAndFitDisks_options, findAndFitDisksVerbose...
        );
    else
        centers_i = cell(n_bands, 1);
        for j = 1:n_bands
            I = imread(filenames{j});
            if i == 1 && j == 1
                image_size = size(I);
            elseif any(image_size ~= size(I))
                error('Not all images have the same dimensions.');
            end
            centers_i{j} = findAndFitDisks(...
                I, mask, bayer_pattern, [], cleanup_radius, k0,...
                findAndFitDisks_options, findAndFitDisksVerbose...
            );
        end
    end
    centers_cell{i} = centers_i;
end
        
%% Fit a dispersion model to the results

if rgb_mode
    % Centers are already associated between colour channels
    centers = cell2mat(centers_cell);
else
    centers = matchByVectors(centers_cell, dispersion_fieldname, reference_index);
end

x_fields = struct(...
    dispersion_fieldname, dispersion_fieldname...
);

disparity = statsToDisparity(...
    centers, reference_index,...
    1, 0, x_fields, bands, bands_to_rgb, statsToDisparityVerbose...
);

if model_from_reference
    centers_for_fitting = repmat(centers(:, reference_index), 1, n_bands);
else
    centers_for_fitting = centers;
end

if strcmp(model_type, 'polynomial')
    if rgb_mode
        [ dispersionfun, dispersion_data ] = xylambdaPolyfit(...
            centers_for_fitting, dispersion_fieldname, max_degree_xy, disparity,...
            dispersion_fieldname, xylambdaFitVerbose...
        );
    else
        [ dispersionfun, dispersion_data ] = xylambdaPolyfit(...
            centers_for_fitting, dispersion_fieldname, max_degree_xy, disparity,...
            dispersion_fieldname, bands, max_degree_lambda, xylambdaFitVerbose...
        );
    end
elseif strcmp(model_type, 'spline')
    if rgb_mode
        [ dispersionfun, dispersion_data ] = xylambdaSplinefit(...
            centers_for_fitting, dispersion_fieldname, disparity,...
            dispersion_fieldname, spline_smoothing_options, xylambdaFitVerbose...
        );
    else
        [ dispersionfun, dispersion_data ] = xylambdaPolyfit(...
            centers_for_fitting, dispersion_fieldname, disparity,...
            dispersion_fieldname, spline_smoothing_options, bands, xylambdaFitVerbose...
        );
    end
else
    error('Unrecognized value of the `model_type` paramter.');
end

%% Visualization

if plot_model
    if rgb_mode
        plotXYLambdaModel(...
            centers_for_fitting, dispersion_fieldname, disparity, dispersion_fieldname,...
            reference_index, dispersionfun...
        );
    else
        plotXYLambdaModel(...
            centers_for_fitting, dispersion_fieldname, disparity, dispersion_fieldname,...
            bands, bands(reference_index), n_lambda_plot, dispersionfun...
        );
    end
end

%% Save results to a file

% Indicate where in the image the model is usable
centers_unpacked = permute(reshape([centers.(dispersion_fieldname)], 2, []), [2 1]);
model_space.corners = [
    min(centers_unpacked(:, 1)), min(centers_unpacked(:, 2));
    max(centers_unpacked(:, 1)), max(centers_unpacked(:, 2))
    ];
model_space.image_size = image_size;
model_space.system = 'image';

save_variables_list = [ parameters_list, {...
        'centers',...
        'disparity',...
        'dispersion_data',...
        'model_space'...
    } ];
uisave(save_variables_list,'RAWDiskDispersionResults');