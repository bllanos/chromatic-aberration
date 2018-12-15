%% Chromatic aberration calibration from captured images
% Obtain dispersion models from RAW images of disk calibration patterns.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% This script can either model dispersion between colour channels, or
% dispersion between images taken under different optical filters.
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
% Four '.mat' files, each containing the following variables:
%
% - 'grouped_filenames': A cell vector of cell vectors of input image
%   filenames retrieved based on the wildcard provided in the parameters
%   section of the script. Each of the inner cell vectors groups the images
%   of the same scene taken under different spectral filters, if spectral
%   dispersion is being modelled. Otherwise, the inner cell vectors are of
%   length one.
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
% - 'model_from_reference': If `true`, dispersion is modelled between bands
%   (colour channels or spectral bands) as a function of positions in the
%   reference band. If `false`, dispersion is modelled as a function of
%   positions in the non-reference bands. The first case is useful for
%   warping the other bands to align with the reference band, such as when
%   correcting chromatic aberration by image warping. The second case is
%   useful for warping an "ideal" image to compare it with an observed
%   aberrated image. In both cases, the dispersion vectors point from the
%   reference band to the other bands.
% - 'model_type': The type of model of dispersion, either 'spline', or
%   'polynomial'. Spline models of dispersion are generated using
%   'xylambdaSplinefit()', whereas polynomial models of dispersion are
%   generated using 'xylambdaPolyfit()'.
%
% One '.mat' file is generated for each possible combination of of
% 'model_from_reference' and 'model_type'. Therefore, only the model of
% dispersion differs between the files. The '.mat' files will be named
% based on the models of dispersion that they contain.
%
% Additionally, the files contain the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.) Note that the set of parameters contains the `bands` variable,
% also output by 'DoubleConvexThickLensDispersion.m'.
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
        'mask_ext',...
        'mask_threshold',...
        'rgb_mode',...
        'bands_regex',...
        'bands',...
        'reference_wavelength',...
        'reference_index',...
        'bands_to_rgb',...
        'bayer_pattern',...
        'cleanup_radius',...
        'k0',...
        'findAndFitDisks_options',...
        'dispersion_fieldname',...
        'max_degree_xy',...
        'max_degree_lambda',...
        'spline_smoothing_options'...
    };

%% Input data and parameters

% ## Input images
%
% Wildcard for 'ls()' to find the images to process. All images are
% expected to be in one directory.
%
% Images are expected to have been preprocessed, such as using
% 'PreprocessRAWImages.m', so that they do not need to be linearized after
% being loaded. For image format files, images will simply be loaded with
% the Image Processing Toolbox 'imread()' function. For '.mat' files, the
% variable to be loaded must be provided in the script parameters.
%
% RAW (non-demosaicked) images are expected.
%
% This script will also search for image masks, using filenames which are
% constructed by appending '_mask' and `mask_ext` (below) to the filepaths
% (after stripping file extensions and wavelength information). Masks are
% used to avoid processing irrelevant portions of images.
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/preprocessed/blended_averaged/pos1_1mmDots_*.mat';
input_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Mask filename extension
mask_ext = 'png';

% Threshold used to binarize mask images, if they are not already binary
% images.
mask_threshold = 0.5; % In a range of intensities from 0 to 1

% ## Spectral information

% Find dispersion between colour channels, as opposed to between images
% taken under different spectral bands
rgb_mode = false;

if rgb_mode
    bands_regex = []; % Not used
    bands = 1:3;
    reference_wavelength = []; % Not used
    reference_index = 2; % Green colour channel
    bands_to_rgb = eye(3);
else
    % Wavelengths will be expected within filenames, extracted using this
    % regular expression
    bands_regex = '_(\d+)nm';
    reference_wavelength = 587.6;
end

% ## Disk fitting
bayer_pattern = 'gbrg'; % Colour-filter pattern
cleanup_radius = 2; % Morphological operations radius for 'findAndFitDisks()'
k0 = 0.1; % `k0` argument of 'findAndFitDisks()'
findAndFitDisks_options.bright_disks = false;
findAndFitDisks_options.mask_as_threshold = false;
findAndFitDisks_options.group_channels = ~rgb_mode;
findAndFitDisks_options.area_outlier_threshold = 2;

% ## Dispersion model generation
dispersion_fieldname = 'center';

% Parameters for polynomial model fitting
max_degree_xy = 12;
max_degree_lambda = 12;

% Parameters for spline model fitting
spline_smoothing_options = struct(...
    'n_iter', [20, 50],...
    'grid_size', [15, 4],...
    'minimum', eps,...
    'maximum', 1e10,...
    'tol', 1e-6 ...
);

% ## Output directory
output_directory = '/home/llanos/Downloads';

% ## Debugging Flags
findAndFitDisksVerbose.verbose_disk_search = true;
findAndFitDisksVerbose.verbose_disk_refinement = false;
findAndFitDisksVerbose.display_final_centers = true;

statsToDisparityVerbose.display_raw_values = true;
statsToDisparityVerbose.display_raw_disparity = true;
statsToDisparityVerbose.filter = struct(...
    dispersion_fieldname, true...
);

xylambdaFitVerbose = true;
plot_model = true;

%% Find the images

if rgb_mode
    [...
        grouped_filenames, path, group_names...
    ] = findAndGroupImages(input_images_wildcard);

else
    [...
        grouped_filenames, path, group_names, bands...
    ] = findAndGroupImages(input_images_wildcard, bands_regex);

    [~, reference_index] = min(abs(bands - reference_wavelength));
    
    % `bands_to_rgb` is used for visualization purposes only, and so does
    % not need to be accurate
    bands_to_rgb = sonyQuantumEfficiency(bands);
    % Normalize, for improved colour saturation
    bands_to_rgb = bands_to_rgb ./ max(max(bands_to_rgb));

    if plot_model
        n_lambda_plot = min(20, length(bands));
    end
end
n_groups = length(grouped_filenames);
n_bands = length(bands);

%% Process the images

centers_cell = cell(n_groups, 1);
image_size = [];
for g = 1:n_groups
    
    % Find any mask
    mask_filename = fullfile(path, [group_names{g}, '_mask.', mask_ext]);
    mask_listing = dir(mask_filename);
    if isempty(mask_listing)
        mask = [];
    else
        mask = imread(mask_filename);
        if size(mask, 3) ~= 1
            error('Expected the mask, "%s", to have only one channel.', mask_filename);
        end
        if ~islogical(mask)
            mask = imbinarize(mask, mask_threshold);
        end
    end
    
    I = loadImage(grouped_filenames{g}{1}, input_images_variable_name);
    if isempty(image_size)
        image_size = size(I);
    elseif any(image_size ~= size(I))
        error('Not all images have the same dimensions.');
    end
    
    if rgb_mode
        centers_cell{g} = findAndFitDisks(...
            I, mask, bayer_pattern, [], cleanup_radius, k0,...
            findAndFitDisks_options, findAndFitDisksVerbose...
        );
    else
        centers_g = cell(n_bands, 1);
        centers_g{1} = findAndFitDisks(...
            I, mask, bayer_pattern, [], cleanup_radius, k0,...
            findAndFitDisks_options, findAndFitDisksVerbose...
        );

        for i = 2:n_bands
            I = loadImage(grouped_filenames{g}{i}, input_images_variable_name);
            if any(image_size ~= size(I))
                error('Not all images have the same dimensions.');
            end
            centers_g{i} = findAndFitDisks(...
                I, mask, bayer_pattern, [], cleanup_radius, k0,...
                findAndFitDisks_options, findAndFitDisksVerbose...
            );
        end
        
        centers_cell{g} = centers_g;
    end     
end
        
%% Fit dispersion models to the results

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

% Indicate where in the image the model is usable
centers_unpacked = permute(reshape([centers.(dispersion_fieldname)], 2, []), [2 1]);
model_space.corners = [
    min(centers_unpacked(:, 1)), min(centers_unpacked(:, 2));
    max(centers_unpacked(:, 1)), max(centers_unpacked(:, 2))
    ];
model_space.corners = max(model_space.corners, 0.5);
model_space.corners(model_space.corners(:, 1) > (image_size(2) - 0.5), 1) = (image_size(2) - 0.5);
model_space.corners(model_space.corners(:, 2) > (image_size(1) - 0.5), 2) = (image_size(1) - 0.5);
model_space.image_size = image_size;
model_space.system = 'image';

save_variables_list = [ parameters_list, {...
    'grouped_filenames',...
    'centers',...
    'disparity',...
    'dispersion_data',...
    'model_space',...
    'model_from_reference',...
    'model_type'...
} ];

for model_from_reference = [true, false]
    if model_from_reference
        centers_for_fitting = repmat(centers(:, reference_index), 1, n_bands);
    else
        centers_for_fitting = centers;
    end
        
    for model_type_cell = {'spline', 'polynomial'}
        model_type = model_type_cell{1};
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
                [ dispersionfun, dispersion_data ] = xylambdaSplinefit(...
                    centers_for_fitting, dispersion_fieldname, disparity,...
                    dispersion_fieldname, spline_smoothing_options, bands, xylambdaFitVerbose...
                );
            end
        else
            error('Unrecognized value of `model_type`.');
        end

        % Visualization
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

        % Save results to a file
        filename = 'RAWDiskDispersionResults';
        if rgb_mode
            filename = [filename, '_RGB_'];
        else
            filename = [filename, '_spectral_'];
        end
        filename = [filename, model_type];
        if model_from_reference
            filename = [filename, '_fromReference'];
        else
            filename = [filename, '_fromNonReference'];
        end
        filename = [filename, '.mat'];
            
        save_data_filename = fullfile(output_directory, filename);
        save(save_data_filename, save_variables_list{:});
    end
end