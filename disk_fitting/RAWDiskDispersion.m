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
% ### Vignetting-corrected images
% 
% Image files with names of the form '*_vignettingCorrected.tif', where '*'
% is the name of the original input image, are saved if
% `save_corrected_images` is `true`, and if masks for calibrating
% vignetting are found with the input images. The images are versions of
% the input images which have been corrected for vignetting, and are used
% to compute dispersion models.
%
% ### Graphical output from 'plotXYLambdaModel()'
% - Displayed if `plot_model` is `true`.
%
% ### Model fitting results
%
% Up to four '.mat' files, each containing the following variables:
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
%   image blob, with k representing either the k-th bandpass filter, or the
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
% One '.mat' file is generated for each possible combination of
% 'model_from_reference' and 'model_type' specified in the script parameters.
% Therefore, only the model of dispersion differs between the files. The '.mat'
% files will be named based on the models of dispersion that they contain.
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
% - Rudakova, V. & Monasse, P. (2014). "Precise correction of lateral
%   chromatic aberration in images" (Guanajuato). 6th Pacific-Rim Symposium
%   on Image and Video Technology, PSIVT 2013. Springer Verlag.
%   doi:10.1007/978-3-642-53842-1_2

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
        'distance_outlier_threshold',...
        'bands_to_rgb',...
        'bayer_pattern',...
        'cleanup_radius',...
        'k0',...
        'findAndFitDisks_options',...
        'dispersion_fieldname',...
        'fill_image',...
        'max_degree_xy_dispersion',...
        'max_degree_lambda',...
        'spline_smoothing_options',...
        'max_degree_xy_vignetting',...
        'quantiles',...
        'model_type_choices',...
        'model_from_reference_choices'...
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
% constructed by appending '_maskDisks' or '_maskVignetting', and
% `mask_ext` (below) to the filepaths (after stripping file extensions and
% wavelength information).
% - Masks with filenames containing '_maskVignetting' define regions for
%   calibrating vignetting corrections to be applied to the image prior to
%   calibrating models of dispersion. If no mask is found for an image, the
%   image will not be corrected for vignetting.
% - Masks with filenames containing '_maskDisks' are used to avoid
%   processing irrelevant portions of images when calibrating models of
%   dispersion. If no mask is found for an image, the entire image will be
%   searched for disks for dispersion calibration.
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/preprocessed_images/exposure_blended/*disks*nm.mat';
input_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Mask filename extension (without the '.')
mask_ext = 'png';

% Threshold used to binarize mask images, if they are not already binary
% images.
mask_threshold = 0.5; % In a range of intensities from 0 to 1

bayer_pattern = 'gbrg'; % Colour-filter pattern

% ## Spectral information

% Find dispersion between colour channels, as opposed to between images
% taken under different spectral bands
rgb_mode = false;

if rgb_mode
    bands_regex = []; % Not used
    n_channels_rgb = 3;
    bands = (1:n_channels_rgb).';
    reference_wavelength = []; % Not used
    reference_index = 2; % Green colour channel
    bands_to_rgb = eye(n_channels_rgb);
else
    % Wavelengths will be expected within filenames, extracted using this
    % regular expression
    bands_regex = '_(\d+)nm';
    reference_wavelength = 587.6;
end
% Threshold number of standard deviations of distance used to reject matches
% between disks
distance_outlier_threshold = 3;

% ## Vignetting correction

% Parameters for polynomial model fitting
max_degree_xy_vignetting = 5;

% Quantiles used for clipping to produce nice output images (for display,
% not for calculation)
quantiles = [0.01, 0.99];

% ## Disk fitting
cleanup_radius = 2; % Morphological operations radius for 'findAndFitDisks()'
k0 = 0.5; % `k0` argument of 'findAndFitDisks()'
findAndFitDisks_options.bright_disks = false;
findAndFitDisks_options.mask_as_threshold = false;
findAndFitDisks_options.group_channels = ~rgb_mode;
findAndFitDisks_options.area_outlier_threshold = 3;

% ## Dispersion model generation
dispersion_fieldname = 'center';

% Force the dispersion model to declare that it is valid over the entire
% image?
fill_image = true;

% Parameters for polynomial model fitting
max_degree_xy_dispersion = 12;
max_degree_lambda = 12;

% Parameters for spline model fitting
spline_smoothing_options = struct(...
    'n_iter', [20, 50],...
    'grid_size', [15, 4],...
    'minimum', eps,...
    'maximum', 1e10,...
    'tol', 1e-6 ...
);

% Which models of dispersion to generate?
model_type_choices = {'spline', 'polynomial'};
model_from_reference_choices = [true, false];

% ## Output directory
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion_newCV_vignettingCorrected/spline';

% ## Debugging Flags
vignettingPolyfitVerbose = false;

findAndFitDisksVerbose.verbose_disk_search = false;
findAndFitDisksVerbose.verbose_disk_refinement = false;
findAndFitDisksVerbose.display_final_centers = true;

statsToDisparityVerbose.display_raw_values = true;
statsToDisparityVerbose.display_raw_disparity = true;
statsToDisparityVerbose.filter = struct(...
    dispersion_fieldname, true...
);

save_corrected_images = true; % Images corrected for vignetting
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
    
    % Find any mask for vignetting calibration
    mask_filename = fullfile(path, [group_names{g}, '_maskVignetting.', mask_ext]);
    mask_listing = dir(mask_filename);
    use_vignetting_correction = ~isempty(mask_listing);
    if use_vignetting_correction
        vignetting_mask = imread(mask_filename);
        if size(vignetting_mask, 3) ~= 1
            error('Expected the vignetting mask, "%s", to have only one channel.', mask_filename);
        end
        if ~islogical(vignetting_mask)
            vignetting_mask = imbinarize(vignetting_mask, mask_threshold);
        end
    end

    % Find any mask for dispersion calibration
    mask_filename = fullfile(path, [group_names{g}, '_maskDisks.', mask_ext]);
    mask_listing = dir(mask_filename);
    if isempty(mask_listing)
        mask = [];
    else
        mask = imread(mask_filename);
        if size(mask, 3) ~= 1
            error('Expected the dispersion calibration mask, "%s", to have only one channel.', mask_filename);
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
    
    % Vignetting correction
    if use_vignetting_correction
        vignettingfun = vignettingPolyfit(...
            I, vignetting_mask, max_degree_xy_vignetting, bayer_pattern, vignettingPolyfitVerbose...
        );
        I = correctVignetting(I, vignettingfun);
        if save_corrected_images
            I_out_debug = clipAndRemap(I, 'uint8', 'quantiles', quantiles);
            [~, I_out_filename] = fileparts(grouped_filenames{g}{1});
            saveImages(...
                'image', output_directory, I_out_filename,...
                I_out_debug, '_vignettingCorrected', []...
                );
        end
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
            
            if use_vignetting_correction
                vignettingfun = vignettingPolyfit(...
                    I, vignetting_mask, max_degree_xy_vignetting, bayer_pattern, vignettingPolyfitVerbose...
                    );
                I = correctVignetting(I, vignettingfun);
                if save_corrected_images
                    I_out_debug = clipAndRemap(I, 'uint8', 'quantiles', quantiles);
                    [~, I_out_filename] = fileparts(grouped_filenames{g}{i});
                    saveImages(...
                        'image', output_directory, I_out_filename,...
                        I_out_debug, '_vignettingCorrected', []...
                        );
                end
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
    % Centers are already matched between colour channels, but we can still
    % filter out outlier matches
    for g = 1:n_groups
        centers_cell{g} = mat2cell(centers_cell{g}, size(centers_cell{g}, 1), ones(1, n_channels_rgb)).';
    end
end
centers = matchByVectors(centers_cell, dispersion_fieldname, reference_index, distance_outlier_threshold);

x_fields = struct(...
    dispersion_fieldname, dispersion_fieldname...
);

disparity = statsToDisparity(...
    centers, reference_index,...
    1, 0, x_fields, bands, bands_to_rgb, statsToDisparityVerbose...
);

% Indicate where in the image the model is usable
if fill_image
    model_space.corners = [
        -Inf, -Inf;
        Inf, Inf
        ];
else
    centers_unpacked = permute(reshape([centers.(dispersion_fieldname)], 2, []), [2 1]);
    model_space.corners = [
        min(centers_unpacked(:, 1)), min(centers_unpacked(:, 2));
        max(centers_unpacked(:, 1)), max(centers_unpacked(:, 2))
        ];
end
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

for model_from_reference = model_from_reference_choices
    if model_from_reference
        centers_for_fitting = repmat(centers(:, reference_index), 1, n_bands);
    else
        centers_for_fitting = centers;
    end

    for model_type_cell = model_type_choices
        model_type = model_type_cell{1};
        if strcmp(model_type, 'polynomial')
            if rgb_mode
                [ dispersionfun, dispersion_data ] = xylambdaPolyfit(...
                    centers_for_fitting, dispersion_fieldname, max_degree_xy_dispersion, disparity,...
                    dispersion_fieldname, xylambdaFitVerbose...
                );
            else
                [ dispersionfun, dispersion_data ] = xylambdaPolyfit(...
                    centers_for_fitting, dispersion_fieldname, max_degree_xy_dispersion, disparity,...
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
