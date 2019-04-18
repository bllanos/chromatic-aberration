%% Chromatic aberration calibration from captured images
% Obtain dispersion models from arbitrary multi-channel (spectral or
% demosaiced-colour) images by registering bands/colour channels.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% This script can create new models of dispersion from scratch, or start
% from existing models of dispersion to speed up convergence.
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
% - 'image_filenames': A cell vector of input image filenames retrieved
%   based on the wildcard provided in the parameters section of the script.
% - 'centers': The independent variables data used for fitting the
%   model of dispersion. `centers` is a structure array, with one field
%   containing the positions of the centres of image patches registered
%   between colour channels or spectral bands. `centers(i, k)` is the
%   centre of the i-th image patch, with k representing either the k-th
%   spectral band, or the k-th colour channel.
% - 'disparity': The dependent variables data used for fitting the
%   model of dispersion. `disparity` is the first output argument of
%   'statsToDisparity()', produced when 'statsToDisparity()' was called
%   with `centers` as one of its input arguments. The format of `disparity`
%   is described in the documentation of 'statsToDisparity()'. `disparity`
%   contains the dispersion vectors between the centres of registered
%   patches for different wavelength bands or colour channels.
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
% - The idea of using mutual information to evaluate image alignment is
%   mentioned in, among other articles,
%
%   Brauers, J., Schulte, B., & Aach, T. (2008). "Multispectral
%     Filter-Wheel Cameras: Geometric Distortion Model and Compensation
%     Algorithms." IEEE Transactions on Image Processing, 17(12),
%     2368-2380. doi:10.1109/TIP.2008.2006605

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 17, 2019

% List of parameters to save with results
parameters_list = {
        'input_images_variable_name',...
        'forward_dispersion_model_filename',...
        'bands_filename',...
        'bands_variable',...
        'bands',...
        'bands_to_rgb',...
        'reference_wavelength',...
        'reference_index',...
        'reg_patch_options',...
        'reg_optimizer',...
        'reg_metric',...
        'reg_pyramid_levels',...
        'distance_outlier_threshold',...
        'dispersion_fieldname',...
        'fill_image',...
        'max_degree_xy_dispersion',...
        'max_degree_lambda',...
        'spline_smoothing_options'...
    };

%% Input data and parameters

% ## Input images
%
% Wildcard for 'ls()' to find the images to process. All images are
% expected to be in one directory.
%
% Images can either be spectral images, or demosaiced colour images,
% although the latter type of images may not give good results because of
% demosaicing artifacts.
input_images_wildcard = '/home/graphicslab/Documents/llanos/Data/20190208_ComputarLens/dataset/channel_scaling/*disks32cm*_dHyper.mat';
input_images_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% Model of dispersion to use as an initial guess
% Can be empty
forward_dispersion_model_filename = '/home/graphicslab/Documents/llanos/Results/Copied elsewhere/20190208_ComputarLens/dispersion/spectral/full_image/RAWDiskDispersionResults_spectral_spline_fromReference.mat';

% ## Spectral information

% Find dispersion between colour channels, as opposed to between spectral
% bands
rgb_mode = false;

if rgb_mode
    n_bands = 3;
    bands = (1:n_bands).';
    reference_wavelength = []; % Not used
    reference_index = 2; % Green colour channel
    bands_to_rgb = eye(n_bands);
    bands_filename = []; % Not used
    bands_variable = []; % Not used
else
    % Path and filename of a '.mat' file containing the wavelengths corresponding to
    % the spectral image.
    bands_filename = '/home/graphicslab/Documents/llanos/Data/20190208_ComputarLens/dataset/channel_scaling/sensor.mat';
    bands_variable = 'bands'; % Variable name in the above file
    reference_wavelength = 587.6;
end

% ## Patch-wise image registration
% (Options for 'registerPatches()')

reg_patch_options = struct('patch_size', [64, 64], 'padding', 16);
% Useful for debugging
% reg_patch_options.target_patch = [1405, 271];

% Use mutual information as an image registration metric
[reg_optimizer, reg_metric] = imregconfig('multimodal');
reg_optimizer.MaximumIterations = 500;

% The `'PyramidLevels'` input argument of 'imregtform()'
reg_pyramid_levels = 3;

% ## Threshold number of standard deviations of distance used to reject
% dispersion vectors
distance_outlier_threshold = 3;

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

% ## Output directory
output_directory = '/home/graphicslab/Documents/llanos/Results/dispersion_registration';

% ## Debugging Flags
registerPatchesVerbose = true;

statsToDisparityVerbose.display_raw_values = true;
statsToDisparityVerbose.display_raw_disparity = true;
statsToDisparityVerbose.filter = struct(...
    dispersion_fieldname, true...
);

xylambdaFitVerbose = true;
plot_model = true;

%% Load calibration data

if ~rgb_mode
    load(bands_filename, bands_variable);
    if exist(bands_variable, 'var')
        bands = eval(bands_variable);
    end
    if ~exist(bands_variable, 'var') || isempty(bands)
        error('No wavelengths loaded.')
    end
    n_bands = length(bands);
    
    [~, reference_index] = min(abs(bands - reference_wavelength));
    
    if plot_model
        n_lambda_plot = min(20, n_bands);
    end
    
    % `bands_to_rgb` is used for visualization purposes only, and so does
    % not need to be accurate
    bands_to_rgb = sonyQuantumEfficiency(bands);
    % Normalize, for improved colour saturation
    bands_to_rgb = bands_to_rgb ./ max(max(bands_to_rgb));
end

has_dispersion = ~isempty(forward_dispersion_model_filename);
if has_dispersion
    [...
        dispersion_data_init, bands_dispersionfun, transform_data...
    ] = loadDispersionModel(forward_dispersion_model_filename, true);

    if rgb_mode && ((n_bands ~= length(bands_dispersionfun)) ||...
           any(n_bands(:) ~= bands_dispersionfun(:)))
            error('Unexpected colour channels used by the model of dispersion.');
    end
end

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Process the images

centers_cell = cell(n_images, 1);
image_size = [];
for g = 1:n_images    

    I = loadImage(image_filenames{g}, input_images_variable_name);
    if isempty(image_size)
        image_size = size(I);
    elseif any(image_size ~= size(I))
        error('Not all images have the same dimensions.');
    end
    
    if has_dispersion
        if g == 1
            dispersionfun_init = makeDispersionForImage(...
                dispersion_data_init, I, transform_data...
            );
        end
        centers_cell{g} = registerPatches(...
            I, reference_index, reg_patch_options, reg_optimizer, reg_metric,...
            reg_pyramid_levels, bands, dispersionfun_init, registerPatchesVerbose...
        );
    else
        centers_cell{g} = registerPatches(...
            I, reference_index, reg_patch_options, reg_optimizer, reg_metric,...
            reg_pyramid_levels, registerPatchesVerbose...
        );
    end
end

%% Fit dispersion models to the results

% Centers are already matched between colour channels/bands, but we can still
% filter out outlier matches
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
    'image_filenames',...
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
        filename = 'RegistrationDispersionResults';
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
