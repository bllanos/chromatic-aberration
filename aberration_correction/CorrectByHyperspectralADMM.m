%% Demosaicing and hyperspectral ADMM-based correction of chromatic aberration
% Convert RAW images to colour images, and correct chromatic aberration by
% estimating a latent hyperspectral or RGB image.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
% A set of RAW images to be demosaiced and corrected for chromatic
% aberration.
%
% Images are expected to have been preprocessed, such as using
% 'AverageRAWImages.m', so that they do not need to be linearized after
% being loaded.  For image format files, images will simply be loaded with
% the Image Processing Toolbox 'imread()' function. For '.mat' files, the
% variable to be loaded must be provided in the script parameters.
%
% All images are expected to have 3 colour channels (Red, Green, Blue)
% (represented in a Bayer pattern as a 2D array). However, the colour
% channels can correspond to narrowband wavelength ranges - This script
% will input a mapping from the colour space of the latent images to the
% colour space of the RAW images.
%
% The images need not have the same pixel dimensions, but they should be
% compatible with the input model of dispersion described below.
%
% ### Model of dispersion
%
% A '.mat' file containing several variables, which is the output of
% 'DoubleConvexThickLensDiskDispersion.m', 'RAWDiskDispersion.m',
% 'DoubleConvexThickLensDispersion.m' or 'BimaterialImages.m', for example.
% The following variables are required:
% - 'dispersion_data': A model of chromatic aberration, modeling the warping
%   from the reference colour channel or wavelength band to the other
%   colour channels or wavelength bands. `dispersion_data` can be converted to
%   a function form using `dispersionfun = makeDispersionfun(dispersion_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `false`.
% The following variables are sometimes required:
% - 'bands': A vector containing the wavelengths or colour channel indices
%   to use as the `lambda` input argument of 'dispersionfunToMatrix()'.
%   `bands` is the wavelength or colour channel information needed to
%   evaluate the dispersion model. This variable is required only if not
%   provided in the colour space conversion data file, or directly in this
%   script (see below).
%
% The following two additional variables are optional. If they are present,
% they will be used for the following purposes:
% - Conversion between the coordinate system in which the model of chromatic
%   aberration was constructed and the image coordinate system.
% - Limiting the correction of chromatic aberration to the region in which
%   the model is valid.
% The first variable, 'model_space' is a structure with same form as the
% `model_space` input argument of 'modelSpaceTransform()'. The second
% variable, `fill`, can be omitted, in which case it defaults to `false`.
% `fill` corresponds to the `fill` input argument of
% 'modelSpaceTransform()'. Refer to the documentation of
% 'modelSpaceTransform.m' for details.
%
% ### Colour space conversion data
% A '.mat' file containing several variables, which is the output of
% 'SonyColorMap.m', for example. The following variables are required:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel or spectral band in the input images to the
%   j-th colour channel or spectral band of the latent images. For example,
%   `sensor_map` is a matrix mapping discretized spectral power
%   distributions to RGB colours.
% - 'channel_mode': A Boolean value indicating whether the latent colour
%   space is a set of colour channels (true) or a set of spectral bands
%   (false).
% The following variables are sometimes required:
% - 'bands': A vector containing the wavelengths or colour channel indices
%   to use as the `lambda` input argument of 'dispersionfunToMatrix()' when
%   modelling the chromatic aberration applied to the latent images, in the
%   forward model of the input images. `bands` is the wavelength or colour
%   channel information needed to evaluate the dispersion model. This
%   variable takes precedence over the same variable provided with the
%   dispersion model (see above), but can be overridden by a version of the
%   variable provided directly in this script (see below).
%
% ### Latent colour space
%
% The 'bands' variable defined in the parameters section of the code below
% can be empty (`[]`), in which case it is loaded from the dispersion model
% data or from the colour space conversion data (see above).
%
% The final value of 'bands' is determined according to the following list,
% in order by decreasing priority:
% - 'bands' defined in this script
% - 'bands' loaded from colour space conversion data
% - 'bands' loaded with the model of dispersion
%
% This script assumes that all of the above sources are compatible with
% respect to the semantics of 'bands'. For instance, they are all using
% colour channel indices, or they are all using wavelength values.
%
% If the final value of 'bands' is not identical to the value of 'bands'
% loaded from the colour space conversion data, the 'sensor_map' variable
% from the colour space conversion data must be resampled along its second
% dimension to correspond to the final value of 'bands'. Note that
% resampling only makes sense when 'bands' is a set of wavelengths, and it
% is not possible for this script to programmatically validate that this is
% the case.
%
% If no variable 'bands' was found in the colour space conversion data, but
% the final value of 'bands' has a length equal to the size of the second
% dimension of 'sensor_map', 'bands' is assumed to be compatible with
% 'sensor_map'. Otherwise, 'sensor_map' is resampled along its second
% dimension, by assuming that its first and last columns correspond to the
% first and last elements of 'bands', and that its remaining columns
% correspond to equally-spaced wavelengths in-between the first and last
% elements of 'bands'.
%
% ## Output
%
% ### Estimated images
%
% One of the following types of images is created for each input image,
% except where specified. The filename of the input image, concatenated
% with a string of parameter information, is represented by '*' below.
% - '*_roi.tif' and '*_roi.mat': A cropped version of the input image
%   (stored in the variable 'I_raw'), containing the portion used as input
%   for ADMM. This region of interest was determined using the
%   `model_space` and `fill` variables saved in the input model of
%   dispersion data file (see above). If these variables were not present,
%   the cropped region is the entire input image. All of the other output
%   images listed below are limited to the region shown in '*_roi.tif'.
% - '*_latent.tif' and '*_latent.mat': The latent image estimated using
%   ADMM (stored in the variable 'I_latent'). The '.tif' image is only
%   output if `save_latent_image_files` is `true`, and an error will be
%   thrown if the latent images are not greyscale or 3-channel images.
% - '*_latent_rgb.tif' and '*_latent_rgb.mat': A colour image (stored in
%   the variable 'I_rgb') created by converting the latent image to the RGB
%   colour space of the camera.
% - '*_latent_warped.tif' and '*_latent_warped.mat': A colour image (stored in
%   the variable 'J_full') created by warping the latent image according to
%   the dispersion model, then converting the image to the RGB colour space
%   of the camera. This output image is, in a sense, a demosaiced version
%   of the input image.
% - '*_reestimated.tif' and '*_reestimated.mat': A simulation (stored in
%   the variable 'J_est') of the input RAW image from the latent image,
%   useful for visually evaluating the convergence of the ADMM algorithm.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'bands': The final value of the 'bands' variable, determined as
%   discussed under the section "Latent colour space" above.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference. 'bands_color' is empty if no
%   such variable was found in the data file, or if the value of the loaded
%   variable was empty.
% - 'bands_dispersionfun': The 'bands' variable loaded from the dispersion
%   model data file, for reference. 'bands_dispersionfun' is empty if no
%   such variable was found in the data file, or if the value of the loaded
%   variable was empty.
% - 'image_bounds': A cell vector, where the i-th cell stores the
%   coordinates of the i-th latent image in the space of the i-th cropped
%   input image. This is the 'image_bounds' output argument of
%   'dispersionfunToMatrix()'.
% - 'image_filenames': A cell vector containing the input image filenames
%   retrieved based on the wildcard provided in the parameters section of
%   the script.
% - 'sensor_map_resampled': The resampled version of the 'sensor_map'
%   variable, determined as discussed under the section "Latent colour
%   space" above.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## Notes
% - The image colour space is not altered by this script; RGB images are
%   produced in the camera's colour space. See 'imreadRAW()' for code to
%   convert an image to sRGB after demosaicing.
% - This script does not distinguish between wavelength bands and colour
%   channels. One can use this script to estimate either a latent
%   hyperspectral image, or a latent aberration-free RGB image (free from
%   lateral chromatic aberration). The latter use case is a baseline that
%   can be compared with the results of 'CorrectByWarping.m'. A latent
%   hyperspectral image can be sharper, in theory, whereas a latent RGB
%   image will retain the within-channel chromatic aberration of the input
%   image. The reason for this difference is the summation of multiple
%   spectral bands into each channel of an RGB image, in contrast to the
%   identity mapping of the colours of a latent RGB image into the colours
%   of the aberrated RGB image. Summation allows multiple sharp bands to
%   form a blurred colour channel.
%
% ## References
% - Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%   single-shot hyperspectral imaging using a prism." ACM Transactions
%   on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%   doi:10.1145/3130800.3130896
% - Boyd, S, et al.. "Distributed Optimization and Statistical Learning via
%   the Alternating Direction Method of Multipliers." Foundations and
%   Trends in Machine Learning, vol. 3, no. 1, pp. 1-122, 2011.
%   doi:10.1561/2200000016

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 25, 2018

% List of parameters to save with results
parameters_list = {
        'bayer_pattern',...
        'dispersion_model_filename',...
        'color_map_filename',...
        'bands_script',...
        'bands_interp_method',...
        'downsampling_factor',...
        'output_directory',...
        'save_latent_image_files',...
        'add_border',...
        'baek2017Algorithm2Options',...
        'rho',...
        'weights',...
        'int_method',...
        'patch_sizes',...
        'paddings',...
        'target_patch',...
        'run_entire_image'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
% '.mat' or image files can be loaded
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180709_TestingSplineModels/ground_truth/splines/swirly_0139_raw_warped.mat';
input_images_variable_name = 'raw_2D'; % Used only when loading '.mat' files

% Colour-filter pattern
bayer_pattern = 'gbrg';

% Model of dispersion
dispersion_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180709_TestingSplineModels/ground_truth/splines/BimaterialImagesData_small.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180709_TestingSplineModels/SonyColorMapData.mat';

% Override the wavelengths or colour channel indices at which to evaluate
% the model of dispersion, if desired.
bands = 430:10:650;
% Interpolation method used when resampling colour space conversion data
bands_interp_method = 'linear';

% Downsampling factor to apply to the estimated latent images relative to
% the input images. If empty (`[]`), downsampling will not occur.
downsampling_factor = 1;

% Output directory for all images and saved parameters
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180709_TestingSplineModels/admm';

% Whether or not to save the latent images to image files, beyond including
% them in the output '.mat' files
save_latent_image_files = false;

% Whether to expand the latent image relative to the input image to cover
% all undistorted coordinates from the input image. This is the
% `add_border` input argument.
add_border = true;

% ## Options for baek2017Algorithm2() (Alternating Direction Method of Multipliers)

% Whether to make the spectral gradient the same size as the image. This is
% the `full_GLambda` input argument.
baek2017Algorithm2Options.full_GLambda = false;

% Penalty parameters in ADMM, the `rho` input argument.
% Sample values seem to be in the range 1-10 (see pages 89, 93, and 95 of
% Boyd et al. 2011)
rho = [ 1, 1 ];

% Weights on the two prior terms, the `weights` input argument.
weights = [ 1, 1 ];

% Convergence tolerances in ADMM, the `tol` input argument.
%
% Reasonable values for the third element are 10^-4 to 10^-3 (page 21 of
% Boyd et al. 2011).
baek2017Algorithm2Options.tol = [ 1e-3, 1e-2, 1e-3 ];

% Maximum number of inner and outer iterations, the `maxit` input argument
baek2017Algorithm2Options.maxit = [ 500, 100 ];

% If the latent space consists of wavelength bands, use this type of
% numerical integration in 'channelConversionMatrix()'. (Otherwise, a value
% of 'none' will automatically be used instead.)
int_method = 'trap';

% ## Options for patch-wise image estimation

% Every combination of rows of `patch_sizes` and elements of `paddings`
% will be tested
patch_sizes = [ % Each row contains a (number of rows, number of columns) pair
    25 25;
]; 
paddings = 11;

% Only estimate a single patch, with its top-left corner at the given (row,
% column) location.
% If empty (`[]`), the entire image will be estimated.
target_patch = [];

% Also compare with whole image estimation
run_entire_image = false;

% ## Debugging Flags
baek2017Algorithm2Verbose = true;

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load calibration data

bands_script = bands;
bands = [];

optional_variable = 'bands';
model_variables_required = { 'dispersion_data', 'model_from_reference' };
model_variables_transform = { 'model_space', 'fill' };
load(...
    dispersion_model_filename,...
    model_variables_required{:}, model_variables_transform{:},...
    optional_variable...
    );
if ~all(ismember(model_variables_required, who))
    error('One or more of the required dispersion model variables is not loaded.')
end
if model_from_reference
    error('Dispersion model is in the wrong frame of reference.')
end

if exist('model_space', 'var')
    crop_image = true;
    if ~exist('fill', 'var')
        fill = false;
    end
else
    crop_image = false;
end

bands_dispersionfun = bands;
bands = [];

model_variables_required = { 'sensor_map', 'channel_mode' };
load(color_map_filename, model_variables_required{:}, optional_variable);
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end

bands_color = bands;

if channel_mode
    baek2017Algorithm2Options.int_method = 'none';
    solvePatchesOptions.int_method = 'none';
else
    baek2017Algorithm2Options.int_method = int_method;
    solvePatchesOptions.int_method = int_method;
end

%% Preprocessing input data

% Select the highest-priority value of `bands`
if ~isempty(bands_script)
    bands = bands_script;
elseif ~isempty(bands_color)
    bands = bands_color;
elseif ~isempty(bands_dispersionfun)
    bands = bands_dispersionfun;
else
    error('The variable `bands` is not defined, or is empty');
end

% Compare with colour space conversion data
n_bands = length(bands);
n_bands_sensor_map = size(sensor_map, 2);
resample_bands = false;
if ~isempty(bands_color)
    if n_bands ~= length(bands_color) || any(bands ~= bands_color)
        % Resampling is needed
        resample_bands = true;
        bands_for_interp = bands_color;
    end
elseif n_bands_sensor_map ~= n_bands
    % Resampling is needed, but will be "blind"
    resample_bands = true;
    bands_for_interp = linspace(bands(1), bands(end), n_bands_sensor_map);
end
n_latent_channels = size(sensor_map, 2);
% Resample colour space conversion data if necessary
if resample_bands
    [sensor_map_resampled, bands] = resampleArrays(...
        bands_for_interp, sensor_map.', bands,...
        bands_interp_method...
        );
    n_bands = length(bands);
    sensor_map_resampled = sensor_map_resampled.';
else
    sensor_map_resampled = sensor_map;
end

if save_latent_image_files && (n_latent_channels ~= 3 && n_latent_channels ~= 1)
    error('Cannot save latent images to image files, because they have %d channels.', n_latent_channels);
end

%% Process the images

img_ext = '.tif';
mat_ext = '.mat';
image_bounds = cell(n_images, 1);
solvePatchesOptions.add_border = add_border;

for i = 1:n_images
    [~, name, ext] = fileparts(image_filenames{i});
    if strcmp(ext, mat_ext)
        if isempty(input_images_variable_name)
            error('A variable name must be given to load input images from %s files.', mat_ext);
        end
        load(image_filenames{i}, input_images_variable_name);
        if exist(input_images_variable_name,'var')
            I_raw = eval(input_images_variable_name);
        else
            error(...
                'The input image variable %s was not loaded from %s.',...
                input_images_variable_name, image_filenames{i}...
                );
        end
    else
        I_raw = imread(image_filenames{i});
    end
    if ~ismatrix(I_raw)
        error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
    end
    
    if crop_image
        [roi, T_roi] = modelSpaceTransform(size(I_raw), model_space, fill);
        dispersionfun = makeDispersionfun(dispersion_data, T_roi);
        I_raw = I_raw(roi(1):roi(2), roi(3):roi(4));
    else
        dispersionfun = makeDispersionfun(dispersion_data);
    end
    
    I_raw = im2double(I_raw);
    image_sampling = size(I_raw);
    if ~isempty(downsampling_factor)
        image_sampling = ceil(image_sampling / downsampling_factor);
    end
    
    for ps = 1:size(patch_sizes, 1)
        for pad = 0:length(paddings)
            if run_entire_image && ps == 1 && pad == 0
                baek2017Algorithm2Options.add_border = add_border;
                [ I_latent, image_bounds{i}, I_rgb, J_full, J_est ] = baek2017Algorithm2(...
                    image_sampling, bayer_pattern, dispersionfun, sensor_map_resampled,...
                    bands, I_raw, rho, weights,...
                    baek2017Algorithm2Options, baek2017Algorithm2Verbose...
                );
                name_params = [name, '_whole'];
            elseif pad > 0
                baek2017Algorithm2Options.add_border = false;
                solvePatchesOptions.patch_size = patch_sizes(ps, :);
                solvePatchesOptions.padding = paddings(pad);
                name_params = [name, sprintf(...
                    '_patch%dx%d_pad%d',...
                    patch_sizes(ps, 1), patch_sizes(ps, 2), paddings(pad)...
                )];
                if isempty(target_patch)
                    [...
                        I_latent, image_bounds{i}, I_rgb, J_full, J_est...
                    ] = solvePatches(...
                        image_sampling, I_raw, bayer_pattern, dispersionfun,...
                        sensor_map_resampled,...
                        bands, solvePatchesOptions, @baek2017Algorithm2,...
                        {...
                            rho, weights,...
                            baek2017Algorithm2Options, baek2017Algorithm2Verbose...
                        }...
                    );
                    
                else
                    [...
                        I_latent, image_bounds{i}, I_rgb, J_full, J_est...
                    ] = solvePatches(...
                        image_sampling, I_raw, bayer_pattern, dispersionfun,...
                        sensor_map_resampled,...
                        bands, solvePatchesOptions, @baek2017Algorithm2,...
                        {...
                            rho, weights,...
                            baek2017Algorithm2Options, baek2017Algorithm2Verbose...
                        }, target_patch...
                    );
                    name_params = [name_params, sprintf(...
                        '_target%dAnd%d',...
                        target_patch(1), target_patch(2)...
                    )];
                end
            else
                continue;
            end
            
            % Save the results
            I_filename = fullfile(output_directory, [name_params '_roi']);
            imwrite(I_raw, [I_filename img_ext]);
            save([I_filename mat_ext], 'I_raw');
            I_filename = fullfile(output_directory, [name_params '_latent']);
            if save_latent_image_files
                imwrite(I_latent, [I_filename img_ext]);
            end
            save([I_filename mat_ext], 'I_latent');
            I_filename = fullfile(output_directory, [name_params '_latent_rgb']);
            imwrite(I_rgb, [I_filename img_ext]);
            save([I_filename mat_ext], 'I_rgb');
            I_filename = fullfile(output_directory, [name_params '_latent_warped']);
            imwrite(J_full, [I_filename img_ext]);
            save([I_filename mat_ext], 'J_full');
            I_filename = fullfile(output_directory, [name_params '_reestimated']);
            imwrite(J_est, [I_filename img_ext]);
            save([I_filename mat_ext], 'J_est');
        end
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'image_filenames',...
        'bands_dispersionfun',...
        'bands_color',...
        'bands',...
        'sensor_map_resampled',...
        'image_bounds'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByHyperspectralADMM.mat');
save(save_data_filename, save_variables_list{:});