%% Demosaicing and hyperspectral ADMM-based correction of chromatic aberration
% Test the L-hypersurface method of Belge et al. 2002 for selecting
% regularization weights
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
%
% #### RAW image
% A RAW image to be demosaiced and corrected for chromatic aberration.
%
% The image is expected to have been preprocessed, such as using
% 'AverageRAWImages.m', so that it does not need to be linearized after
% being loaded.  For image format files, the image will simply be loaded
% with the Image Processing Toolbox 'imread()' function. For '.mat' files,
% the variable to be loaded must be provided in the script parameters.
%
% The image is expected to have 3 colour channels (Red, Green, Blue)
% (represented in a Bayer pattern as a 2D array). However, the colour
% channels can correspond to narrowband wavelength ranges - This script
% will input a mapping from the colour space of the latent image to the
% colour space of the RAW image.
%
% #### True image
% A spectral or colour image serving as the ground truth for image
% estimation. The true image is needed to compare the weights selected
% using the L-hypersurface with the weights giving the lowest error
% relative to the true image.
%
% The true image must be associated with a '.mat' file containing a
% vector with the variable 'bands'. 'bands' must have the same length as
% the third dimension of the true image, and must contain the colour
% channel indices or wavelengths corresponding to the true image. 'bands'
% is used to evaluate the dispersion model. Note that 'bands' takes
% precedence over the same variable defined in 'SetFixedParameters.m'.
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
% - 'bands': A vector containing the wavelengths or colour channel indices
%   corresponding to the second dimension of 'sensor_map'. 'bands' is
%   required to resample 'sensor_map' so that it maps from the colour space
%   of the latent image to the colour space of the input RAW image.
%
% ## Output
%
% ### Graphical output
%
% Figures are opened showing the search path taken by the fixed-point
% iterative method of Belge et al. 2002 for selecting regularization
% weights. The search path can be shown on plots of the L-hypersurface, and
% of the data fitting error hypersurface, depending on the amount of
% graphical output requested. (Sampling these surfaces is
% computationally-expensive.) Lastly, an additional figure will show the
% location of the image patch used for selecting the regularization
% weights.
%
% Graphical output will not be produced if there are more than three
% regularization weights to be chosen.
%
% ### Estimated images
%
% One of each of the following types of images is created, depending on the
% type of latent image (spectral or colour). The images are produced under
% the regularization weights chosen by the method of Belge et al. 2002. The
% filename of the input image, concatenated with a string of parameter
% information, is represented by '*' below.
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
% - '*_warped.tif' and '*_warped.mat': A version of the latent image
%   (stored in the variable 'I_warped') created by warping the latent image
%   according to the dispersion model. The '.tif' image is only output if
%   `save_latent_image_files` is `true`, and an error will be thrown if the
%   latent images are not greyscale or 3-channel images.
% - '*_rgb.tif' and '*_rgb.mat': A colour image (stored in the variable
%   'I_rgb') created by converting the latent image to the RGB colour space
%   of the camera.
% - '*_rgb_warped.tif' and '*_rgb_warped.mat': A colour image (stored in
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
% - 'bands': The value of the 'bands' variable loaded with the true latent
%   image.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference. 'bands_color' is empty if no
%   such variable was found in the data file, or if the value of the loaded
%   variable was empty.
% - 'image_bounds': The coordinates of the latent image in the space of the
%   cropped input image. This is the 'image_bounds' output argument of
%   'dispersionfunToMatrix()'.
% - 'input_image_filename': The input image filename found using the
%   wildcard provided in the parameters section of the script.
% - 'true_image_filename': The latent image filename found using the
%   wildcard provided in the parameters section of the script.
% - 'sensor_map_resampled': The resampled version of the 'sensor_map'
%   variable, generated for compatibility with the true latent image.
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
% - This script uses the first row of `weights` defined in
%   'SetFixedParameters.m' to initialize the fixed-point algorithm for
%   selecting regularization weights based on the image content
%   (implemented in 'selectWeights()'). Elements of `weights(1, :)` can be
%   set to zero to disable the corresponding regularization terms. Note
%   that the number of nonzero elements of `weights(1, :)` determines the
%   dimensionality of the visualizations output by this script.
% - This script could estimate downsampled images (configured by adjusting
%   `downsampling_factor` in 'SetFixedParameters.m'), if it were to use
%   'solvePatches()' instead of 'solvePatchesAligned()' for patch-based
%   image estimation. In that case, however, this script could not process
%   very large images, except at higher downsampling factors.
%   'solvePatchesAligned()' can process large images, but cannot downsample
%   images.
% - This script only uses the first row of `patch_sizes`, and the first
%   element of `paddings`, defined in 'SetFixedParameters.m'.
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
% - Belge, M, Kilmer, M. E., & Miller, E. L.. "Efficient determination of
%   multiple regularization parameters in a generalized L-curve
%   framework." Inverse Problems, vol. 18, pp. 1161-1183, 2002.
%   doi:10.1088/0266-5611/18/4/314

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 24, 2018

% List of parameters to save with results
parameters_list = {
        'reverse_dispersion_model_filename',...
        'color_map_filename',...
        'output_directory',...
        'target_patch'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the image to process.
% '.mat' or image files can be loaded
input_image_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/lacelike*raw.tif';
input_image_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the true image.
% '.mat' or image files can be loaded
true_image_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/lacelike*hyper.mat';
true_image_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% Data file containing the colour channels or wavelengths associated with
% the true image
true_image_bands_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/BimaterialImagesData.mat';

% Model of dispersion
reverse_dispersion_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/BimaterialImagesData.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/NikonD5100ColorMapData.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% ## Options for the L-hypersurface method of Belge et al. 2002

% The top-left corner (row, column) of the image patch to use for
% regularization weights selection. If empty (`[]`), the patch will be
% selected by the user.
target_patch = [];

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Load the images

input_image_filename = listFiles(input_image_wildcard);
[I_raw, name] = loadImage(input_image_filename{1}, input_image_variable_name);

if ~ismatrix(I_raw)
    error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
end

true_image_filename = listFiles(true_image_wildcard);
I_gt = loadImage(true_image_filename{1}, true_image_variable_name);

bands = [];
bands_variable_name = 'bands';
load(true_image_bands_filename, bands_variable_name);
if isempty(bands)
    error('No wavelength band or colour channel information is associated with the true image.')
end
bands_gt = bands;

%% Load calibration data

[...
    dispersion_data, ~, transform_data...
] = loadDispersionModel(reverse_dispersion_model_filename, false, false);

bands = [];
model_variables_required = { 'sensor_map', 'channel_mode' };
load(color_map_filename, model_variables_required{:}, bands_variable_name);
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end
if isempty(bands)
    error('No (non-empty) variable `bands` loaded from colour space conversion data.');
end

bands_color = bands;
bands = bands_gt;

if channel_mode
    baek2017Algorithm2Options.int_method = 'none';
    solvePatchesOptions.int_method = 'none';
    selectWeightsOptions.int_method = 'none';
else
    baek2017Algorithm2Options.int_method = int_method;
    solvePatchesOptions.int_method = int_method;
    selectWeightsOptions.int_method = int_method;
end

%% Preprocessing input data

n_bands = length(bands);
% Resample colour space conversion data if necessary
if n_bands ~= length(bands_color) || any(bands ~= bands_color)
    [sensor_map_resampled, bands] = resampleArrays(...
        bands_color, sensor_map.', bands,...
        bands_interp_method...
        );
    if n_bands ~= length(bands)
        error('The colour space conversion data does not cover a sufficiently large range of wavelengths.');
    end
    sensor_map_resampled = sensor_map_resampled.';
else
    sensor_map_resampled = sensor_map;
end

% Crop images to the region of valid dispersion
[dispersionfun, I_raw] = makeDispersionForImage(...
    dispersion_data, I_raw, transform_data...
);
image_sampling = size(I_raw);

roi = modelSpaceTransform(...
    [size(I_gt, 1), size(I_gt, 2)],...
    transform_data.model_space, transform_data.fill...
);
if ~isempty(roi)
    I_gt = I_gt(roi(1):roi(2), roi(3):roi(4), :);
end
if any([size(I_gt, 1), size(I_gt, 2)] ~= image_sampling)
    error([
        'The RAW version of the image has different spatial dimensions fro',...
        'm the true latent image.'...
    ]);
end

%% Process the image

solvePatchesOptions.add_border = add_border; % Not used by solvePatchesAligned()
baek2017Algorithm2Options.add_border = false;
patch_size = patch_sizes(1, :);
padding = paddings(1);
solvePatchesOptions.patch_size = patch_size;
solvePatchesOptions.padding = padding;
weights = weights(1, :);

if ~isempty(downsampling_factor)
    if downsampling_factor ~= 1
        warning([...
            '`downsampling_factor` is ignored, because solvePatchesAligned(',...
            ') will be used instead of solvePatches().'...
        ]);
    end
    % image_sampling = ceil(image_sampling / downsampling_factor);
end

name_params = [name, sprintf(...
    '_patch%dx%d_pad%d_weights%ew%ew%e_',...
    patch_size(1), patch_size(2), padding,...
    weights(1), weights(2), weights(3)...
)];
[...
    I_latent, image_bounds, I_rgb, J_full, J_est, I_warped...
] = solvePatchesAligned(...
    I_raw, bayer_pattern, dispersionfun,...
    sensor_map_resampled,...
    bands, solvePatchesOptions, @baek2017Algorithm2,...
    {...
        weights, rho,...
        baek2017Algorithm2Options, baek2017Algorithm2Verbose...
    }...
);

% Save the results
saveImages(...
    output_directory, name_params,...
    I_raw, 'roi', 'I_raw',...
    I_latent, 'latent', 'I_latent',...
    I_rgb, 'rgb', 'I_rgb',...
    J_full, 'rgb_warped', 'J_full',...
    J_est, 'reestimated', 'J_est',...
    I_warped, 'warped', 'I_warped'...
);

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'input_image_filename',...
        'true_image_filename',...
        'bands_color',...
        'bands',...
        'sensor_map_resampled',...
        'image_bounds'...
    } ];
save_data_filename = fullfile(output_directory, 'ValidateLHypersurface.mat');
save(save_data_filename, save_variables_list{:});