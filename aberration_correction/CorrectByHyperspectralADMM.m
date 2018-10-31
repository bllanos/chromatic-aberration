%% Demosaicing and hyperspectral ADMM-based correction of chromatic aberration
% Convert RAW images to colour images, and simultaneously correct chromatic
% aberration, by estimating a latent hyperspectral or RGB image.
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
% 'RAWDiskDispersion.m', 'DoubleConvexThickLensDispersion.m' or
% 'BimaterialImages.m', for example. The following variables are required:
% - 'dispersion_data': A model of chromatic aberration, modeling the warping
%   from the reference colour channel or wavelength band to the other
%   colour channels or wavelength bands. `dispersion_data` can be converted to
%   a function form using `dispersionfun = makeDispersionfun(dispersion_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `false`.
% - 'bands': A vector containing the wavelengths or colour channel indices
%   to use as the `lambda` input argument of 'dispersionfunToMatrix()'.
%   `bands` is the wavelength or colour channel information needed to
%   evaluate the dispersion model.
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
%   corresponding to the second dimension of 'sensor_map'.
%
% ## Output
%
% ### Estimated images
%
% One of each of the following types of images is created for each input
% image, except where specified. The filename of the input image, concatenated
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
%   output if the latent images are greyscale or 3-channel images.
% - '*_warped.tif' and '*_warped.mat': A version of the latent image
%   (stored in the variable 'I_warped') created by warping the latent image
%   according to the dispersion model. The '.tif' image is only output if
%   the latent images are greyscale or 3-channel images.
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
% ### Regularization weights images
%
% If `use_fixed_weights` in 'SetFixedParameters.m' is `false`, then the
% image estimation algorithm will automatically choose weights on the
% regularization terms in the ADMM optimization problem. For the i-th
% enabled regularization term in the ADMM optimization problem, an image
% will be output, as the variable 'I_weights', in the file
% '*_weight${i}Image.mat', where '*' represents the filename of the input
% image concatenated with a string of parameter information. A pixel in the
% image will contain the weight on the i-th regularization term used when
% estimating the pixel.
%
% Furthermore, base-10 logarithmic compressions of the images will be shown
% in figures, and the figures will be saved as '*_weight${i}Image.fig'.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'bands': A vector containing the wavelengths or colour channel
%   indices at which the latent images are sampled.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file.
% - 'color_weights': A matrix for converting pixels in the latent image to
%   colour, as determined by the 'sensor_map' variable loaded from the
%   colour space conversion data file, and by the type of numerical
%   intergration to perform.
% - 'spectral_weights': A matrix for converting pixels in the latent image
%   to the spectral, or colour channel space, described by 'bands_color'.
% - 'image_filenames': A cell vector containing the input image filenames
%   retrieved based on the wildcard provided in the parameters section of
%   the script.
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
% - One can use this script to estimate either a latent hyperspectral
%   image, or a latent aberration-free RGB image (free from lateral
%   chromatic aberration). A latent hyperspectral image can be sharper, in
%   theory, whereas a latent RGB image will retain the within-channel
%   chromatic aberration of the input image. The reason for this difference
%   is the summation of multiple spectral bands into each channel of an RGB
%   image, in contrast to the identity mapping of the colours of a latent
%   RGB image into the colours of the aberrated RGB image. Summation allows
%   multiple sharp bands to form a blurred colour channel.
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
        'reverse_dispersion_model_filename',...
        'color_map_filename',...
        'output_directory',...
        'target_patch',...
        'run_entire_image'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
% '.mat' or image files can be loaded
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/lacelike*raw.tif';
input_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Model of dispersion
% Can be empty
reverse_dispersion_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/BimaterialImagesData.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/Downloads/SonyColorMapData.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% ## Options for patch-wise image estimation

% Only estimate a single patch, with its top-left corner at the given (row,
% column) location.
% If empty (`[]`), the entire image will be estimated.
target_patch = [];

% Also compare with (or only run) whole image estimation. Only enable this
% for small images.
run_entire_image = false;

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load calibration data

has_dispersion = ~isempty(reverse_dispersion_model_filename);
if has_dispersion
    [...
        dispersion_data, bands_dispersionfun, transform_data...
    ] = loadDispersionModel(reverse_dispersion_model_filename, false, false);
end

model_variables_required = { 'sensor_map', 'channel_mode', 'bands' };
load(color_map_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end

bands_color = bands;

if channel_mode
    if has_dispersion && ...
       ((length(bands_color) ~= length(bands_dispersionfun)) ||...
       any(bands_color ~= bands_dispersionfun))
        error('When estimating a colour image, the same colour channels must be used by the model of dispersion.');
    end
    color_weights = sensor_map;
    spectral_weights = eye(length(bands_color));
else
    [color_weights, spectral_weights, bands] = samplingWeights(...
      sensor_map, bands_color, bands_color, samplingWeightsOptions, samplingWeightsVerbose...
    );
end

%% Process the images

if isempty(patch_sizes) && ~run_entire_image
    error('Neither patch-based image estimation, nor whole image estimation, were requested.');
end

if ~isempty(target_patch)
    solvePatchesADMMOptions.patch_options.target_patch = target_patch;
    if run_entire_image
        error('The entire image will not be estimated if a target patch is set.');
    end
end

for i = 1:n_images
    [I_raw, name] = loadImage(image_filenames{i}, input_images_variable_name);

    if ~ismatrix(I_raw)
        error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
    end
    
    if has_dispersion
        [dispersionfun, I_raw] = makeDispersionForImage(...
            dispersion_data, I_raw, transform_data...
        );
    else
        dispersionfun = [];
    end
    
    image_sampling = size(I_raw);
    
    for ps = 0:size(patch_sizes, 1)
        for pad = 0:length(paddings)
            if run_entire_image && ps == 0 && pad == 0
                solvePatchesADMMOptions.patch_options.patch_size = image_sampling;
                solvePatchesADMMOptions.patch_options.padding = 0;
                name_params = [name, '_whole'];
            elseif ps > 0 && pad > 0
                solvePatchesADMMOptions.patch_options.patch_size = patch_sizes(ps, :);
                solvePatchesADMMOptions.patch_options.padding = paddings(pad);
                name_params = [name, sprintf(...
                    '_patch%dx%d_pad%d',...
                    patch_sizes(ps, 1), patch_sizes(ps, 2), paddings(pad)...
                )];
            else
                continue;
            end
            
            if ~isempty(target_patch)
                name_params = [name_params, sprintf(...
                    '_target%dAnd%d',...
                    target_patch(1), target_patch(2)...
                )];
            end
            [...
                I_latent,...
                I_rgb,...
                weights_images,...
                J_full,...
                J_est,...
                I_warped...
            ] = solvePatchesADMM(...
              [], I_raw, bayer_pattern, dispersionfun,...
              color_weights, bands,...
              solvePatchesADMMOptions.admm_options,...
              solvePatchesADMMOptions.reg_options,...
              solvePatchesADMMOptions.patch_options,...
              solvePatchesADMMVerbose...
            );

            if ~use_fixed_weights
                enabled_weights = solvePatchesADMMOptions.reg_options.enabled;
                n_active_weights = sum(enabled_weights);
                to_all_weights = find(enabled_weights);
                for w = 1:n_active_weights
                    aw = to_all_weights(w);
                    saveImages(...
                        'data', output_directory, name_params,...
                        weights_images(:, :, w), sprintf('weight%dImage', aw), 'I_weights'...
                    );

                    fg = figure;
                    imagesc(log10(weights_images(:, :, w)));
                    c = colorbar;
                    c.Label.String = sprintf('log_{10}(weight %d)', aw);
                    xlabel('Image x-coordinate')
                    ylabel('Image y-coordinate')
                    title(sprintf('Per-patch regularization weight %d', aw));
                    savefig(...
                        fg,...
                        fullfile(output_directory, [name_params  sprintf('weight%dImage.fig', aw)]),...
                        'compact'...
                        );
                    close(fg);
                end
            end
            
            % Save the results
            saveImages(...
                output_directory, name_params,...
                I_raw, '_roi', 'I_raw',...
                I_latent, '_latent', 'I_latent',...
                I_rgb, '_rgb', 'I_rgb',...
                J_full, '_rgb_warped', 'J_full',...
                J_est, '_reestimated', 'J_est',...
                I_warped, '_warped', 'I_warped'...
            );
        end
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'bands',...
        'bands_color',...
        'color_weights',...
        'spectral_weights',...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByHyperspectralADMM.mat');
save(save_data_filename, save_variables_list{:});