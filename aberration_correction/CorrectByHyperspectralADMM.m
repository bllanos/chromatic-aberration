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
% 'PreprocessRAWImages.m', so that they do not need to be linearized after
% being loaded. For image format files, images will simply be loaded with
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
% 'SonyColorMap.m', for example, or other scripts in 'sensor/'. The following
% variables are required:
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
% Refer to the section 'Output images' in SetFixedParameters.m'.
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
% estimating the pixel. If `target_patch_weights` in the script parameters is
% not empty, then weights are chosen based only on the given target patch, and
% 'I_weights' will have the dimensions defined by the patch size.
%
% Also, base-10 logarithmic compressions of the regularization weights images
% will be shown in figures, and the figures will be saved as
% '*_weight${i}Image.fig'.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'bands': A vector containing the wavelengths or colour channel
%   indices at which the latent images are sampled.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file.
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
%   produced in the camera's colour space.
% - One can use this script to estimate either a latent hyperspectral
%   image, or a latent aberration-free RGB image (free from lateral
%   chromatic aberration). A latent hyperspectral image can be sharper, in
%   theory, whereas a latent RGB image will retain the within-channel
%   chromatic aberration of the input image. The reason for this difference
%   is the summation of multiple spectral bands into each channel of an RGB
%   image, in contrast to the identity mapping of the colours of a latent
%   RGB image into the colours of the aberrated RGB image. Summation allows
%   multiple sharp bands to form a blurred colour channel.
% - If `solvePatchesMultiADMMOptions.sampling_options.show_steps` is
%   `true`, and the images are being estimated in a spectral space, as
%   opposed to a colour space, then saved value of 'bands' will be a cell
%   vector describing the multiple spectral resolutions used for image
%   estimation. One set of output images, including regularization weights
%   images, will be saved for each spectral resolution. Refer to the
%   documentation of 'solvePatchesMultiADMM.m' for more information.
% - Presently, this script ignores `criteria` in 'SetFixedParameters.m', and
%   only selects regularization weights using either the minimum distance
%   criteria, or based on similarity with a demosaicing result, depending on the
%   value of `solvePatchesADMMOptions.reg_options.demosaic` in
%   'SetFixedParameters.m'.
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
    'target_patch_weights',...
    'run_entire_image'...
};

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
% '.mat' or image files can be loaded
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/exposure_blending/d2_book*_unfiltered.mat';
input_images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Model of dispersion
% Can be empty
reverse_dispersion_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/spectral/full_image/RAWDiskDispersionResults_spectral_polynomial_fromNonReference.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/SonyColorMapData.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads/new_patch';

% ## Options for patch-wise image estimation

% Only estimate a single patch, with its top-left corner at the given (row,
% column) location. If empty (`[]`), the entire image will be estimated.
target_patch = [];

% Only select regularization weights for a single patch, with its top-left
% corner at the given (row, column) location. If empty (`[]`), regularization
% weights will be selected for each patch separately. (THIS IS SLOW)
target_patch_weights = [141, 1115];

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
end

if channel_mode
    options = solvePatchesADMMOptions;
    multi_step = false;
else
    options = solvePatchesMultiADMMOptions;
    multi_step = options.sampling_options.show_steps;
end

%% Process the images

if isempty(patch_sizes) && ~run_entire_image
    error('Neither patch-based image estimation, nor whole image estimation, were requested.');
end
has_target_patch = ~isempty(target_patch);
use_target_patch_weights = ~isempty(target_patch_weights) && ~use_fixed_weights;
if has_target_patch && run_entire_image
    error('The entire image will not be estimated if a target patch is set.');
end

if ~use_fixed_weights
    n_weights = length(options.reg_options.enabled);
    enabled_weights = options.reg_options.enabled;
    n_active_weights = sum(enabled_weights);
end

if save_all_images
    extra_images = cell(3, 1);
else
    extra_images = cell(0, 1);
end

for i = 1:n_images
    [I_raw, name] = loadImage(image_filenames{i}, input_images_variable_name);

    if ~ismatrix(I_raw)
        error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
    end
    
    if has_dispersion
        [dispersionfun, I_raw] = makeDispersionForImage(...
            dispersion_data, I_raw, transform_data, true...
        );
    else
        dispersionfun = [];
    end
    
    image_sampling = size(I_raw);
    
    for ps = 0:size(patch_sizes, 1)
        for pad = 0:length(paddings)
            if run_entire_image && ps == 0 && pad == 0
                options.patch_options.patch_size = image_sampling;
                options.patch_options.padding = 0;
                name_params = [name, '_whole'];
            elseif ps > 0 && pad > 0
                options.patch_options.patch_size = patch_sizes(ps, :);
                options.patch_options.padding = paddings(pad);
                name_params = [name, sprintf(...
                    '_patch%dx%d_pad%d',...
                    patch_sizes(ps, 1), patch_sizes(ps, 2), paddings(pad)...
                )];
            else
                continue;
            end
            
            if has_target_patch
                name_params = [name_params, sprintf(...
                    '_target%dAnd%d',...
                    target_patch(1), target_patch(2)...
                )];
            end
            if use_target_patch_weights
                name_params = [name_params, sprintf(...
                    '_weightsTarget%dAnd%d',...
                    target_patch_weights(1), target_patch_weights(2)...
                )];
                options.patch_options.target_patch = target_patch_weights;
            elseif has_target_patch
                options.patch_options.target_patch = target_patch;
            end
            if use_fixed_weights
                name_params = [name_params, '_fixedWeights'];
            end
            if ~use_target_patch_weights ||...
                    (has_target_patch && use_target_patch_weights && all(target_patch == target_patch_weights))
                if channel_mode
                    [...
                        I_latent,...
                        I_rgb,...
                        weights_images,...
                        extra_images{:}...
                    ] = solvePatchesADMM(...
                      [], I_raw, bayer_pattern, dispersionfun,...
                      sensor_map, bands_color,...
                      options.admm_options,...
                      options.reg_options,...
                      options.patch_options,...
                      solvePatchesADMMVerbose...
                    );
                else
                    [...
                        bands,...
                        I_latent,...
                        I_rgb,...
                        weights_images,...
                        extra_images{:}...
                    ] = solvePatchesMultiADMM(...
                        [], I_raw, bayer_pattern, dispersionfun, sensor_map, bands_color,...
                        options.sampling_options,...
                        options.admm_options,...
                        options.reg_options,...
                        options.patch_options,...
                        solvePatchesMultiADMMVerbose...
                    );
                end
            elseif use_target_patch_weights
                if channel_mode
                    [ ~, ~, weights_images ] = solvePatchesADMM(...
                      [], I_raw, bayer_pattern, dispersionfun,...
                      sensor_map, bands_color,...
                      options.admm_options,...
                      options.reg_options,...
                      options.patch_options,...
                      solvePatchesADMMVerbose...
                    );
                else
                    show_steps_original = options.sampling_options.show_steps;
                    options.sampling_options.show_steps = true;
                    [~ , ~, ~, weights_images ] = solvePatchesMultiADMM(...
                        [], I_raw, bayer_pattern, dispersionfun, sensor_map, bands_color,...
                        options.sampling_options,...
                        options.admm_options,...
                        options.reg_options,...
                        options.patch_options,...
                        solvePatchesMultiADMMVerbose...
                    );
                    options.sampling_options.show_steps = show_steps_original;
                end
                
                reg_options_original = options.reg_options;
                if channel_mode
                    weights = reshape(weights_images(1, 1, :), 1, []);
                    options.reg_options.minimum_weights = zeros(1, n_weights);
                    options.reg_options.maximum_weights = zeros(1, n_weights);
                    options.reg_options.minimum_weights(enabled_weights) = weights;
                    options.reg_options.maximum_weights(enabled_weights) = weights;
                else
                    options.reg_options.multi_weights = reshape(weights_images(1, 1, :), n_active_weights, []).';
                end
                if has_target_patch
                    options.patch_options.target_patch = target_patch;
                else
                    options.patch_options = rmfield(options.patch_options, 'target_patch');
                end
                
                if channel_mode
                    [...
                        I_latent,...
                        I_rgb,...
                        ~,...
                        extra_images{:}...
                    ] = solvePatchesADMM(...
                      [], I_raw, bayer_pattern, dispersionfun,...
                      sensor_map, bands_color,...
                      options.admm_options,...
                      options.reg_options,...
                      options.patch_options,...
                      solvePatchesADMMVerbose...
                    );
                else
                    [...
                        bands,...
                        I_latent,...
                        I_rgb,...
                        ~,...
                        extra_images{:}...
                    ] = solvePatchesMultiADMM(...
                        [], I_raw, bayer_pattern, dispersionfun, sensor_map, bands_color,...
                        options.sampling_options,...
                        options.admm_options,...
                        options.reg_options,...
                        options.patch_options,...
                        solvePatchesMultiADMMVerbose...
                    );
                end
                options.reg_options = reg_options_original;
            end

            if multi_step
                bands_all = bands;
            else
                bands_all = {bands};
            end
            
            spectral_inc = 0;
            color_inc = 0;
            raw_inc = 0;
            weights_inc = 0;
            for t = 1:length(bands_all)
                if multi_step
                    name_params_t = [name_params, sprintf('_step%d', t)];
                else
                    name_params_t = name_params;
                end
                n_bands = length(bands_all{t});
                
                if ~use_fixed_weights
                    to_all_weights = find(enabled_weights);
                    for w = 1:n_active_weights
                        aw = to_all_weights(w);
                        saveImages(...
                            'data', output_directory, name_params_t,...
                            weights_images(:, :, weights_inc + w), sprintf('weight%dImage', aw), 'I_weights'...
                        );

                        fg = figure;
                        imagesc(log10(weights_images(:, :, weights_inc + w)));
                        c = colorbar;
                        c.Label.String = sprintf('log_{10}(weight %d)', aw);
                        xlabel('Image x-coordinate')
                        ylabel('Image y-coordinate')
                        title(sprintf('Per-patch regularization weight %d', aw));
                        savefig(...
                            fg,...
                            fullfile(output_directory, [name_params_t  sprintf('weight%dImage.fig', aw)]),...
                            'compact'...
                            );
                        close(fg);
                    end
                    weights_inc = weights_inc + n_active_weights;
                end

                % Save the results
                saveImages(...
                    output_directory, name_params_t,...
                    I_latent(:, :, (spectral_inc + 1):(spectral_inc + n_bands)), '_latent', 'I_latent'...
                );
                saveImages(...
                    'image', output_directory, name_params_t,...
                    I_rgb(:, :, (color_inc + 1):(color_inc + size(sensor_map, 1))), '_rgb', 'I_rgb'...
                );
                if save_all_images
                    saveImages(...
                        output_directory, name_params_t,...
                        I_raw, '_roi', 'I_raw',...
                        extra_images{1}(:, :, (color_inc + 1):(color_inc + size(sensor_map, 1))), '_rgb_warped', 'J_full',...
                        extra_images{2}(:, :, (raw_inc + 1):(raw_inc + size(I_raw, 3))), '_reestimated', 'J_est',...
                        extra_images{3}(:, :, (spectral_inc + 1):(spectral_inc + n_bands)), '_warped', 'I_warped'...
                    );
                end
                
                spectral_inc = spectral_inc + n_bands;
                color_inc = color_inc + size(sensor_map, 1);
                raw_inc = raw_inc + size(I_raw, 3);
            end
        end
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'bands',...
        'bands_color',...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByHyperspectralADMM.mat');
save(save_data_filename, save_variables_list{:});