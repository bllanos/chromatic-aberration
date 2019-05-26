%% Correction of chromatic aberration by image warping
% Correct chromatic aberration by warping colour channels or spectral bands
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
% A set of full-colour or spectral images to be corrected for chromatic
% aberration.
%
% Spectral images must be accompanied by a '.mat' file containing the
% wavelengths at which the spectral images are defined, or describing the
% spectral basis representation more generally, as discussed in the
% parameters section below.
%
% Images may have different pixel dimensions, provided that they are
% compatible with the input model of dispersion described below.
%
% ### Model of dispersion
%
% A '.mat' file containing several variables, which is the output of
% 'RAWDiskDispersion.m', for example. The following variables are required:
% - 'dispersion_data': A model of chromatic aberration, modeling the warping
%   from the reference colour channel or wavelength band to the other
%   colour channels or wavelength bands. `dispersion_data` can be converted to
%   a function form using `dispersionfun = makeDispersionfun(dispersion_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `true`.
% - 'bands': A vector containing the wavelengths or colour channel indices
%   at which the dispersion model was originally fit. In the case of colour
%   channels, 'bands' is needed to check if the dispersion model is
%   compatible with the colour space conversion data (see below).
%   Otherwise, it is not used.
%
% The following two additional variables are optional. If they are present,
% they will be used for the following purposes:
% - Conversion between the coordinate system in which the model of
%   chromatic aberration was constructed and the image coordinate system.
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
%   j-th colour channel or spectral band of the input images. For example,
%   `sensor_map` is a matrix mapping discretized spectral power
%   distributions to RGB colours.
% - 'channel_mode': A Boolean value indicating whether the latent colour
%   space is a set of colour channels (true) or a set of spectral bands
%   (false).
% - 'bands': A vector containing the wavelengths or colour channel indices
%   corresponding to the second dimension of 'sensor_map'.
%
% Colour space conversion data is only required if images are to be converted
% from spectra to colour, or between colour spaces.
%
% ## Output
%
% ### Warped images
%
% For each input image, a '.mat' file with a name ending in '_unwarped' is
% created. If colour conversion is requested, an 'I_unwarped' variable in the
% file contains the corrected image, converted to the output colour space
% according to the colour space conversion data described above. Otherwise, if
% spectral images are being processed, but are not being converted to colour,
% the 'I_unwarped' variable in the file contains the corrected spectral image.
% Images are cropped to the region in which the dispersion model is valid.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'image_filenames': A cell vector containing the input image filenames
%   retrieved based on the wildcard provided in the parameters section of
%   the script.
% - 'time': The time, in seconds, taken to produce the each output image.
%   'time' is a vector with a length equal to the number of images.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## Notes
% - Saving colour images is recommended, even if the input images are spectral
%   images (as also mentioned in the documentation of 'warpImageSpectral.m'):
%   - Creating warped spectral images requires a large amount of memory,
%     especially if `disable_interpolation` below is `false`.
%   - Warping changes the spectral sampling of the spectral images, because
%     they need to be converted from a basis representation to discretized
%     spectral signatures. Unless the conversion is an identity mapping, it
%     is not possible to compare them with the original spectral images.
%
% ## References
% - Rudakova, V. & Monasse, P. (2014). "Precise correction of lateral
%   chromatic aberration in images" (Guanajuato). 6th Pacific-Rim Symposium
%   on Image and Video Technology, PSIVT 2013. Springer Verlag.
%   doi:10.1007/978-3-642-53842-1_2
% - J. Brauers and T. Aach. "Geometric Calibration of Lens and Filter
%   Distortions for Multispectral Filter-Wheel Cameras," IEEE Transactions
%   on Image Processing, vol. 20, no. 2, pp. 496-505, 2011.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 23, 2019

% List of parameters to save with results
parameters_list = {
    'input_images_variable_name',...
    'is_spectral',...
    'convert_to_color',...
    'color_map_filename',...
    'sampling_filename',...
    'forward_dispersion_model_filename',...
    'disable_interpolation'...
};

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
% '.mat' or image files can be loaded
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/*_dHyper.mat';
input_images_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% Whether the input images are spectral images or colour images
is_spectral = true;

% Colour space conversion data
convert_to_color = false; % Whether or not to convert images to a different colour space
if convert_to_color
    color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/sensor.mat';
else
    color_map_filename = []; % Unused
end

% Spectral sampling
if is_spectral
    % Path and filename of a '.mat' file containing the information needed to
    % convert control bands to raw colour channels, or to physical spectral bands.
    % The file must contain a vector, `bands`, of wavelengths corresponding to the
    % spectral bands of the spectral images. It must also contain a variable,
    % `findSamplingOptions`, such as defined in 'SetFixedParameters.m'.
    sampling_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dHyper_dispersion/RunOnDataset_20190421_ComputarLens_dHyper_dispersion.mat';
    
    % Whether or not to prevent spectral resampling and interpolation during
    % warping. When `disable_interpolation` is `true`, the `options` input
    % argument of `dispersionfunToMatrix` only has a 'bands_in' field. It only
    % makes sense for this parameter to be `true` when the spectral images were
    % obtained using narrowband optical filters.
    disable_interpolation = true;
else
    sampling_filename = []; % Unused
    disable_interpolation = []; % Unused
end

% Model of dispersion
forward_dispersion_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/spectral/full_image/RAWDiskDispersionResults_spectral_polynomial_fromReference.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/warping_filtered_images/bandpassFiltered_warpCorrected_asSpectralImages';

% Parameters which do not usually need to be changed. Some of these
% parameters will be overridden by the input data.
run('SetFixedParameters.m')

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load calibration data

[...
    dispersion_data, bands_dispersionfun, transform_data...
] = loadDispersionModel(forward_dispersion_model_filename, true, true);

if convert_to_color
    [sensor_map, channel_mode, bands_color] = loadColorMap(color_map_filename);
    if is_spectral ~= channel_mode
        error('`is_spectral` must match the `channel_mode` flag in the colour space conversion data.');
    end
    if channel_mode
        if ((length(bands) ~= length(bands_dispersionfun)) ||...
            any(bands(:) ~= bands_dispersionfun(:)))
            error('When warping colour images, the same colour channels must be used by the model of dispersion.');
        end
        color_weights = sensor_map;
    end
end

if is_spectral
    [findSamplingOptions, bands] = loadVariables(sampling_filename, {'findSamplingOptions', 'bands'});
    n_bands = length(bands);
    if convert_to_color
        color_weights = colorWeights(...
            sensor_map, bands_color, bands, findSamplingOptions...
        );
    end
end

if is_spectral && convert_to_color
    dispersion_options = struct(...
        'int_method', findSamplingOptions.int_method,...
        'support_threshold', findSamplingOptions.support_threshold,...
        'bands_padding', findSamplingOptions.bands_padding,...
        'interpolant', findSamplingOptions.interpolant,...
        'interpolant_ref', findSamplingOptions.interpolant_ref...
    );
    if disable_interpolation
        dispersion_options.resolution = 0;
    else
        dispersion_options.resolution = dispersionfunToMatrixOptions.resolution;
    end
elseif is_spectral && ~disable_interpolation
    dispersion_options = struct(...
        'resolution', dispersionfunToMatrixOptions.resolution,...
        'int_method', findSamplingOptions.int_method,...
        'support_threshold', findSamplingOptions.support_threshold,...
        'bands_padding', findSamplingOptions.bands_padding,...
        'interpolant', findSamplingOptions.interpolant,...
        'interpolant_ref', findSamplingOptions.interpolant_ref,...
        'bands_in', bands...
    );
else
    dispersion_options = struct('bands_in', bands);
end

%% Process the images

time = zeros(n_images, 1);
for i = 1:n_images
    [I, name] = loadImage(image_filenames{i}, input_images_variable_name);
    
    [dispersionfun, I] = makeDispersionForImage(...
        dispersion_data, I, transform_data, true...
    );
    
    time_start = tic;
    if is_spectral && convert_to_color
        I = warpImageSpectral(...
            I, bands, sensor_map, bands_color, dispersion_options,...
            imageFormationPatchOptions, dispersionfun...
        );        
    else
        I = dispersionfunToMatrix(...
            dispersionfun, dispersion_options, I, false...
        );
        if convert_to_color
            I = channelConversion(I, color_weights);
        end
    end
    time(i) = toc(time_start);
    saveImages(...
        'data', output_directory, name,...
        I, '_unwarped', 'I_unwarped'...
    );
end

%% Save parameters to a file
save_variables_list = [ parameters_list, {...
        'image_filenames', 'time'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByWarpingData.mat');
save(save_data_filename, save_variables_list{:});