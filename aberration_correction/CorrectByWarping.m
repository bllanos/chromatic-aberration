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
% wavelengths at which the spectral images are defined.
%
% Images may have different pixel dimensions, provided that they are
% compatible with the input model of chromatic aberration described below.
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
% ## Output
%
% ### Warped images
%
% For each input image, a '.mat' file with a name ending in '_unwarped' is
% created. The 'I_unwarped' variable in the file contains the corrected
% image. Images are cropped to the region in which the dispersion model is
% valid.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'image_filenames': A cell vector containing the input image filenames
%   retrieved based on the wildcard provided in the parameters section of
%   the script.
% - 'time': The time, in seconds, taken to produce the output images.
%   'time' is a vector with a length equal to the number of images.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
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
    'rgb_mode',...
    'bands',...
    'bands_filename',...
    'bands_variable',...
    'forward_dispersion_model_filename'...
};

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
% '.mat' or image files can be loaded
input_images_wildcard = 'C:\Users\GraphicsLab\Documents\llanos\Results\run_on_dataset_ignoreDispersion\*_latent.mat';
input_images_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Indicate whether the images are to be treated as colour images or
% spectral images
rgb_mode = false;

if rgb_mode
    n_bands = 3;
    bands = (1:n_bands).';
    bands_filename = []; % Not used
    bands_variable = []; % Not used
else
    % Path and filename of a '.mat' file containing the wavelengths corresponding to
    % the spectral image.
    bands_filename = 'C:\Users\GraphicsLab\Documents\llanos\Results\run_on_dataset_ignoreDispersion\RunOnDataset_20190208_ComputarLens_rawCaptured_ignoreDispersion.mat';
    bands_variable = 'bands'; % Variable name in the above file
end

% Model of dispersion
forward_dispersion_model_filename = 'C:\Users\GraphicsLab\Documents\llanos\Data\20190208_ComputarLens\dispersion\spectral\full_image\RAWDiskDispersionResults_spectral_polynomial_fromReference.mat';

% Output directory for all images and saved parameters
output_directory = 'C:\Users\GraphicsLab\Documents\llanos\Results\run_on_dataset_ignoreDispersion_warpCorrected';

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load dispersion model

[...
    dispersion_data, bands_dispersionfun, transform_data...
] = loadDispersionModel(forward_dispersion_model_filename, true, true);

if rgb_mode
    if ((length(bands) ~= length(bands_dispersionfun)) ||...
       any(bands(:) ~= bands_dispersionfun(:)))
        error('When warping colour images, the same colour channels must be used by the model of dispersion.');
    end
    
    dispersion_options = struct('bands_in', bands);
else
    load(bands_filename, bands_variable);
    if exist(bands_variable, 'var')
        bands = eval(bands_variable);
    end
    if ~exist(bands_variable, 'var') || isempty(bands)
        error('No wavelengths loaded.')
    end
    n_bands = length(bands);
    
    dispersion_options = dispersionfunToMatrixOptions;
    dispersion_options.bands_in = bands;
end

%% Process the images

time = zeros(n_images, 1);
for i = 1:n_images
    [I, name] = loadImage(image_filenames{i}, input_images_variable_name);
    
    [dispersionfun, I] = makeDispersionForImage(...
        dispersion_data, I, transform_data, true...
    );
    
    time_start = tic; 
    I = dispersionfunToMatrix(...
        dispersionfun, dispersion_options, I, false...
    );
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