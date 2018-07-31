%% Demosaicing and basic correction of chromatic aberration
% Convert RAW images to colour images, and correct chromatic aberration by
% warping colour channels.
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
% being loaded. For image format files, images will simply be loaded with
% the Image Processing Toolbox 'imread()' function. For '.mat' files, the
% variable to be loaded must be provided in the script parameters.
%
% All images are expected to have 3 colour channels (Red, Green, Blue)
% (represented in a Bayer pattern as a 2D array). However, the colour
% channels can correspond to narrowband wavelength ranges - This script
% will input the wavelengths corresponding to the colour channels.
%
% Images may have different pixel dimensions, provided that they are
% compatible with the input model of chromatic aberration described below.
%
% ### Model of chromatic aberration
% A '.mat' file containing several variables, which is the output of
% 'DoubleConvexThickLensDiskDispersion.m', 'RAWDiskDispersion.m',
% 'DoubleConvexThickLensDispersion.m' or 'BimaterialImages.m', for example.
% The following variables are required:
% - 'dispersion_data': A model of chromatic aberration, modelling the
%   warping from the reference colour channel to the other colour channels.
%   `dispersion_data` can be converted to a function form using
%   `dispersionfun = makeDispersionfun(dispersion_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `true`.
% - 'bands': A vector containing the colour channel indices to use as the
%   `lambda` input argument of 'dispersionfunToMatrix()'. `bands` is the
%   colour channel information needed to evaluate the dispersion model.
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
% ## Output
%
% ### Colour images
% One of each of the following types of images is created for each input
% image. The filename of the input image is represented by '*' below.
% - '*_roi.tif' and '*_roi.mat': A cropped version of the input image,
%   (stored in the variable 'I_raw') containing the portion to be
%   corrected. This region of interest was determined using the
%   `model_space` and `fill` variables saved in the input model of
%   chromatic aberration data file (see above). If these variables were not
%   present, the cropped region is the entire input image. All of the other
%   output images listed below are limited to the region shown in
%   '*_roi.tif'.
% - '*_color_bilinear.tif' and '*_color_bilinear.mat': A colour image
%   (stored in the variable 'I_color_bilinear') created by bilinear
%   interpolation of each colour channel using 'bilinearDemosaic()'.
% - '*_color_demosaic.tif' and '*_color_demosaic.mat': A colour image
%   (stored in the variable 'I_color_demosaic') created by demosaicing using
%   MATLAB's 'demosaic()' function.
% - '*_color_bilinear_warped.tif' and '*_color_bilinear_warped.mat': A
%   colour image (stored in the variable 'I_color_bilinear_warped') created
%   by bilinear interpolation of each colour channel using
%   'bilinearDemosaic()', followed by warping to correct for chromatic
%   aberration.
% - '*_color_demosaic_warped.tif' and '*_color_demosaic_warped.mat': A
%   colour image (stored in the variable 'I_color_demosaic_warped') created
%   by demosaicing using MATLAB's 'demosaic()' function, followed by
%   warping to correct for chromatic aberration.
%
% Both '.mat' and '.tif' files are output, where applicable, for
% monochromatic or three-channel images, to provide both easy display
% ('.tif' files) and lossless storage ('.mat' files).
%
% ### Parameters
%
% A '.mat' file containing the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## Notes
% - The image colour space is not altered by this script. See 'imreadRAW()'
%   for code to convert an image to sRGB after demosaicing.
% - There is some quantization error in the images demosaiced with MATLAB's
%   'demosaic()' function, because this function only operates on integer
%   datatypes.
%
% ## References
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 22, 2018

% List of parameters to save with results
parameters_list = {
        'forward_dispersion_model_filename',...
        'output_directory'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
% '.mat' or image files can be loaded
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180709_TestingSplineModels/ground_truth/splines/swirly_0138_raw_warped.mat';
input_images_variable_name = 'raw_2D'; % Can be empty unless loading images from '.mat' files

% Model of dispersion
forward_dispersion_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180709_TestingSplineModels/DoubleConvexThickLensDispersionResults_spline_modelFromReference_true_fill.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/Downloads';

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load dispersion model

[...
    dispersion_data, bands, transform_data...
] = loadDispersionModel(forward_dispersion_model_filename, true);

%% Process the images

n_demosaicing_methods = 2;
for i = 1:n_images
    [I_raw, name] = loadImage(image_filenames{i}, input_images_variable_name);
    
    if ~ismatrix(I_raw)
        error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
    end
    
    [dispersionfun, I_raw] = makeDispersionForImage(...
        dispersion_data, I_raw, transform_data...
    );    
    saveImages(...
        output_directory, name,...
        I_raw, '_roi', 'I_raw'...
    );
    
    image_sampling_in = size(I_raw);
    W = dispersionfunToMatrix(...
        dispersionfun, bands, image_sampling_in, image_sampling_in,...
        [0, 0, image_sampling_in(2),  image_sampling_in(1)], false...
    );

    for j = 1:n_demosaicing_methods
        if j == 1
            I_color_warped = bilinearDemosaic(I_raw, bayer_pattern);
        elseif j == 2
            if isa(I_raw, 'double')
                I_raw_int = im2uint16(I_raw);
                I_color_warped = im2double(demosaic(I_raw_int, bayer_pattern));
            else
                I_color_warped = im2double(demosaic(I_raw, bayer_pattern));
            end
        else
            error('No demosaicing method associated with index %d.', j);
        end
        I_color = warpImage(I_color_warped, W, image_sampling_in);
        
        % Save the results
        if j == 1
            saveImages(...
                output_directory, name,...
                I_color, '_color_bilinear', 'I_color_bilinear',...
                I_color_warped, '_color_bilinear_warped', 'I_color_bilinear_warped'...
            );
        elseif j == 2
            saveImages(...
                output_directory, name,...
                I_color, '_color_demosaic', 'I_color_demosaic',...
                I_color_warped, '_color_demosaic_warped', 'I_color_demosaic_warped'...
            );
        else
            error('No demosaicing method associated with index %d.', j);
        end
    end
end

%% Save parameters to a file
save_variables_list = [ parameters_list, {...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByWarpingParameters.mat');
save(save_data_filename, save_variables_list{:});