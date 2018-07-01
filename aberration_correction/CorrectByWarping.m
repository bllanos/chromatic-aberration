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
% being loaded. Images will simply be loaded with the Image Processing
% Toolbox 'imread()' function. All images are expected to have 3 colour
% channels (Red, Green, Blue) (represented in a Bayer pattern as a 2D
% array). However, the colour channels can correspond to narrowband
% wavelength ranges - This script will input the wavelengths corresponding
% to the colour channels.
%
% Images may have different pixel dimensions, provided that they are
% compatible with the input model of chromatic aberration described below.
%
% ### Polynomial model of chromatic aberration
% A '.mat' file containing several variables, which is the output of
% 'DoubleConvexThickLensDiskDispersion.m', 'RAWDiskDispersion.m',
% 'DoubleConvexThickLensDispersion.m' or 'BimaterialImages.m', for example.
% The following variables are required:
% - 'polyfun_data': A polynomial model of chromatic aberration, modeling the
%   warping from the reference colour channel to the other colour channels.
%   `polyfun_data` can be converted to a function form using `polyfun =
%   makePolyfun(polyfun_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `true`.
% - 'bands': A vector containing the wavelengths or colour channel indices
%   to use as the `lambda` input argument of 'polyfunToMatrix()'. `bands`
%   is the wavelength or colour channel information needed to evaluate the
%   dispersion model.
%
% The following two additional variables are optional. If they are present,
% they will be used for the following purposes:
% - Conversion between the coordinate system in which the polynomial model
%   of chromatic aberration was constructed and the image coordinate
%   system.
% - Limiting the correction of chromatic aberration to the region in which
%   the polynomial model is valid.
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
% One of the following types of images is created for each input image. The
% filename of the input image is represented by '*' below.
% - '*_roi.tif': A cropped version of the input image, containing the
%   portion to be corrected. This region of interest was determined using
%   the `model_space` and `fill` variables saved in the input polynomial
%   model of chromatic aberration data file (see above). If these variables
%   were not present, the cropped region is the entire input image. All of
%   the other output images listed below are limited to the region shown in
%   '*_roi.tif'.
% - '*_color_bilinear.tif': A colour image created by bilinear
%   interpolation of each colour channel using 'bilinearDemosaic()'.
% - '*_color_demosaic.tif': A colour image created by demosaicing using
%   MATLAB's 'demosaic()' function.
% - '*_color_bilinear_warped.tif': A colour image created by bilinear
%   interpolation of each colour channel using 'bilinearDemosaic()',
%   followed by warping to correct for chromatic aberration.
% - '*_color_demosaic_warped.tif': A colour image created by demosaicing using
%   MATLAB's 'demosaic()' function, followed by warping to correct for
%   chromatic aberration.
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
        'bayer_pattern',...
        'polynomial_model_filename',...
        'output_directory'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process.
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180530_CorrectionMethodsBasicComparison/input_images/*raw*';

% Colour-filter pattern
bayer_pattern = 'gbrg';

% Polynomial model of dispersion
polynomial_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180530_CorrectionMethodsBasicComparison/RAWDiskDispersionResults_true.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180530_CorrectionMethodsBasicComparison/output_images_warping';

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load dispersion model

model_variables_required = { 'polyfun_data', 'model_from_reference', 'bands' };
model_variables_transform = { 'model_space', 'fill' };
load(polynomial_model_filename, model_variables_required{:}, model_variables_transform{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the dispersion model variables is not loaded.')
end
if ~model_from_reference
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

%% Process the images

ext = '.tif';

n_demosaicing_methods = 2;
for i = 1:n_images
    [~, name] = fileparts(image_filenames{i});
    I_raw = imread(image_filenames{i});
    if ~ismatrix(I_raw)
        error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
    end
    
    if crop_image
        [roi, T_roi] = modelSpaceTransform(size(I_raw), model_space, fill);
        polyfun = makePolyfun(polyfun_data, T_roi);
        I_raw = I_raw(roi(1):roi(2), roi(3):roi(4));
    else
        polyfun = makePolyfun(polyfun_data);
    end
    
    I_raw_filename = fullfile(output_directory, [name '_roi' ext]);
    imwrite(I_raw, I_raw_filename);
    
    image_sampling_in = size(I_raw);
    W = polyfunToMatrix(...
        polyfun, bands, image_sampling_in, image_sampling_in,...
        [0, 0, image_sampling_in(2),  image_sampling_in(1)], false...
    );

    for j = 1:n_demosaicing_methods
        if j == 1
            I_color = bilinearDemosaic(I_raw, bayer_pattern);
        elseif j == 2
            I_color = im2double(demosaic(I_raw, bayer_pattern));
        else
            error('No demosaicing method associated with index %d.', j);
        end
        I_color_warped = warpImage(I_color, W, image_sampling_in);
        
        % Save the results
        if j == 1
            I_color_filename = [name '_color_bilinear' ext];
            I_color_warped_filename = [name '_color_bilinear_warped' ext];
        elseif j == 2
            I_color_filename = [name '_color_demosaic' ext];
            I_color_warped_filename = [name '_color_demosaic_warped' ext];
        else
            error('No demosaicing method associated with index %d.', j);
        end
        I_color_filename = fullfile(output_directory, I_color_filename);
        I_color_warped_filename = fullfile(output_directory, I_color_warped_filename);
        imwrite(I_color, I_color_filename);
        imwrite(I_color_warped, I_color_warped_filename);
    end
end

%% Save parameters to a file
save_variables_list = [ parameters_list, {...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByWarpingParameters.mat');
save(save_data_filename, save_variables_list{:});