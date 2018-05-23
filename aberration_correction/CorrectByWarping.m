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
% Toolbox 'imread()' function. All images are expected to have the same
% pixel dimensions and 3 colour channels (Red, Green, Blue) (represented in
% a Bayer pattern as a 2D array).
%
% ### Polynomial model of chromatic aberration
% A '.mat' file containing several variables, which is the output of
% 'DoubleConvexThickLensDiskDispersion.m', 'RAWDiskDispersion.m' or
% 'DoubleConvexThickLensDispersion.m'. The following variables are
% required:
% - 'polyfun_data': A polynomial model of chromatic aberration, modeling the
%   warping from the reference colour channel to the other colour channels.
%   `polyfun_data` can be converted to a function form using `polyfun =
%   makePolyfun(polyfun_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `true`.
%
% ## Output
%
% ### Colour images
% One of the following types of images is created for each input image. The
% filename of the input image is represented by '*' below.
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
% Additionally, the file contains a 'bands' variable used as the 'lambda'
% input argument of 'polyfunToMatrix()', for reference. 'bands' is equal to
% `[1 2 3]`, representing the Red, Green, and Blue colour channels.
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
input_images_wildcard = '';

% Colour-filter pattern
bayer_pattern = 'gbrg';

% Polynomial model of dispersion
polynomial_model_filename = '';

% Output directory for all images and saved parameters
output_directory = '';

%% Find the images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

%% Load dispersion model

model_variables_required = { 'polyfun_data', 'model_from_reference' };
load(polynomial_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the dispersion model variables is not loaded.')
end
if ~model_from_reference
    error('Dispersion model is in the wrong frame of reference.')
end

%% Process the images

polyfun = makePolyfun(polyfun_data);
bands = 1:3;
ext = '.tif';

n_demosaicing_methods = 2;
for i = 1:n_images
    [~, name] = fileparts(image_filenames{i});
    I_raw = imread(image_filenames{i});
    if i == 1
        image_sampling_in = size(I_raw);
        if length(image_sampling_in) > 2
            error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
        end
        W = polyfunToMatrix(...
            polyfun, bands, image_sampling_in, image_sampling_in,...
            [0, 0, image_sampling_in], false...
        );
    else
        if any(image_sampling_in ~= size(I_raw))
            error('Not all images have the same dimensions.');
        end
    end
    for j = 1:n_demosaicing_methods
        if j == 1
            I_color = bilinearDemosaic(I_raw, bayer_pattern);
        elseif j == 2
            I_color = demosaic(I_raw, bayer_pattern);
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
        imwrite(I_color, I_color_filename);
        imwrite(I_color_warped, I_color_warped_filename);
    end
end

%% Save parameters to a file
save_variables_list = [ parameters_list, {...
        'bands',...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectByWarpingParameters.mat');
save(save_data_filename, save_variables_list{:});