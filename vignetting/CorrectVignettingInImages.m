%% Vignetting correction in a set of images
% Create a model of vignetting from a reference image, and then use it to
% correct both the reference image, and a set of other images.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
%
% Two filepath wildcards, pointing to either monochromatic, or multi-channel,
% reference, and other images, respectively.
%
% #### Reference image
%
% A single reference image is expected. If the reference image is a
% multi-channel image, a single channel will be used to calibrate the
% correction. For both monochromatic and multi-channel reference images, the
% vignetting correction will be calibrated only from pixels in the Green filter
% locations of the Bayer pattern (as implemented in 'vignettingPolyfit()'. Even
% if the sensor lacks a Bayer filter, the smaller number of pixels used for
% calibration is useful for reducing computation time. If the reference image
% was produced by demosaicing, then only using pixels at the original Green
% filter locations may provide some robustness to demosaicing artifacts (which
% should not have a significant impact, in any case).
%
% #### Other images
%
% The other images will be corrected for vignetting according to the model of
% vignetting created from the reference image. As 'correctVignetting()' applies
% the same correction to all image channels, the other images need not have the
% same numbers of channels as the reference image.
%
% ### Vignetting calibration mask
%
% A filepath wildcard pointing to a single-channel image where bright pixels
% indicate regions that should all have the same intensity in the reference
% image. The vignetting model will fit such that these pixels will have
% intensities close to one after the image is corrected for vignetting.
%
% ## Output
%
% ### Corrected images
%
% A corrected image is output for each input image, including the reference
% image. The output files are given the same names as the input files, postfixed
% with '_vc'. For input images with two or more than three channels, the output
% images are in '.mat' format. Both '.mat' and '.tif' files are output for
% monochromatic or three-channel images to provide both easy display ('.tif'
% files) and lossless storage ('.mat' files). In this case, the '.tif' files do
% not contain the original data, but contain versions of the output images
% post-processed with 'clipAndRemap()' for better viewing.
%
% The variable names under which the images are stored in '.mat' files are
% specified in the parameters section of this script.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'polyfun_data': A structure containing the model of vignetting used to
%   correct the images. 'polyfun_data' is the output argument of the same name
%   of 'vignettingPolyfit()'. Refer to the documentation of
%   'vignettingPolyfit.m' for details.
% - 'reference_filename': A cell vector containing the reference image filename
%   retrieved based on the wildcard provided in the parameters section of the
%   script.
% - 'mask_filename': A cell vector containing the mask image filename
%   retrieved based on the wildcard provided in the parameters section of the
%   script.
% - 'other_filenames': A cell vector containing the other images' filenames
%   retrieved based on the wildcard provided in the parameters section of the
%   script.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 27, 2019

% List of parameters to save with results
parameters_list = {
        'mask_threshold',...
        'bayer_pattern',...
        'reference_channel',...
        'max_degree',...
        'quantiles',...
        'output_directory'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the reference image
reference_wildcard = fullfile('.', 'demo_data', 'hdr_averaged_images', 'd1_disks32cmV2_550nm.mat');
reference_variable_name = 'I_raw'; % Used for input and output '.mat' file

% Wildcard for 'ls()' to find the mask
mask_wildcard = fullfile('.', 'demo_data', 'hdr_averaged_images', 'd1_disks32cmV2_maskVignetting.png');
mask_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the other images (can be empty)
other_wildcard = [];
other_variable_name = 'I_raw'; % Used for input and output '.mat' files

% Threshold used to binarize the mask image, if it is not already a binary
% image.
mask_threshold = 0.5; % In a range of intensities from 0 to 1

bayer_pattern = 'gbrg'; % Colour-filter pattern

% Channel to use to calibrate the vignetting model. If the reference image is a
% single channel image, the single channel will be used regardless of the value
% of this parameter.
reference_channel = 2;

% Maximum degree of the polynomial model of vignetting
max_degree = 5;

% Quantiles used for clipping to produce nice output images in image-format
% files
quantiles = [0.01, 0.99];

% Output directory
output_directory = '${DIRPATH}';

% ## Debugging Flags
vignettingPolyfitVerbose = true;

%% Calibrate vignetting

reference_filename = listFiles(reference_wildcard);
if length(reference_filename) ~= 1
    error('A single reference image is expected.');
end
I_reference = loadImage(reference_filename{1}, reference_variable_name);

mask_filename = listFiles(mask_wildcard);
if length(mask_filename) ~= 1
    error('A single mask image is expected.');
end
mask = loadImage(mask_filename{1}, mask_variable_name);
if size(mask, 3) ~= 1
    error('Expected the vignetting calibration mask, "%s", to have only one channel.', mask_filename);
end
if ~islogical(mask)
    mask = imbinarize(mask, mask_threshold);
end

if size(I_reference, 3) == 1
    I_reference_c = I_reference;
else
    I_reference_c = I_reference(:, :, reference_channel);
end
[polyfun, polyfun_data] = vignettingPolyfit(...
    I_reference_c, mask, max_degree, bayer_pattern, vignettingPolyfitVerbose...
);

%% Correct all images for vignetting

other_filenames = listFiles(other_wildcard);
all_filenames = [reference_filename, other_filenames];
for i = 1:length(all_filenames)
    if i == 1
        I = I_reference;
        variable_name = reference_variable_name;
    else
        I = loadImage(all_filenames{i}, other_variable_name);
        variable_name = other_variable_name;
    end
    I = correctVignetting(I, polyfun);
    [~, I_out_filename] = fileparts(all_filenames{i});
    saveImages(...
        'data', output_directory, I_out_filename,...
        I, '_vc', variable_name...
    );
    if any(size(I, 3) == [1 3])
        I_out_remapped = clipAndRemap(I, 'uint8', 'quantiles', quantiles);
        saveImages(...
            'image', output_directory, I_out_filename,...
            I_out_remapped, '_vc', []...
        );
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'polyfun_data',...
        'reference_filename',...
        'mask_filename',...
        'other_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'CorrectVignettingInImagesData.mat');
save(save_data_filename, save_variables_list{:});