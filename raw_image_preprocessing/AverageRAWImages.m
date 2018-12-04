%% Average RAW images
% Use the `dirreadRAW()` function to average together raw images, and
% optionally perform demosaicing.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
%
% None, but image files are saved by `dirreadRAW()`.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2017

%% Input data and parameters

% Generate full-colour images (`true`) or colour-filter array data
% (`false`)
full_color = false;

% Image processing directions
% Refer to the documentation of `imreadRAW()` for details.
ops.linearize = true;
if full_color
    ops.demosaic = true;
    ops.convertColor = true;
    ops.wb = true;
    align = 'gbrg'; % Colour-filter pattern
    wb = [1 1 1]; % White-balance multipliers
else
    ops.demosaic = false;
    ops.convertColor = false;
    ops.wb = false;
end
% Display intermediate results during processing, for each image
verbose = false;

% Directory containing the input images
in_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/images_withProjector';

% Directory in which to save the output images
out_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/averaged_images';

% Output filename extension (no dot)
ext = 'tif';

% Input filename wildcard
wildcard = '*.tif';

% Portion of filename not used to determine if the image is unique
regex = '_\d{4}-\d{2}-\d{2}-\d{6}-\d{4}';

%% Process the images
if full_color
    output_files = dirreadRAW(...
        in_directory, out_directory, ext, wildcard, regex, ops, align, wb, verbose...
        );
else
    output_files = dirreadRAW(...
        in_directory, out_directory, ext, wildcard, regex, ops, verbose...
        );
end