%% Real image dataset preparation
% Preprocess RAW images by subtracting dark frames, averaging to reduce
% noise, and reducing clipping by merging multiple exposures.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
% All input images are RAW images to be loaded with 'imreadRAW()':
% - Images to be used as experimental data. For each scene, multiple
%   exposure settings may have been used to improve dynamic range.
% - Dark frame images, taken under each of the exposure settings used to
%   capture the data images.
%
% Several replicates of each image can be loaded and averaged to reduce
% noise.
%
% ## Output
%
% ### Output images
%
% Depending on the parameters below, the following two types of images are
% saved to '.mat' files in double precision format:
% - Images corrected by dark frame subtraction, and blended between
%   exposures to reduce clipping.
% - Images averaged across replicates before being processed as above.
%
% The first type of images is recommended as input to image restoration
% algorithms. The second type of images is recommended as input to
% calibration algorithms and as ground truth. For either type of images,
% dark frame subtraction is performed by subtracting dark frames that were
% averaged across replicates. Images may have small negative values because
% of dark frame subtraction, but are normalized to have maximum values no
% greater than one.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'peaks': A vector containing the peak values used to normalize images,
%   output by 'blendExposures()'. Refer to the documentation of
%   'blendExposures.m' for details.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## References
% The approach taken in this script is inspired by:
%
%   Darrodi, M. M., Finlayson, G., Goodman, T., & Mackiewicz, M. (2015).
%   Reference data set for camera spectral sensitivity estimation. Journal
%   of the Optical Society of America A: Optics and Image Science, and
%   Vision, 32(3), 381-391. doi:10.1364/JOSAA.32.000381

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 4, 2018

% List of parameters to save with results
parameters_list = {
        'darkSubtract_dir',...
        'var_name',...
        'wildcards',...
        'darkSubtract_regex',...
        'blendExposures_dir',...
        'blendExposures_regex',...
        'range',...
        'align'...
    };

%% Input data and parameters

% ## Input arguments for 'darkSubtract()'

% Directories in which to store dark-subtracted images
darkSubtract_dir.out_averaged = '/home/llanos/Downloads/data/dark_subtracted_averaged'; % Averaged images
darkSubtract_dir.out_single = '/home/llanos/Downloads/data/dark_subtracted_original'; % Non-averaged images

% Image variable name
var_name = 'I_raw';

% Wildcards for 'ls()' to find the RAW images to process
% Data images
wildcards.in = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/images_withProjector/*1mmDots*.tif';
% Dark frame images
wildcards.dark = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/images_withProjector/*dark*.tif';

% Regular expression for removing the portion of a filename that differs
% between replicates of an image
darkSubtract_regex.dedup = '_\d{4}-\d{2}-\d{2}-\d{6}-\d{4}';
% Regular expression for extracting the portion of a filename that must
% match between an image and the corresponding dark frame
darkSubtract_regex.dedup = '(pos\d+).*(\d+ms)';

% ## Input arguments for 'blendExposures()'

% Directories in which to store final images
blendExposures_dir.out_reference = '/home/llanos/Downloads/data/blended_averaged'; % Averaged images
blendExposures_dir.other_paths = '/home/llanos/Downloads/data/blended_original'; % Non-averaged images

% Regular expressions identifying exposure settings. Each cell vector is a
% group of exposures that can be blended together. Within each group, the
% exposures must be ordered from lowest to highest.
blendExposures_regex = {...
    {'50ms', '200ms', '400ms', '600ms'},...
    {'0030ms', '0325ms'}...
};

% Range of pixel values used to calibrate scaling factors between exposures
range = [0, 1];

% Colour-filter pattern code
align = 'gbrg';

%% Load and preprocess the images

darkSubtract_output_files = darkSubtract(...
    darkSubtract_dir, var_name, wildcards, darkSubtract_regex...
);

[output_files, peaks] = blendExposures(...
    blendExposures_dir, var_name,...
    darkSubtract_output_files.out_averaged,...
    darkSubtract_output_files.out_single,...
    blendExposures_regex, range, align...
);

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'peaks'...
    } ];
save_data_filename = fullfile(output_directory, 'PreprocessRAWImages.mat');
save(save_data_filename, save_variables_list{:});