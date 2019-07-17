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
% - 'scaling_factors': A cell vector containing the multipliers used by
%   'blendExposures()' to convert between exposures. Refer to the
%   documentation of 'blendExposures.m' for details.
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
        'radius',...
        'bayer_pattern'...
    };

%% Input data and parameters

% ## Input arguments for 'darkSubtract()'

% Directories in which to store dark-subtracted images
darkSubtract_dir.out_averaged = fullfile('.', 'demo_data', 'averaged_images'); % Averaged images
%darkSubtract_dir.out_single = '${DIRPATH}'; % Non-averaged images

% Image variable name to use
var_name = 'I_raw';

% Wildcards for 'ls()' to find the RAW images to process
% Data images
wildcards.in = fullfile('.', 'demo_data', 'captured_images', 'images_filtered_light', 'disks', '*.tif');
% Dark frame images
wildcards.dark = fullfile('.', 'demo_data', 'captured_images', 'dark_frames', '*.tif');

% Regular expression for removing the portion of a filename that differs
% between replicates of an image
darkSubtract_regex.dedup = '_\d{4}-\d{2}-\d{2}-\d{6}-\d{4}';
% Regular expression for extracting the portion of a filename that must
% match between an image and the corresponding dark frame
darkSubtract_regex.dark_match = 'd(\d).+(_[\d.]+ms)';

% ## Input arguments for 'blendExposures()'

% Directories in which to store final images
blendExposures_dir.out_reference = fullfile('.', 'demo_data', 'hdr_averaged_images'); % Averaged images
% Non-averaged images would need to have identical filenames other than the
% exposures in order, to be processed correctly by 'blendExposures()'
% (otherwise they are identified as from different "scenes"). This would
% mean selecting one replicant at random per exposure, to avoid filename
% conflicts after removing the timestamps and sequence numbers from
% filenames. But, in any case, merging images across exposures changes
% image noise characteristics, whereas the point of having non-averaged
% images is to preserve noise characteristics. Therefore, there is no reason to
% process the non-averaged images.
%blendExposures_dir.other_paths = '${DIRPATH}';

% Regular expressions identifying exposure settings. Each cell vector is a
% group of exposures that can be blended together. Within each group, the
% exposures must be ordered from lowest to highest.
blendExposures_regex = {...
    {'_25ms', '_50ms', '_100ms', '_250ms', '_500ms', '_2000ms', '_3916ms'}...
};

blendExposures_dir.out_reference = repmat({blendExposures_dir.out_reference}, length(blendExposures_regex), 1);
%blendExposures_dir.other_paths = repmat({blendExposures_dir.other_paths}, length(blendExposures_regex), 1);

% Range of pixel values used to calibrate scaling factors between exposures
range = [0.02, 0.95];

% Radius of the greyscale erosion disk structuring element used to mitigate blooming
% Set to zero for sensors which do not exhibit blooming (e.g. CMOS sensors)
radius = 1;

% Colour-filter pattern code
bayer_pattern = 'gbrg';

% ## Other parameters

% Directory in which to save the final output '.mat' file containing
% parameters and saved variables
output_directory = fullfile('.', 'demo_data', 'hdr_averaged_images');

% ## Debugging Flags

blendExposuresVerbose = true;

%% Load and preprocess the images

darkSubtract_output_files = darkSubtract(...
    darkSubtract_dir, var_name, wildcards, darkSubtract_regex...
);

[output_files, scaling_factors] = blendExposures(...
    blendExposures_dir, var_name,...
    darkSubtract_output_files.out_averaged,...
    {},... %darkSubtract_output_files.out_single,...
    blendExposures_regex, range, radius, bayer_pattern, blendExposuresVerbose...
);

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'scaling_factors'...
    } ];
save_data_filename = fullfile(output_directory, 'PreprocessRAWImages.mat');
save(save_data_filename, save_variables_list{:});
