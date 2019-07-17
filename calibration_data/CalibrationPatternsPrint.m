%% Generate Chromatic Aberration Calibration Patterns
% Save images to serve as printed targets during disparity and defocus
% calibration
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
% - Uniformly black image
% - Uniformly white image
% - Pattern of black disks on a white background
% - Pattern of white disks on a black background
%
% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 25, 2018

%% Input data and parameters

% Disk radius in millimetres
% Unfortunately, 'diskImagePrint()' seems to produce disks which are
% slightly smaller than specified here.
disk_radius = 2;

% Disk separation in millimetres
disk_separation = 8;

% Output Directory
output_directory = '${DIRPATH}';

% ## Printing parameters

% Paper size (width, height), in inches
paper_size_inches = [18, 28];

% Paper orientation
paper_orientation = 'landscape';

% Image size on page (width, height), in inches, independent of how the
% paper is oriented (portrait or landscape). (Allow at least 0.25 inches on
% all sides)
image_size = [28, 18] - 0.25 * 2;

% Filename Suffix
filename_suffix = '_print';

%% Uniform calibration patterns

fg_black = uniformImagePrint( image_size, 0 );
fg_white = uniformImagePrint( image_size, 1 );

%% Disk calibration patterns

fg_disks_black = diskImagePrint(...
    image_size, disk_radius, disk_separation,...
    1, 0 ...
);
fg_disks_white = diskImagePrint(...
    image_size, disk_radius, disk_separation,...
    0, 1 ...
);

%% Save images

extension = [filename_suffix '.eps'];

param_str = [
    '_', num2str(image_size(1)), 'x', num2str(image_size(2))
];

disks_param_str = [
    param_str, 'r', num2str(disk_radius), 'sep', num2str(disk_separation)...
    ];

% Uniform black
printFigure(...
    fullfile(output_directory, ['black' param_str extension]),...
    paper_size_inches, paper_orientation, image_size, fg_black...
)

% Uniform white
printFigure(...
    fullfile(output_directory, ['white' param_str extension]),...
    paper_size_inches, paper_orientation, image_size, fg_white...
)

% Black disks
printFigure(...
    fullfile(output_directory, ['disksBlack' disks_param_str extension]),...
    paper_size_inches, paper_orientation, image_size, fg_disks_black...
)

% White disks
printFigure(...
    fullfile(output_directory, ['disksWhite' disks_param_str extension]),...
    paper_size_inches, paper_orientation, image_size, fg_disks_white...
)

close all