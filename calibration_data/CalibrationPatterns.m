%% Generate Chromatic Aberration Calibration Patterns
% Save images for display during disparity and defocus calibration
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
% Two versions of each image will be output: One for a digital display, and
% one for printing. Additionally, the individual colour channels of each
% image for digital display will be saved as separate files.
%
% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 31, 2017

%% Input data and parameters

% ## Digital display parameters

% Image pixel dimensions (height, width)
% - [1080 1920] for the LG 27MP37HQ monitor
% - [1050 1680] for the Acer AL2216W monitor
image_size = [1050 1680];

% Image resolution (pixels per inch)
resolution_ppi = [];
% or pixel pitch (mm per pixel)
% - 0.311 mm for the LG 27MP37HQ monitor
% - 0.282 mm for the Acer AL2216W monitor
resolution_pitch = 0.282;

% Filename Suffix
filename_suffix_digital = '_Acer';

% ## Printing parameters

% Paper size (width, height), in inches
paper_size_inches = [8.5, 11];

% Paper orientation
paper_orientation = 'landscape';

% Figure size on page (width, height), in inches (Allow at least 0.25
% inches on all sides)
figure_size_printed = paper_size_inches - 0.25 * 2;

% Filename Suffix
filename_suffix_print = '_print';

% ## Common parameters

% Disk radius in millimetres
disk_radius = 4;

% Disk separation in millimetres
disk_separation = 11;

% Output Directory
output_directory = '/home/llanos/Downloads/patterns';

%% Unit conversions
inchPerMM = unitsratio('inch', 'mm');
MMPerInch = unitsratio('mm', 'inch');
if ~isempty(resolution_ppi) && ~isempty(resolution_pitch)
    error('Specify only one of `resolution_ppi` or `resolution_pitch`.');
elseif ~isempty(resolution_ppi)
    use_ppi = true;
    resolution_mm = resolution_ppi * inchPerMM;
elseif ~isempty(resolution_pitch)
    use_ppi = false;
    resolution_mm = 1 / resolution_pitch;
    resolution_ppi = 1 / (resolution_pitch * inchPerMM);
else
    error('Specify one of `resolution_ppi` or `resolution_pitch`.');
end

%% Uniform calibration patterns

[I_black, I_black_rgb, fg_black] = uniformImage( image_size, 0 );
[I_white, I_white_rgb, fg_white] = uniformImage( image_size, 1 );

%% Disk calibration patterns

[I_disks_black, I_disks_black_rgb] = diskImage(...
    image_size, resolution_mm, disk_radius, disk_separation,...
    1, 0 ...
);
[I_disks_white, I_disks_white_rgb] = diskImage(...
    image_size, resolution_mm, disk_radius, disk_separation,...
    0, 1 ...
);

%% Save images

tiff_extension = [filename_suffix_digital '.tiff'];
print_extension = [filename_suffix_print '.pdf'];

param_str = [
    '_', num2str(image_size(1)), 'x', num2str(image_size(2))
];

if use_ppi
    param_str = [
        param_str, 'resPPI',...
        num2str(resolution_ppi)...
    ];
else
    param_str = [
        param_str, 'resPitch',...
        num2str(resolution_pitch)...
    ];
end

disks_param_str = [
    param_str, 'r', num2str(disk_radius), 'sep', num2str(disk_separation)...
    ];

general_tiff_options = {
    'Resolution', resolution_ppi
};

n_channels = 3;
channel_labels = {'R', 'G', 'B'};

% Uniform black
tiff_options = [general_tiff_options {
    'Description', 'Uniform black'
}];

imwrite(...
    I_black,...
    fullfile(output_directory, ['black' param_str tiff_extension]),...
    tiff_options{:}...
    );

printFigure(...
    fullfile(output_directory, ['black' param_str print_extension]),...
    paper_size_inches, paper_orientation, figure_size_printed, fg_black...
)

% Uniform white
tiff_options = [general_tiff_options {
    'Description', 'Uniform white'
}];

imwrite(...
    I_white,...
    fullfile(output_directory, ['white' param_str tiff_extension]),...
    tiff_options{:}...
    );

printFigure(...
    fullfile(output_directory, ['white' param_str print_extension]),...
    paper_size_inches, paper_orientation, figure_size_printed, fg_white...
)

for i = 1:n_channels
    imwrite(...
        I_white_rgb{i},...
        fullfile(output_directory, ['white', channel_labels{i}, param_str, tiff_extension]),...
        tiff_options{:}...
    );
end

% Black disks
tiff_options = [general_tiff_options {
    'Description', 'Black disks'
}];

imwrite(...
    I_disks_black,...
    fullfile(output_directory, ['disksBlack' disks_param_str tiff_extension]),...
    tiff_options{:}...
    );

for i = 1:n_channels
    imwrite(...
        I_disks_black_rgb{i},...
        fullfile(...
            output_directory,...
            [...
                'disksBlack', channel_labels{i}, disks_param_str,...
                tiff_extension...
            ]...
        ),...
        tiff_options{:}...
    );
end

% White disks
tiff_options = [general_tiff_options {
    'Description', 'White disks'
}];

imwrite(...
    I_disks_white,...
    fullfile(output_directory, ['disksWhite' disks_param_str tiff_extension]),...
    tiff_options{:}...
    );

for i = 1:n_channels
    imwrite(...
        I_disks_white_rgb{i},...
        fullfile(...
            output_directory,...
            [...
                'disksWhite', channel_labels{i}, disks_param_str,...
                tiff_extension...
            ]...
        ),...
        tiff_options{:}...
    );
end