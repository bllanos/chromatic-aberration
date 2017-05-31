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
% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 31, 2017

%% Input data and parameters

% Image pixel dimensions (height, width)
image_size = [1050 1680];

% Image resolution (pixels per inch)
resolution_ppi = [];
% or pixel pitch (mm per pixel)
resolution_pitch = 0.282;

% Disk radius in millimetres
disk_radius = 4;

% Disk separation in millimetres
disk_separation = 11;

% Output Directory
output_directory = 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170531_CalibrationPatterns';

%% Unit conversions
inchPerMM = unitsratio('inch', 'mm');
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

I_black = uniformImage( image_size, 0 );
I_white = uniformImage( image_size, 1 );

%% Disk calibration patterns

I_disks_black = diskImage(...
    image_size, resolution_mm, disk_radius, disk_separation,...
    1, 0 ...
);
I_disks_white = diskImage(...
    image_size, resolution_mm, disk_radius, disk_separation,...
    0, 1 ...
);

%% Save images

tiff_extension = '.tiff';

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

filename_base = [param_str tiff_extension];

disks_filename_base = [
    param_str, 'r', num2str(disk_radius), 'sep', num2str(disk_separation),...
    tiff_extension
    ];

general_tiff_options = {
    'Resolution', resolution_ppi
};

% Uniform black
tiff_options = [general_tiff_options {
    'Description', 'Uniform black'
}];

imwrite(...
    I_black,...
    fullfile(output_directory, ['black' filename_base]),...
    tiff_options{:}...
    );

% Uniform white
tiff_options = [general_tiff_options {
    'Description', 'Uniform white'
}];

imwrite(...
    I_white,...
    fullfile(output_directory, ['white' filename_base]),...
    tiff_options{:}...
    );

% Black disks
tiff_options = [general_tiff_options {
    'Description', 'Black disks'
}];

imwrite(...
    I_disks_black,...
    fullfile(output_directory, ['disksBlack' disks_filename_base]),...
    tiff_options{:}...
    );

% White disks
tiff_options = [general_tiff_options {
    'Description', 'White disks'
}];

imwrite(...
    I_disks_white,...
    fullfile(output_directory, ['disksWhite' disks_filename_base]),...
    tiff_options{:}...
    );