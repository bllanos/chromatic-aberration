%% Visualize models of dispersion as point spread functions
% Synthesize image patches which show the chromatic dispersion pattern generated
% by a single-pixel spectral impulse, for both spectral and colour models of
% dispersion.
%
% ## References
%
% Colour correction calibration method:
%
%   Finlayson, G., Gong, H., and Fisher, R.B. "Color Homography: Theory and
%   Applications." IEEE Transactions on Pattern Analysis and Machine
%   Intelligence, vol. 41, no. 1, pp. 20-33, 2019.
%   doi:10.1109/TPAMI.2017.2760833

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 3, 2019

% Presently, this script just saves some key information, not a complete set of
% parameters
parameters_list = {};

%% Input data and parameters

% Model of spectral dispersion
spectral_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/spectral/full_image/RAWDiskDispersionResults_spectral_polynomial_fromNonReference.mat';

% Model of colour dispersion
color_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/rgb/full_image/RAWDiskDispersionResults_RGB_polynomial_fromNonReference.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/SonyColorMapData.mat';

% Upper bound on the image position shift caused by dispersion. This will become
% the patch half-side length (rounded down to the nearest integer).
dispersion_size = 10; % pixels

% Magnification factor at which to save the results
magnification = 4;

% (row, column) pixel indices in the image at which to evaluate dispersion
location = [1526, 1963];

% The size of the whole image (height, width)
image_sampling = [2048, 2448];

% Sample spectral reflectances
reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180626_SpectralCharacterizationOfSetup/spectra_averaged.csv';

% Which samples to use from the file of reflectances
reflectance_columns = 13:36; % All samples, used for colour calibration
% Indices within the samples of reflectances to use for the graphical output
selected_reflectance_indices = [10, 13, 14, 15, 16, 19];

% CIE D-illuminant
illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
illuminant_temperature = 6504; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D

% The index of the white patch on the ColorChecker
white_index = 19;

% Whitepoint to use for XYZ to sRGB conversion
whitepoint = [1, 1, 1];

% CIE tristimulus functions
xyzbar_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180614_ASTM_E308/Table1_CIE1931_2DegStandardObserver.csv';

% ## Output parameters

% Relative side length of the patch to replace with an inset of the sample colour
inset_fraction = 0.2;

% Class of the output images
class_images = 'uint8';

% Refer to the MATLAB documentation on "Figure Properties"
% Output figure paper size, in inches
output_size_page = [15, 15];
% Output figure width, in inches
output_width = 9.5;
% Horizontal and vertical offsets from the lower left corner of the page, in
% inches
output_margin = [0.25, 0.25];

% Output directory
output_directory = '/home/llanos/Downloads';

% Parameters which do not usually need to be changed
run('SetFixedParameters.m')

%% Load calibration data

[...
    dispersion_data, ~, transform_data...
] = loadDispersionModel(spectral_model_filename, false, false);
[~, T_roi] = modelSpaceTransform(...
    image_sampling, transform_data.model_space, transform_data.fill, false...
);
df_spectral = makeDispersionfun(dispersion_data, T_roi);

[...
    dispersion_data, ~, transform_data...
] = loadDispersionModel(color_model_filename, false, false);
[~, T_roi] = modelSpaceTransform(...
    image_sampling, transform_data.model_space, transform_data.fill, false...
);
df_color = makeDispersionfun(dispersion_data, T_roi);

model_variables_required = { 'sensor_map', 'channel_mode', 'bands' };
load(color_map_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end
if channel_mode
    error('The input space of the colour conversion data must be a spectral space, not a space of colour channels.')
end
bands_color = bands;

%% Load reflectance data

sample_table = readtable(reflectances_filename);
variable_names = sample_table.Properties.VariableNames;
bands_measured = sample_table.(variable_names{1});
reflectances = sample_table{:, reflectance_columns};

n_reflectances = length(selected_reflectance_indices);
disp('Names of sample reflectances:')
for pc = 1:n_reflectances
    disp(variable_names{reflectance_columns(selected_reflectance_indices(pc))});
end
fprintf('\n');

%% Load the illuminant and convert reflectances to radiances

illuminant_data = csvread(illuminant_filename);
bands_illuminant = illuminant_data(:, 1);
S_illuminant = illuminant_data(:, 2:end);
spd_illuminant = ciedIlluminant(...
    illuminant_temperature, bands_illuminant, S_illuminant, bands_illuminant...
);

[bands_radiance, radiances] = reflectanceToRadiance(...
    bands_illuminant, spd_illuminant,...
    bands_measured, reflectances...
);
n_bands_radiance = length(bands_radiance);

%% Colour calibration

% Find theoretical sample colours
xyzbar_table = readtable(xyzbar_filename);
lambda_xyzbar = xyzbar_table{:, 1};
xyzbar = xyzbar_table{:, 2:end};

[~, reflectances_xyz] = reflectanceToColor(...
    bands_measured, ones(size(bands_measured)),... % Equal energy radiator
    bands_measured, reflectances,...
    lambda_xyzbar, xyzbar...
);

% Simulate actual colours
color_weights_options = findSamplingOptions;
color_weights_options.interpolant = findSamplingOptions.interpolant_ref;
color_weights = colorWeights(...
    sensor_map, bands_color, bands_radiance, color_weights_options...
);

rgb = channelConversion(radiances.', color_weights, 2);

% Colour homography
M_homog = ransachomocal_luv(...
    rgb,...
    reflectances_xyz,...
    reflectances_xyz(white_index, :),...
    rgb...
);
xyz_est = uea_homocvt(rgb, M_homog);
% Normalize by the white patch's green intensity
M_homog = M_homog ./ xyz_est(white_index, 2);

%% Preprocessing before image patch montage creation

output_postfix = '.pdf';
patch_size = repelem(2 * dispersion_size + 1, 1, 2);
font_size = max(12, floor(0.05 * max(patch_size) * magnification));
text_offset = [font_size * 3, font_size];

inset_size = floor(inset_fraction * patch_size);
roi_inset = [
    patch_size(2) - inset_size(2) + 1, patch_size(2),...
    1, inset_size(1)
];
inset_half_width = floor((inset_size - 1) / 2);

n_rows = 2; % Two models of dispersion
n_cols = n_reflectances;
n_boxes = n_rows * n_cols;
n_channels_rgb = 3;
bands_rgb = (1:n_channels_rgb).';

%% Image patch montage creation

I_montage = zeros([n_rows * patch_size(2), n_cols * patch_size(1), n_channels_rgb]);
I_montage_normalized = zeros(size(I_montage));

for col = 1:n_cols
    patch_spectral = zeros([patch_size, n_bands_radiance]);
    patch_spectral(dispersion_size + 1, dispersion_size + 1, :) = reshape(...
        radiances(:, selected_reflectance_indices(col)), 1, 1, []...
    );
    patch_color = channelConversion(patch_spectral, color_weights);
    
    inset = repmat(patch_color(dispersion_size + 1, dispersion_size + 1, :), [inset_size 1]);
    
    for row = 1:n_rows
        roi_montage = [
            ((row - 1) * patch_size(2)) + 1, row * patch_size(2),...
            ((col - 1) * patch_size(1)) + 1, col * patch_size(1)...
        ];
    
        if row == 1
            options = struct('bands_in', bands_rgb);
            df = df_color;
            patch_to_warp = patch_color;
        elseif row == 2
            options = dispersionfunToMatrixOptions;
            options.bands_in = bands_radiance;
            options.bands_out = bands_color;
            options.color_map = sensor_map;
            options.interpolant = options.interpolant_ref;
            df = df_spectral;
            patch_to_warp = patch_spectral;
        else
            error('There are only two models of dispersion.');
        end
        
        [ patch_warped, bands_dispersion ] = dispersionfunToMatrix(...
            df, options, patch_to_warp, true, flip(location) - 1 ...
        );
        if row == 2 && col == 1
            disp('Wavelengths at which dispersion is being sampled:');
            disp(bands_dispersion)
        end
    
        for c = 1:n_channels_rgb
            I_montage_normalized(roi_montage(1):roi_montage(2), roi_montage(3):roi_montage(4), c) = ...
                patch_warped(:, :, c) ./ max(max(patch_warped(:, :, c)));
        end

        patch_warped(roi_inset(1):roi_inset(2), roi_inset(3):roi_inset(4), :) = inset;

        I_montage(roi_montage(1):roi_montage(2), roi_montage(3):roi_montage(4), :) = patch_warped;
    end
end

I_montage = channelConversion(I_montage, M_homog);
I_montage = xyz2rgb(...
    I_montage, 'ColorSpace', 'srgb',...
    'WhitePoint', whitepoint,...
    'OutputType', class_images...
);

I_montage = imresize(I_montage, magnification, 'nearest');
I_montage_normalized = imresize(I_montage_normalized, magnification, 'nearest');

%% Show and save results

output_filename_prefix = fullfile(...
    output_directory,...
    sprintf(...
        'locationX%dY%d_dispersionResolution%g',...
        location(2), location(1), dispersionfunToMatrixOptions.resolution...
    )...
);

% Prepare annotations
b = 1;
labels = cell(n_boxes, 1);
reflectance_str = cell(n_boxes, 1);
type_str = cell(n_boxes, 1);
for col = 1:n_cols
    for row = 1:n_rows
        labels{b} = sprintf('(%d)', selected_reflectance_indices(col));
        reflectance_str{b} = variable_names{reflectance_columns(selected_reflectance_indices(col))};
        if row == 1
            type_str{b} = 'Color';
        elseif row == 2
            type_str{b} = 'Spectral';
        else
            error('There are only two models of dispersion.');
        end
        b = b + 1;
    end
end

% Create a figure for all channels, and two for each channel
for i = 1:(n_channels_rgb * 2 + 1)
    fg = figure;
    if i == 1
        output_filename_prefix_i = output_filename_prefix;
        imshow(I_montage);
    else
        c = ceil((i - 1) / 2);
        if mod(i, 2)
            imagesc(I_montage_normalized(:, :, c));
            output_filename_prefix_i = sprintf('%s_channel%d_normalized', output_filename_prefix, c);
        else
            imagesc(I_montage(:, :, c));
            output_filename_prefix_i = sprintf('%s_channel%d', output_filename_prefix, c);
        end
        colormap gray
        colorbar;
    end
    axis off
    axis image

    set(fg, 'PaperUnits', 'inches');
    set(fg, 'PaperSize', output_size_page);
    output_height = output_width * size(I_montage, 1) / size(I_montage, 2);
    set(fg, 'PaperPosition', [output_margin(1) output_margin(2) output_width, output_height]);
    print(...
        fg, sprintf('%s%s', output_filename_prefix_i, output_postfix),...
        '-dpdf', '-r300'...
    );

    % Annotate and output
    b = 1;
    for col = 1:n_cols
        for row = 1:n_rows
            x = col * patch_size(1) * magnification - text_offset(1);
            y = row * patch_size(2) * magnification - text_offset(2);
            text(x, y, labels{b}, 'FontSize', font_size, 'Color', 'w');
            b = b + 1;
        end
    end
    % For some reason, the annotations are not printing in the correct colour,
    % but remain black.
    print(...
        fg, sprintf('%s_annotated%s', output_filename_prefix_i, output_postfix),...
        '-dpdf', '-r300'...
    );
    close(fg);
end

% Output a legend
labels_table = reshape(reshape(labels, n_rows, n_cols).', [], 1);
reflectance_str_table = reshape(reshape(reflectance_str, n_rows, n_cols).', [], 1);
type_str_table = reshape(reshape(type_str, n_rows, n_cols).', [], 1);
t = table(labels_table, reflectance_str_table, type_str_table, 'VariableNames', {'Label_OverColumnsThenRows', 'Reflectance', 'ModelOfdispersion'});
writetable(...
    t, sprintf('%s_legend.csv', output_filename_prefix)...
);