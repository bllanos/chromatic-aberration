%% Visual evaluation of image estimation along lines
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Reference spectral images
%
% One or more images to be used to plot comparison lines in figures of spectral
% intensity along a line segment in the image plane. The images must have the
% same spectral sampling.
%
% Reference spectral images can be omitted if there are no spectral images to
% evaluate.
%
% ### Estimated spectral images
%
% Spectral images to be evaluated, all having a common spectral sampling,
% although not necessarily the same sampling as the reference images.
%
% ### Spectral sampling spaces
% - A '.mat' file containing the wavelengths at which the reference spectral
%   images are defined.
% - A '.mat' file containing the wavelengths at which the estimated spectral
%   images are defined, as well as the spectral resampling parameters governing
%   their conversion to other spectral sampling spaces.
% - A '.mat' file containing the wavelengths of the spectral bands which are to
%   be plotted in the output figures. Usually, this will be the same set of
%   wavelengths as the spectral bands of the reference spectral images, if any
%   are provided. Otherwise, there will be some error in the plotlines for the
%   reference images arising from interpolating their spectral information.
%
% ### Colour space conversion data
% A '.mat' file containing several variables, which is the output of
% 'SonyColorMap.m', for example. The following variables are required:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel of the output sensor response images to the
%   j-th element of 'bands' (below).
% - 'channel_mode': A Boolean value indicating whether the input colour
%   space is a set of colour channels (true) or a set of spectral bands
%   (false). A value of `false` is required.
% - 'bands': A vector containing the wavelengths corresponding to the
%   second dimension of 'sensor_map'.
%
% ### Reference colour image
%
% One or more images to be used to plot comparison lines in figures of colour
% channel intensities along a line segment in the image plane.
%
% ### Nice (i.e. colour-corrected) colour image
%
% A colour image used to show where the line segment evaluated is in the image.
%
% ### Estimated colour images
%
% Colour images to be evaluated.
%
% ## Output
%
% In the output figure legends, and the output files, the names of the
% algorithms being evaluated will be the unique portions of the estimated image
% filenames.
%
% ### Intermediate data and parameters
% A '.mat' file containing the following variables, as appropriate:
% - 'bands_reference': A vector containing the wavelengths of the spectral
%   bands in the reference spectral image.
% - 'bands_estimated': A vector containing the wavelengths of the spectral
%   bands used in the estimated images.
% - 'bands_plot': A vector containing the wavelengths of the spectral
%   bands used in the figures.
% - 'spectral_weights': A matrix for converting pixels in the spectral
%   space of the estimated spectral images to the spectral space of
%   'bands_plot'.
% - 'spectral_weights_reference': A matrix for converting pixels in the spectral
%   space of the reference spectral image to the spectral space of
%   'bands_plot'.
% - 'algorithms': A structure containing the names of the algorithms which were
%   evaluated, extracted from the names of the input estimated images.
%   'algorithms' has the following fields:
%   - 'spectral': A cell vector of character vectors containing the names of the
%     algorithms which produced the estimated spectral images.
%   - 'color': A cell vector of character vectors containing the names of the
%     algorithms which produced the estimated colour images.
% - 'algorithm_filenames': A structure of the same form as 'algorithms'
%   containing the names and paths of the input estimated images.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`.
%
% The file is saved as 'ShowImageLineData.mat'.
%
% ### Graphical output
%
% For each possible pairing of a reference image and an estimated image, one
% figure is generated. For estimated spectral images, figures are generated for
% spectral image evaluation, and also for colour image evaluation. Figures are
% saved to files whose names contain the endpoints of the image line segment
% being analyzed, and the indices of the reference and estimated images.
%
% Additionally, a patch of the "nice" colour image (see above) is saved with the
% line segment annotated.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 23, 2019

% List of parameters to save with results
parameters_list = {
    'reference_spectral_variable_name',...
    'names_spectral',...
    'reference_color_variable_name',...
    'names_color',...
    'bands_filename_reference',...
    'bands_variable_reference',...
    'annotation_image_filename',...
    'annotation_image_variable_name',...
    'sampling_filename',...
    'bands_filename_plot',...
    'bands_variable_plot',...
    'spectral_variable_name',...
    'color_map_filename',...
    'color_variable_name',...
    'line_ends',...
    'line_resolution'...
};

%% Input data and parameters

% Wildcard for 'ls()' to find the reference spectral images. At least one image
% must be found, unless there are no spectral images to evaluate (in which case
% the wildcard can be replaced with `[]`). '.mat' or image files can be loaded.
reference_spectral_wildcard = '/home/llanos/Downloads/testing/reference/spectral/*.mat';
reference_spectral_variable_name = 'I_hyper'; % Used only when loading '.mat' files

% Names for the reference spectral images, used in figures
names_spectral = {'Captured', 'Unwarped'};

% Wildcard for 'ls()' to find the reference colour images. At least one image
% must be found. '.mat' or image files can be loaded.
reference_color_wildcard = '/home/llanos/Downloads/testing/reference/rgb/*.mat';
reference_color_variable_name = 'I_rgb'; % Used only when loading '.mat' files

% Names for the reference colour images, used in figures
names_color = {'Captured', 'RGB unwarped', 'Spectral unwarped'};

% Path and filename of a '.mat' file containing the wavelengths corresponding to
% the reference spectral images
bands_filename_reference = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/sensor.mat';
bands_variable_reference = 'bands'; % Variable name in the above file

% Path and filename of a colour-corrected image to be used to produce a figure
% annotated with the line segment. ('.mat' or image files can be loaded)
annotation_image_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/warping_filtered_images/bandpassFiltered_correctedRGB/d2_book_dHyper_correctedRGB.tif';
annotation_image_variable_name = []; % Used only when loading '.mat' files

% Path and filename of a '.mat' file containing the wavelengths corresponding to
% the estimated spectral images
sampling_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/run_on_dataset_dHyper_dispersion/RunOnDataset_20190421_ComputarLens_dHyper_dispersion.mat';

% Path and filename of a '.mat' file containing the wavelengths to be used for
% plotting
bands_filename_plot = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/sensor.mat';
bands_variable_plot = 'bands'; % Variable name in the above file

% Wildcard for 'ls()' to find the estimated spectral images to process (can be
% empty). '.mat' or image files can be loaded.
spectral_wildcard = '/home/llanos/Downloads/testing/estimated/spectral/*.mat';
spectral_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/channel_scaling/sensor.mat';

% Wildcard for 'ls()' to find the estimated colour images to process (can be
% empty). '.mat' or image files can be loaded.
color_wildcard = '/home/llanos/Downloads/testing/estimated/rgb/*.mat';
color_variable_name = 'I_rgb'; % Used only when loading '.mat' files

% The image line segment to analyze (start pixel column, start pixel row, end
% pixel column, and end pixel row)
line_ends = [380, 104, 480, 109];

% Pixel resolution at which to sample the line
line_resolution = 0.1;

% ## Printing parameters for the annotated line segment image

% Refer to the MATLAB documentation on "Figure Properties"
% Output figure paper size, in inches
output_size_page = [15, 15];
% Output figure widths, in inches
output_width = 9.5;
% Horizontal and vertical offsets from the lower left corner of the page, in
% inches
output_margin = [0.25, 0.25];

% Output directory
output_directory = '/home/llanos/Downloads/testing';

%% Load wavelengths and spectral to colour conversion information

has_spectral = ~isempty(spectral_wildcard);
has_color = ~isempty(color_wildcard);

bands_estimated = [];
bands_reference = [];
bands_plot = [];
spectral_weights = [];
spectral_weights_reference = [];
if has_spectral
    bands_plot = loadVariables(bands_filename_plot, bands_variable_plot);
    
    [findSamplingOptions, bands_estimated] = loadVariables(sampling_filename, {'findSamplingOptions', 'bands'});
    n_bands_estimated = length(bands_estimated);

    [sensor_map, ~, bands_color] = loadColorMap(color_map_filename, false);    
    color_weights = colorWeights(...
        sensor_map, bands_color, bands_estimated, findSamplingOptions...
    );

    spectral_weights = resamplingWeights(...
        bands_plot, bands_estimated, findSamplingOptions.interpolant, findSamplingOptions.bands_padding...
    );
    
    bands_reference = loadVariables(bands_filename_reference, bands_variable_reference);
    n_bands_reference = length(bands_reference);

    spectral_weights_reference = resamplingWeights(...
        bands_plot, bands_reference, findSamplingOptions.interpolant_ref, findSamplingOptions.bands_padding...
    );
end
n_bands_plot = length(bands_plot);

%% Show the location of the line segment in the image

I_annotation = loadImage(annotation_image_filename, annotation_image_variable_name);
image_sampling = [size(I_annotation, 1), size(I_annotation, 2)];
if any(line_ends([1, 3]) > image_sampling(2)) || any(line_ends([2, 4]) > image_sampling(1)) || any(line_ends < 1)
    error('The line segment crosses the image boundaries.');
end

fg = figure;
line_length = sqrt((line_ends(3) - line_ends(1)) ^ 2 + (line_ends(4) - line_ends(2)) ^ 2);
pad = line_length * 0.1;
roi = [
    max(1, floor(min(line_ends([2, 4])) - pad));
    min(image_sampling(1), ceil(max(line_ends([2, 4])) + pad));
    max(1, floor(min(line_ends([1, 3])) - pad));
    min(image_sampling(2), ceil(max(line_ends([1, 3])) + pad))
];

I_annotation_roi = I_annotation(roi(1):roi(2), roi(3):roi(4), :);
imshow(I_annotation_roi);
hold on
plot(...
    line_ends([1, 3]) - roi(3) + 0.5, line_ends([2, 4]) - roi(1) + 0.5,...
    'Color', 'r', 'LineWidth', 2,...
    'Marker', 'none', 'LineStyle', '-'...
);
hold off

output_postfix = '.pdf';
output_filename_prefix = fullfile(output_directory, sprintf(...
    'lineX%dY%dX%dY%d',...
    line_ends(1), line_ends(2), line_ends(3), line_ends(4)...
));

set(fg, 'PaperUnits', 'inches');
set(fg, 'PaperSize', output_size_page);
output_height = output_width * size(I_annotation_roi, 1) / size(I_annotation_roi, 2);
set(fg, 'PaperPosition', [output_margin(1) output_margin(2) output_width, output_height]);
print(...
    fg, sprintf('%s_annotated%s', output_filename_prefix, output_postfix),...
    '-dpdf', '-r300'...
);
close(fg);

%% Prepare for image comparison

% Images to be evaluated, and evaluation parameters
algorithm_filenames.spectral = {};
n_spectral_algorithms = 0;
if has_spectral
    algorithm_filenames.spectral = listFiles(spectral_wildcard);
    n_spectral_algorithms = length(algorithm_filenames.spectral);
    if n_spectral_algorithms == 0
        has_spectral = false;
        warning('No estimated spectral images found.');
    end
end

algorithm_filenames.color = {};
if has_color
    algorithm_filenames.color = listFiles(color_wildcard);
end
if has_color || has_spectral
    algorithm_names = trimCommon([algorithm_filenames.spectral; algorithm_filenames.color]);
else
    algorithm_names = {};
end
algorithms.spectral = algorithm_names(1:n_spectral_algorithms);
algorithms.color = algorithm_names((n_spectral_algorithms + 1):end);
algorithms_escaped.spectral = strrep(algorithms.spectral, '_', '\_');
algorithms_escaped.color = strrep(algorithms.color, '_', '\_');
n_color_algorithms = length(algorithms.color);
if n_color_algorithms == 0 && has_color
    has_color = false;
    warning('No estimated colour images found.');
end

% Reference images
n_spectral_references = 0;
if has_spectral
    reference_filenames.spectral = listFiles(reference_spectral_wildcard);
    n_spectral_references = length(reference_filenames.spectral);
    if n_spectral_references ~= length(names_spectral)
        error('The number of names for the spectral reference images does not match the number of reference images.');
    end
    if n_spectral_references == 0
        error('No reference spectral images found.');
    end
end

reference_filenames.color = listFiles(reference_color_wildcard);
n_color_references = length(reference_filenames.color);
if n_color_references ~= length(names_color)
    error('The number of names for the colour reference images does not match the number of reference images.');
end
if n_color_references == 0
    error('No reference colour images found.');
end

evaluation_plot_colors_spectral = jet(n_bands_plot);
n_channels_rgb = 3;
evaluation_plot_colors_color = eye(n_channels_rgb);
evaluation_plot_styles = {'--', ':', '-.'};

% Sampling positions along the line
n_line_samples = ceil(line_length / line_resolution);
line_coords = line_ends([1, 2]) + (...
    repmat(linspace(0, 1, n_line_samples).', 1, 2) .*...
    repmat(line_ends([3, 4]) - line_ends([1, 2]), n_line_samples, 1)...
) - 0.5; % Convert from pixel indices to pixel centre coordinates
% Plots will be with respect to the coordinate that changes the most
if abs(line_ends(3) - line_ends(1)) >= abs(line_ends(4) - line_ends(2))
    line_x_plot = line_coords(:, 1);
    plot_label_x = 'Image x-coordinate';
else
    line_x_plot = line_coords(:, 2);
    plot_label_x = 'Image y-coordinate';
end

% Enumerate the positions of all pixels
[X, Y] = meshgrid(1:image_sampling(2), 1:image_sampling(1));
X = X - 0.5; % Place coordinates at pixel centres
Y = Y - 0.5;

%% Spectral image comparison

for r = 1:n_spectral_references
    for i = 1:n_spectral_algorithms
        fg = figure;
        legend_str = cell(2 * n_bands_plot, 1);
        hold on
        I = loadImage(reference_filenames.spectral{r}, reference_spectral_variable_name);
        line_samples = zeros(n_line_samples, n_bands_reference);
        for c = 1:n_bands_reference
            line_samples(:, c) = interp2(X, Y, I(:, :, c), line_coords(:, 1), line_coords(:, 2));
        end
        % Spectral sampling conversion
        line_samples_plot = channelConversion(line_samples, spectral_weights_reference, 2);

        for c = 1:n_bands_plot
            plot(...
                line_x_plot, line_samples_plot(:, c),...
                'Color', evaluation_plot_colors_spectral(c, :), 'LineWidth', 2,...
                'Marker', 'none', 'LineStyle', '-'...
            );
            legend_str{c} = sprintf(...
                '%g nm, %s', bands_plot(c), names_spectral{r}...
            );
        end

        I = loadImage(algorithm_filenames.spectral{i}, spectral_variable_name);
        line_samples = zeros(n_line_samples, n_bands_estimated);
        for c = 1:n_bands_estimated
            line_samples(:, c) = interp2(X, Y, I(:, :, c), line_coords(:, 1), line_coords(:, 2));
        end
        % Spectral sampling conversion
        line_samples_plot = channelConversion(line_samples, spectral_weights, 2);

        for c = 1:n_bands_plot
            plot(...
                line_x_plot, line_samples_plot(:, c),...
                'Color', evaluation_plot_colors_spectral(c, :), 'LineWidth', 2,...
                'Marker', 'none', 'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
            );
            legend_str{n_bands_plot + c} = sprintf(...
                '%g nm, %s', bands_plot(c), algorithms_escaped.spectral{i}...
            );
        end
        hold off

        legend(legend_str{:});
        xlabel(plot_label_x);
        ylabel('Radiance');
        ax = gca;
        ax.Color = [0.5 0.5 0.5];

        savefig(...
            fg,...
            sprintf('%s_spectralPlot_ref%d_test%d.fig', output_filename_prefix, r, i), 'compact'...
        );
        close(fg);
    end
end

%% Colour image comparison

for r = 1:n_color_references
    for i_all = 1:(n_spectral_algorithms + n_color_algorithms)
        fg = figure;
        legend_str = cell(2 * n_channels_rgb, 1);
        hold on
        I = loadImage(reference_filenames.color{r}, reference_color_variable_name);
        line_samples_plot = zeros(n_line_samples, n_channels_rgb);
        for c = 1:n_channels_rgb
            line_samples_plot(:, c) = interp2(X, Y, I(:, :, c), line_coords(:, 1), line_coords(:, 2));
        end

        for c = 1:n_channels_rgb
            plot(...
                line_x_plot, line_samples_plot(:, c),...
                'Color', evaluation_plot_colors_color(c, :), 'LineWidth', 2,...
                'Marker', 'none', 'LineStyle', '-'...
            );
            legend_str{c} = sprintf(...
                'Channel %d, %s', c, names_color{r}...
            );
        end

        if i_all <= n_spectral_algorithms
            i = i_all;
            I = loadImage(algorithm_filenames.spectral{i}, spectral_variable_name);
            line_samples = zeros(n_line_samples, n_bands_estimated);
            for c = 1:n_bands_estimated
                line_samples(:, c) = interp2(X, Y, I(:, :, c), line_coords(:, 1), line_coords(:, 2));
            end
            % Convert to colour
            line_samples_plot = channelConversion(line_samples, color_weights, 2);

            for c = 1:n_channels_rgb
                plot(...
                    line_x_plot, line_samples_plot(:, c),...
                    'Color', evaluation_plot_colors_color(c, :), 'LineWidth', 2,...
                    'Marker', 'none', 'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
                );
                legend_str{n_channels_rgb + c} = sprintf(...
                    'Channel %d, %s', c, algorithms_escaped.spectral{i}...
                );
            end
            
            output_filename = sprintf('%s_colorPlot_ref%d_testSpectral%d.fig', output_filename_prefix, r, i);
        else
            i = i_all - n_spectral_algorithms;
            I = loadImage(algorithm_filenames.color{i}, color_variable_name);
            line_samples_plot = zeros(n_line_samples, n_channels_rgb);
            for c = 1:n_channels_rgb
                line_samples_plot(:, c) = interp2(X, Y, I(:, :, c), line_coords(:, 1), line_coords(:, 2));
            end

            for c = 1:n_channels_rgb
                plot(...
                    line_x_plot, line_samples_plot(:, c),...
                    'Color', evaluation_plot_colors_color(c, :), 'LineWidth', 2,...
                    'Marker', 'none', 'LineStyle', evaluation_plot_styles{mod(i - 1, length(evaluation_plot_styles)) + 1}...
                );
                legend_str{n_channels_rgb + c} = sprintf(...
                    'Channel %d, %s', c, algorithms_escaped.color{i}...
                );
            end
            
            output_filename = sprintf('%s_colorPlot_ref%d_testColor%d.fig', output_filename_prefix, r, i);
        end
        hold off

        legend(legend_str{:});
        xlabel(plot_label_x);
        ylabel('Intensity');
        ax = gca;
        ax.Color = [0.5 0.5 0.5];

        savefig(fg, output_filename, 'compact');
        close(fg);
    end
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {
    'bands_reference', 'bands_estimated', 'bands_plot', 'spectral_weights',...
    'spectral_weights_reference', 'algorithms', 'algorithm_filenames'...
} ];
save_data_filename = fullfile(output_directory, 'ShowImageLineData.mat');
save(save_data_filename, save_variables_list{:});