%% Visualization of relative spectral goodness-of-fit
%
% Plot changes in goodness-of-fit between hyperspectral images

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 17, 2018

%% Input data and parameters

% Wildcard for 'ls()' to find the true image.
% '.mat' or image files can be loaded
true_image_wildcard = 'dispersion0.000000e+00_latent.mat';
true_image_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the baseline estimated image.
% '.mat' or image files can be loaded
base_image_wildcard = 'noise0.000000e+00_dispersion0.000000e+00_patch24x24_pad16_mdc_latent.mat';
base_image_variable_name = 'I_latent'; % Used only when loading '.mat' files

% Wildcards for 'ls()' to find the comparison images. Use one wildcard per
% image.
% '.mat' or image files can be loaded
other_images_wildcard = {
    'noise0.000000e+00_dispersion1.000000e-01_patch24x24_pad16_mdc_latent.mat',...
    'noise0.000000e+00_dispersion3.000000e-01_patch24x24_pad16_mdc_latent.mat',...
    'noise0.000000e+00_dispersion1.000000e+00_patch24x24_pad16_mdc_latent.mat',...
    'noise0.000000e+00_dispersion2.000000e+00_patch24x24_pad16_mdc_latent.mat',...
    'noise0.000000e+00_dispersion3.000000e+00_patch24x24_pad16_mdc_latent.mat'...
    };
% The same variable name is used as for the base image.

% Names for the images
true_image_name = 'Ground truth';
base_image_name = '0';
other_images_names = { '0.1', '0.3', '1', '2', '3' };

% Number of horizontal lines along which to plot goodness-of-fit
n_lines_x = 5;

% Number of vertical lines along which to plot goodness-of-fit
n_lines_y = 5;

% Output directory for all figures
output_directory = '/home/llanos/Downloads';

%% Compute goodness-of-fit

n_other_images = length(other_images_wildcard);
if n_other_images ~= length(other_images_names)
    error('There must be as many other images as there are names for those images.');
end

true_image_filename = listFiles(true_image_wildcard);
I_gt = loadImage(true_image_filename{1}, true_image_variable_name);

n_all_images = n_other_images + 1;
image_width = size(I_gt, 2);
image_height = size(I_gt, 1);
gof = zeros(image_height, image_width, n_all_images);
I_gt_dot = dot(I_gt, I_gt, 3);
for i = 1:n_all_images
    if i == 1
        image_filename = listFiles(base_image_wildcard);
    else
        image_filename = listFiles(other_images_wildcard{i - 1});
    end
    I = loadImage(image_filename{1}, base_image_variable_name);
    gof(:, :, i) = abs(dot(I_gt, I, 3)) ./...
            sqrt(I_gt_dot .* dot(I, I, 3));
end

%% Visualize relative goodness-of-fit

for i = 1:n_other_images
    fg = figure;
    imagesc(gof(:, :, i + 1) - gof(:, :, 1));
    xlabel('X');
    ylabel('Y');
    c = colorbar;
    c.Label.String = 'Change in goodness-of-fit';
    title(sprintf(...
        'Goodness of fit of image %s relative to image %s',...
        other_images_names{i}, base_image_name...
    ))
    savefig(...
        fg,...
        fullfile(output_directory, sprintf('relativeGOF%d.fig', i)), 'compact'...
    );
    close(fg);
end

%% Visualize goodness-of-fit along scanlines

line_x = round(linspace(1, image_width, n_lines_x + 2));
line_x = line_x(2:(end - 1));
line_y = round(linspace(1, image_height, n_lines_y + 2));
line_y = line_y(2:(end - 1));

plot_colors = jet(n_all_images);
plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
plot_styles = {'--', ':', '-.'};
legend_str = [base_image_name, other_images_names];

line_coord = [line_x, line_y];
for ln = 1:(n_lines_x + n_lines_y)
    fg = figure;
    if ln > n_lines_x
        coord = 'y';
        gof_line = gof(line_coord(ln), :, :);
        plot_x = 1:image_width;
        other_coord = 'x';
    else
        coord = 'x';
        gof_line = gof(:, line_coord(ln), :);
        plot_x = 1:image_height;
        other_coord = 'y';
    end
    xlabel(sprintf('Image %s-coordinate', other_coord));
    ylabel('Spectral radiance goodness-of-fit');
    ylim([0, 1]);
    title(sprintf(...
        'Spectral radiance goodness-of-fit along line %s = %d',...
        coord, line_coord(ln)...
    ));
    hold on
    for i = 1:n_all_images
        plot(...
            plot_x, squeeze(gof_line(:, :, i)),...
            'Color', plot_colors(i, :), 'LineWidth', 2,...
            'Marker', plot_markers{1 + mod(i - 1, length(plot_markers))},...
            'LineStyle', plot_styles{1 + mod(i - 1, length(plot_styles))}...
        );
    end
    hold off
    legend(legend_str);

    savefig(...
        fg,...
        fullfile(output_directory, sprintf('lineGOF_%s%d.fig', coord, line_coord(ln))),...
        'compact'...
    );
    close(fg);
end