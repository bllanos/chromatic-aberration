%% Visualization of the differences in dispersion within colour channels
%
% Plot dispersion magnitudes along a line from the image centre for the
% wavelengths at which each colour channel's sensitivity rises above or
% drops below a given threshold.
%
% This script assumes that each colour channel only has a single interval
% in its domain where it rises above the given proportion of its maximum
% value over its entire domain.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 21, 2018

%% Input data and parameters

% Spectral model of dispersion
spectral_model_filename = fullfile('.', 'demo_data', 'dispersion_models', 'registration', 'RegistrationDispersionResults_spectral_polynomial_fromReference.mat');

% Camera spectral sensitivity data
color_map_filename = fullfile('.', 'demo_data', 'multispectral_images', 'sensor.mat');

% Threshold in relative sensitivty within a colour channel at which to mark
% the "edges" of the colour channel in the spectrum
channel_threshold = 0.5;

% Number of points to plot in a line joining the bottom right corner of the
% dispersion model's domain with the origin of the dispersion model's
% coordinate system
n_points = 500;

%% Processing

% Find the wavelengths bounding the high-sensitivity regions of the camera's
% colour channels
[sensor_map, ~, bands] = loadColorMap(color_map_filename, false);

cutoff_wavelength_sensitivites = max(sensor_map, [], 2);
sensor_map_relative = sensor_map ./ repmat(cutoff_wavelength_sensitivites, 1, size(sensor_map, 2));
cutoff_wavelength_sensitivites = cutoff_wavelength_sensitivites * channel_threshold;
n_channels = size(sensor_map, 1);
cutoff_wavelengths = zeros(n_channels, 2);
for c = 1:n_channels
    channel_thresholded = sensor_map_relative(c, :) > channel_threshold;
    index = max(find(channel_thresholded, 1), 2);
    cutoff_wavelengths(c, 1) = (cutoff_wavelength_sensitivites(c) - sensor_map(c, index - 1))...
        * (bands(index) - bands(index - 1))...
        / (sensor_map(c, index) - sensor_map(c, index - 1))...
        + bands(index - 1);
    index = min(find(channel_thresholded, 1, 'last') + 1, length(bands));
    cutoff_wavelengths(c, 2) = (cutoff_wavelength_sensitivites(c) - sensor_map(c, index - 1))...
        * (bands(index) - bands(index - 1))...
        / (sensor_map(c, index) - sensor_map(c, index - 1))...
        + bands(index - 1);
end

figure;
plot_colors = eye(3);
hold on
for c = 1:n_channels
    plot(bands, sensor_map(c, :), 'Color', plot_colors(c, :), 'LineWidth', 2);
    scatter(cutoff_wavelengths(c, :), repelem(cutoff_wavelength_sensitivites(c), 2), [], [0, 0, 0], 'filled');
end
hold off
xlabel('Wavelength [nm]')
ylabel('Relative quantum efficiency')

% Evaluate the spectral model of dispersion
load(spectral_model_filename);
dispersionfun = makeDispersionfun(dispersion_data);

if strcmp(model_space.system, 'geometric')
    origin = [0, 0];
elseif strcmp(model_space.system, 'image')
    origin = flip(model_space.image_size / 2);
end

x = linspace(origin(1), model_space.corners(3), n_points + 1).';
y = linspace(origin(2), model_space.corners(4), n_points + 1).';
x = x(2:end);
y = y(2:end);
lambda = repelem(reshape(cutoff_wavelengths.', [], 1), n_points);
x_rep = repmat(x, numel(cutoff_wavelengths), 1);
y_rep = repmat(y, numel(cutoff_wavelengths), 1);

dispersion = dispersionfun([x_rep, y_rep, lambda]);
% Project the dispersion onto the vector pointing away from the image center
from_center = [x_rep, y_rep] - repmat(origin, length(x_rep), 1);
from_center = from_center ./ repmat(sqrt(dot(from_center, from_center, 2)), 1, 2);
dispersion_mag = dot(dispersion, from_center, 2);
dispersion_mag = reshape(dispersion_mag, n_points, size(cutoff_wavelengths, 2), size(cutoff_wavelengths, 1));

% Visualize dispersion magnitudes
distance_to_origin = [x - origin(1), y - origin(2)];
distance_to_origin = sqrt(distance_to_origin(:, 1) .^ 2 + distance_to_origin(:, 2) .^ 2);

figure;
hold on
legend_str = cell(numel(cutoff_wavelengths), 1);
for c = 1:n_channels
    plot(distance_to_origin, dispersion_mag(:, 1, c), '-', 'Color', plot_colors(c, :), 'LineWidth', 2);
    legend_str{c * 2 - 1} = sprintf('\\lambda = %g', cutoff_wavelengths(c, 1));
    plot(distance_to_origin, dispersion_mag(:, 2, c), ':', 'Color', plot_colors(c, :), 'LineWidth', 2);
    legend_str{c * 2} = sprintf('\\lambda = %g', cutoff_wavelengths(c, 2));
end
hold off
legend(legend_str);
xlabel('Distance to image centre [mm]')
ylabel('Dispersion along vector to image centre [mm]')
title(sprintf('Dispersion at the %g percent sensitivity thresholds for each colour channel', channel_threshold * 100))