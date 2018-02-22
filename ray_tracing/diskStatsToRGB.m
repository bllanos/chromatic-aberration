function [ stats_rgb ] = diskStatsToRGB( stats_wavelengths, wavelengths_to_rgb, threshold )
% DISKSTATSTORGB  Create the RGB equivalent of monochromatic PSF statistics
%
% ## Syntax
% stats_rgb = diskStatsToRGB( stats_wavelengths, wavelengths_to_rgb, threshold )
%
% ## Description
% stats_rgb = diskStatsToRGB( stats_wavelengths, wavelengths_to_rgb, threshold )
%   Aggregate PSF statistics for individual wavelengths into PSF statistics
%   for colour channels.
%
% ## Input Arguments
%
% stats_wavelengths -- PSF statistics for individual wavelengths
%   A 3D structure array, where the dimensions represent the following:
%   - Light sources
%   - Wavelengths of light
%   - Depths of the light sources (The same light source is presented at
%     different depths)
%
%   Each element of `stats_wavelengths` is a structure of the form output
%   by 'analyzePSF()', describing the blur circle, predicted by geometric
%   optics. For instance, elements were generated using the `psfFn` output
%   argument of 'opticsToPSF()'.
%
% wavelengths_to_rgb -- RGB quantum efficiencies
%   RGB quantum efficiencies for the wavelengths corresponding to the
%   second dimension of `stats_wavelengths`. The i-th row of this k x 3
%   matrix represents the RGB sensitivities corresponding to the i-th
%   wavelength, the wavelength for which `stats_wavelengths(:, i, :)` was
%   produced.
%
% threshold -- Quantum efficiency threshold
%   When generating the `radius` field of the output argument, `stats_rgb`,
%   only point spread functions for a given light source at a given depth
%   with intensities greater than or equal to `threshold` times the
%   intensity of the point spread function with the largest intensity for
%   this light source and depth will be considered. `threshold` accounts
%   for the fact that quantum efficiencies are never zero for any
%   wavelength, but only a certain range of wavelengths will have an
%   appreciable contribution to the image.
%
% ## Output Arguments
%
% stats_rgb -- PSF statistics for color channels
%   A version of `stats_wavelengths`, where the second dimension now has
%   size 3, and represents RGB colour channels instead of single
%   wavelengths. The fields of each element of the structure array are
%   produced by combining the data for the individual wavelengths, as
%   described in the comments in the body of this function. As presently
%   implemented, however, not all fields in `stats_wavelengths` are
%   processed.
%
% ## Notes
% - This function is primarily intended as a helper function of
%   'doubleSphericalLensPSF2()'.
%
% See also opticsToPSF, doubleSphericalLensPSF2

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 20, 2018

nargoutchk(1,1);
narginchk(3,3);

n_points = size(stats_wavelengths, 1);
n_wavelengths = size(wavelengths_to_rgb, 1);
if n_wavelengths ~= size(stats_wavelengths, 2)
    error('Mismatch between the number of wavelengths in `stats_wavelengths`, and the number of rows in `wavelengths_to_rgb`.')
end
n_depths = size(stats_wavelengths, 3);

names = fieldnames(stats_wavelengths);
n_names = length(names);

% Preliminary processing
stats_cell = struct2cell(stats_wavelengths);
for i = 1:n_names
    stats_cell_i = squeeze(stats_cell(i, :, :, :));
    stats_cell_i = reshape(stats_cell_i, n_points, 1, n_wavelengths, n_depths);
    stats_mat_i = cell2mat(stats_cell_i);
    stats_mat.(names{i}) = stats_mat_i;
end

n_channels = size(wavelengths_to_rgb, 2);
wavelengths_to_rgb_weights = sum(wavelengths_to_rgb, 1);
wavelengths_to_rgb_weights = reshape(wavelengths_to_rgb_weights, 1, 1, 1, 1, n_channels);
wavelengths_to_rgb_reshape = reshape(wavelengths_to_rgb, 1, 1, n_wavelengths, 1, n_channels);
% Normalize
wavelengths_to_rgb_reshape = wavelengths_to_rgb_reshape ./ repmat(...
    wavelengths_to_rgb_weights, 1, 1, n_wavelengths, 1, 1 ...
    );

% Mean position
% Perform an average over wavelengths, weighted by quantum efficiencies
mean_position_rgb = repmat(stats_mat.mean_position, 1, 1, 1, 1, n_channels);
wavelengths_to_rgb_rep = repmat(wavelengths_to_rgb_reshape, n_points, size(mean_position_rgb, 2), 1, n_depths, 1);
mean_position_rgb = mean_position_rgb .* wavelengths_to_rgb_rep;
mean_position_rgb = sum(mean_position_rgb, 3);
mean_position_rgb = squeeze(permute(mean_position_rgb, [1, 2, 5, 4, 3]));
stats_rgb_mat.mean_position = mean_position_rgb;

% Mean value
%
% Sum the intensities (accounting for quantum efficiencies) from PSF disks
% overlapping the mean position. This is an approximation - In reality, the
% intensity inside each disk is not constant.
mean_value_rep = repmat(stats_mat.mean_value, 1, 1, 1, 1, n_channels);
wavelengths_to_rgb_rep = repmat(wavelengths_to_rgb_reshape, n_points, size(mean_value_rep, 2), 1, n_depths, 1);
mean_value_rep = mean_value_rep .* wavelengths_to_rgb_rep;
r_sq = stats_mat.radius .^ 2;
distance_filter = zeros(n_points, size(mean_value_rep, 2), n_wavelengths, n_depths, n_channels);
for c = 1:n_channels
    mean_position_c_rep = repmat(mean_position_rgb(:, :, c, :), 1, 1, n_wavelengths, 1);
    distance_sq = mean_position_c_rep - stats_mat.mean_position;
    distance_sq = dot(distance_sq, distance_sq, 2);
    distance_filter(:, :, :, :, c) = (distance_sq <= r_sq);
end
mean_value_rgb = sum(mean_value_rep .* distance_filter, 3);
mean_value_rgb = squeeze(permute(mean_value_rgb, [1, 2, 5, 4, 3]));
stats_rgb_mat.mean_value = mean_value_rgb;

% Max position
%
% This is a combinatorial problem: Consider all regions where the PSF disks
% overlap, and find the region with the highest resulting intensity. I
% don't think that solving this problem is warranted, especially
% considering that I would be approximating the intensity inside each disk
% as constant.

% Max value
% Same concerns as for max position.

% Radius
%
% Take the PSFs for the wavelengths that have intensities above `threshold`
% times the intensity of the most intense PSF, and find the largest
% distances between points on their edges.
radii_mat = zeros(n_points, size(stats_mat.radius, 2), 1, n_depths, n_channels);
max_mean_value_rep = repmat(max(mean_value_rep, [], 3), 1, 1, n_wavelengths, 1, 1);
separation_filter = (mean_value_rep >= threshold * max_mean_value_rep);
for w = 1:n_wavelengths
    position_w_rep = repmat(stats_mat.mean_position(:, :, w, :), 1, 1, n_wavelengths, 1);
    separation = position_w_rep - stats_mat.mean_position;
    separation_distance = dot(separation, separation, 2);
    separation_distance = separation_distance +...
        repmat(stats_mat.radius(:, :, w, :), 1, 1, n_wavelengths, 1) +...
        stats_mat.radius;
    separation_distance = repmat(separation_distance, 1, 1, 1, 1, n_channels);
    separation_distance = separation_distance .* separation_filter;
    separation_distance = max(separation_distance, [], 3);
    radii_mat = max(radii_mat, separation_distance);
end
radius_rgb = squeeze(permute(radii_mat, [1, 2, 5, 4, 3]));
stats_rgb_mat.radius = radius_rgb;

% Split back into a structure array
stats_rgb_cell = cell(n_names, n_points, n_channels, n_depths);
for i = 1:n_names
    stats_mat_i = stats_rgb_mat.(names{i});
    stats_cell_i = mat2cell(...
        stats_mat_i, ones(n_points, 1), size(stats_mat_i, 2), ones(n_channels, 1), ones(n_depths, 1)...
    );
    stats_cell_i = reshape(stats_cell_i, 1, n_points, n_channels, n_depths);
    stats_rgb_cell(i, :, :, :) = stats_cell_i;
end
stats_rgb = cell2struct(stats_rgb_cell, names, 1);

end
