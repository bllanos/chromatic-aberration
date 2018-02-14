function [ image_bounds ] = imageBoundaries( image_bounds, stats )
% IMAGEBOUNDARIES  Pad a grid of image positions to estimate image boundaries
%
% ## Syntax
% image_bounds = imageBoundaries( image_bounds, stats )
%
% ## Description
% image_bounds = imageBoundaries( image_bounds, stats )
%   Find image boundaries that extend beyond the given image positions by
%   the approximate spacing between image positions.
%
% ## Input Arguments
%
% image_bounds -- Image domain
%   The rectangular domain of the image. `image_bounds` is a vector
%   containing the following elements:
%   1 - The x-coordinate of the bottom left corner of the image
%   2 - The y-coordinate of the bottom left corner of the image
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
%   If `image_bounds` is empty, and `stats` describes more than one light
%   source, this function will automatically estimate the image boundaries.
%
% stats -- Image point statistics
%   Statistics computed for the images produced by a grid of light sources,
%   where the first dimension of `stats` is assumed to index the light
%   sources. Only `stats.mean_position` is used by this function. The
%   'mean_position' field is assumed to contain a two element vector
%   consisting of an image x and y-coordinate, respectively. The image
%   positions are assumed to approximately form an axis-aligned grid - This
%   function assumes the grid is relatively undistorted, allowing for
%   automatic estimation of the numbers of positions along each axis.
%
% ## Output Arguments
%
% image_bounds -- Image domain
%   If the input argument `image_bounds` is non-empty, the output argument
%   is a copy of the input argument. Otherwise, `image_bounds` has the same
%   form as described in the documentation of the input argument, but is
%   calculated from `stats.mean_position`. Specifically, the image
%   boundaries are set so that the image is at least large enough to
%   contain two extra rows, and two extra columns, of image positions in
%   the grid.
%
% ## Notes
% - This function is primarily intended as a helper function of
%   'doubleSphericalLensPSF()'.
%
% See also doubleSphericalLensPSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 21, 2017

nargoutchk(1,1);
narginchk(2,2);

n_lights = size(stats, 1);
if isempty(image_bounds) && n_lights > 1
    X_image_ideal = [stats.mean_position];
    X_image_ideal = reshape(X_image_ideal, 2, []);
    X_image_ideal = X_image_ideal.';
    min_x = min(X_image_ideal(:, 1));
    min_y = min(X_image_ideal(:, 2));
    max_x = max(X_image_ideal(:, 1));
    max_y = max(X_image_ideal(:, 2));
    if (min_x == max_x) || (min_y == max_y)
        error('Lights have zero spread in the x and/or y directions. A grid cannot be inferred.')
    end

    % Find a best-fit grid spacing
    fit = Inf(n_lights, 2);
    n_points = size(X_image_ideal, 1);
    X_image_ideal_rel = X_image_ideal - repmat(mean(X_image_ideal, 1), n_points, 1);
    for i = 2:n_lights
        d_x = (max_x - min_x) / (i - 1);
        d_y = (max_y - min_y) / (i - 1);
        if mod(i, 2) == 0
            X_image_ideal_rel_i = X_image_ideal_rel - repmat([d_x, d_y] / 2, n_points, 1);
        else
            X_image_ideal_rel_i = X_image_ideal_rel;
        end
        X_image_ideal_residual = X_image_ideal_rel_i ./ repmat([d_x, d_y], n_points, 1);
        X_image_ideal_residual = abs(X_image_ideal_residual - round(X_image_ideal_residual));
        fit(i, :) = sum(X_image_ideal_residual, 1);
    end
    [~, n_lights_x] = min(fit(:, 1));
    [~, n_lights_y] = min(fit(:, 2));

    left_buffer = (max_x - min_x) / (n_lights_x - 1);
    right_buffer = left_buffer;

    bottom_buffer = (max_y - min_y) / (n_lights_y - 1);
    top_buffer = bottom_buffer;

    image_bounds = [
        min_x - left_buffer,...
        min_y - bottom_buffer
    ];
    image_width = max_x - image_bounds(1) + right_buffer;
    image_height = max_y - image_bounds(2) + top_buffer;
    image_bounds = [
        image_bounds,...
        image_width image_height...
    ];
elseif isempty(image_bounds)
    error('Image boundaries cannot be automatically estimated for single light sources.')
end

end
