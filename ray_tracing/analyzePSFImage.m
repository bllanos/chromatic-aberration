function [ stats ] = analyzePSFImage( psf_image, image_bounds, mask, varargin )
% ANALYZEPSFIMAGE  Describe a point spread function represented by an image
%
% ## Syntax
% stats = analyzePSFImage( psf_image, image_bounds, mask [, verbose] )
%
% ## Description
% stats = analyzePSFImage( psf_image, image_bounds, mask [, verbose] )
%   Returns statistics describing the density function represented by the
%   grid of samples.
%
% ## Input Arguments
%
% psf_image -- Image samples
%   An image_height x image_width array containing samples of a density
%   function on a rectangular grid. Presently, the interpretation of the
%   image is unimportant, but this function was written in a context where
%   `psf_image` represented a point spread function.
%
% image_bounds -- Image domain
%   The rectangular domain, in world coordinates of the image `psf_image`.
%   `image_bounds` is a vector containing the following elements:
%   1 - The x-coordinate of the bottom left corner of the image
%   2 - The y-coordinate of the bottom left corner of the image
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
%   The pixel at pixel indices `(i, j)` has its center at coordinates
%   ```
%   [...
%     image_bounds(1) + (i - 0.5) * image_bounds(3) / size(psf_image, 2),...
%     image_bounds(2) + (size(psf_image, 1) - j + 0.5) * image_bounds(4) / size(psf_image, 1)...
%   ]
%   ```
%
% mask -- Region of interest
%   A binary image, with the same dimensions as `psf_image`, indicating
%   which pixels to use. `mask` can be the region in which the image has
%   non-negligible values, for example.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% stats -- Distribution statistics
%   A structure describing `psf_image`, using the following fields:
%   - mean_position: The centroid of `psf_image`, computed by weighting
%     each pixel position (in world coordinates) by the pixel's value.
%   - mean_value: The evaluation of `psf_image` at `mean_position`, using
%     bilinear interpolation of neighbouring pixels.
%   - max_position: A two-element row vector containing the x and
%     y-world coordinates of the peak pixel value in `psf_image`.
%
%     If there are multiple local maxima in `psf_image`, and their average
%     location has a value which is lower than more than one of them, then
%     `max_position` is an NaN 1 x 2 array.
%   - max_value: A scalar containing the peak value corresponding to
%     `max_position`. `max_value` is the bilinear interpolation of
%     `psf_image` at `max_position`, and is NaN if `max_position` is NaN.
%   - radius: The weighted mean of the distances of pixels in `psf_image`
%     from `mean_position`, where the weights are the values of the pixels.
%     `radius` is related to the second moments of `psf_image`, and is
%     expressed in world units.
%
% ## Notes
% - To preallocate a structure array, to be filled with subsequent calls to
%   'analyzePSFImage()', use 'preallocateStats()'.
% - The analysis is delegated to 'analyzePSF()'.
%
% See also preallocateStats, analyzePSF, densifyRaysImage

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 17, 2017

nargoutchk(1, 1);
narginchk(3, 4);

image_width_px = size(psf_image, 2);
image_height_px = size(psf_image, 1);

image_position_ind = find(mask);
n_points = length(image_position_ind);
[image_position_row, image_position_col] = ind2sub(...
    [image_height_px, image_width_px], image_position_ind...
    );

% Create an image used to determine adjacency between pixels
index_image = nan(image_height_px, image_width_px);
convex_hull_indices = convhull([image_position_row, image_position_col]);
convhull_mask = roipoly(...
    index_image, image_position_col(convex_hull_indices), image_position_row(convex_hull_indices)...
);
index_image(convhull_mask) = 0;
index_image(image_position_ind) = 1:n_points;

image_position = [...
    image_bounds(1) + ((image_position_col - 0.5) * image_bounds(3) / image_width_px),...
    image_bounds(2) + ((image_height_px - image_position_row + 0.5) * image_bounds(4) / image_height_px)...
  ];

psf_values = psf_image(image_position_ind);

% Use an 8-neighbour adjacency configuration
n_neighbors = 8;
v_adj_ind = [
    image_position_row - 1, image_position_col    ;
    image_position_row    , image_position_col - 1;
    image_position_row - 1, image_position_col - 1;
    image_position_row + 1, image_position_col - 1;
    image_position_row - 1, image_position_col + 1;
    image_position_row + 1, image_position_col    ;
    image_position_row    , image_position_col + 1;
    image_position_row + 1, image_position_col + 1
    ];
v_adj_ind(...
    v_adj_ind(:, 1) < 1 | v_adj_ind(:, 1) > image_height_px, 1 ...
    ) = nan;
v_adj_ind(...
    v_adj_ind(:, 2) < 1 | v_adj_ind(:, 2) > image_width_px, 2 ...
    ) = nan;
v_adj_ind = sub2ind([image_height_px, image_width_px], v_adj_ind(:, 1), v_adj_ind(:, 2));
v_adj_ind = reshape(v_adj_ind, [], n_neighbors);

v_adj_mapped = nan(size(v_adj_ind));
v_adj_ind_filter = isfinite(v_adj_ind);
v_adj_mapped(v_adj_ind_filter) = index_image(v_adj_ind(v_adj_ind_filter));

v_adj = cell(n_points, 1);
n_neighbors_max = min(n_points - 1, n_neighbors);
for i = 1:n_points
    v_adj_mapped_i = v_adj_mapped(i, :);
    if any(v_adj_mapped_i == 0)
        % Just find the 8 nearest-neighbours
        distances_sq_i = image_position - repmat(image_position(i, :), n_points, 1);
        distances_sq_i = dot(distances_sq_i, distances_sq_i, 2);
        [~, sorting_map] = sort(distances_sq_i);
        v_adj{i} = sorting_map(2:(n_neighbors_max + 1));
    else
        v_adj{i} = v_adj_mapped_i(isfinite(v_adj_mapped_i));
    end
end

stats = analyzePSF( psf_values, image_position, v_adj, varargin{:} );

end

