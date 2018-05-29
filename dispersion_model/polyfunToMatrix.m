function [ W, image_bounds_out ] = polyfunToMatrix(polyfun, lambda, image_sampling_in, varargin)
% POLYFUNTOMATRIX  Convert a polynomial warp model to a warp matrix
%
% ## Syntax
% W = polyfunToMatrix(...
%   polyfun, lambda, image_sampling_in [, image_sampling_out, image_bounds, negate]...
% )
% [ W, image_bounds ] = polyfunToMatrix(...
%   polyfun, lambda, image_sampling_in [, image_sampling_out, image_bounds, negate]...
% )
%
% ## Description
% W = polyfunToMatrix(...
%   polyfun, lambda, image_sampling_in [, image_sampling_out, image_bounds, negate]...
% )
%   Returns a matrix for distorting a vectorized image.
%
% [ W, image_bounds ] = polyfunToMatrix(...
%   polyfun, lambda, image_sampling_in [, image_sampling_out, image_bounds, negate]...
% )
%   Additionally returns the boundaries of the undistorted image in the
%   coordinate system of the distorted image.
%
% ## Input Arguments
%
% polyfun -- Polynomial model of image warp
%   A function handle, produced by 'makePolyfun()'. `polyfun(X_1)`, where
%   `X_1` is a three-element row vector (x_1, y_1, lambda), returns the
%   displacement vector from the distorted position `X_1` to its
%   undistorted position (x_2, y_2, lambda).
%
% lambda -- Wavelength bands
%   A vector of length 'n' containing the wavelengths or colour channel
%   indices at which to evaluate the polynomial warp model
%
% image_sampling_in -- Distorted image dimensions
%   A two-element vector containing the image height, and width,
%   respectively, of the distorted image (in units of pixels)
%
% image_sampling_out -- Undistorted image dimensions
%   A two-element vector containing the image height, and width,
%   respectively, of the undistorted image (in units of pixels). If
%   `image_sampling_out` is empty (`[]`) or not passed, it will be set
%   equal to `image_sampling_in`.
%
% image_bounds -- Undistorted image coordinate frame
%   A 4-element vector, containing the following elements:
%   1 - The pixel x coordinate of the top left corner of the image
%   2 - The pixel y coordinate of the top left corner of the image (using
%       the convention that the image y-axis increases downwards)
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
%   `image_bounds` describes the position and size of the undistorted
%   image in the space of the distorted image. Units are pixels in the
%   distorted image. Note that the last two elements of `image_bounds` do
%   not need to be equal to the elements of `image_sampling_out`, as the
%   two images need not have equal pixel sizes.
%
%   If `image_bounds` is empty or not passed, it will be set to the
%   bounding box of the undistorted coordinates of the distorted image.
%   More precisely, in this case, all pixels in the distorted image will be
%   generated, as implemented in `W`, by bilinear interpolation of four
%   different pixels in the undistorted image. Otherwise, there may be pixels
%   in the distorted image whose undistorted coordinates are on or outside the
%   border of the undistorted image. In this case, the border pixels of the
%   undistorted image will be replicated, such that the four pixels used for
%   bilinear interpolation will not be all unique.
%
% negate -- Warp inversion flag
%   If `true`, the warp vectors calculated by `polyfun` will be negated
%   before use. (Defaults to `false` if not passed.)
%
% ## Output Arguments
%
% W -- Warp matrix
%   A (n_px_in x n)-by-(n_px_out x n) sparse array (n = length(lambda))
%   which warps the undistorted image to the distorted image, according to
%   the equation:
%     `I_distorted = W * I_undistorted`
%   `I_undistorted` is a vectorized form of an image where all pixels have
%   been rearranged from columnwise order into a column vector.
%   `I_undistorted(i + n_px_out * (k - 1))` is the value of the k-th
%   wavelength band, or colour channel, at the i-th pixel of the
%   undistorted image. Similarly, `I_distorted(j + n_px_in * (k - 1))` is
%   the value of the k-th wavelength band or colour channel at the j-th
%   pixel of the distorted image.
%
% image_bounds -- Undistorted image coordinate frame
%   A copy of the `image_bounds` input argument, or its calculated version,
%   if the `image_bounds` input argument was empty or not passed. Refer to
%   the documentation of the `image_bounds` input argument above.
%
% ## Algorithm
%
% `W` is created by first finding the undistorted position of each pixel in
% the distorted image, using `polyfun`. Second, the function determines
% bilinear interpolation weights on the pixels in the undistorted image
% corresponding to the undistorted position, and enters these weights into
% `W`. Therefore, `W(j + n_px_in * (k - 1), i + n_px_out * (k - 1))` is the
% weight of the i-th pixel in the undistorted image in the bilinear
% interpolation calculation yielding the j-th pixel in the distorted image,
% for the k-th wavelength band.
%
% ## Notes
% - This function is presently not suitable for downsampling an image
%   during warping, as proper downsampling would involve blurring to avoid
%   aliasing.
%
% ## References
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.
%   - See the end of section 3 on computing the values of a warped image.
% - J. Brauers and T. Aach. "Geometric Calibration of Lens and Filter
%   Distortions for Multispectral Filter-Wheel Cameras," IEEE Transactions
%   on Image Processing, vol. 20, no. 2, pp. 496-505, 2011.
%   - Bilinear interpolation is used to compensate for chromatic
%     aberration (page 501).
% - Bilinear interpolation formulae:
%   https://en.wikipedia.org/wiki/Bilinear_interpolation
%
% See also xylambdaPolyfit, makePolyfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 6, 2018

nargoutchk(1, 2);
narginchk(3, 6);

if ~isempty(varargin)
    image_sampling_out = varargin{1};
else
    image_sampling_out = image_sampling_in;
end
set_bounds = (length(varargin) > 1) && ~isempty(varargin{2});
if set_bounds
    image_bounds_out = varargin{2};
end
if length(varargin) > 2
    negate = varargin{3};
else
    negate = false;
end

if length(image_sampling_in) ~= 2
    error('The `image_sampling_in` input argument must contain an image height and width only.');
end
if length(image_sampling_out) ~= 2
    error('The `image_sampling_out` input argument must contain an image height and width only.');
end

% Enumerate the positions of all pixels in the distorted image
[X, Y] = meshgrid(1:image_sampling_in(2), 1:image_sampling_in(1));
x = X(:) - 0.5; % Place coordinates at pixel centres
y = Y(:) - 0.5;
n_px_in = length(x);

% Find the undistorted positions of all pixels
n_lambda = length(lambda);
if size(lambda, 2) > size(lambda, 1)
    lambda = lambda.';
end
x_all_bands = repmat(x, n_lambda, 1);
y_all_bands = repmat(y, n_lambda, 1);
lambda_all_bands = repelem(lambda, n_px_in);
disparity = polyfun([ x_all_bands, y_all_bands, lambda_all_bands ]);
if negate
    disparity = -disparity;
end
x_out_all_bands = x_all_bands + disparity(:, 1);
y_out_all_bands = y_all_bands + disparity(:, 2);
n_px_lambda_in = length(x_out_all_bands);

% Find the undistorted image boundaries
if ~set_bounds
    max_x = max(x_out_all_bands);
    max_y = max(y_out_all_bands);
    min_x = min(x_out_all_bands);
    min_y = min(y_out_all_bands);

    image_width_minus1 = max_x - min_x;
    pixel_width = image_width_minus1 / (image_sampling_out(2) - 1);
    image_height_minus1 = max_y - min_y;
    pixel_height = image_height_minus1 / (image_sampling_out(1) - 1);
    image_bounds_out = [
        min_x - 0.5, min_y - 0.5,...
        pixel_width * image_sampling_out(2),...
        pixel_height * image_sampling_out(1)...
    ];
end
pixel_scale_x = image_sampling_out(2) / image_bounds_out(3);
pixel_scale_y = image_sampling_out(1) / image_bounds_out(4);

% Convert undistorted coordinates to pixel coordinates in the undistorted
% image
x_out_all_bands = pixel_scale_x * (x_out_all_bands - image_bounds_out(1));
y_out_all_bands = pixel_scale_y * (y_out_all_bands - image_bounds_out(2));

% Find the positions of neighbouring pixels in the undistorted image
offsets = [
     0.5 -0.5; % Top right (Q21)
    -0.5 -0.5; % Top left (Q11)
    -0.5  0.5; % Bottom left (Q12)
     0.5  0.5; % Bottom right (Q22)
     ];
n_offsets = size(offsets, 1);

neighbour_x = zeros(n_px_lambda_in, n_offsets);
neighbour_y = zeros(n_px_lambda_in, n_offsets);
for i = 1:n_offsets
    neighbour_x(:, i) = x_out_all_bands + offsets(i, 1);
    neighbour_y(:, i) = y_out_all_bands + offsets(i, 2);
end
neighbour_index_x = ceil(reshape(neighbour_x, [], 1));
neighbour_index_y = ceil(reshape(neighbour_y, [], 1));

% Replicate pixels outside the image boundaries
neighbour_index_x(neighbour_index_x < 1) = 1;
neighbour_index_x(neighbour_index_x > image_sampling_out(2)) = image_sampling_out(2);
neighbour_index_y(neighbour_index_y < 1) = 1;
neighbour_index_y(neighbour_index_y > image_sampling_out(1)) = image_sampling_out(1);

% Convert back to pixel coordinates
neighbour_x = reshape(neighbour_index_x - 0.5, [], n_offsets);
neighbour_y = reshape(neighbour_index_y - 0.5, [], n_offsets);

% Find bilinear interpolation weights
neighbour_weights = zeros(n_px_lambda_in, n_offsets);
x1 = neighbour_x(:, 2);
x2 = neighbour_x(:, 1);
y1 = neighbour_y(:, 2);
y2 = neighbour_y(:, 3);
dx = x2 - x1;
dy = y2 - y1;
tox2 = x2 - x_out_all_bands;
tox1 = x_out_all_bands - x1;
toy2 = y2 - y_out_all_bands;
toy1 = y_out_all_bands - y1;

neighbour_weights(:, 1) = tox1 .* toy2; % Q21
neighbour_weights(:, 2) = tox2 .* toy2; % Q11
neighbour_weights(:, 3) = tox2 .* toy1; % Q12
neighbour_weights(:, 4) = tox1 .* toy1; % Q22
neighbour_weights = neighbour_weights ./ (dx .* dy);

neighbour_weights = reshape(neighbour_weights, [], 1);

% Find linear indices for pixels
neighbour_index_linear = sub2ind(...
    [image_sampling_out, n_lambda],...
    neighbour_index_y,...
    neighbour_index_x,...
    repmat(repelem((1:n_lambda).', n_px_in), n_offsets, 1)...
);

% Assemble the sparse matrix
indices_in = repmat((1:n_px_lambda_in).', n_offsets, 1);
W = sparse(...
    indices_in,...
    neighbour_index_linear,...
    neighbour_weights,...
    n_px_lambda_in, prod(image_sampling_out) * n_lambda...
    );
end
