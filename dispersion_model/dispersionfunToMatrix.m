function [ W, image_bounds_out ] = dispersionfunToMatrix(dispersionfun, lambda, image_sampling_in, varargin)
% DISPERSIONFUNTOMATRIX  Convert a warp model to a warp matrix, or use it
% directly to warp an image
%
% ## Syntax
% W = dispersionfunToMatrix(...
%   dispersionfun, lambda, image_sampling_in [, image_sampling_out,...
%   image_bounds, negate, offset]...
% )
% I_in = dispersionfunToMatrix(...
%   dispersionfun, lambda, image_sampling_in, I_out [,...
%   image_bounds, negate, offset]...
% )
% [ ____, image_bounds ] = dispersionfunToMatrix(____)
%
% ## Description
% W = dispersionfunToMatrix(...
%   dispersionfun, lambda, image_sampling_in [, image_sampling_out,...
%   image_bounds, negate, offset]...
% )
%   Returns a matrix for distorting a vectorized image.
%
% I_in = dispersionfunToMatrix(...
%   dispersionfun, lambda, image_sampling_in, I_out [,...
%   image_bounds, negate, offset]...
% )
%   Returns the distorted version of the input image.
%
% [ ____, image_bounds ] = dispersionfunToMatrix(____)
%   Additionally returns the boundaries of the undistorted image in the
%   coordinate system of the distorted image.
%
% ## Input Arguments
%
% dispersionfun -- Model of image warp
%   A function handle, produced by 'makeDispersionfun()'.
%   `dispersionfun(X_1)`, where `X_1` is a three-element row vector (x_1,
%   y_1, lambda), returns the displacement vector from the distorted
%   position `X_1` to its undistorted position (x_2, y_2, lambda).
%
% lambda -- Wavelength bands
%   A vector of length 'n' containing the wavelengths or colour channel
%   indices at which to evaluate the warp model
%
% image_sampling_in -- Distorted image dimensions
%   A two-element vector containing the image height, and width,
%   respectively, of the distorted image (in units of pixels)
%
% I_out -- Undistorted image
%   A 2D or 3D array containing an image to be warped. If `I_out` is a 2D array,
%   it must have a size other than 2 x 1 or 1 x 2 in order for it to be
%   distinguishable from `image_sampling_out`.
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
%   `image_bounds` describes the position and size of the undistorted image
%   in the space of the distorted image. Units are pixels in the distorted
%   image. Note that the last two elements of `image_bounds` do not need to
%   be equal to the elements of `image_sampling_out`, as the two images
%   need not have equal pixel sizes. If the two images are constrained to
%   overlap exactly, then `image_bounds` should be set as follows:
%     `image_bounds = [0, 0, image_sampling_in(2), image_sampling_in(1)]`
%
%   If `image_bounds` is empty or not passed, it will be set to the
%   bounding box of the undistorted coordinates of the distorted image.
%   More precisely, in this case, all pixels in the distorted image will be
%   generated, as implemented in `W`, by bilinear interpolation of four
%   different pixels in the undistorted image.
%
%   Otherwise, if `image_bounds` is not empty, and describes a space which
%   is smaller than the bounding box of the undistorted coordinates of the
%   distorted image, there may be pixels in the distorted image whose
%   undistorted coordinates are on or outside the border of the undistorted
%   image. If so, the border pixels of the undistorted image will be
%   replicated, such that the four pixels used for bilinear interpolation
%   will not be all unique.
%
% negate -- Warp inversion flag
%   If `true`, the warp vectors calculated by `dispersionfun` will be
%   negated before use. (Defaults to `false` if not passed.)
%
% offset -- Coordinate offset
%   A two element vector whose elements are offsets added to the pixel x-
%   and y-coordinates, respectively, prior to calling `dispersionfun` to
%   evaluate the undistorted position of the pixel. `offset` is useful for
%   calculating the warp matrix for a sub-image, without needing to modify
%   `dispersionfun` to use the local coordinates of the sub-image. Defaults
%   to `[0, 0]` if not passed.
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
% I_in -- Output image
%   An image_sampling_in(1) x image_sampling_in(2) x length(lambda) array
%   containing the version of `I_out` distorted according to `dispersionfun`.
%   The call syntax where `I_in` is returned instead of `W` is more memory
%   efficient, because `W` consumes a lot of memory. However, calling this
%   function for different values of `I_out`, to apply the same warp to
%   different images, is equivalent to recomputing `W`, and is therefore slower
%   than computing `W` and then applying it to the images.
%
% image_bounds -- Undistorted image coordinate frame
%   A copy of the `image_bounds` input argument, or its calculated version,
%   if the `image_bounds` input argument was empty, or was not passed.
%   Refer to the documentation of the `image_bounds` input argument above
%   for details.
%
% ## Algorithm
%
% `W` is created by first finding the undistorted position of each pixel in
% the distorted image, using `dispersionfun`. Second, the function determines
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
% See also xylambdaPolyfit, xylambdaSplinefit, makeDispersionfun, warpImage

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 6, 2018

nargoutchk(1, 2);
narginchk(3, 7);

image_passed = false;

n_lambda = length(lambda);
if ~isempty(varargin)
    image_sampling_out = varargin{1};
    image_passed = ndims(image_sampling_out) == 3 ||...
    ~(all(size(image_sampling_out) == [2 1]) ||...
      all(size(image_sampling_out) == [1 2]));
    if image_passed
        I_out = image_sampling_out;
        image_sampling_out = [size(I_out, 1), size(I_out, 2)];
        if size(I_out, 3) ~= n_lambda
            error('`I_out` must have a size of `length(lambda)` in its third dimension.');
        end
    end
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
if length(varargin) > 3
    offset = varargin{4};
else
    offset = [0, 0];
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
if size(lambda, 2) > size(lambda, 1)
    lambda = lambda.';
end
x_all_bands = repmat(x, n_lambda, 1);
y_all_bands = repmat(y, n_lambda, 1);
lambda_all_bands = reshape(repelem(lambda, n_px_in), [], 1);
disparity = dispersionfun([...
    x_all_bands + offset(1),...
    y_all_bands + offset(2),...
    lambda_all_bands...
    ]);
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
offsets = (0.5 + eps) * [
     1 -1; % Top right (Q21)
    -1 -1; % Top left (Q11)
    -1  1; % Bottom left (Q12)
     1  1; % Bottom right (Q22)
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

% Convert back to pixel coordinates
neighbour_x = reshape(neighbour_index_x - 0.5, [], n_offsets);
neighbour_y = reshape(neighbour_index_y - 0.5, [], n_offsets);

% Replicate pixels outside the image boundaries
neighbour_index_x(neighbour_index_x < 1) = 1;
neighbour_index_x(neighbour_index_x > image_sampling_out(2)) = image_sampling_out(2);
neighbour_index_y(neighbour_index_y < 1) = 1;
neighbour_index_y(neighbour_index_y > image_sampling_out(1)) = image_sampling_out(1);

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
    repmat(reshape(repelem((1:n_lambda).', n_px_in), [], 1), n_offsets, 1)...
);

if image_passed
    W = zeros(n_px_lambda_in, 1);
    for i = 1:n_offsets
        W = W + (...
            neighbour_weights(((i - 1) * n_px_lambda_in + 1):(i * n_px_lambda_in)) .* ...
            I_out(neighbour_index_linear(((i - 1) * n_px_lambda_in + 1):(i * n_px_lambda_in)))...
        );
    end
    W = reshape(W, [image_sampling_in, n_lambda]);
else
    % Assemble the sparse matrix
    indices_in = repmat((1:n_px_lambda_in).', n_offsets, 1);
    W = sparse(...
        indices_in,...
        neighbour_index_linear,...
        neighbour_weights,...
        n_px_lambda_in, prod(image_sampling_out) * n_lambda...
        );
end
end
