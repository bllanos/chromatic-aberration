function [ W, lambda_dispersion ] = dispersionfunToMatrix(dispersionfun, spectral_options, image_sampling, varargin)
% DISPERSIONFUNTOMATRIX  Convert a warp model to a warp matrix, or use it
% directly to warp an image
%
% ## Syntax
% W = dispersionfunToMatrix(...
%   dispersionfun, spectral_options, image_sampling [, negate, offset]...
% )
% I_out = dispersionfunToMatrix(...
%   dispersionfun, spectral_options, I_in [, negate, offset]...
% )
% [ ____, lambda_dispersion ] = dispersionfunToMatrix(____)
%
% ## Description
% W = dispersionfunToMatrix(...
%   dispersionfun, spectral_options, image_sampling [, negate, offset]...
% )
%   Returns a matrix for distorting a vectorized image.
%
% I_out = dispersionfunToMatrix(...
%   dispersionfun, spectral_options, I_in [, negate, offset]...
% )
%   Returns the distorted version of the input image.
%
% [ ____, lambda_dispersion ] = dispersionfunToMatrix(____)
%   Additionally returns the spectral bands at which the dispersion model was
%   sampled.
%
% ## Input Arguments
%
% dispersionfun -- Model of image warp
%   A function handle, produced by 'makeDispersionfun()'. `dispersionfun(X_1)`,
%   where `X_1` is a three-element row vector (x_1, y_1, lambda), returns the
%   displacement vector from the distorted position `X_1` to its undistorted
%   position (x_2, y_2, lambda).
%
% spectral_options -- Spectral sampling options
%   `spectral_options` is a structure of options controlling the spectral
%   sampling of dispersion, and allowing for simultaneous warping, and spectral
%   resampling or conversion to colour.
%
%   The following fields are mandatory:
%   - 'bands_in': A vector of length 'n' containing the wavelengths or colour
%     channel indices at which the images to be warped are sampled. 'bands_in'
%     must have the same length as the size of `I_in` in its third dimension.
%     `bands_in` must contain evenly-spaced values if spectral resampling is to
%     occur (depending on the options below).
%
%   The following fields are optional:
%   - 'resolution': A non-negative scalar providing the desired approximate
%     spacing, in pixels, between the images for consecutive wavelengths at
%     which the dispersion is to be sampled. If 'resolution' is zero or is
%     missing, the dispersion function will be evaluated at the elements of
%     'bands_in'. Otherwise, it will be evaluated at a sequence of wavelengths
%     which would, if dispersion varied linearly with wavelength, have images
%     which are shifted relative to those for adjacent wavelengths in the
%     sequence by approximately 'resolution' pixels.
%   - 'bands_out': A vector of length 'm' containing the wavelengths at which
%     the warped images must be sampled. If 'bands_out' is missing, it will be
%     set equal to 'bands_in'. Either 'bands_out', or 'color_map' and
%     'color_bands' can be present, not both.
%   - 'color_map': A 2D array with 'm' rows, where `color_map(i, j)` is the
%     sensitivity of the i-th colour channel of the warped image to the j-th
%     spectral band in `color_bands`. In contrast with the `sensitivity`
%     argument of 'solvePatchesADMM()', `color_map` is not a colour conversion
%     matrix, as it does not perform the desired numerical integration, over the
%     spectrum, that is part of colour conversion. Either 'bands_out', or
%     'color_map' and 'color_bands' can be present, not both.
%   - 'color_bands': A vector, of length equal to the size of the second
%     dimension of 'color_map', containing the wavelengths at which the
%     sensitivity functions in `color_map` have been sampled. `color_bands(j)`
%     is the wavelength corresponding to `color_map(:, j)`. The values in
%     'color_bands' are expected to be evenly-spaced. Either 'bands_out', or
%     'color_map' and 'color_bands' can be present, not both.
%   - 'int_method': The numerical integration method to use when
%     integrating over the responses of colour channels to compute colour
%     values. `int_method` is passed to `integrationWeights()` as its `method`
%     input argument. 'int_method' is not used if 'color_map' and 'color_bands'
%     are not present.
%   - 'bands_padding': Padding to use when computing spectral interpolation
%     matrices. Refer to the documentation of the 'padding' input argument of
%     'resamplingWeights()' for details, in 'resamplingWeights.m'.
%   - 'interpolant': The function convolved with spectral signals to interpolate
%     them during resampling. 'interpolant' is passed to 'resamplingWeights()'
%     as its `f` input argument. Refer to the documentation of
%     'resamplingWeights.m' for more details.
%
%   The above optional fields should not be present if the dispersion matrix is
%   to operate on colour channels instead of spectral bands, because it is
%   nonsensical to resample colour channels. It is not correct to instead set
%   'resolution' to zero, and 'bands_out' to 'bands_in', because 'interpolant'
%   may still not result in an identity mapping.
%
% image_sampling -- Distorted image dimensions
%   A two-element vector containing the image height, and width, respectively,
%   of the distorted image (in units of pixels)
%
% I_in -- Undistorted image
%   A 2D or 3D array containing an image to be warped. If `I_in` is a 2D array,
%   it must have a size other than 2 x 1 or 1 x 2 in order for it to be
%   distinguishable from `image_sampling`.
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
%   A (n_px x m)-by-(n_px x n) sparse array.
%   - n is `length(spectral_options.bands_in)`, and
%   - m is `length(spectral_options.bands_out)` or
%     `size(spectral_options.color_map, 1)`, depending on which of
%     `spectral_options.bands_out` or `spectral_options.color_map` exists.
%
%   `W` warps the undistorted image to the distorted image, and may
%   simutaneously change the spectral sampling of the image or convert it to
%   colour, depending on the settings in `spectral_options`. `W` operates
%   according to the equation:
%     `I_distorted = W * I_undistorted`
%   `I_undistorted` is a vectorized form of an image where all pixels have been
%   rearranged from columnwise order into a column vector. `I_undistorted(i +
%   n_px * (k - 1))` is the value of the k-th wavelength band, or colour
%   channel, at the i-th pixel of the undistorted image. Similarly,
%   `I_distorted(j + n_px * (p - 1))` is the value of the p-th wavelength band
%   or colour channel at the j-th pixel of the distorted image.
%
% I_out -- Output image
%   An image_sampling(1) x image_sampling(2) x m array containing the
%   version of `I_in` distorted according to `dispersionfun`. The call syntax
%   where `I_out` is returned instead of `W` is more memory efficient, because
%   `W` consumes a lot of memory. However, calling this function for different
%   values of `I_in`, to apply the same warp to different images, is equivalent
%   to recomputing `W`, and is therefore slower than computing `W` and then
%   applying it to the images.
%
% lambda_dispersion -- Dispersion evaluation wavelengths
%   A vector containing the wavelengths at which `dispersionfun` was calculated.
%   If 'spectral_options.resolution' is zero or is missing, `lambda_dispersion`
%   is equal to 'spectral_options.bands_in'. Otherwise, it is calculated as
%   described in the following section.
%
% ## Algorithm
%
% In the case where spectral resampling or conversion to colour is not required,
% `W` is created by first finding the undistorted position of each pixel in the
% distorted image, using `dispersionfun`. Second, the function determines
% bilinear interpolation weights on the pixels in the undistorted image
% corresponding to the undistorted position, and enters these weights into `W`.
% Therefore, `W(j + n_px * (k - 1), i + n_px * (k - 1))` is the weight of
% the i-th pixel in the undistorted image in the bilinear interpolation
% calculation yielding the j-th pixel in the distorted image, for the k-th
% colour channel or spectral band.
%
% If spectral resampling or conversion to colour is required, `W` is the product
% of three matrices, `W = C_1 * W_0 * C_2`. `C_2` converts the undistorted image
% to the spectral sampling space used to calculate dispersion. `W_0` warps the
% image, and `C_1` converts the result to the output spectral bands or colour
% channels.
%
% When `spectral_options.resolution` is greater than zero, the function
% determines the wavelengths at which to evaluate `dispersionfun` as follows:
% - Select the four pixels that are the centers of the four quadrants of the
%   distorted image. (The center of the image is not selected, because
%   dispersion is often negligible in the image center. The corners of the image
%   are not selected because dispersion is not modelled with as much accuracy at
%   the edges of the image.)
% - Select lower and upper wavelengths, `lambda_0` and `lambda_1`, at which to
%   evaluate dispersion by finding the smallest and largest wavelengths,
%   respectively, across `spectral_options.bands_in`,
%   `spectral_options.bands_out`, and `spectral_options.color_bands`.
% - Evaluate `dispersionfun` at `lambda_0` and `lambda_1`, at each of the four
%   pixels previously selected.
% - Find the longest displacement between the dispersion vectors for `lambda_0`
%   and `lambda_1` across the four pixels.
% - Divide the longest displacement by `spectral_options.resolution` to obtain
%   the number of sub-intervals, `n_intervals`, into which to divide the
%   wavelength interval from lambda_0 to lambda_1.
% - The wavelengths at which to evaluate the model of dispersion are
%   `lambda_dispersion = linspace(lambda_0, lambda_1, n_intervals + 1)`.
%
% ## Notes
% - Previous versions of this function were able to change the image size,
%   position, and resolution (upsampling only) during warping. This extra
%   functionality was eliminated to simplify the code.
% - There may be pixels in the distorted image whose undistorted coordinates are
%   on or outside the border of the undistorted image. If so, the border pixels
%   of the undistorted image will be replicated, such that the four pixels used
%   for bilinear interpolation will not be all unique.
%
% ## References
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.
%   - See the end of Section 3 on computing the values of a warped image.
% - J. Brauers and T. Aach. "Geometric Calibration of Lens and Filter
%   Distortions for Multispectral Filter-Wheel Cameras," IEEE Transactions
%   on Image Processing, vol. 20, no. 2, pp. 496-505, 2011.
%   - Bilinear interpolation is used to compensate for chromatic aberration
%     (page 501).
% - Bilinear interpolation formulae:
%   https://en.wikipedia.org/wiki/Bilinear_interpolation
%
% See also xylambdaPolyfit, xylambdaSplinefit, makeDispersionfun, warpImage

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 6, 2018

nargoutchk(1, 2);
narginchk(3, 5);

% Parse input arguments
if size(spectral_options.bands_in, 2) > size(spectral_options.bands_in, 1)
    spectral_options.bands_in = spectral_options.bands_in.';
end
n_lambda = length(spectral_options.bands_in);

image_passed = (ndims(image_sampling) == 3) ||...
    ~(all(size(image_sampling) == [2 1]) ||...
      all(size(image_sampling) == [1 2]));
if image_passed
    I_in = image_sampling;
    image_sampling = [size(I_in, 1), size(I_in, 2)];
    if size(I_in, 3) ~= n_lambda
        error('`I_in` must have a size of `length(spectral_options.bands_in)` in its third dimension.');
    end
else
    if length(image_sampling) ~= 2
        error('The `image_sampling` input argument must contain an image height and width only.');
    end
end

negate = false;
offset = [0, 0];
if ~isempty(varargin)
    negate = varargin{1};
    if length(varargin) > 1
        offset = varargin{2};
    end
end

% Define neighbours for bilinear interpolation
offsets = (0.5 + eps) * [
     1 -1; % Top right (Q21)
    -1 -1; % Top left (Q11)
    -1  1; % Bottom left (Q12)
     1  1; % Bottom right (Q22)
     ];
n_offsets = size(offsets, 1);

% Build the solution one image column at a time
n_px = prod(image_sampling);
n_px_lambda = n_px * n_lambda;
if image_passed
    W = zeros([image_sampling, n_lambda]);
else
    W = spalloc(...
        n_px_lambda, n_px_lambda, n_px_lambda * n_offsets...
    );
end

% Preallocation
n_col_lambda = image_sampling(1) * n_lambda;
neighbour_x = zeros(n_col_lambda, n_offsets);
neighbour_y = zeros(n_col_lambda, n_offsets);
neighbour_weights = zeros(n_col_lambda, n_offsets);
index_out_lambda = reshape(repelem(1:n_lambda, image_sampling(1)), [], 1);
neighbour_index_lambda = repmat(index_out_lambda, n_offsets, 1);

for col = 1:image_sampling(2)
    % Enumerate the positions of pixels in the distorted image
    [X, Y] = meshgrid(col, 1:image_sampling(1));
    x_all_bands = repmat(X(:) - 0.5, n_lambda, 1); % Place coordinates at pixel centres
    y_all_bands = repmat(Y(:) - 0.5, n_lambda, 1);
    
    % Find the undistorted positions of all pixels
    lambda_all_bands = reshape(repelem(spectral_options.bands_in, image_sampling(1)), [], 1);
    disparity = dispersionfun([...
        x_all_bands + offset(1),...
        y_all_bands + offset(2),...
        lambda_all_bands...
        ]);
    if negate
        disparity = -disparity;
    end
    x_in_all_bands = x_all_bands + disparity(:, 1);
    y_in_all_bands = y_all_bands + disparity(:, 2);

    % Find the positions of neighbouring pixels in the undistorted image
    for i = 1:n_offsets
        neighbour_x(:, i) = x_in_all_bands + offsets(i, 1);
        neighbour_y(:, i) = y_in_all_bands + offsets(i, 2);
    end
    neighbour_index_x = ceil(reshape(neighbour_x, [], 1));
    neighbour_index_y = ceil(reshape(neighbour_y, [], 1));

    % Convert back to pixel coordinates
    neighbour_x = reshape(neighbour_index_x - 0.5, [], n_offsets);
    neighbour_y = reshape(neighbour_index_y - 0.5, [], n_offsets);

    % Replicate pixels outside the image boundaries
    neighbour_index_x(neighbour_index_x < 1) = 1;
    neighbour_index_x(neighbour_index_x > image_sampling(2)) = image_sampling(2);
    neighbour_index_y(neighbour_index_y < 1) = 1;
    neighbour_index_y(neighbour_index_y > image_sampling(1)) = image_sampling(1);

    % Find bilinear interpolation weights
    x1 = neighbour_x(:, 2);
    x2 = neighbour_x(:, 1);
    y1 = neighbour_y(:, 2);
    y2 = neighbour_y(:, 3);
    dx = x2 - x1;
    dy = y2 - y1;
    tox2 = x2 - x_in_all_bands;
    tox1 = x_in_all_bands - x1;
    toy2 = y2 - y_in_all_bands;
    toy1 = y_in_all_bands - y1;

    neighbour_weights(:, 1) = tox1 .* toy2; % Q21
    neighbour_weights(:, 2) = tox2 .* toy2; % Q11
    neighbour_weights(:, 3) = tox2 .* toy1; % Q12
    neighbour_weights(:, 4) = tox1 .* toy1; % Q22
    neighbour_weights = neighbour_weights ./ (dx .* dy);

    % Find linear indices for pixels
    neighbour_index_linear = sub2ind(...
        [image_sampling, n_lambda],...
        neighbour_index_y,...
        neighbour_index_x,...
        neighbour_index_lambda...
    );
    index_out_linear = (1:image_sampling(1)).' + (col - 1) * image_sampling(1) +...
        n_px * (index_out_lambda - 1);

    if image_passed
        for i = 1:n_offsets
            W(index_out_linear) = W(index_out_linear) + (...
                neighbour_weights(:, i) .* ...
                I_in(neighbour_index_linear(((i - 1) * n_col_lambda + 1):(i * n_col_lambda)))...
            );
        end
    else
        W(index_out_linear + (neighbour_index_linear - 1) * n_px_lambda) = neighbour_weights(:);
    end
end

end
