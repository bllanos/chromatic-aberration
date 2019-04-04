function [ W, bands_dispersion ] = dispersionfunToMatrix(dispersionfun, options, image_sampling, varargin)
% DISPERSIONFUNTOMATRIX  Convert a warp model to a warp matrix, or use it
% directly to warp an image
%
% ## Syntax
% W = dispersionfunToMatrix(...
%   dispersionfun, options, image_sampling [, negate, offset]...
% )
% I_out = dispersionfunToMatrix(...
%   dispersionfun, options, I_in [, negate, offset]...
% )
% [ ____, bands_dispersion ] = dispersionfunToMatrix(____)
%
% ## Description
% W = dispersionfunToMatrix(...
%   dispersionfun, options, image_sampling [, negate, offset]...
% )
%   Returns a matrix for distorting a vectorized image.
%
% I_out = dispersionfunToMatrix(...
%   dispersionfun, options, I_in [, negate, offset]...
% )
%   Returns the distorted version of the input image.
%
% [ ____, bands_dispersion ] = dispersionfunToMatrix(____)
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
% options -- Spectral sampling options
%   `options` is a structure of options controlling the spectral
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
%     the warped images must be sampled (prior to conversion to colour, if
%     colour conversion is enabled). If 'bands_out' is missing, it will be set
%     equal to 'bands_in'. The values in 'bands_out' are expected to be
%     evenly-spaced.
%   - 'color_map': A 2D array with 'm' rows, where `color_map(i, j)` is the
%     sensitivity of the i-th colour channel of the warped image to the j-th
%     spectral band in `bands_out`. `color_map` is not a colour conversion
%     matrix, as it does not perform the desired numerical integration, over the
%     spectrum, that is part of colour conversion.
%   - 'support_threshold': A fraction indicating what proportion of the
%     peak magnitude of a colour channel's sensitivity function the user
%     considers to be effectively zero (i.e. no sensitivity).
%     'support_threshold' is used to find the endpoints of the sequence of
%     wavelengths at which dispersion is evaluated. 'support_threshold' is not
%     used if 'color_map' is not present, and defaults to zero if not passed.
%   - 'int_method': The numerical integration method to use when
%     integrating over the responses of colour channels to compute colour
%     values. `int_method` is used by `colorWeights()`, but is not used if
%     'color_map' is not present.
%   - 'bands_padding': Padding to use when computing spectral interpolation
%     matrices. Refer to the documentation of the 'padding' input argument of
%     'resamplingWeights()' for details, in 'resamplingWeights.m'.
%   - 'interpolant': The function convolved with spectral signals to interpolate
%     them from the sampling space of 'bands_in' to the sampling space of
%     'bands_dispersion' or 'bands_out'. 'interpolant' is passed to
%     'resamplingWeights()' as its `f` input argument. Refer to the
%     documentation of 'resamplingWeights.m' for more details.
%   - 'interpolant_ref': Similar to 'interpolant', but for interpolation from
%     the sampling space of 'bands_dispersion' to the sampling space of
%     'bands_out'.
%
%   The above optional fields should not be present if the dispersion matrix is
%   to operate on colour channels instead of spectral bands, because colour
%   channels cannot be resampled. It is not correct to instead set 'resolution'
%   to zero, and 'bands_out' to 'bands_in' when working with colour channels,
%   because 'interpolant' may still not result in an identity mapping.
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
%   - n is `length(options.bands_in)`, and
%   - m is `length(options.bands_out)`
%
%   `W` warps the undistorted image to the distorted image, and may
%   simutaneously change the spectral sampling of the image or convert it to
%   colour, depending on the settings in `options`. `W` operates
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
% bands_dispersion -- Dispersion evaluation wavelengths
%   A vector containing the wavelengths at which `dispersionfun` was evaluated.
%   If 'options.resolution' is zero or is missing, `bands_dispersion`
%   is equal to 'options.bands_in'. Otherwise, it is calculated as
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
% When `options.resolution` is greater than zero, the function
% determines the wavelengths at which to evaluate `dispersionfun` as follows:
% - Select the four pixels that are the centers of the four quadrants of the
%   distorted image. (The center of the image is not selected, because
%   dispersion is often negligible in the image center. The corners of the image
%   are not selected because dispersion is not modelled with as much accuracy at
%   the edges of the image.)
% - Select lower and upper wavelengths, `lambda0` and `lambda1`, at which to
%   evaluate dispersion, by finding the smallest and largest wavelengths,
%   respectively, across `options.bands_in`, and
%   `options.bands_out`. If `options.support_threshold` and
%   `options.color_map` are both present, then `lambda0` and `lambda1` are
%   constrained to the region of the spectrum in which `options.color_map` has
%   values which are fractions of at least `options.support_threshold` of the
%   peak values of their respective channels.
% - Evaluate `dispersionfun` at `lambda0` and `lambda1`, at the four pixels
%   previously selected.
% - Find the mean displacement between the dispersion vectors for `lambda0`
%   and `lambda1` across the four pixels.
% - Divide the mean displacement by `options.resolution` to obtain
%   the number of sub-intervals, `n_intervals`, into which to divide the
%   wavelength interval from `lambda0` to `lambda1`.
% - The wavelengths at which to evaluate the model of dispersion are
%   `bands_dispersion = linspace(lambda_0, lambda_1, n_intervals + 1)`.
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
% See also xylambdaPolyfit, xylambdaSplinefit, makeDispersionfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 6, 2018

nargoutchk(1, 2);
narginchk(3, 5);

% Parse input arguments
bands_in = options.bands_in;
if size(bands_in, 2) > size(bands_in, 1)
    bands_in = bands_in.';
end
n_bands_in = length(bands_in);

image_passed = (ndims(image_sampling) == 3) ||...
    ~(all(size(image_sampling) == [2 1]) ||...
      all(size(image_sampling) == [1 2]));
if image_passed
    I_in = image_sampling;
    image_sampling = [size(I_in, 1), size(I_in, 2)];
    if size(I_in, 3) ~= n_bands_in
        error('`I_in` must have a size of `length(options.bands_in)` in its third dimension.');
    end
    I_in = reshape(I_in, [], 1);
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

% Process spectral sampling options
do_pre_resampling = isfield(options, 'resolution') && (options.resolution ~= 0);
has_bands_out = isfield(options, 'bands_out');
do_post_resampling = has_bands_out || do_pre_resampling;
if has_bands_out
    bands_out = options.bands_out;
else
    bands_out = bands_in;
end
if size(bands_out, 2) > size(bands_out, 1)
    bands_out = bands_out.';
end
n_bands_out = length(bands_out);
if has_bands_out && (n_bands_out > 1)
    diff_bands = diff(bands_out);
    if max(abs(diff_bands - diff_bands(1))) > 1e-6
        error('`options.bands_out` must contain equally-spaced values.')
    end
end

has_color_map = isfield(options, 'color_map');
if has_color_map
    if ~has_bands_out
        error('If `options.color_map` exists, then `options.bands_out` must also exist.');
    end
    color_map = options.color_map;
    if size(color_map, 2) ~= n_bands_out
        error('`options.color_map` must have as many columns as the length of `options.bands_out`.');
    end
end

if do_pre_resampling
    if options.resolution < 0
        error('`options.resolution` must be non-negative.');
    end
    
    % Choose a sequence of wavelengths at which to evaluate dispersion
    lambda0 = min(min(bands_in), min(bands_out));
    lambda1 = max(max(bands_in), max(bands_out));
    if has_color_map && isfield(options, 'support_threshold')
        if options.support_threshold < 0
            error('`options.support_threshold` must be non-negative.');
        end
        color_map_relative = abs(color_map) ./ repmat(max(abs(color_map), [], 2), 1, n_bands_out);
        color_map_relative_logical = any(color_map_relative >= options.support_threshold, 1);
        ind = find(color_map_relative_logical, 1, 'first');
        lambda0 = max(bands_out(ind), lambda0);
        ind = find(color_map_relative_logical, 1, 'last');
        lambda1 = min(bands_out(ind), lambda1);
    end
    sample_points = [
        (image_sampling(2) / 4 * [1; 3; 1; 3]) + offset(1),...
        (image_sampling(1) / 4 * [1; 1; 3; 3]) + offset(2)
    ] - 0.5;
    n_sample_points = size(sample_points, 1);
    sample_points = [repelem(sample_points, 2, 1), repmat([lambda0; lambda1], n_sample_points, 1)]; 
    distance_samples = dispersionfun(sample_points);
    distance_samples = diff(distance_samples, 1, 1);
    distance_samples = distance_samples(1:2:end, :);
    distance_samples = dot(distance_samples, distance_samples, 2);
    distance_samples = mean(sqrt(distance_samples));
    n_bands_dispersion = max(ceil(distance_samples / options.resolution), 1) + 1;
    bands_dispersion = linspace(lambda0, lambda1, n_bands_dispersion);
else
    bands_dispersion = bands_in;
    n_bands_dispersion = n_bands_in;
end

if n_bands_in > 1 && (do_pre_resampling || do_post_resampling)
    diff_bands = diff(bands_in);
    if max(abs(diff_bands - diff_bands(1))) > 1e-6
        error('`options.bands_in` must contain equally-spaced values.')
    end
end

% Create spectral resampling matrices
if do_pre_resampling
    pre_map = resamplingWeights(...
        bands_dispersion, bands_in, options.interpolant, options.bands_padding...
    );
end
if do_post_resampling
    color_weights_options = options;
    if do_pre_resampling
        post_interpolant = options.interpolant;
    else
        post_interpolant = options.interpolant_ref;
        color_weights_options.interpolant = options.interpolant_ref;
    end
    if has_color_map
        post_map = colorWeights(...
            color_map, bands_out, bands_dispersion, color_weights_options...
        );
        n_bands_out = size(color_map, 1);
    else
        post_map = resamplingWeights(...
            bands_out, bands_dispersion, post_interpolant, options.bands_padding...
        );
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

% Allocate output array
n_px = prod(image_sampling);
n_px_lambda_in = n_px * n_bands_in;
n_px_lambda_out = n_px * n_bands_out;
n_px_lambda_dispersion = n_px * n_bands_dispersion;
if image_passed
    W = zeros([image_sampling, n_bands_out]);
else
    n_values_max = n_px_lambda_out * n_offsets * n_bands_in;
    n_values = 0;
    W_rows = zeros(n_values_max, 1);
    W_cols = zeros(n_values_max, 1);
    W_vals = zeros(n_values_max, 1);
end

% Convert per-pixel spectral sampling conversion matrices into sparse matrices
% which can operate on all pixels in an image column
n_col_lambda = image_sampling(1) * n_bands_dispersion;
if do_post_resampling
    post_map = sparse(...
        repmat(repmat((0:(n_bands_out - 1)).', n_bands_dispersion, 1), image_sampling(1), 1) * image_sampling(1) + repelem((1:image_sampling(1)).', n_bands_out * n_bands_dispersion, 1),...
        repmat(repelem((0:(n_bands_dispersion - 1)).', n_bands_out, 1), image_sampling(1), 1) * image_sampling(1) + repelem((1:image_sampling(1)).', n_bands_out * n_bands_dispersion, 1),...
        repmat(reshape(post_map, [], 1), image_sampling(1), 1),...
        image_sampling(1) * n_bands_out, n_col_lambda...
    );
end
if do_pre_resampling
    % This will consume a lot of memory for large images
    pre_map = sparse(...
        repmat(repmat((0:(n_bands_dispersion - 1)).', n_bands_in, 1), n_px, 1) * n_px + repelem((1:n_px).', n_bands_in * n_bands_dispersion, 1),...
        repmat(repelem((0:(n_bands_in - 1)).', n_bands_dispersion, 1), n_px, 1) * n_px + repelem((1:n_px).', n_bands_in * n_bands_dispersion, 1),...
        repmat(reshape(pre_map, [], 1), n_px, 1),...
        n_px_lambda_dispersion, n_px_lambda_in...
    );
end

% Preallocation before the loop
indices_col = repmat((1:n_col_lambda).', n_offsets, 1);
neighbour_x = zeros(n_col_lambda, n_offsets);
neighbour_y = zeros(n_col_lambda, n_offsets);
neighbour_weights = zeros(n_col_lambda, n_offsets);
index_out_lambda = repelem((1:n_bands_dispersion).', image_sampling(1), 1);
neighbour_index_lambda = repmat(index_out_lambda, n_offsets, 1);
lambda_all_bands = reshape(repelem(bands_dispersion, image_sampling(1)), [], 1);

% Build the output array one image column at a time
for col = 1:image_sampling(2)
    % Enumerate the positions of pixels in the distorted image
    [X, Y] = meshgrid(col, 1:image_sampling(1));
    x_all_bands = repmat(X(:) - 0.5, n_bands_dispersion, 1); % Place coordinates at pixel centres
    y_all_bands = repmat(Y(:) - 0.5, n_bands_dispersion, 1);
    
    % Find the undistorted positions of all pixels
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
    %
    % MATLAB's sparse matrix constructor will add values for coefficients given
    % the same index in the sparse matrix.
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
        [image_sampling, n_bands_dispersion],...
        neighbour_index_y,...
        neighbour_index_x,...
        neighbour_index_lambda...
    );
    
    % Build the dispersion and resampling matrix for this column
    W_col = sparse(...
        indices_col,...
        neighbour_index_linear,...
        reshape(neighbour_weights, [], 1),...
        n_col_lambda, n_px_lambda_dispersion...
        );
    if do_post_resampling
        W_col_post = post_map * W_col;
    else
        W_col_post = W_col;
    end
    if do_pre_resampling
        % Notice that this is done after applying `post_map`, as `post_map`
        % usually reduces the size of the matrix, whereas `pre_map` usually
        % increases the size of the matrix.
        W_col_all = W_col_post * pre_map;
    else
        W_col_all = W_col_post;
    end

    if image_passed
        W(:, col, :) = reshape(W_col_all * I_in, image_sampling(1), 1, n_bands_out);
    else
        W_rows_col = repmat(...
                repmat(...
                    ((image_sampling(1) * (col - 1) + 1):(image_sampling(1) * col)).',...
                    n_bands_out, 1 ...
                ) + repelem(...
                    (0:n_px:(n_px * (n_bands_out - 1))).',...
                    image_sampling(1), 1 ...
                ),...
            n_px_lambda_in, 1 ...
        );
        W_cols_col = repelem((1:n_px_lambda_in).', image_sampling(1) * n_bands_out, 1);
        [row_filter, col_filter, W_vals_col] = find(W_col_all);
        linear_filter = sub2ind(size(W_col_all), row_filter, col_filter);
        n_values_new = n_values + length(W_vals_col);
        W_rows((n_values + 1):n_values_new) = W_rows_col(linear_filter);
        W_cols((n_values + 1):n_values_new) = W_cols_col(linear_filter);
        W_vals((n_values + 1):n_values_new) = W_vals_col;
        n_values = n_values_new;
    end
end

if ~image_passed
    W = sparse(...
        W_rows(1:n_values),...
        W_cols(1:n_values),...
        W_vals(1:n_values),...
        n_px_lambda_out, n_px_lambda_in...
    );
end

end
