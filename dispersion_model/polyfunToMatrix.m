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
%   2 - The pixel y coordinate of the top left corner of the image
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
%   `image_bounds` describes the position and size of the undistorted
%   image in the space of the distorted image. Units are pixels in the
%   distorted image. Note that the last two elements of `image_bounds` do
%   not need to be equal to the elements of `image_sampling_out`, as the
%   two need not have equal pixel sizes.
%
%   If `image_bounds` is empty or not passed, it will be set to the
%   bounding box of the undistorted coordinates of the distorted image.
%   More precisely, all pixels in the distorted image will be generated, as
%   implemented in `W`, by bilinear interpolation of four pixels in the
%   undistorted image; There will not be any pixels in the distorted image
%   which originate from less than four pixels in the undistorted image.
%
% negate -- Warp inversion flag
%   If `true`, the warp vectors calculated by `polyfun` will be negated
%   before use. (Defaults to `false` if not passed.)
%
% ## Output Arguments
%
% W -- Warp matrix
%   A (n_px_in x n)-by-(n_px_out x n) array (n = length(lambda)) which
%   warps the undistorted image to the distorted image, according to the
%   equation:
%     `I_distorted = W * I_undistorted`
%   `I_undistorted` is a vectorized form of an image where all pixels have
%   been rearranged from columnwise order into matrix rows.
%   `I_undistorted(i, k)` is the value of the k-th wavelength band or
%   colour channel at the i-th pixel of the undistorted image. Similarly,
%   `I_distorted(j, k)` is the value of the k-th wavelength band or colour
%   channel at the j-th pixel of the distorted image.
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
% corresponding to the undistored position, and enters these weights into
% `W`. Therefore, `W(j,i)` is the weight of the i-th pixel in the
% undistorted image in the bilinear interpolation calculation yielding the
% j-th pixel in the distorted image.
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
%
% See also xylambdaPolyfit, makePolyfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 6, 2018

nargoutchk(1, 2);
narginchk(3, 6);

end

