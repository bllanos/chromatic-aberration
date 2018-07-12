function [ I, image_bounds, varargout ] = solvePatches(...
    image_sampling, add_border, J_2D, align, dispersionfun, sensitivity,...
    lambda, patch_size, padding, f, f_args, varargin...
    )
% SOLVEPATCHES  Run an image estimation algorithm on image patches
%
% ## Syntax
% I = solvePatches(...
%   image_sampling, add_border, J, align, dispersionfun, sensitivity,...
%   lambda, patch_size, padding, f, f_args...
% )
% [ I, image_bounds ] = solvePatches(___)
% [ I, image_bounds, I_rgb ] = solvePatches(___)
% [ I, image_bounds, I_rgb, J_full ] = solvePatches(___)
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatches(___)
% [...
%   I, image_bounds, I_rgb, J_full, J_est, f_output...
% ] = solvePatches(___, target_patch [, n_output])
%
% ## Description
% I = solvePatches(...
%   image_sampling, add_border, J, align, dispersionfun, sensitivity,...
%   lambda, patch_size, padding, f, f_args...
% )
%   Run an image estimation algorithm on image patches and stitch together
%   the results from all patches.
%
% [ I, image_bounds ] = solvePatches(___)
%   Additionally returns the boundaries of the output image in the
%   coordinate system of the input image.
%
% [ I, image_bounds, I_rgb ] = solvePatches(___)
%   Additionally returns the RGB equivalent of the output image.
%
% [ I, image_bounds, I_rgb, J_full ] = solvePatches(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped by the dispersion model.
%
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatches(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% [...
%   I, image_bounds, I_rgb, J_full, J_est, f_output...
% ] = solvePatches(___, target_patch, n_output)
%   Run an image estimation algorithm on a single image patch, and return
%   its additional output arguments.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the height and width, respectively, of
%   the output image `I`, in pixels.
%
% add_border -- Image boundary conditions
%   A Boolean value indicating whether or not the `image_bounds` input
%   argument of 'dispersionfunToMatrix()' should be empty. If `true`, the
%   output image `I` will be large enough to contain the un-warped
%   coordinates of all pixels in `J`. If `false`, the output image `I` will
%   be clipped to the region occupied by `J`.
%
% J -- Input RAW image
%   A 2D array containing the raw colour-filter pattern data of an image.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% dispersionfun -- Model of dispersion
%   A function handle, produced by 'makeDispersionfun()'.
%   `dispersionfun(X)`, where `X` is a three-element row vector (x, y,
%   lambda), returns the dispersion vector for the position (x, y) in `J`
%   corresponding to light with wavelength or colour channel index
%   `lambda`. The dispersion vector corrects for lateral chromatic
%   aberration by pointing from the corresponding position in the reference
%   spectral band or colour channel to position (x, y). This function will
%   negate the dispersion vectors produced by `dispersionfun()` in order to
%   create a warp matrix from `I` to `J`.
%
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th input colour channel or spectral
%   band of `I`. `sensitivity` is a matrix mapping colours in `I` to
%   colours in `J`.
%
% lambda -- Wavelength bands
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersionfun`. 'c' is the desired number of spectral bands or colour
%   channels in `I`.
%
% patch_size -- Patch size
%   A two-element vector containing the height and width, respectively, of
%   the image patches to be estimated. `patch_size` does not include
%   padding used to eliminate artifacts from the patch-wise estimation.
%
% padding -- Patch padding
%   A scalar containing the pixel width of the border surrounding each
%   image patch. When running the image estimation algorithm, this function
%   passes the algorithm a patch of the input image which is large enough
%   to estimate the patch of dimensions `patch_size` in the output image as
%   well as a border of width `padding` around the patch. (Note that the
%   patch of the input image accommodates for the warping given by
%   `dispersionfun`.) The final output image, `I`, is assembled from only
%   the central regions of patches, so this function discards the border
%   region estimated by the algorithm.
%
% f -- Image estimation algorithm
%   The handle to a function implementing an image estimation algorithm.
%   The first output argument of the function must be the output latent
%   image to be estimated. The first few input arguments of the function
%   must be the following:
%   - image_sampling: A two-element vector containing the height and width,
%     respectively, of the output image patch.
%   - align: A four-character character vector, specifying the Bayer tile
%     pattern of the input image patch from `J`.
%   - dispersion_matrix: A warp matrix from coordinates in the output image
%     patch to coordinates in the input image patch.
%   - sensitivity: A copy of this function's argument of the same name.
%   - lambda: A copy of this function's argument of the same name.
%   - J_patch: A patch of the input image `J`.
%
% f_args -- Additional image estimation algorithm parameters
%   A cell vector of input arguments to `f`, beyond those listed above.
%
% target_patch -- Single patch coordinates
%   A two-element vector containing the image coordinates of the top-left
%   corner of the image patch to be estimated. When `target_patch` is
%   passed, all output arguments are calculated for a single image patch,
%   rather than for the entire image. While a border around `I` will have
%   been estimated, with a width given by `padding`, it will not be
%   included in the output. The border region will not be used when
%   calculating `I_rgb`, `J_full`, and `J_est`.
%
% n_output -- Number of additional output arguments
%   The number of additional output arguments to request from `f` beyond
%   the output image patch `I`. These additional output arguments are
%   returned in `f_output`. Defaults to zero if not passed.
%
% ## Output Arguments
%
% I -- Latent image
%   An image_sampling(1) x image_sampling(2) x length(lambda) array,
%   storing the latent image estimated on a patch-wise basis using the
%   function `f`.
%
%   If `target_patch` is passed, then `I` is an patch_size(1) x
%   patch_size(2) x length(lambda) array containing the patch of the latent
%   image with its top left corner at position `target_patch`.
%
% image_bounds -- Latent image coordinate frame
%   The boundaries of `I` expressed in the coordinate frame of `J`.
%   `image_bounds` has the form of the output argument of the same name of
%   'dispersionfunToMatrix()'. If `add_border` is `false`, and
%   `target_patch` is not passed, `image_bounds` will be equal to `[0, 0,
%   size(J, 2), size(J, 1)]`.
%
% I_rgb -- Latent RGB image
%   The RGB equivalent of the latent image, generated using the colour
%   space conversion data in `sensitivity`. An image_sampling(1) x
%   image_sampling(2) x 3 array.
%
%   If `target_patch` is passed, then `I_rgb` is an patch_size(1) x
%   patch_size(2) x length(lambda) array.
%
% J_full -- Warped latent RGB image
%   An RGB image produced by warping `I` according to the dispersion model,
%   followed by conversion to the RGB colour space of `J`. An size(J, 1) x
%   size(J, 2) x 3 array.
%
%   If `target_patch` is passed, then `J_full` is close in size to a
%   patch_size(1) x patch_size(2) x 3 array, but its exact spatial
%   dimensions depend on the nature of the warp between `I` and `J`.
%
% J_est -- Re-estimated input RAW image
%   The mosaiced version of `J_full`, representing the result of passing
%   `I` through the forward model of dispersion and image capture. An array
%   with the same dimensions as `J`.
%
%   If `target_patch` is passed, then `J_est` is a 2D array with the same
%   sizes in its first two dimensions as `J_full`.
%
% f_output -- Additional output arguments of the image estimation algorithm
%   When `target_patch` is passed, this function can return additional
%   output arguments from `f`, as there is no concern over how to combine
%   them from multiple image patches. `f_output` is a cell vector of length
%   `n_output`.
%
% ## References
%
% This function combines image patches by placing them side-by-side and
% discarding the overlapping regions. There are other methods for combining
% overlapping image patches, such as by averaging them. The following
% article contains a brief discussion of patch-wise solutions, in the
% context of chromatic aberration correction:
%
% Sun, T., Peng, Y., and Heidrich, W. (2017). "Revisiting cross-channel
%   information transfer for chromatic aberration correction." In 2017 IEEE
%   International Conference on Computer Vision (ICCV) (pp. 3268–3276).
%   doi:10.1109/ICCV.2017.352
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 12, 2018

narginchk(11, 13);

end