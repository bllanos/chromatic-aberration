function [ weights, patch_lim, I_patch ] = selectWeights(...
    J, align, dispersionfun, sensitivity,...
    lambda, options, f, f_args, varargin...
    )
% SELECTWEIGHTS  Use the L-hypersurface method to select regularization weights
%
% ## Syntax
% weights = selectWeights(...
%   J, align, dispersionfun, sensitivity, lambda, options,...
%   f, f_args [, target_patch, verbose]...
% )
% [ weights, patch_lim ] = selectWeights(___)
% [ weights, patch_lim, I_patch ] = selectWeights(___)
% [ weights, patch_lim, I_patch, search ] = selectWeights(___)
%
% ## Description
% weights = selectWeights(...
%   J, align, dispersionfun, sensitivity, lambda, options,...
%   f, f_args [, target_patch, verbose]...
% )
%   Returns the regularization weights selected using the L-hypersurface
%   method of Belge et al. 2002.
%
% [ weights, patch_lim ] = selectWeights(___)
%   Additionally returns the boundaries of the image patch used for weight
%   selection
%
% [ weights, patch_lim, I_patch ] = selectWeights(___)
%   Additionally returns the latent image patch estimated using the
%   selected weights.
%
% [ weights, patch_lim, I_patch, search ] = selectWeights(___)
%   Additionally returns the path taken by the fixed-point iterative
%   algorithm used to select the weights.
%
% ## Input Arguments
%
% J -- Input image
%   A 2D or 3D array containing the input image for the latent image
%   estimation problem.
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
%   create a warp matrix from `I_patch` to `J`.
%
%   `dispersionfun` can be empty (`[]`), if there is no model of
%   dispersion.
%
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th input colour channel or spectral band
%   of `I_patch`. `sensitivity` is a matrix mapping colours in `I_patch` to
%   colours in `J`.
%
% lambda -- Wavelength bands
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersionfun`. 'c' is the desired number of spectral bands or colour
%   channels in `I_patch`.
%
% f -- Image estimation algorithm
%   The handle to a function implementing an image estimation algorithm.
%   The first few input arguments of the function must be the following:
%   - image_sampling: A two-element vector containing the height and width,
%     respectively, of the output image patch.
%   - align: A four-character character vector, specifying the Bayer tile
%     pattern of the input image patch from `J`.
%   - dispersion_matrix: Either a warp matrix from coordinates in the
%     output image patch to coordinates in the input image patch, or an
%     empty array (`[]`), if there is no model of dispersion.
%   - sensitivity: A copy of this function's argument of the same name.
%   - lambda: A copy of this function's argument of the same name.
%   - J_patch: A patch of the input image `J`.
%   - weights: A vector the same length as `options.enabled_weights`,
%     containing the weights on the regularization terms in the image
%     estimation problem.
%
%   The first output argument of the function must be the latent image
%   being estimated, `I_patch`.
%
%   The second output argument of the function must be a vector of
%   `length(options.enabled_weights) + 1` elements. The first element is
%   the data fitting error, the mean squared error between `J_patch` and
%   the version of `J_patch` re-estimated from `I_patch`. The remaining
%   elements are the regularization errors. In other words, they are the
%   norms of the vectors obtained from `I_patch` by the regularization
%   operators used in the image estimation problem. Preferably, the
%   regularization errors are normalized by the numbers of elements in the
%   vectors. Also, it may be desirable for the elements of the second
%   output argument to be calculated after discarding a border from the
%   image, to eliminate border artifacts.
%
% f_args -- Additional image estimation algorithm parameters
%   A cell vector of input arguments to `f`, beyond those listed above.
%   `f_args` can be an empty cell array, `{}`.
%
% target_patch -- Patch coordinates
%   A two-element vector containing the row and column, respectively, of
%   the top-left corner of the image patch to be estimated. The
%   L-hypersurface method will only be applied to a patch of the image, for
%   efficiency. Note that the patch has the same coordinates in both the
%   input image and the latent image, as the two images are aligned. (Refer
%   to notes below for the rationale.)
%
%   When `target_patch` empty, or is not passed, a figure will
%   be opened for the user to select a pixel as the centre of the target
%   patch.
%
% options -- Options and small parameters
%   A structure with the following fields:
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of the image patch to be estimated. `patch_size`
%     includes any padding used to eliminate border artifacts, although
%     such padding is not represented in this function, but must be
%     enforced in the manner in which `f` calculates its second output
%     argument. (Refer to the documentation of `f` above for details.)
%   - 'enabled_weights': A logical vector with a length equal to the number
%     of regularization terms in the image estimation problem.
%     `enabled_weights(i)` is `true` if the i-th regularization term will
%     be included in the L-hypersurface method. If `enabled_weights(i)` is
%     `false`, the weight for the i-th regularization term will be set to
%     zero, rather than set using the L-hypersurface method.
%   - 'minimum_weights': A vector the same length as 'enabled_weights'
%     specifying suggested minimum values for the regularization weights.
%     The elements of 'minimum_weights' will only be used if the matrix
%     used in the data-fitting term is singular. Otherwise, the smallest
%     singular value of the data fitting matrix operator will be used.
%     Refer to section 3.4 of Belge et al. 2002 for details.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show the progress of the
%   iterative search for `weights`.
%
% ## Output Arguments
%
% weights -- Selected regularization weights
%   The value of the `weights` input argument to `f`, chosen by the
%   L-hypersurface method of Belge et al. 2002. `weights` is a vector of
%   the same length as `options.enabled_weights`.
%
% patch_lim -- Patch boundaries
%   A 2 x 2 array, with the following elements:
%   - `patch_lim(1, 1)` is the row index of the top left corner of the
%     patch in the image.
%   - `patch_lim(2, 1)` is the row index of the bottom right corner of the
%     patch in the image.
%   - `patch_lim(1, 2)` is the column index of the top left corner of the
%     patch in the image.
%   - `patch_lim(2, 2)` is the column index of the bottom right corner of
%     the patch in the image.
%
%   The patch defined by `patch_lim` is the patch of the latent image
%   (`I_patch`) estimated when searching for good regularization weights.
%   It may not have the size defined by `options.patch_size`, depending on
%   the location of the patch relative to the image boundaries.
%
% I_patch -- Latent image patch
%   The patch of the latent image estimated under the regularization
%   weights output in `weights`. The first two dimensions of `I_patch` are
%   defined by `patch_lim`, while its third dimension has a size equal to
%   the length of `lambda`.
%
% search -- L-hypersurface method search path
%   A structure with the following fields:
%   - 'weights': An array of dimensions n x length(options.enabled_weights),
%     where 'n' is the number of iterations used by the L-hypersurface
%     fixed-point algorithm to select the final regularization weights.
%     `weights` is equal to `search.weights(n, :)`. `search.weights(i, :)`
%     is the estimate of the regularization weights at the i-th iteration.
%   - 'err': An array of dimensions n x (length(options.enabled_weights) + 1)
%     containing the second output argument of `f` at each iteration.
%
% ## References
%
% The fixed-point algorithm for selecting regularization weights, a proxy
% for an exhaustive search of the L-hypersurface, is described in:
%
%   Belge, M, Kilmer, M. E., & Miller, E. L.. "Efficient determination of
%     multiple regularization parameters in a generalized L-curve
%     framework." Inverse Problems, vol. 18, pp. 1161-1183, 2002.
%     doi:10.1088/0266-5611/18/4/314
%
% ## Notes
% - The dispersion matrix mapping `I_patch` to the `J_patch` input argument
%   of `f`, given as the `dispersion_matrix` input argument of `f`, will be
%   constructed such that `I_patch` and `J_patch` have the same
%   resolutions, and coincide in space. As such, some pixels in `I_patch`
%   may not map to positions in `J_patch`, and therefore will be
%   under-constrained. However, determining the boundaries of `I_patch`
%   such that all pixels map to the region of `J_patch` is a difficult
%   problem (without calculating the dispersion matrix for, at worst case,
%   the entire image).

end