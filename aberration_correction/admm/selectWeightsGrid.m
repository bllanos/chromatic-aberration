function [ weights, patch_lim, I_patch, varargout ] = selectWeightsGrid(...
    J, align, dispersionfun, sensitivity,...
    lambda, rho, baek2017Algorithm2Options, options, varargin...
    )
% SELECTWEIGHTSGRID  Use the grid search method to select regularization weights
%
% ## Syntax
% weights = selectWeightsGrid(...
%   J, align, dispersionfun, sensitivity, lambda, rho,...
%   baek2017Algorithm2Options, options [, target_patch, verbose]...
% )
% [ weights, patch_lim ] = selectWeightsGrid(___)
% [ weights, patch_lim, I_patch ] = selectWeightsGrid(___)
% [ weights, patch_lim, I_patch, search ] = selectWeightsGrid(___)
%
% ## Description
% weights = selectWeightsGrid(...
%   J, align, dispersionfun, sensitivity, lambda, rho,...
%   baek2017Algorithm2Options, options [, target_patch, verbose]...
% )
%   Returns the regularization weights selected using the minimum distance
%   criterion grid search method of Song et al. 2016
%
% [ weights, patch_lim ] = selectWeightsGrid(___)
%   Additionally returns the boundaries of the image patch used for weight
%   selection
%
% [ weights, patch_lim, I_patch ] = selectWeightsGrid(___)
%   Additionally returns the latent image patch estimated using the
%   selected weights.
%
% [ weights, patch_lim, I_patch, search ] = selectWeightsGrid(___)
%   Additionally returns the path taken by the grid search algorithm used
%   to select the weights.
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
%   as the `sensorAlignment` input argument of `demosaic()`. `align` can
%   also be empty, indicating that the input image is not mosaiced.
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
% rho -- Penalty parameters
%   The `rho` input argument to use with 'baek2017Algorithm2()'.
%
% baek2017Algorithm2Options -- ADMM image estimation options
%   The `options` input argument to use with 'baek2017Algorithm2()'.
%
% target_patch -- Patch coordinates
%   A two-element vector containing the row and column, respectively, of
%   the top-left corner of the image patch to be estimated. The grid search
%   method will only be applied to a patch of the image, for efficiency,
%   and also to customize the results to image features of interest. Note
%   that the patch has the same coordinates in both the input image and the
%   latent image, as the two images are aligned. (Refer to the 'Notes'
%   section below for the rationale.)
%
%   When `target_patch` empty, or is not passed, a figure will
%   be opened for the user to select a pixel as the centre of the target
%   patch.
%
% options -- Options and small parameters
%   A structure with the following fields:
%   - 'int_method': The numerical integration method used for spectral to
%     colour space conversion. `int_method` is passed to
%     'projectionMatrix()' as its `int_method` argument. Refer to the
%     documentation of 'projectionMatrix.m' for details. If 'int_method' is
%     'none', as should be the case when colour conversion is from a set of
%     colour channels, not a spectral space, numerical integration will not
%     be performed.
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of the image patch to be estimated. `patch_size`
%     includes any padding used to eliminate border artifacts, although
%     such padding is not represented in this function, but must be
%     enforced in the manner in which `f` calculates its second output
%     argument. (Refer to the documentation of `f` above for details.)
%   - 'enabled_weights': A logical vector with a length equal to the number
%     of regularization terms in the image estimation problem.
%     `enabled_weights(i)` is `true` if the i-th regularization term will
%     be included in the grid search method. If `enabled_weights(i)` is
%     `false`, the weight for the i-th regularization term will be set to
%     zero, rather than set using the grid search method.
%   - 'maxit': The maximum number of grid refinement iterations to perform
%     in the method of Song et al. 2016.
%   - 'clip_weights': If `true`, the origin of the minimum distance
%     criterion will be set based on 'minimum_weights' and
%     'maximum_weights', rather than chosen semi-automatically, and the
%     output weights will be constrained to lie between the minimum and
%     maximum weights.
%   - 'minimum_weights': A vector, of the same length as 'enabled_weights',
%     specifying minimum values for the regularization weights. The
%     elements of 'minimum_weights' will be used if the matrix operator of
%     the data-fitting term in the image estimation problem is singular.
%     They will also be used if 'clip_weights' is `true`. Refer to section
%     3.4 of Belge et al. 2002 for details: If the matrix is singular, "an
%     appropriate multiple of machine precision" is suggested, or a value
%     "just large enough so that the regularization problem is well
%     defined."
%   - 'maximum_weights': A vector, of the same length as 'enabled_weights',
%     specifying maximum values for the regularization weights. These
%     values will be used if 'clip_weights' is `true`, and otherwise this
%     field is optional.
%   - 'tol': The first element is the threshold value of the relative
%     change in any regularization weight from one iteration to the next.
%     When all weights change less than this threshold, iteration
%     terminates. Refer to equation 31 of Belge et al. 2002. (Belge et al.
%     set a value of 10^-4, but suggest that larger tolerances should be
%     acceptable.)
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show the progress of the
%   iterative search for `weights`, as well as warnings when the search
%   behaves in unexpected ways.
%
% ## Output Arguments
%
% weights -- Selected regularization weights
%   The value of the `weights` input argument to `f`, chosen by the minimum
%   distance criterion grid search method of Song et al. 2016. `weights` is
%   a vector of the same length as `options.enabled_weights`.
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
%   the location of the patch relative to the image boundaries. Note that
%   the patch of the latent image corresponds to a patch in the input image
%   `J` at the same coordinates. (Refer to the 'Notes' section below for
%   details.)
%
% I_patch -- Latent image patch
%   The patch of the latent image estimated under the regularization
%   weights output in `weights`. The first two dimensions of `I_patch` are
%   defined by `patch_lim`, while its third dimension has a size equal to
%   the length of `lambda`.
%
% search -- Grid search method search path
%   A structure with the following fields:
%   - 'weights': An array of dimensions n x length(options.enabled_weights),
%     where 'n' is the number of iterations used by the grid search
%     algorithm to select the final regularization weights. The `weights`
%     output argument is copied to `search.weights(n, :)`.
%     `search.weights(i, :)` is the estimate of the regularization weights
%     at the i-th iteration.
%   - 'err': An array of dimensions n x (length(options.enabled_weights) + 1)
%     containing the second output argument of 'baek2017Algorithm2()' at
%     each iteration. `err(i)` is the second output argument of
%     'baek2017Algorithm2()' when called using `search.weights(i, :)` as
%     the `weights` input argument 'baek2017Algorithm2()'.
%   - 'origin': The reference point of the minimum distance criterion used
%     for the grid search algorithm. A vector of length
%     (length(options.enabled_weights) + 1), where the first element
%     corresponds to the residual, and the remaining elements correspond to
%     the regularization terms. The elements corresponding to disabled
%     regularization terms are set to zero.
%   - 'origin_min_weights': The minimum regularization weight values,
%     corresponding to the first element of 'origin'. A vector with the
%     same length as `options.enabled_weights`, with the elements
%     corresponding to disabled regularization terms set to zero.
%   - 'origin_max_weights': The maximum regularization weight values,
%     corresponding to the subsequent elements of 'origin'. A vector with
%     the same length as `options.enabled_weights`, with the elements
%     corresponding to disabled regularization terms set to zero.
%
% ## References
%
% A fixed-point algorithm for selecting regularization weights, a proxy
% for an exhaustive search of the L-hypersurface, is described in:
%
%   Belge, M, Kilmer, M. E., & Miller, E. L.. "Efficient determination of
%     multiple regularization parameters in a generalized L-curve
%     framework." Inverse Problems, vol. 18, pp. 1161-1183, 2002.
%     doi:10.1088/0266-5611/18/4/314
%
% Note that this method is derived for an unconstrained optimization
% problem. It is likely not appropriate to use it in combination with the
% non-negativity constraint supported by 'baek2017Algorithm2()', for
% example. However, I use the same method of setting the reference point of
% the minimum distance function described in this article.
%
% The following article discusses the grid-search method for minimizing the
% minimum distance function, which forms the basis for most of the code in
% this file:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
%
% ## Notes
% - In contrast to 'selectWeights()', this function is more specific to
%   'baek2017Algorithm2()', because it needs to know how to turn off the
%   non-negativity constraint during image estimation, in order to find the
%   reference point of the minimum distance criterion.
% - The dispersion matrix mapping `I_patch` to the `J` input argument
%   of 'baek2017Algorithm2()', given as the `dispersion` input argument of
%   'baek2017Algorithm2()', will be constructed such that `I_patch` and `J`
%   have the same resolutions, and coincide in space. As such, some pixels
%   in `I_patch` may not map to positions in `J`, and therefore will be
%   under-constrained during image estimation. However, determining the
%   boundaries of `I_patch` such that all pixels map to the region of `J`
%   is a difficult problem (when trying to do so without calculating the
%   dispersion matrix for, at worst case, the entire image).
%
% See also selectWeights, projectionMatrix, baek2017Algorithm2,
% solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 4, 2018

end
