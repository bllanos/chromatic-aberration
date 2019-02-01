function [ I_3D, varargout ] = solvePatchesADMM(...
    I_in, J_2D, align, dispersionfun, sensitivity, lambda,...
    admm_options, reg_options, patch_options, varargin...
)
% SOLVEPATCHESADMM  Run ADMM (loosely) as in Algorithm 2 of Baek et al. 2017, with weight selection and patch-wise decomposition
%
% ## Usage
%
% This is a super function which combines most of the functionality of
% 'baek2017Algorithm2()`, 'selectWeightsGrid()', 'trainWeights()', and
% 'solvePatchesAligned()' for improved efficiency.
%
% ## Syntax
% I = solvePatchesADMM(...
%   I_in, J, align, dispersionfun, sensitivity, lambda,...
%   admm_options, reg_options, patch_options [, verbose]...
% )
% [ I, I_rgb ] = solvePatchesADMM(___)
% [ I, I_rgb, weights_images ] = solvePatchesADMM(___)
% [ I, I_rgb, weights_images, J_full ] = solvePatchesADMM(___)
% [ I, I_rgb, weights_images, J_full, J_est ] = solvePatchesADMM(___)
% [ I, I_rgb, weights_images, J_full, J_est, I_warped ] = solvePatchesADMM(___)
% [ I, I_rgb, weights_images, J_full, J_est, I_warped, search ] = solvePatchesADMM(...)
%
% ## Description
% I = solvePatchesADMM(...
%   I_in, J, align, dispersionfun, sensitivity, lambda,...
%   admm_options, reg_options, patch_options [, verbose]...
% )
%   Estimate a latent RGB or hyperspectral image `I` from dispersion in
%   the input RAW image `J`.
%
% [ I, I_rgb ] = solvePatchesADMM(___)
%   Additionally returns the RGB equivalent of the latent image.
%
% [ I, I_rgb, weights_images ] = solvePatchesADMM(___)
%   Additionally returns images illustrating the weights selected for the
%   regularization terms in the optimization problem.
%
% [ I, I_rgb, weights_images, J_full ] = solvePatchesADMM(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped according to the model of dispersion.
%
% [ I, I_rgb, weights_images, J_full, J_est ] = solvePatchesADMM(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% [ I, I_rgb, weights_images, J_full, J_est, I_warped ] = solvePatchesADMM(___)
%   Additionally returns the version of the latent image warped according
%   to the model of dispersion.
%
% [ I, I_rgb, weights_images, J_full, J_est, I_warped, search ] = solvePatchesADMM(...)
%   Additionally returns the search path taken to select regularizaion
%   weights for a single image patch. This call syntax is available only
%   when `patch_options.target_patch` exists.
%
% ## Input Arguments
%
% I_in -- True image structure
%   `I_in` is used to select regularization weights based on similarity
%   with another image (such as the true image). When `I_in` is empty
%   (`[]`), regularization weights will be selected using the minimum
%   distance criterion described in Song et al. 2016, or, if
%   `reg_options.demosaic` is `true`, based on similarity with a demosaiced
%   version of `J`.
%
%   When `I_in` is not empty, it must be a structure with the following
%   fields:
%   - 'I': A 2D or 3D array containing the image, against which the
%     estimated latent image will be compared. 'I' and `J` must have the
%     same sizes in their first two dimensions. (i.e. They must have the
%     same image resolutions.)
%   - 'spectral_weights': A 2D array, where `spectral_weights(i, j)` is the
%     sensitivity of the i-th colour channel or spectral band of 'I' to the
%     j-th colour channel or spectral band of the estimated latent image.
%     `spectral_weights` is a matrix mapping colours in the estimated
%     latent image to colours in the representation of 'I'.
%     `spectral_weights` must account for any numerical intergration that
%     is part of colour/spectral conversion.
%
% J -- Input RAW image
%   A 2D array containing the raw colour-filter pattern data of an image,
%   from which the latent image `I` is to be reconstructed.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% dispersionfun -- Model of dispersion
%   `dispersionfun` can be empty (`[]`), if there is no model of
%   dispersion. Otherwise, `dispersionfun` must be a function handle, such
%   as produced by 'makeDispersionfun()'. `dispersionfun(X)`, where `X` is
%   a three-element row vector (x, y, l), returns the dispersion vector for
%   the position (x, y) in `J` corresponding to light with wavelength or
%   colour channel index `l`. The dispersion vector points from the
%   corresponding position in the reference spectral band or colour channel
%   to position (x, y). This function will negate dispersion vectors in
%   order to create a warp matrix from `I` to `J`.
%
% sensitivity -- Colour space basis vectors
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th colour channel or spectral band of
%   `I`. `sensitivity` is a matrix mapping colours/spectra in `I` to
%   colours in `J`.
%
% lambda -- Spectral bands or colour channel indices
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersionfun`. 'c' is the desired number of spectral bands or colour
%   channels in `I`, and will be the size of `I` in its third dimension.
%
% admm_options -- Image estimation algorithm options
%   `admm_options` is a structure with the following fields, containing
%   settings for the ADMM-based image estimation algorithm:
%   - 'rho': A three or four-element vector containing penalty parameters
%     used in the ADMM framework. The first three elements correspond to
%     the regularization terms described in the documentation of
%     `reg_options` below. The fourth element is a penalty parameter for a
%     non-negativity constraint on the solution, and is only required if
%     the 'nonneg' field is `true`.
%   - 'full_GLambda': A Boolean value used as the `replicate` input
%     argument of 'spectralGradient()' when creating the spectral gradient
%     matrix for regularizing the spectral dimension of the latent image.
%     Refer to the documentation of 'spectralGradient.m' for details.
%     'full_GLambda' is not used if spectral regularization is disabled by
%     `reg_options.enabled`.
%   - 'maxit': A two-element vector. The first element contains the maximum
%     number of iterations to use with MATLAB's 'pcg()' function during the
%     I-minimization step of the ADMM algorithm. The second element of
%     `maxit` contains the maximum number of ADMM iterations to perform.
%   - 'norms': A three-element logical vector, corresponding to the
%     regularization terms listed in the documentation of `reg_options`
%     below. Each element specifies whether to use the L1 norm (`true`) or
%     an L2 norm (`false`) of the corresponding regularization penalty
%     vector. If some elements of 'norms' are `false`, the ADMM iterations
%     are simplified by eliminating slack variables. If all elements are
%     `false`, and 'nonneg' is `false`, then ADMM reduces to a
%     least-squares solution.
%   - 'nonneg': A Boolean scalar specifying whether or not to enable a
%     non-negativity constraint on the estimated image. If `true`, 'rho'
%     must have four elements.
%   - 'tol': A three-element vector containing convergence tolerances. The
%     first element is the tolerance value to use with MATLAB's 'pcg()'
%     function, such as when solving the I-minimization step of the ADMM
%     algorithm. The second and third elements are the absolute and
%     relative tolerance values for the ADMM algorithm, as explained in
%     Section 3.3.1 of Boyd et al. 2011.
%   - 'varying_penalty_params': If empty (`[]`), the penalty parameters
%     passed in 'rho' will be fixed across all ADMM iterations. Otherwise,
%     'varying_penalty_params' is a three-element vector containing the
%     parameters 'tau_incr', 'tau_decr', and 'mu', respectively, in
%     equation 3.13 of Boyd et al. 2011. In this case, the penalty
%     parameters in the ADMM iterations will vary, so as to speed up
%     convergence. Refer to Section 3.4.1 of Boyd et al. 2011 for a full
%     explanation.
%   - 'init': A character vector specifying how to initialize the latent
%     image `out.I`:
%     - 'zero': `out.I` will be a zero vector.
%     - 'uniform': `out.I` will be a pattern of uniform spectral
%       intensities that best fits the input image `J`.
%
% reg_options -- Regularization weight selection options
%   There are three regularization terms which can be enabled:
%   1 - Regularization of the spatial gradient of the image, as in Equation
%       6 of Baek et al. 2017.
%   2 - Regularization of the spectral gradient of the spatial gradient of
%       the image, as in Equation 6 in Baek et al. 2017.
%   3 - A second-order gradient prior designed to penalize colour-filter
%       array artifacts, implemented in antiMosaicMatrix.m.
%
%   `reg_options` is a structure with the following fields, controlling
%   regularization, and regularization weight selection:
%   - 'demosaic': A logical scalar indicating whether or not to select
%     regularization weights based on similarity with a demosaicking
%     result. Presently, the demosaicing result is the bilinear
%     interpolation of the Green channel (the second channel) of `J`. If
%     `I_in` is not empty, and 'demosaic' is `true`, regularization weights
%     will still be selected based on similarity with `I_in.I`. If `I_in`
%     is empty, and 'demosaic' is `false`, then regularization weights will
%     be selected using the minimum distance criterion.
%   - 'enabled': A logical vector, where each element indicates whether the
%     corresponding regularization term listed above is enabled.
%   - 'n_iter': A two-element vector, where the first element is the
%     maximum number of grid refinement iterations to perform in the method
%     of Song et al. 2016. The second is the minimum number of iterations
%     to perform (even if the process has converged, as specified by
%     'tol'). Note that an iterative grid search is used regardless of the
%     criterion being optimized (either the minimum distance criterion, or
%     similarity with an image).
%   - 'minimum_weights': A vector containing minimum values for the
%     regularization weights. If `I_in` is not empty or 'demosaic' is
%     `true`, 'minimum_weights' can contain any values thought to be
%     reasonable values for producing a good image. In contrast, if `I_in`
%     is empty, and 'demosaic' is `false`, the minimum distance criterion
%     will be used to select regularization weights, and the origin of the
%     minimum distance criterion will be set using 'minimum_weights'.
%     Therefore, in this case, 'minimum_weights' should contain the
%     smallest values that make the image estimation problem solvable.
%   - 'maximum_weights': A vector containing maximum values for the
%     regularization weights. 'maximum_weights' may be used to set the
%     origin of the minimum distance criterion, like 'minimum_weights', and
%     so needs to contain very large values if `I_in` is empty, and
%     'demosaic' is `false`. Otherwise, 'maximum_weights' can contain any
%     values thought to be reasonable values for producing a good image.
%   - 'low_guess': A vector containing predicted lower bounds for the
%     regularization weights. The search algorithm will examine
%     regularization weights smaller than these values if necessary (as
%     small as 'minimum_weights'), but otherwise, 'low_guess' may reduce
%     the number of iterations required for convergence.
%   - 'high-guess': A vector containing predicted upper bounds for the
%     regularization weights, used in the same way as 'low_guess'.
%   - 'tol': The threshold value of the relative change in the objective
%     criterion for regularization weights selection from one iteration to
%     the next. When the change is less than this threshold, and the
%     minimum number of iterations has been reached, iteration terminates.
%
%   If 'minimum_weights' and 'maximum_weights' are identical,
%   regularization weight selection is disabled, and these values are used
%   as the final regularization weights.
%
% patch_options -- Options for patch-wise image estimation
%   A structure containing the following fields:
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of image patches to be estimated individually.
%     `patch_size` does not include padding used to eliminate patch border
%     artifacts. Patches along the bottom and right edges of the image may
%     be made smaller to fit within the image's borders.
%   - 'padding': A scalar containing the pixel width of the border
%     surrounding each image patch. The image patches actually estimated
%     are of size `patch_size + padding`, but a border of width 'padding'
%     is stripped when combining patches to form the output image. Note
%     that patches along the edges of the image are not padded to extend
%     outside the image's borders, and so will only have padding towards
%     the interior of the image. 'padding' helps mitigate artifacts from
%     patch-wise image estimation, and should be at least as large as the
%     amount of dispersion in the image formation model.
%   - 'target_patch': An optional field. If it exists, 'target_patch' is a
%     two-element vector containing the row and column, respectively, of
%     the top-left corner of the single image patch to be estimated. For
%     simplicity, so that it is not necessary to handle edge cases, wherein
%     the patch may not be a valid colour-filter array image, the elements
%     of 'target_patch' must be odd integers. When `target_patch` is
%     passed, all output arguments are calculated for a single image patch,
%     rather than for the entire image. While a border around the patch
%     will have been estimated, with a width given by 'padding', it will
%     not be included in the output. Prior to its removal, the border
%     region will be used to calculate output images aside from `I`, to
%     limit artifacts from image warping.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show progress
%   information. This flag is also passed internally to
%   'weightsLowMemory()'.
%
% ## Output Arguments
%
% I -- Latent image
%   A size(J, 1) x size(J, 2) x length(lambda) array, storing the estimated
%   latent image corresponding to `J`.
%
% I_rgb -- Latent colour image
%   The colour equivalent of the latent image, generated using the colour
%   space conversion data in `sensitivity`. An size(J, 1) x size(J, 2) x
%   size(sensitivity, 1) array.
%
% weights_images -- Selected regularization weights
%   A size(J, 1) x size(J, 2) x w array, where 'w' is the number of `true`
%   values in `reg_options.enabled`. `weights_images(p, q, i)` contains the
%   weight, on the regularization term corresponding to the i-th `true`
%   value in `reg_options.enabled`, used when estimating the pixel at
%   location (p, q).
%
%   If regularization weights are fixed, i.e. if
%   `all(reg_options.minimum_weights == reg_options.maximum_weights)` is
%   `true`, then `weights_images` is empty (`[]`);
%
% J_full -- Warped latent colour image
%   A colour image produced by warping `I` according to the dispersion
%   model, followed by conversion to the colour space of `J`. An size(J, 1)
%   x size(J, 2) x size(sensitivity, 1) array.
%
% J_est -- Re-estimated input RAW image
%   The mosaiced version of `J_full`, representing the result of passing
%   `I` through the forward model of dispersion and image capture. An array
%   with the same dimensions as `J`.
%
% I_warped -- Warped latent image
%   An size(J, 1) x size(J, 2) x length(lambda) array, storing the latent
%   image warped according to the dispersion model.
%
% search -- Grid search method for regularization weight selection search path
%   The `search` output argument of 'weightsLowMemory()'. `search` is
%   specific to a single image patch, and therefore can only be output when
%   `patch_options.target_patch` exists. When `I_in` is empty, `search`
%   will have the additional fields that it has in the case where
%   `in_penalties` is an input argument of 'weightsLowMemory()'.
%
% ## Notes
%
% - If all elements of `reg_options.enabled` are `false`, and
%   `admm_options.nonneg` is `false`, then this function will find a
%   solution to the problem:
%     argmin_i (||M * Omega * Phi * i - j||_2) ^ 2
%   where `M` performs mosaicing, `Omega` converts colours to the colour
%   space of `J`, and `Phi` warps the image according to the dispersion
%   model. `i` and `j` are the vectorized versions of the latent and input
%   images, respectively. If the linear system is underdetermined, the
%   function will find the minimum-norm least squares solution. If the
%   problem is sufficiently determined, as it may be in cases where `i` has
%   fewer spectral bands or colour channels than `j` has colour channels,
%   then the function will find the iterative approximation to the exact
%   solution (using MATLAB's 'pcg()' function), or will find a
%   least-squares solution (in the overdetermined case).
% - In contrast to 'selectWeightsGrid()', the origin of the minimum
%   distance criterion is always set using `reg_options.minimum_weights`
%   and `reg_options.maximum_weights`, rather than chosen
%   semi-automatically.
% - This function imitates the behaviour of 'selectWeightsGrid()' with
%   `'normalized'` as the value of its `options.scaling` argument.
% - During regularization weight selection, image patches are not stripped
%   of a border before calculating their data-fitting errors,
%   regularization errors, or similarity with the true image, for
%   simplicity. (Contrast with the behaviour of 'trainWeights()' and
%   'selectWeightsGrid()' in combination with 'baek2017Algorithm2()'.)
% - The estimated image is spatially-registered with the input image,
%   conforming to the behaviour of 'baek2017Algorithm2()' with its
%   `options.add_border` input argument set to `false`. This simplifies
%   patch-wise image estimation, as explained in the documentation of
%   'solvePatchesAligned()'.
% - `patch_options.patch_size`, `patch_options.padding`, and the dimensions
%   of `I_in.I` and `J` must all be even integers for the images and image
%   patches to map to valid Bayer colour-filter array patterns.
%
% ## References
%
% This function implements Algorithm 2 in the first set of supplemental
% material of the following article:
%
%   Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%     single-shot hyperspectral imaging using a prism." ACM Transactions
%     on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%     doi:10.1145/3130800.3130896
%
% Depending on the options passed, this function implements variants of the
% algorithm: L2-norm priors instead of L1-norm priors, an extra prior
% designed to remove colour-filter array artifacts, and a non-negativity
% constraint. I implemented the non-negativity constraint by adding an
% extra term to the ADMM x-minimization step, and an additional
% z-minimization and dual update step. This is different from the
% constrained optimization examples in Boyd et. al. 2011, sections 4.2.5
% and 5.2, but I think it matches the format given at the start of Chapter
% 5.
%
% A non-negativity constraint was used in (among other works):
%
%   Park, J.-I., Lee, M.-H., Grossberg, M. D., & Nayar, S. K. (2007).
%     "Multispectral Imaging Using Multiplexed Illumination." In 2007 IEEE
%     International Conference on Computer Vision (ICCV).
%     doi:10.1109/ICCV.2007.4409090
%
% The method for initializing the latent image is based on Equation 6 of
% the following article. Note that this article also briefly discusses
% methods for combining overlapping image patches, for future reference.
% This function simply discards the overlapping regions of image patches,
% which may lead to visible boundaries between patches, but which prevents
% any poor behaviour of the image estimation algorithm near the edges of
% patches from affecting the result.
%
%   Sun, T., Peng, Y., & Heidrich, W. (2017). "Revisiting cross-channel
%     information transfer for chromatic aberration correction." In 2017
%     IEEE International Conference on Computer Vision (ICCV) (pp.
%     3268–3276). doi:10.1109/ICCV.2017.352
%
% For more information on ADMM (Alternating Direction Method of
% Multipliers), read:
%
%   Boyd, S, et al.. "Distributed Optimization and Statistical Learning via
%     the Alternating Direction Method of Multipliers." Foundations and
%     Trends in Machine Learning, vol. 3, no. 1, pp. 1-122, 2011.
%     doi:10.1561/2200000016
%
% The following article discusses the grid-search method for minimizing the
% minimum distance criterion used to select regularization weights. To
% select regularization weights which maximize the similarity of the
% estimated and true/demosaiced images, I also use the grid-search method.
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
%
% ## Future Work
%
% There are several modifications and expansions which may improve the
% performance of ADMM:
% - Section 3.4.1 of Boyd et al. 2011, "Varying Penalty Parameter"
%   - Now activated by a non-empty `admm_options.varying_penalty_params`
%     vector.
%   - Further refinement may be possible by "taking into account
%     the relative magnitudes of [the primal and dual convergence
%     thresholds]."
% - Section 4.3.2 of Boyd et al. 2011, "Early Termination"
%
% See also baek2017Algorithm2LowMemory, weightsLowMemory,
% solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 5, 2018

narginchk(9, 10);
nargoutchk(1, 7);

verbose = false;
if ~isempty(varargin)
    verbose = varargin{1};
end

if verbose
    tic
end

% Input argument parsing

output_search = nargout > 6;
do_single_patch = isfield(patch_options, 'target_patch');
if output_search && ~do_single_patch
    error('`search` can only be output when `patch_options.target_patch` exists.');
end

has_dispersion = ~isempty(dispersionfun);
if has_dispersion && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

nargout_before_weights = 2;
output_weights = (nargout > nargout_before_weights) && ~all(reg_options.minimum_weights == reg_options.maximum_weights);
input_I_in = ~isempty(I_in);
if nargout > nargout_before_weights
    n_auxiliary_images = nargout - 2;
else
    n_auxiliary_images = 1;
end
image_sampling = [size(J_2D, 1), size(J_2D, 2)];
patch_size = patch_options.patch_size;
padding = patch_options.padding;
n_bands = length(lambda);
n_channels_rgb = size(sensitivity, 1);
enabled_weights = reg_options.enabled;
n_active_weights = sum(enabled_weights);
if all(~enabled_weights) && output_weights
    error('Cannot output selected regularization weights, because all regularization terms are disabled.');
end
if enabled_weights(2) && n_bands < 2
    error('Cannot enable spectral regularization when `lambda` has length one.');
end
use_min_norm = all(~enabled_weights) && ~admm_options.nonneg;

if any(mod(image_sampling, 2) ~= 0)
    error('The input image `J` must have dimensions which are even integers, or it cannot represent a colour-filter array image.');
end
if any(mod(patch_size, 2) ~= 0)
    error('`patch_options.patch_size` must be even integers to produce patches which are valid colour-filter array images.');
end
if mod(padding, 2) ~= 0
    error('`patch_options.padding` must be an even integer to produce patches which are valid colour-filter array images.');
end
if input_I_in
    if any(image_sampling ~= [size(I_in.I, 1), size(I_in.I, 2)])
        error('The spatial dimensions of `I_in.I` must match those of `J`.')
    end
    if n_bands ~= size(I_in.spectral_weights, 2)
        error('The number of wavelengths in `lambda` must equal the size of `I_in.spectral_weights` in its second dimension.');
    end
    if size(I_in.spectral_weights, 1) ~= size(I_in.I, 3)
        error('The number of rows of `I_in.spectral_weights` must equal the size of `I_in.I` in its third dimension.');
    end
    I_in_j.spectral_weights = I_in.spectral_weights;
elseif reg_options.demosaic
    I_in_j.spectral_weights = sensitivity(2, :);
else
    I_in_j = struct;
end

if do_single_patch
    target_patch = patch_options.target_patch;
    if target_patch(1) < 1 || target_patch(1) > image_sampling(1) ||...
       target_patch(2) < 1 || target_patch(2) > image_sampling(2)
        error('The target patch is outside of the image bounds.');
    end
    if any(mod(target_patch, 2) ~= 1)
        error(['To prevent having an odd-sized patch because of clippin',...
            'g at the image boundaries, the target patch corner must hav',...
            'e odd integer coordinates.']);
    end
end

if n_bands ~= size(sensitivity, 2)
    error('The number of wavelengths in `lambda` must equal the size of `sensitivity` in its second dimension.');
end

if verbose
    disp('Splitting the input image into columns...');
end

% Channel indices in the input concatenation of images
n_channels_in = 0;
if input_I_in
    channels_in.I_in = [ 1, size(I_in.I, 3) ]; % True image
    n_channels_in = n_channels_in + channels_in.I_in(2);
end
channels_in.J = n_channels_in + [1, size(J_2D, 3)]; % Input image
n_channels_in = n_channels_in + channels_in.J(2);

% Channel indices in the output concatenation of images
n_channels_out = 0;
channels_out.I = n_channels_out + [1, n_bands]; % Estimated latent image
n_channels_out = n_channels_out + channels_out.I(2);
if n_auxiliary_images > 0
    channels_out.I_rgb = n_channels_out + [1, n_channels_rgb];
    n_channels_out = n_channels_out + channels_out.I_rgb(2);
    if output_weights
        channels_out.I_weights = n_channels_out + [1, n_active_weights];
        n_channels_out = n_channels_out + channels_out.I_weights(2);
    end
    if n_auxiliary_images > 1
        channels_out.J_full = n_channels_out + [1, n_channels_rgb];
        n_channels_out = n_channels_out + channels_out.J_full(2);
        if n_auxiliary_images > 2
            channels_out.J_est = n_channels_out + [1, size(J_2D, 3)];
            n_channels_out = n_channels_out + channels_out.J_est(2);
            if n_auxiliary_images > 3
                channels_out.I_warped = n_channels_out + [1, n_bands];
                n_channels_out = n_channels_out + channels_out.I_warped(2);
            end
        end
    end
end

% Divide the input images into columns which will be sent to individual
% parallel workers
if do_single_patch
    n_i = 1;
    n_j = 1;
    patch_offset = target_patch - 1;
else
    n_i = ceil(image_sampling(1) / patch_size(1));
    n_j = ceil(image_sampling(2) / patch_size(2));
    patch_offset = [0, 0];
end
n_patches = n_i * n_j;
columns_in = cell(1, n_j);
for j = 1:n_j
    cols_ind_in = [
        max((j - 1) * patch_size(2) + 1 - padding + patch_offset(2), 1);
        min(j * patch_size(2) + padding + patch_offset(2), image_sampling(2))
    ];
    columns_in{j} = zeros(image_sampling(1), diff(cols_ind_in) + 1, n_channels_in);
    if input_I_in
        columns_in{j}(:, :, channels_in.I_in(1):channels_in.I_in(2)) = I_in.I(:, cols_ind_in(1):cols_ind_in(2), :);
    end
    columns_in{j}(:, :, channels_in.J(1):channels_in.J(2)) = J_2D(:, cols_ind_in(1):cols_ind_in(2), :);
end

if verbose
    fprintf('\tDone.\n');
    disp('Parallel processing of columns...');
end

% Process each column
columns_out = cell(1, n_j);
if output_search
    search_out = cell(1, n_j);
end
parfor j = 1:n_j
    
    column_in_j = columns_in{j};
    image_sampling_p = [0, size(column_in_j, 2)];
    corner = [0, (j - 1) * patch_size(2) + 1 + patch_offset(2)];
    cols_trim_out = [ min(padding + 1, corner(2)), 0 ];
    cols_trim_out(2) = cols_trim_out(1) + min(patch_size(2) - 1, image_sampling(2) - corner(2));
    col_out_width = diff(cols_trim_out) + 1;
    column_out_j = zeros(size(column_in_j, 1), col_out_width, n_channels_out);
    
    % Process each patch within the column
    for i = 1:n_i
        corner(1) = (i - 1) * patch_size(1) + 1 + patch_offset(1);
        rows_trim_out = [ min(padding + 1, corner(1)), 0 ];
        rows_trim_out(2) = rows_trim_out(1) + min(patch_size(1) - 1, image_sampling(1) - corner(1));
        patch_lim_rows = [
            max(corner(1) - padding, 1);
            min(corner(1) + patch_size(1) + padding - 1, image_sampling(1));
        ];
        patch_end_row = min(corner(1) + patch_size(1) - 1, image_sampling(1));
        image_sampling_p(1) = diff(patch_lim_rows) + 1;
        numel_p = prod(image_sampling_p) * n_bands;
        
        if has_dispersion
            dispersion_matrix_p = dispersionfunToMatrix(...
                dispersionfun, lambda, image_sampling_p, image_sampling_p,...
                [0, 0, image_sampling_p(2), image_sampling_p(1)], true,...
                flip(corner) - 1 ...
            );
        else
            dispersion_matrix_p = [];
        end

        % Solve for the output patch
        in_admm = initBaek2017Algorithm2LowMemory(...
            column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.J(1):channels_in.J(2)),...
            align, dispersion_matrix_p,...
            sensitivity, enabled_weights, admm_options...
        );
        if use_min_norm
            % Minimum-norm least squares solution to the non-regularized problem
            if(verbose)
                fprintf('Computing the minimum-norm least squares solution...\n');
            end
            in_admm.J = reshape(column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.J(1):channels_in.J(2)), [], 1);
            patches_I_ij = lsqminnorm(in_admm.M_Omega_Phi, in_admm.J);
            if(verbose)
                fprintf('\t...done.\n');
            end
                        
        elseif input_I_in || reg_options.demosaic
            I_in_p = I_in_j;
            if input_I_in
                dispersion_matrix_p_weights = [];
                I_in_p.I = column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.I_in(1):channels_in.I_in(2));
            else
                dispersion_matrix_p_weights = dispersion_matrix_p;
                I_in_p.I = reshape(bilinearDemosaic(...
                    column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.J(1):channels_in.J(2)),...
                    align, [false, true, false]...
                ), [], 1);
            end
            in_weightsLowMemory = initWeightsLowMemory(I_in_p, dispersion_matrix_p_weights, numel_p);
            if output_search
                [...
                    patches_I_ij, weights, ~, ~, search_out{j}...
                ] = weightsLowMemory(...
                    admm_options, reg_options, in_weightsLowMemory,...
                    in_admm, I_in_p, verbose...
                );
            else
                [...
                    patches_I_ij, weights...
                ] = weightsLowMemory(...
                    admm_options, reg_options, in_weightsLowMemory,...
                    in_admm, I_in_p, verbose...
                );
            end

        else
            in_penalties = initPenalties(in_admm.M_Omega_Phi, in_admm.G);
            in_weightsLowMemory = initWeightsLowMemory([], [], numel_p);
            if output_search
                [...
                    patches_I_ij, weights, ~, ~, ~, search_out{j}...
                ] = weightsLowMemory(...
                    admm_options, reg_options, in_weightsLowMemory,...
                    in_admm, in_penalties, verbose...
                );
            else
                [...
                    patches_I_ij, weights...
                ] = weightsLowMemory(...
                    admm_options, reg_options, in_weightsLowMemory,...
                    in_admm, in_penalties, verbose...
                );
            end
        end
        
        patches_I_ij_3D = reshape(patches_I_ij, [image_sampling_p n_bands]);
        column_out_j(...
            corner(1):patch_end_row, :, channels_out.I(1):channels_out.I(2)...
        ) = patches_I_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
    
        if n_auxiliary_images > 0
            patches_I_rgb_ij = in_admm.Omega * patches_I_ij;
            patches_I_rgb_ij_3D = reshape(patches_I_rgb_ij, [image_sampling_p n_channels_rgb]);
            column_out_j(...
                corner(1):patch_end_row, :, channels_out.I_rgb(1):channels_out.I_rgb(2)...
            ) = patches_I_rgb_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
        
            if n_auxiliary_images > 1
                if has_dispersion
                    patches_I_warped_ij = dispersion_matrix_p * patches_I_ij;
                    patches_J_full_ij = in_admm.Omega * patches_I_warped_ij;
                    patches_J_full_ij_3D = reshape(patches_J_full_ij, [image_sampling_p n_channels_rgb]);
                    column_out_j(...
                        corner(1):patch_end_row, :, channels_out.J_full(1):channels_out.J_full(2)...
                    ) = patches_J_full_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                else
                    patches_I_warped_ij = patches_I_ij_3D;
                    patches_J_full_ij = patches_I_rgb_ij;
                    column_out_j(...
                        corner(1):patch_end_row, :, channels_out.J_full(1):channels_out.J_full(2)...
                    ) = column_out_j(...
                        corner(1):patch_end_row, :, channels_out.I_rgb(1):channels_out.I_rgb(2)...
                    );
                end
                
                if n_auxiliary_images > 2
                    patches_J_est_ij_3D = reshape(in_admm.M * patches_J_full_ij, [image_sampling_p size(J_2D, 3)]);
                    column_out_j(...
                        corner(1):patch_end_row, :, channels_out.J_est(1):channels_out.J_est(2)...
                    ) = patches_J_est_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                    
                    if n_auxiliary_images > 3
                        patches_I_warped_ij_3D = reshape(patches_I_warped_ij, [image_sampling_p n_bands]);
                        column_out_j(...
                            corner(1):patch_end_row, :, channels_out.I_warped(1):channels_out.I_warped(2)...
                        ) = patches_I_warped_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                    end
                end
            end
        end
        
        if output_weights
            column_out_j(...
                corner(1):patch_end_row, :, channels_out.I_weights(1):channels_out.I_weights(2)...
            ) = repmat(reshape(weights(enabled_weights), 1, 1, n_active_weights), diff(rows_trim_out) + 1, col_out_width, 1);
        end
        
        if verbose
            fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i, n_patches);
        end
    end
    columns_out{j} = column_out_j;
end

if verbose
    fprintf('\tDone.\n');
    disp('Recombining results from patches...');
end

% Recombine patches
if do_single_patch
    images_out = columns_out{1}(target_patch(1):min(target_patch(1) + patch_size(1) - 1, image_sampling(1)), :, :);
else
    images_out = zeros(image_sampling(1), image_sampling(2), n_channels_out);
    offset = 0;
    for j = 1:n_j
        width_j = size(columns_out{j}, 2);
        images_out(:, (offset + 1):(offset + width_j), :) = columns_out{j};
        offset = offset + width_j;
    end
end
I_3D = images_out(:, :, channels_out.I(1):channels_out.I(2));

if n_auxiliary_images > 0
    varargout = cell(1, nargout - 1);
    varargout{1} = images_out(:, :, channels_out.I_rgb(1):channels_out.I_rgb(2));
    if n_auxiliary_images > 1
        varargout{3} = images_out(:, :, channels_out.J_full(1):channels_out.J_full(2));
        if n_auxiliary_images > 2
            varargout{4} = images_out(:, :, channels_out.J_est(1):channels_out.J_est(2));
            if n_auxiliary_images > 3
                varargout{5} = images_out(:, :, channels_out.I_warped(1):channels_out.I_warped(2));
            end
        end
    end
end
if output_weights
    varargout{2} = images_out(:, :, channels_out.I_weights(1):channels_out.I_weights(2));
elseif nargout > nargout_before_weights
    varargout{2} = [];
end
if output_search
    varargout(end) = search_out(1);
end

if verbose
    fprintf('\tDone.\n');
    toc
end

end