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
% [ I, I_rgb, J_full ] = solvePatchesADMM(___)
% [ I, I_rgb, J_full, J_est ] = solvePatchesADMM(___)
% [ I, I_rgb, J_full, J_est, I_warped ] = solvePatchesADMM(___)
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
% [ I, I_rgb, J_full ] = solvePatchesADMM(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped according to the model of dispersion.
%
% [ I, I_rgb, J_full, J_est ] = solvePatchesADMM(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% [ I, I_rgb, J_full, J_est, I_warped ] = solvePatchesADMM(___)
%   Additionally returns the version of the latent image warped according
%   to the model of dispersion.
%
% ## Input Arguments
%
% I_in -- True image
%   A 2D or 3D array containing the ground truth latent image, against
%   which the estimated latent image will be compared. When `I_in` is not
%   empty, the regularization weights used in image estimation will be
%   selected to optimize the similarity of the estimated and true images.
%   When `I_in` is empty (`[]`), regularization weights will be selected
%   using the minimum distance criterion described in Song et al. 2016.
%
%   `I_in` and `J` must have the same sizes in their first two dimensions.
%   (i.e. They must have the same image resolutions.)
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
%   Note that the `size(I_in, 3)` must equal 'c' if `I_in` is not empty.
%
% admm_options -- Image estimation algorithm options
%   `admm_options` is a structure with the following fields, containing
%   settings for the ADMM-based image estimation algorithm:
%   - 'rho': A three or four-element vector containing penalty parameters
%     used in the ADMM framework. The first three elements correspond to
%     the regularization terms described in the documentation of
%     `reg_options` below. The fourth element is a penalty parameter for a
%     non-negativity constraint on the solution, and is only required if
%     `nonneg` is `true`.
%   - 'full_GLambda': A Boolean value used as the `replicate` input
%     argument of 'spectralGradient()' when creating the spectral gradient
%     matrix for regularizing the spectral dimension of the latent image.
%     Refer to the documentation of 'spectralGradient.m' for details.
%     'full_GLambda' is not used if spectral regularization is disabled by
%     `reg_options.enabled`.
%   - 'int_method': The numerical integration method used for spectral to
%     colour space conversion. `int_method` is passed to
%     'channelConversionMatrix()' as its `int_method` argument. Refer to
%     the documentation of 'channelConversionMatrix.m' for details. If
%     'int_method' is 'none', numerical integration will not be performed.
%     'int_method' should be 'none' when `I` contains colour channels as
%     opposed to spectral bands.
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
%   - 'enabled': A logical vector, where each element indicates whether the
%     corresponding regularization term listed above is enabled.
%   - 'n_iter': A two-element vector, where the first element is the
%     maximum number of grid refinement iterations to perform in the method
%     of Song et al. 2016. The second is the minimum number of iterations
%     to perform (even if the process has converged, as specified by
%     'tol'). Note that an iterative grid search is used regardless of the
%     criterion being optimized (either the minimum distance criterion, or
%     similarity with the true image).
%   - 'minimum_weights': A vector containing minimum values for the
%     regularization weights. If `I_in` is not empty, 'minimum_weights' can
%     contain any values thought to be reasonable values for producing a
%     good image. In contrast, if `I_in` is empty, the minimum distance
%     criterion will be used to select regularization weights, and the
%     origin of the minimum distance criterion will be set using
%     'minimum_weights'. Therefore, in this case, 'minimum_weights' should
%     contain the smallest values that make the image estimation problem
%     solvable.
%   - 'maximum_weights': A vector containing maximum values for the
%     regularization weights. 'maximum_weights' is also used to set the
%     origin of the minimum distance criterion, like 'minimum_weights', and
%     so needs to contain very large values if `I_in` is empty. If `I_in`
%     is not empty, 'maximum_weights' can contain any values thought to be
%     reasonable values for producing a good image.
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
%     artifacts. Patches along the bottom and right edges of the image will
%     be made smaller to fit within the image's borders.
%   - 'padding': A scalar containing the pixel width of the border
%     surrounding each image patch. The image patches actually estimated
%     are of size `patch_size + padding`, but a border of width 'padding'
%     is stripped when combining patches to form the output image. Note
%     that patches along the edges of the image are not padded to extend
%     outside the image's borders, and so will only have padding towards
%     the interior of the image. 'padding' helps mitigate artifacts from
%     patch-wise image estimation, and should be large enough to
%     accommodate for dispersion.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show progress
%   information.
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
% J_full -- Warped latent RGB image
%   An RGB image produced by warping `I` according to the dispersion model,
%   followed by conversion to the RGB colour space of `J`. An size(J, 1) x
%   size(J, 2) x 3 array.
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
% ## Notes
%
% - If all elements of `reg_options.enabled` are `false`, and
%   `admm_options.nonneg` is `false`, then this function will find a
%   solution to the problem:
%     argmin_i (||M * Omega * Phi * i - j||_2) ^ 2
%   where `M` performs mosaicing, `Omega` converts colours to the colour
%   space of `J`, and `Phi` warps the image according to the dispersion
%   model. `i` and `j` are the vecctorized versions of the latent and input
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
%   simplicity.
% - The estimated image is spatially-registered with the input image,
%   conforming to the behaviour of 'baek2017Algorithm2()' with its
%   `options.add_border` input argument set to `false`. This simplifies
%   patch-wise image estimation, as explained in the documentation of
%   'solvePatchesAligned()'.
% - `patch_options.patch_size`, `patch_options.padding`, and the dimensions
%   of `I_in` and `J` must all be even integers for the images and image
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
% estimated and true images, I also use the grid-search method.
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
% See also baek2017Algorithm2, selectWeightsGrid, trainWeights,
% solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 5, 2018

narginchk(9, 10);
nargoutchk(1, 5);

end