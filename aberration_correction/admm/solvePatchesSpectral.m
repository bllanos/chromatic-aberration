function [ bands, I_3D, varargout ] = solvePatchesSpectral(...
    I_in, J_2D, align, dispersionfun, color_map, color_bands, sampling_options,...
    admm_options, reg_options, patch_options, varargin...
)
% SOLVEPATCHESSPECTRAL  Run ADMM for spectral image estimation
%
% ## Usage
%
% This function is the equivalent of 'solvePatchesColor()' for spectral image
% estimation. It assumes that the image can be estimated at arbitrary spectral
% resolutions. As such, it will choose the appropriate spectral resolution, and
% can estimate the image at successively higher spectral resolutions if desired.
% 'solvePatchesColor()' is appropriate when the output image is composed of
% colour channels, which are categorical, not sampled, representations of
% spectral information.
%
% ## Syntax
% [ bands, I ] = solvePatchesSpectral(...
%   I_in, J, align, dispersionfun, color_map, color_bands, sampling_options,...
%   admm_options, reg_options, patch_options [, verbose]...
% )
% [ bands, I, I_rgb ] = solvePatchesSpectral(___)
% [ bands, I, I_rgb, weights_images ] = solvePatchesSpectral(___)
% [ bands, I, I_rgb, weights_images, J_full ] = solvePatchesSpectral(___)
% [ bands, I, I_rgb, weights_images, J_full, J_est ] = solvePatchesSpectral(___)
% [ bands, I, I_rgb, weights_images, J_full, J_est, I_warped ] = solvePatchesSpectral(___)
% [ bands, I, I_rgb, weights_images, J_full, J_est, I_warped, search ] = solvePatchesSpectral(...)
%
% ## Description
% [ bands, I ] = solvePatchesSpectral(...
%   I_in, J, align, dispersionfun, color_map, color_bands, sampling_options,...
%   admm_options, reg_options, patch_options [, verbose]...
% )
%   Estimate a latent RGB or hyperspectral image `I` from dispersion in
%   the input RAW image `J`, and also return the wavelengths corresponding
%   to the spectral bands of `I`.
%
% [ bands, I, I_rgb ] = solvePatchesSpectral(___)
%   Additionally returns the RGB equivalent of the latent image.
%
% [ bands, I, I_rgb, weights_images ] = solvePatchesSpectral(___)
%   Additionally returns images illustrating the weights selected for the
%   regularization terms in the optimization problem.
%
% [ bands, I, I_rgb, weights_images, J_full ] = solvePatchesSpectral(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped according to the model of dispersion.
%
% [ bands, I, I_rgb, weights_images, J_full, J_est ] = solvePatchesSpectral(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% [ bands, I, I_rgb, weights_images, J_full, J_est, I_warped ] = solvePatchesSpectral(___)
%   Additionally returns the version of the latent image warped according
%   to the model of dispersion.
%
% [ bands, I, I_rgb, weights_images, J_full, J_est, I_warped, search ] = solvePatchesSpectral(...)
%   Additionally returns the search path taken to select regularization
%   weights for a single image patch, at the highest spectral resolution.
%   This call syntax is available only when `patch_options.target_patch`
%   exists.
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
%   - 'spectral_bands': A vector containing the wavelengths corresponding
%     to the third dimension of 'I'.
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
%   `dispersionfun` can be empty (`[]`), if there is no model of dispersion.
%   Otherwise, `dispersionfun` must be a function handle, such as produced by
%   'makeDispersionfun()'. `dispersionfun(X)`, where `X` is a three-element row
%   vector (x, y, l), returns the dispersion vector for the position (x, y) in
%   `J` corresponding to light with wavelength `l`. The dispersion vector points
%   from the corresponding position in the reference spectral band to position
%   (x, y). This function will negate dispersion vectors in order to create a
%   warp matrix from `I` to `J`.
%
% color_map -- Colour channel spectral sensitivities
%   A 2D array, where `color_map(i, j)` is the sensitivity of the i-th colour
%   channel of `J` to the j-th spectral band in `color_bands`. `color_map` is
%   not a colour conversion matrix, as it does not perform the desired numerical
%   integration, over the spectrum, that is part of colour conversion.
%
% color_bands -- Wavelength bands for colour channel sensitivities
%   A vector, of length equal to the size of the second dimension of
%   `color_map`, containing the wavelengths at which the sensitivity
%   functions in `color_map` have been sampled. `color_bands(j)` is the
%   wavelength corresponding to `color_map(:, j)`. The values in
%   `color_bands` are expected to be evenly-spaced.
%
% sampling_options -- Spectral sampling options
%   `sampling_options` is a structure with the following fields used by
%   'findSampling()' and 'dispersionfunToMatrix()':
%   - 'resolution': A non-negative scalar providing the desired approximate
%     spacing, in pixels, between the images for consecutive wavelengths at
%     which the dispersion is to be sampled. If 'resolution' is zero or is
%     missing, the dispersion function will be evaluated at the same spectral
%     sampling as the estimated image.
%   - 'int_method': The numerical integration method to use when
%     integrating over the responses of colour channels to compute colour
%     values. `int_method` is used by `colorWeights()`.
%   - 'power_threshold': An option used by 'findSampling()' to select
%     spectral sampling points.
%   - 'n_bands': An option used by 'findSampling()' to select
%     spectral sampling points.
%   - 'support_threshold': A fraction indicating what proportion of the
%     peak magnitude of a colour channel's sensitivity function the user
%     considers to be effectively zero (i.e. no sensitivity).
%   - 'bands_padding': An option used by 'findSampling()' to control how
%     spectral signals are extrapolated.
%   - 'interpolant': The function convolved with spectral signals to interpolate
%     them from the sampling space of the estimated image to another sampling
%     space. 'interpolant' is passed to 'resamplingWeights()' as its `f` input
%     argument. Refer to the documentation of 'resamplingWeights.m' for more
%     details.
%   - 'interpolant_ref': Similar to 'interpolant', but for interpolating between
%     the spectral sampling spaces of quantities other than the estimated image.
%
%   The following fields are in addition to those used by
%   'findSampling()' and 'dispersionfunToMatrix()':
%   - 'progression': A character vector describing how to select the
%     intermediate numbers of bands at which to estimate the latent image. The
%     final number of bands, denoted by 'S' below, is determined by calling
%     'findSampling()' with the options given above.
%     - 'sequential': Use the following numbers of wavelengths: 3, 4, 5, ..., S.
%       If 'S' is less than 3, then the sequence is just 'S'.
%     - 'doubling': Use the following numbers of wavelengths: 4, 8, ...,
%       2 ^ floor(log_2(S)), S. (If 'S' is equal to `2 ^ floor(log_2(S))`, then
%       the sequence ends with `2 ^ floor(log_2(S))`.) If 'S' is less than 4,
%       then the sequence is just 'S'.
%     - 'last': Estimate the latent image at 'S' bands directly.
%   - 'show_steps': A logical scalar which, if `true`, causes the function
%     to return the images for all intermediate numbers of wavelengths, not
%     only for 'S' wavelengths. When 'show_steps' is `true`:
%     - `bands` and `search` are cell vectors, where the i-th cell contains
%       the data for the i-th number of wavelengths.
%     - `I`, `I_rgb`, `weights_images`, `J_full`, `J_est`, and `I_warped`
%       are extended in the third dimension by stacking results for
%       individual numbers of wavelengths.
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
%   - 'tol': A two-element vector containing convergence tolerances. The
%     first element is the tolerance value to use with MATLAB's 'pcg()'
%     function, such as when solving the I-minimization step of the ADMM
%     algorithm. The second element is the relative tolerance for the ADMM
%     algorithm, as explained in Section 3.3.1 of Boyd et al. 2011.
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
%   3 - Regularization of the spatial Laplacian of the image, as in Song et al.
%       2016.
%
%   `reg_options` is a structure with the following fields, controlling
%   regularization, and regularization weight selection:
%   - 'demosaic': A logical scalar indicating whether or not to select
%     regularization weights based on similarity with a demosaicking result.
%     Presently, the demosaicing result is the bilinearly interpolated version
%     of the colour channels of `J`. If `I_in` is not empty, and 'demosaic' is
%     `true`, regularization weights will still be selected based on similarity
%     with `I_in.I`. If `I_in` is empty, and 'demosaic' is `false`, then
%     regularization weights will be selected using the minimum distance
%     criterion.
%   - 'demosaic_channels': In the context where 'demosaic' is `true` (see
%     above), a logical vector indicating which channels of the demosaicing
%     result should be used to evaluate similarity.
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
%   With 'solvePatchesColor()', if 'minimum_weights' and 'maximum_weights'
%   were identical, regularization weight selection was disabled, and these
%   values were used as the final regularization weights.
%
%   With this function, it is not recommended to make 'minimum_weights' and
%   'maximum_weights' identical, as it will fix the regularization weights
%   for all spectral resolutions to the same values. Instead, `reg_options`
%   can be given an additional field:
%   - 'multi_weights': An optional field, containing a matrix with the same
%     number of columns as the length of 'enabled'. 'multi_weights' must
%     have the same number of rows as the length of the sequence of
%     spectral resolutions used in image estimation. If present,
%     'multi_weights' is used to fix regularization weights specific to
%     each spectral resolution. It takes precedence over 'minimum_weights'
%     and 'maximum_weights'.
%
%     As it may be hard to know the number of spectral resolutions in
%     advance, the recommended approach is to call this function on one or
%     more target patches (see `patch_options.target_patch` below) and set
%     `sampling_options.show_steps` to `true`. The `weights_images` output
%     argument for these calls will then give both the number of spectral
%     resolutions, and the regularization weights selected for each
%     spectral resolution.
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
% In the following, 'S' is the final number of spectral bands, described in
% the documentation of `sampling_options` above. The sizes and datatypes of
% the output arguments are given below for the case where
% `sampling_options.show_steps` is `false`.
%
% bands -- Spectral bands
%   A vector of length 'S' containing the wavelengths corresponding to the
%   spectral bands in `I`.
%
% I -- Latent image
%   A size(J, 1) x size(J, 2) x S array, storing the estimated latent image
%   corresponding to `J`.
%
% I_rgb -- Latent colour image
%   The colour equivalent of the latent image, generated using the colour
%   space conversion data in `color_map`. An size(J, 1) x size(J, 2) x
%   size(color_map, 1) array.
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
%   `true`, or `reg_options.multi_weights` exists, then `weights_images` is
%   empty (`[]`);
%
% J_full -- Warped latent colour image
%   A colour image produced by warping `I` according to the dispersion
%   model, followed by conversion to the colour space of `J`. An size(J, 1)
%   x size(J, 2) x size(color_map, 1) array.
%
% J_est -- Re-estimated input RAW image
%   The mosaiced version of `J_full`, representing the result of passing
%   `I` through the forward model of dispersion and image capture. An array
%   with the same dimensions as `J`.
%
% I_warped -- Warped latent image
%   An size(J, 1) x size(J, 2) x S array, storing the latent image warped
%   according to the dispersion model.
%
% search -- Grid search method for regularization weight selection search path
%   The `search` output argument of 'weightsLowMemory()'. `search` is
%   specific to a single image patch, and therefore can only be output when
%   `patch_options.target_patch` exists. When `I_in` is empty, `search`
%   will have the additional fields that it has in the case where
%   `in_penalties` is an input argument of 'weightsLowMemory()'. `search`
%   describes the search path taken under spectral resolution 'S'. If
%   `sampling_options.show_steps` is `true`, then `search` is a cell
%   vector, where the i-th cell stores the search path for the i-th
%   spectral resolution.
%
% ## Notes
%
% - All notes in the documentation of 'solvePatchesColor()' apply to this
%   function as well.
%
% ## References
% - See references of 'solvePatchesColor()'.
%
% ## Future Work
%
% - See future work for 'solvePatchesColor()'.
% - The multi-resolution image estimation approach taken by this function
%   could be leveraged for more sophisticated regularization methods:
%   - Anisotropic smoothing: The image spatial gradient could be subject to
%     a penalty that varies with its direction and with its location in the
%     image, as a function of the image gradient in a lower spectral
%     resolution result. For instance, at strong edges in the previous
%     result, a larger spatial gradient should be allowed, provided that
%     it is aligned perpendicular to the edges.
%   - Damping: There could be an additional regularization term
%     penalizing the difference from the previous, lower spectral
%     resolution result.
% - The patch size could be varied as a function of the spectral
%   resolution.
%
% See also solvePatchesColor, baek2017Algorithm2LowMemory, weightsLowMemory

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 6, 2018

narginchk(10, 11);
nargoutchk(2, 8);

verbose = false;
if ~isempty(varargin)
    verbose = varargin{1};
end

if verbose
    tic
end

% Input argument parsing

output_search = nargout > 7;
do_single_patch = isfield(patch_options, 'target_patch');
if output_search && ~do_single_patch
    error('`search` can only be output when `patch_options.target_patch` exists.');
end

has_dispersion = ~isempty(dispersionfun);
if has_dispersion && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

nargout_before_weights = 3;
has_multi_weights = isfield(reg_options, 'multi_weights');
output_weights = (nargout > nargout_before_weights)...
    && ~all(reg_options.minimum_weights == reg_options.maximum_weights)...
    && ~has_multi_weights;
if has_multi_weights
    multi_weights = reg_options.multi_weights;
    reg_options = rmfield(reg_options, 'multi_weights');
else
    multi_weights = [];
end

input_I_in = ~isempty(I_in);
if nargout > nargout_before_weights
    n_auxiliary_images = nargout - nargout_before_weights;
else
    n_auxiliary_images = 1;
end
image_sampling = [size(J_2D, 1), size(J_2D, 2)];
patch_size = patch_options.patch_size;
padding = patch_options.padding;
n_channels_rgb = size(color_map, 1);
n_channels_rgb1 = n_channels_rgb - 1;
enabled_weights = reg_options.enabled;
n_active_weights = sum(enabled_weights);
if all(~enabled_weights) && output_weights
    error('Cannot output selected regularization weights, because all regularization terms are disabled.');
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
    if length(I_in.spectral_bands) ~= size(I_in.I, 3)
        error('The length of `I_in.spectral_bands` must equal the size of `I_in.I` in its third dimension.');
    end
end
if length(color_bands) ~= size(color_map, 2)
    error('The number of columns of `color_map` must equal the length of `color_bands`.');
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

% Choose the desired spectral resolutions
if input_I_in
    [...
      color_weights, spectral_weights_final, bands_final...
    ] = findSampling(...
      color_map, color_bands, I_in.spectral_bands, sampling_options...
    );
else
    [...
      color_weights, ~, bands_final...
    ] = findSampling(color_map, color_bands, {}, sampling_options);
end
n_bands_final = length(bands_final);
if enabled_weights(2) && n_bands_final < 2
    error('Cannot enable spectral regularization because the output image will have a single channel or spectral band.');
end
if strcmp(sampling_options.progression, 'sequential')
    n_bands_all = min(3, n_bands_final):n_bands_final;
elseif strcmp(sampling_options.progression, 'doubling')
    if n_bands_final > 3
        np2 = nextpow2(n_bands_final);
        if 2 ^ np2 == n_bands_final
            n_bands_all = pow2(0:np2);
        else
            n_bands_all = [pow2(0:(np2 - 1)) n_bands_final];
        end
    else
        n_bands_all = n_bands_final;
    end
elseif strcmp(sampling_options.progression, 'last')
    n_bands_all = n_bands_final;
else
    error('Unrecognized value of `sampling_options.progression`.');
end
n_steps = length(n_bands_all);
if has_multi_weights
    if size(multi_weights, 1) ~= n_steps
        error(['`reg_options.multi_weights` is expected to have as many r'...
            'ows as there are spectral resolutions, %d.'], n_steps);
    end
end
color_weights_all = cell(n_steps, 1);
color_weights_all{end} = color_weights;
upsampling_weights = cell(n_steps - 1, 1);
if input_I_in
    spectral_weights_all = cell(n_steps, 1);
    spectral_weights_all{end} = spectral_weights_final;
else
    spectral_weights_all = [];
end
bands_all = cell(n_steps, 1);
bands_all{end} = bands_final;
for t = 1:(n_steps - 1)
    sampling_options.n_bands = n_bands_all(t);
    if input_I_in
        [...
          color_weights_all{t}, spectral_weights_all{t}, bands_all{t}...
        ] = findSampling(...
          color_map, color_bands, I_in.spectral_bands, sampling_options...
        );
    else
        [...
          color_weights_all{t}, ~, bands_all{t}...
        ] = findSampling(color_map, color_bands, {}, sampling_options);
    end
end
for t = 1:(n_steps - 1)
    upsampling_weights{t} = resamplingWeights(...
        bands_all{t + 1}, bands_all{t},...
        sampling_options.interpolant, sampling_options.bands_padding...
    );
end
show_steps = sampling_options.show_steps;
if show_steps
    n_bands_total = sum(n_bands_all);
    n_channels_rgb_total = n_channels_rgb * n_steps;
else
    n_bands_total = n_bands_final;
    n_channels_rgb_total = n_channels_rgb;
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
sz_J3 = size(J_2D, 3);
channels_in.J = n_channels_in + [1, sz_J3]; % Input image
n_channels_in = n_channels_in + channels_in.J(2);

% Channel indices in the output concatenation of images
n_channels_out = 0;
channels_out.I = n_channels_out + [1, n_bands_total]; % Estimated latent image
n_channels_out = n_channels_out + channels_out.I(2);
if n_auxiliary_images > 0
    channels_out.I_rgb = n_channels_out + [1, n_channels_rgb_total];
    n_channels_out = n_channels_out + channels_out.I_rgb(2);
    if output_weights
        if show_steps
            channels_out.I_weights = n_channels_out + [1, n_active_weights * n_steps];
        else
            channels_out.I_weights = n_channels_out + [1, n_active_weights];
        end
        n_channels_out = n_channels_out + channels_out.I_weights(2);
    end
    if n_auxiliary_images > 1
        channels_out.J_full = n_channels_out + [1, n_channels_rgb_total];
        n_channels_out = n_channels_out + channels_out.J_full(2);
        if n_auxiliary_images > 2
            if show_steps
                channels_out.J_est = n_channels_out + [1, sz_J3 * n_steps];
            else
                channels_out.J_est = n_channels_out + [1, sz_J3];
            end
            n_channels_out = n_channels_out + channels_out.J_est(2);
            if n_auxiliary_images > 3
                channels_out.I_warped = n_channels_out + [1, n_bands_total];
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
    patches_I_ij = [];
    column_in_j = columns_in{j};
    image_sampling_p = [0, size(column_in_j, 2)];
    corner = [0, (j - 1) * patch_size(2) + 1 + patch_offset(2)];
    cols_trim_out = [ min(padding + 1, corner(2)), 0 ];
    cols_trim_out(2) = cols_trim_out(1) + min(patch_size(2) - 1, image_sampling(2) - corner(2));
    col_out_width = diff(cols_trim_out) + 1;
    column_out_j = zeros(size(column_in_j, 1), col_out_width, n_channels_out);
    
    if output_search
        if show_steps
            search_out{j} = cell(n_steps, 1);
        else
            search_out{j} = cell(1, 1);
        end
    end
    
    if has_dispersion
        dispersion_options_I_warped = struct(...
            'resolution', sampling_options.resolution,...
            'support_threshold', sampling_options.support_threshold,...
            'bands_padding', sampling_options.bands_padding,...
            'interpolant', sampling_options.interpolant,...
            'interpolant_ref', sampling_options.interpolant_ref...
        );
        dispersion_options = dispersion_options_I_warped;
        dispersion_options.bands_out = color_bands;
        dispersion_options.color_map = color_map;
        dispersion_options.int_method = sampling_options.int_method;
    end
    
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
        
        spectral_inc = 0;
        color_inc = 0;
        raw_inc = 0;
        weights_inc = 0;
        output_step = 1;
        for t = 1:n_steps
            if has_dispersion
                dispersion_options.bands_in = bands_all{t};
                dispersion_matrix_p = dispersionfunToMatrix(...
                    dispersionfun, dispersion_options, image_sampling_p, true,...
                    flip(corner) - 1 ...
                );
            else
                dispersion_matrix_p = [];
            end

            % Solve for the output patch
            n_bands_t = n_bands_all(t);
            n_px_p = prod(image_sampling_p);
            numel_p = n_px_p * n_bands_t;
            n_bands_t1 = n_bands_t - 1;
            color_weights_t = color_weights_all{t};
            enabled_weights_t = enabled_weights;
            if n_bands_t == 1
                enabled_weights_t(2) = false;
            end
            in_admm = initBaek2017Algorithm2LowMemory(...
                column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.J(1):channels_in.J(2)),...
                align, dispersion_matrix_p,...
                color_weights_t, enabled_weights_t, admm_options...
            );
            if t > 1
                % Initialize with the result of the previous step
                in_admm.I = channelConversion(patches_I_ij, upsampling_weights{t - 1}, 1);
            end
            if show_steps
                output_step = t;
                if t > 1
                    spectral_inc = spectral_inc + n_bands_all(t - 1);
                    color_inc = color_inc + n_channels_rgb;
                    raw_inc = raw_inc + sz_J3;
                    weights_inc = weights_inc + n_active_weights;
                end
            end
            
            if ~use_min_norm
                reg_options_p = reg_options;
                reg_options_p.enabled = enabled_weights_t;
                if has_multi_weights
                    reg_options_p.minimum_weights = multi_weights(t, :);
                    reg_options_p.maximum_weights = multi_weights(t, :);
                end
            end
            
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

            elseif input_I_in  || reg_options_p.demosaic
                if input_I_in
                    dispersion_matrix_p_weights = [];
                    I_in_p = struct(...
                        'I', column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.I_in(1):channels_in.I_in(2)),...
                        'spectral_weights', spectral_weights_all{t}...
                    );
                else
                    n_demosaic_channels = sum(reg_options.demosaic_channels);
                    demosaic_row_indices = zeros(n_demosaic_channels * n_px_p, 1);
                    dc_output_ind = 0;
                    for dc = 1:length(reg_options.demosaic_channels)
                        if reg_options.demosaic_channels(dc)
                            dc_output_ind = dc_output_ind + 1;
                            demosaic_row_indices(((dc_output_ind - 1) * n_px_p + 1):(dc_output_ind * n_px_p)) = (((dc - 1) * n_px_p + 1):(dc * n_px_p)).';
                        end
                    end
                    dispersion_matrix_p_weights = dispersion_matrix_p(demosaic_row_indices, :);
                    I_in_p = struct(...
                        'I', reshape(bilinearDemosaic(...
                                column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.J(1):channels_in.J(2)),...
                                align, reg_options_p.demosaic_channels...
                             ), [], 1, n_demosaic_channels),...
                        'spectral_weights', color_weights_all{t}(reg_options_p.demosaic_channels, :)...
                    );
                end
                in_weightsLowMemory = initWeightsLowMemory(I_in_p, dispersion_matrix_p_weights, numel_p);
                if output_search && (show_steps || t == n_steps)
                    [...
                        patches_I_ij, weights, ~, ~, search_out{j}{output_step}...
                    ] = weightsLowMemory(...
                        admm_options, reg_options_p, in_weightsLowMemory,...
                        in_admm, I_in_p, verbose...
                    );
                else
                    [...
                        patches_I_ij, weights...
                    ] = weightsLowMemory(...
                        admm_options, reg_options_p, in_weightsLowMemory,...
                        in_admm, I_in_p, verbose...
                    );
                end

            else
                in_penalties = initPenalties(in_admm.M_Omega_Phi, in_admm.G);
                in_weightsLowMemory = initWeightsLowMemory([], [], numel_p);
                if output_search && (show_steps || t == n_steps)
                    [...
                        patches_I_ij, weights, ~, ~, ~, search_out{j}{output_step}...
                    ] = weightsLowMemory(...
                        admm_options, reg_options_p, in_weightsLowMemory,...
                        in_admm, in_penalties, verbose...
                    );
                else
                    [...
                        patches_I_ij, weights...
                    ] = weightsLowMemory(...
                        admm_options, reg_options_p, in_weightsLowMemory,...
                        in_admm, in_penalties, verbose...
                    );
                end
            end

            if show_steps || t == n_steps
                patches_I_ij_3D = reshape(patches_I_ij, [image_sampling_p n_bands_t]);
                column_out_j(...
                    corner(1):patch_end_row, :, (spectral_inc + channels_out.I(1)):(spectral_inc + channels_out.I(1) + n_bands_t1)...
                ) = patches_I_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);

                if n_auxiliary_images > 0
                    patches_I_rgb_ij = in_admm.Omega * patches_I_ij;
                    patches_I_rgb_ij_3D = reshape(patches_I_rgb_ij, [image_sampling_p n_channels_rgb]);
                    column_out_j(...
                        corner(1):patch_end_row, :, (color_inc + channels_out.I_rgb(1)):(color_inc + channels_out.I_rgb(1) + n_channels_rgb1)...
                    ) = patches_I_rgb_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);

                    if n_auxiliary_images > 1
                        if has_dispersion
                            patches_J_full_ij = dispersion_matrix_p * patches_I_ij;
                            patches_J_full_ij_3D = reshape(patches_J_full_ij, [image_sampling_p n_channels_rgb]);
                            column_out_j(...
                                corner(1):patch_end_row, :, (color_inc + channels_out.J_full(1)):(color_inc + channels_out.J_full(1) + n_channels_rgb1)...
                            ) = patches_J_full_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                        else
                            patches_J_full_ij = patches_I_rgb_ij;
                            column_out_j(...
                                corner(1):patch_end_row, :, (color_inc + channels_out.J_full(1)):(color_inc + channels_out.J_full(1) + n_channels_rgb1)...
                            ) = column_out_j(...
                                corner(1):patch_end_row, :, (color_inc + channels_out.I_rgb(1)):(color_inc + channels_out.I_rgb(1) + n_channels_rgb1)...
                            );
                        end

                        if n_auxiliary_images > 2
                            patches_J_est_ij_3D = reshape(in_admm.M * patches_J_full_ij, [image_sampling_p sz_J3]);
                            column_out_j(...
                                corner(1):patch_end_row, :, (raw_inc + channels_out.J_est(1)):(raw_inc + channels_out.J_est(1) + sz_J3 - 1)...
                            ) = patches_J_est_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);

                            if n_auxiliary_images > 3
                                if has_dispersion
                                    dispersion_options_I_warped.bands_in = bands_all{t};
                                    patches_I_warped_ij_3D = dispersionfunToMatrix(...
                                        dispersionfun, dispersion_options_I_warped, patches_I_ij_3D, true,...
                                        flip(corner) - 1 ...
                                    );
                                    column_out_j(...
                                        corner(1):patch_end_row, :, (spectral_inc + channels_out.I_warped(1)):(spectral_inc + channels_out.I_warped(1) + n_bands_t1)...
                                    ) = patches_I_warped_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                                else
                                    column_out_j(...
                                        corner(1):patch_end_row, :, (spectral_inc + channels_out.I_warped(1)):(spectral_inc + channels_out.I_warped(1) + n_bands_t1)...
                                    ) = column_out_j(...
                                        corner(1):patch_end_row, :, (spectral_inc + channels_out.I(1)):(spectral_inc + channels_out.I(1) + n_bands_t1)...
                                    );
                                end
                            end
                        end
                    end
                end

                if output_weights
                    column_out_j(...
                        corner(1):patch_end_row, :, (weights_inc + channels_out.I_weights(1)):(weights_inc + channels_out.I_weights(1) + n_active_weights - 1)...
                    ) = repmat(reshape(weights(enabled_weights), 1, 1, n_active_weights), diff(rows_trim_out) + 1, col_out_width, 1);
                end
            end
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

if show_steps
    bands = bands_all;
else
    bands = bands_final;
end

if n_auxiliary_images > 0
    varargout = cell(1, nargout - 2);
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
    if show_steps
        varargout(end) = search_out{1};
    else
        varargout(end) = search_out{1}(1);
    end
end

if verbose
    fprintf('\tDone.\n');
    toc
end

end