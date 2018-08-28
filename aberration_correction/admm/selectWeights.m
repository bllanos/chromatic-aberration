function [ weights, patch_lim, I_patch, varargout ] = selectWeights(...
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
%   method of Belge et al. 2002
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
%   magnitudes (e.g. squared L2-norm) of the vectors obtained by applying
%   the regularization operators used in the image estimation problem to
%   `I_patch`. Preferably, the regularization errors are normalized by the
%   numbers of elements in the vectors. Also, it may be desirable for the
%   elements of the second output argument to be calculated after
%   discarding a border from the image, to eliminate border artifacts.
%
% f_args -- Additional image estimation algorithm parameters
%   A cell vector of input arguments to `f`, beyond those listed above.
%   `f_args` can be an empty cell array, `{}`.
%
% target_patch -- Patch coordinates
%   A two-element vector containing the row and column, respectively, of
%   the top-left corner of the image patch to be estimated. The
%   L-hypersurface method will only be applied to a patch of the image, for
%   efficiency, and also to customize the results to image features of
%   interest. Note that the patch has the same coordinates in both the
%   input image and the latent image, as the two images are aligned. (Refer
%   to the 'Notes' section below for the rationale.)
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
%     be included in the L-hypersurface method. If `enabled_weights(i)` is
%     `false`, the weight for the i-th regularization term will be set to
%     zero, rather than set using the L-hypersurface method.
%   - 'maxit': Maximum number of fixed-point iterations to perform in the
%     L-hypersurface method of Belge et al. 2002.
%   - 'minimum_weights': A vector of the same length as 'enabled_weights'
%     specifying minimum values for the regularization weights. The
%     elements of 'minimum_weights' will only be used if the matrix
%     operator of the data-fitting term in the image estimation problem is
%     singular. Otherwise, the smallest singular value of the data fitting
%     matrix operator, "H", will be used. Refer to section 3.4 of Belge et
%     al. 2002 for details: "An appropriate multiple of machine precision"
%     is suggested, or a value "just large enough so that the
%     regularization problem is well defined."
%   - 'initial_weights': A vector the same length as 'enabled_weights',
%     giving an initial guesses for the regularization weights.
%     'initial_weights' is used to initialize the fixed-point algorithm.
%   - 'tol': The threshold value of the relative change in any
%     regularization weight from one iteration to the next. When all
%     weights change less than this threshold, iteration terminates. Refer
%     to equation 31 of Belge et al. 2002. (Belge et al. set a value of
%     10^-4, but suggest that larger tolerances should be acceptable.)
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
% search -- L-hypersurface method search path
%   A structure with the following fields:
%   - 'weights': An array of dimensions n x length(options.enabled_weights),
%     where 'n' is one greater than the number of iterations used by the
%     L-hypersurface fixed-point algorithm to select the final
%     regularization weights. The `weights` output argument is copied to
%     `search.weights(n, :)`. `search.weights(i + 1, :)` is the estimate of
%     the regularization weights at the i-th iteration, since
%     `search.weights(1, :)` is a copy of `options.initial_weights`.
%   - 'err': An array of dimensions n x (length(options.enabled_weights) + 1)
%     containing the second output argument of `f` at each iteration.
%     `err(i)` is the second output argument of `f` when called using
%     `search.weights(i, :)` as the `weights` input argument of `f`.
%   - 'origin': The origin point of the minimum distance function used for
%     the fixed-point algorithm. A vector of length
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
%   under-constrained during image estimation. However, determining the
%   boundaries of `I_patch` such that all pixels map to the region of
%   `J_patch` is a difficult problem (when trying to do so without
%   calculating the dispersion matrix for, at worst case, the entire
%   image).
%
% See also projectionMatrix, baek2017Algorithm2, solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 22, 2018

nargoutchk(1, 4);
narginchk(8, 10);

if any(options.initial_weights(options.enabled_weights) <= 0)
    error('All initial guesses for regularization weights must be greater than zero.');
end
if any(options.minimum_weights <= 0)
    error('All minimum values for regularization weights must be greater than zero.');
end

target_patch = [];
verbose = false;
if ~isempty(varargin)
    target_patch = varargin{1};
    if length(varargin) > 1
        verbose = varargin{2};
    end
end

% Select the image patch
if isempty(target_patch)
    fg = figure;
    imshow(J);
    title('Choose the center of the image patch')
    [x,y] = ginput(1);
    target_patch = [
        max(1, round(y) - floor(options.patch_size(1) / 2)),...
        max(1, round(x) - floor(options.patch_size(2) / 2))...
    ];
    close(fg);
end
image_sampling = [size(J, 1), size(J, 2)];
[ patch_lim, ~ ] = patchBoundaries(...
  image_sampling, options.patch_size, 0, target_patch...
);

% Construct arguments for the image estimation algorithm
if isempty(align)
    align_f = [];
else
    align_f = offsetBayerPattern(patch_lim(1, :), align);
end
image_sampling_f = diff(patch_lim, 1, 1) + 1;
if ~isempty(dispersionfun)
    dispersion_f = dispersionfunToMatrix(...
        dispersionfun, lambda, image_sampling_f, image_sampling_f,...
        [0, 0, image_sampling_f(2), image_sampling_f(1)], true, flip(target_patch) - 1 ...
    );
else
    dispersion_f = [];
end
J_f = J(patch_lim(1, 1):patch_lim(2, 1), patch_lim(1, 2):patch_lim(2, 2), :);

% Select the origin of the minimum distance function
% See Section 3.4 of Belge et al. 2002.
enabled_weights = options.enabled_weights;
n_weights = length(enabled_weights);
n_active_weights = sum(enabled_weights);
n_err = n_weights + 1;
to_active_weights = double(enabled_weights);
to_active_weights(enabled_weights) = 1:n_active_weights;
to_all_weights = find(enabled_weights);
H = projectionMatrix(...
    image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
    image_sampling_f, options.int_method, []...
);
s_min = svds(H, 1, 'smallest');
H_is_singular = size(H, 1) < size(H, 2) || s_min < eps;
if H_is_singular
    % Matrix is singular
    origin_min_weights = reshape(options.minimum_weights, 1, n_weights);
else
    % Select the smallest singular value
    origin_min_weights = repmat(s_min ^ 2, 1, n_weights);
end
origin_min_weights(~enabled_weights) = 0;
origin_max_weights = repmat(svds(H, 1) ^ 2, 1, n_active_weights);

[~, err] = f(...
    image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
    J_f, origin_min_weights, f_args{:}...
);
origin = [log(err(1)), zeros(1, n_active_weights)];
for w = 1:n_active_weights
    point = zeros(1, n_weights);
    point(to_all_weights(w)) = origin_max_weights(w);
    [~, err] = f(...
        image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
        J_f, point, f_args{:}...
    );
    origin(w + 1) = log(err(to_all_weights(w) + 1));
end
if verbose
    fprintf('The origin is (%d', origin(1));
    for aw = 1:n_weights
        if enabled_weights(aw)
            fprintf(', %d', origin(to_active_weights(aw)));
        else
            fprintf(', _');
        end
    end
    fprintf('), corresponding to:\n\tminimum weights (');
    for aw = 1:n_weights
        if enabled_weights(aw)
            fprintf('%d', origin_min_weights(aw));
        else
            fprintf('_');
        end
        if aw < n_weights
            fprintf(', ');
        end
    end
    fprintf(')\n\t');
    if H_is_singular
        fprintf('(from `options.minimum_weights`)');
    else
        fprintf('(from smallest matrix singular value)');
    end
    fprintf('\n\tmaximum weights (');
    for aw = 1:n_weights
        if enabled_weights(aw)
            fprintf('%d', origin_max_weights(to_active_weights(aw)));
        else
            fprintf('_');
        end
        if aw < n_weights
            fprintf(', ');
        end
    end
    fprintf(')\n');
end

output_path = (nargout > 3);
if output_path
    search.origin = zeros(1, n_weights + 1);
    search.origin([false, enabled_weights]) = origin(2:end);
    search.origin_min_weights = origin_min_weights;
    search.origin_max_weights = zeros(1, n_weights);
    search.origin_max_weights(enabled_weights) = origin_max_weights;
end

% Fixed-point iteration
if any(options.initial_weights < origin_min_weights) ||...
   any(options.initial_weights(enabled_weights) > origin_max_weights)
    error('The initial weights are outside the bounds defining the origin of the minimum distance function.');
end

if output_path
    search.weights = zeros(options.maxit + 1, n_weights);
    search.err = zeros(options.maxit + 1, n_err);
end
weights_prev = options.initial_weights;
weights_prev(~enabled_weights) = 0;
weights = zeros(1, n_weights);
changes = Inf(1, n_weights);
converged = false;
for iter = 1:options.maxit
    [~, err] = f(...
        image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
        J_f, weights_prev, f_args{:}...
    );
    if output_path
        search.weights(iter, :) = weights_prev;
        search.err(iter, :) = err;
    end
    
    % Equation 21 of Belge et al. 2002.
    for w = 1:n_active_weights
        aw = to_all_weights(w);
        weights(aw) = (err(1) * (log(err(aw + 1)) - origin(w + 1))) / ...
            (err(aw + 1) * (log(err(1)) - origin(1)));
        if weights(aw) < 0
            warning('Negative iterate encountered in selectWeights(): weights(%d) = %g', aw, weights(aw));
        end
    end
    
    % Check for convergence (Equation 31 of Belge et al. 2002.)
    for w = 1:n_active_weights
        aw = to_all_weights(w);
        changes(aw) = abs(weights(aw) - weights_prev(aw)) / abs(weights_prev(aw));
    end
    if verbose
        fprintf('%d: err = ( %g', iter, err(1));
        for e = 2:n_err
            if enabled_weights(e - 1)
                fprintf(', %g', err(e));
            else
                fprintf(', _');
            end
        end
        fprintf(')\n%d:   weights = (', iter);
        for aw = 1:n_weights
            if enabled_weights(aw)
                fprintf('%d', weights(aw));
            else
                fprintf('_');
            end
            if aw < n_weights
                fprintf(', ');
            end
        end
        fprintf(')\n%d:   changes = (', iter);
        for aw = 1:n_weights
            if enabled_weights(aw)
                fprintf('%d', changes(aw));
            else
                fprintf('_');
            end
            if aw < n_weights
                fprintf(', ');
            end
        end
        fprintf(')\n');
    end
    
    converged = all(changes <= options.tol);
    if converged
        break;
    end
    weights_prev = weights;
end


if verbose
    if converged
        fprintf('Convergence after %d iterations.\n', iter);
    else
        fprintf('Maximum number of iterations, %d, reached without convergence.\n', iter);
    end
end

[I_patch, err] = f(...
    image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
    J_f, weights, f_args{:}...
);
if output_path
    search.weights(iter + 1, :) = weights;
    search.weights = search.weights(1:(iter + 1), :);
    search.err(iter + 1, :) = err;
    search.err = search.err(1:(iter + 1), :);
    varargout{1} = search;
end

end
