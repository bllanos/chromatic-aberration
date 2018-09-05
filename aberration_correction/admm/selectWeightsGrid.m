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
%     change in the minimum distance criterion from one iteration to the
%     next. When the change is less than this threshold, iteration
%     terminates.
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
%   - 'weights': An array of dimensions
%     (n + 1) x length(options.enabled_weights), where 'n' is the number of
%     iterations used by the grid search algorithm to select the final
%     regularization weights. The `weights` output argument is copied to
%     `search.weights(n + 1, :)`. `search.weights(i, :)` is the estimate of
%     the regularization weights at the i-th iteration.
%   - 'err': An array of dimensions
%     (n + 1) x (length(options.enabled_weights) + 1) containing the second
%     output argument of 'baek2017Algorithm2()' at each iteration. `err(i)`
%     is the second output argument of 'baek2017Algorithm2()' when called
%     using `search.weights(i, :)` as the `weights` input argument
%     'baek2017Algorithm2()'.
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

    function [err_raw, I] = getErr(weights)
        [I, err_raw] = baek2017Algorithm2(...
            image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
            J_f, weights, rho, baek2017Algorithm2Options...
        );
        err_raw = err_raw.';
    end

nargoutchk(1, 4);
narginchk(8, 10);

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
% See Section 3.4 of Belge et al. 2002 and Section IV-B of Song et al. 2016
enabled_weights = options.enabled_weights;
err_filter = [true, enabled_weights];
n_weights = length(enabled_weights);
n_active_weights = sum(enabled_weights);
n_err = n_weights + 1;
to_active_weights = double(enabled_weights);
to_active_weights(enabled_weights) = 1:n_active_weights;
to_all_weights = find(enabled_weights);

origin_min_weights = reshape(options.minimum_weights, 1, n_weights);
if options.clip_weights
    origin_max_weights = reshape(options.maximum_weights(enabled_weights), 1, n_active_weights);
else
    H = projectionMatrix(...
        image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
        image_sampling_f, options.int_method, []...
    );
    s_min = svds(H, 1, 'smallest');
    H_is_singular = size(H, 1) < size(H, 2) || s_min < eps;
    if ~H_is_singular
        % Select the smallest singular value
        origin_min_weights = repmat(s_min ^ 2, 1, n_weights);
    end
    origin_max_weights = repmat(svds(H, 1) ^ 2, 1, n_active_weights);
end
origin_min_weights(~enabled_weights) = 0;

nonneg = baek2017Algorithm2Options.nonneg;
% Select the origin in the unconstrained version of the problem
baek2017Algorithm2Options.nonneg = false;

err = getErr(origin_min_weights);
origin = [err(1), zeros(1, n_active_weights)];
for w = 1:n_active_weights
    point = zeros(1, n_weights);
    point(to_all_weights(w)) = origin_max_weights(w);
    err = getErr(point);
    origin(w + 1) = err(to_all_weights(w) + 1);
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
    if options.clip_weights || H_is_singular
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
    fprintf(')\n\t');
    if options.clip_weights
        fprintf('(from `options.maximum_weights`)');
    else
        fprintf('(from largest matrix singular value)');
    end
    fprintf('\n');
end

baek2017Algorithm2Options.nonneg = nonneg;

output_path = (nargout > 3);
if output_path
    search.origin = [origin(1), zeros(1, n_weights)];
    search.origin([false, enabled_weights]) = origin(2:end);
    search.origin_min_weights = origin_min_weights;
    search.origin_max_weights = zeros(1, n_weights);
    search.origin_max_weights(enabled_weights) = origin_max_weights;
end

% Grid search iteration

if output_path
    search.weights = zeros(options.maxit + 1, n_weights);
    search.err = zeros(options.maxit + 1, n_err);
end

grid_side_length = 4;
grid_vectors = zeros(grid_side_length, n_active_weights);
for w = 1:n_active_weights
    grid_vectors(:, w) = logspace(...
                log10(origin_min_weights(to_all_weights(w))),...
                log10(origin_max_weights(w)),...
                grid_side_length...
    ).';
end
eval_side_length = grid_side_length - 2;
n_samples = eval_side_length ^ n_active_weights;
weights_samples = zeros(n_samples, n_weights);
err_samples = zeros(n_samples, n_err);

origin_rep = repmat(origin, n_samples, 1);

converged = false;
subs = cell(n_active_weights, 1);
eval_side_length_rep = repmat(eval_side_length, 1, n_active_weights);
for iter = 1:options.maxit
    
    % Generate the grid of weights
    for w = 1:n_active_weights
        weights_samples(:, to_all_weights(w)) = repmat(...
            repelem(grid_vectors(2:(end - 1), w), eval_side_length ^ (w - 1)),...
            eval_side_length ^ (n_active_weights - w), 1 ...
        );
    end
    
    % Evaluate the response surface at the grid points
    parfor s = 1:n_samples
        [~, err_samples(s, :)] = baek2017Algorithm2(...
            image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
            J_f, weights_samples(s, :), rho, baek2017Algorithm2Options...
        );
    end
    
    % Evaluate the minimum distance criterion
    distance_sq = sum((err_samples(:, err_filter) - origin_rep).^2, 2);
    [distance_sq_current, min_ind] = min(distance_sq);
    
    if output_path
        search.weights(iter, :) = weights_samples(min_ind, :);
        search.err(iter, :) = err_samples(min_ind, :);
    end
    
    % Check for convergence
    change = NaN;
    if iter == 1 || distance_sq_current <= distance_sq_prev
        if iter > 1
            change = abs(distance_sq_current - distance_sq_prev) ./ abs(distance_sq_prev);
        end
        converged = (change < options.tol);
        distance_sq_prev = distance_sq_current;
        weights = weights_samples(min_ind, :);
    end
    
    if verbose
        fprintf('%d:   weights = (', iter);
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
        fprintf('), distance = %g (change: %g)\n', distance_sq_prev, change);
    end
    
    if converged
        break;
    end
    
    % Set up the grid sampling for the next iteration
    [subs{:}] = ind2sub(eval_side_length_rep, min_ind);
    subs_vector = cell2mat(subs) + 1;
    for w = 1:n_active_weights
        grid_vectors(:, w) = logspace(...
                    log10(grid_vectors(subs_vector(w) - 1, w)),...
                    log10(grid_vectors(subs_vector(w) + 1, w)),...
                    grid_side_length...
        ).';
    end
end

if verbose
    if converged
        fprintf('Convergence after %d iterations.\n', iter);
    else
        fprintf('Maximum number of iterations, %d, reached without convergence.\n', iter);
    end
end

% Clip the final weights
if options.clip_weights
    for aw = 1:n_weights
        if enabled_weights(aw)
            if weights(aw) < options.minimum_weights(aw)
                if verbose
                    warning('Clipped weight %d from %g to the minimum value %g.', aw, weights(aw), options.minimum_weights(aw));
                end
                weights(aw) = options.minimum_weights(aw);
            elseif weights(aw) > options.maximum_weights(aw)
                if verbose
                    warning('Clipped weight %d from %g to the maximum value %g.', aw, weights(aw), options.maximum_weights(aw));
                end
                weights(aw) = options.maximum_weights(aw);
            end
        end
    end
end

[err, I_patch] = getErr(weights);
if output_path
    output_index = iter + 1;
    search.weights(output_index, :) = weights;
    search.weights = search.weights(1:output_index, :);
    search.err(output_index, :) = err;
    search.err = search.err(1:output_index, :);
    varargout{1} = search;
end

end
