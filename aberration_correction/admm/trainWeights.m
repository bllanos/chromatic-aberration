function [ weights, patch_lim, I_patch, varargout ] = trainWeights(...
    I, J, align, dispersionfun, sensitivity,...
    lambda, options, f, f_args, varargin...
    )
% TRAINWEIGHTS  Use a grid search method to select ideal regularization weights
%
% ## Syntax
% weights = trainWeights(...
%   I, J, align, dispersionfun, sensitivity, lambda, options,...
%   f, f_args [, target_patch, verbose]...
% )
% [ weights, patch_lim ] = trainWeights(___)
% [ weights, patch_lim, I_patch ] = trainWeights(___)
% [ weights, patch_lim, I_patch, search ] = trainWeights(___)
%
% ## Description
% weights = trainWeights(...
%   I, J, align, dispersionfun, sensitivity, lambda, options,...
%   f, f_args [, target_patch, verbose]...
% )
%   Returns the regularization weights selected by minimizing the mean
%   square error
%
% [ weights, patch_lim ] = trainWeights(___)
%   Additionally returns the boundaries of the image patch used for weight
%   selection
%
% [ weights, patch_lim, I_patch ] = trainWeights(___)
%   Additionally returns the latent image patch estimated using the
%   selected weights.
%
% [ weights, patch_lim, I_patch, search ] = trainWeights(___)
%   Additionally returns the path taken by the grid search algorithm used
%   to select the weights.
%
% ## Input Arguments
%
% I -- True image
%   A 2D or 3D array containing the ground truth latent image, against
%   which the estimated latent image will be compared.
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
% f_args -- Additional image estimation algorithm parameters
%   A cell vector of input arguments to `f`, beyond those listed above.
%   `f_args` can be an empty cell array, `{}`.
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
%   - 'border': A scalar containing the pixel width of the border region to
%     exclude when calculating the mean square error between the true and
%     estimated image patches. 'border' mitigates the impact of edge
%     artifacts on the results.
%   - 'enabled_weights': A logical vector with a length equal to the number
%     of regularization terms in the image estimation problem.
%     `enabled_weights(i)` is `true` if the i-th regularization term will
%     be included in the grid search method. If `enabled_weights(i)` is
%     `false`, the weight for the i-th regularization term will be set to
%     zero, rather than set using the grid search method.
%   - 'n_iter': The first element is the maximum number of grid refinement
%     iterations to perform. The second is the minimum number of iterations
%     to perform (which takes priority over termination based on 'tol').
%   - 'minimum_weights': A vector, of the same length as 'enabled_weights',
%     specifying minimum values for the regularization weights.
%   - 'maximum_weights': A vector, of the same length as 'enabled_weights',
%     specifying maximum values for the regularization weights.
%   - 'tol': The threshold value of the relative change in the mean square
%     error from one iteration to the next. When the change is less than
%     this threshold, iteration terminates.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show the progress of the
%   iterative search for `weights`, as well as warnings when the search
%   behaves in unexpected ways.
%
% ## Output Arguments
%
% weights -- Selected regularization weights
%   The value of the `weights` input argument to `f`, chosen by minimizing
%   the mean square error between the estimated image, `I_patch`, and `I`,
%   using a grid search method (as in Song et al. 2016). `weights` is a
%   vector of the same length as `options.enabled_weights`.
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
%     containing the mean square error between the estimated image patch
%     and the patch of `I` at each iteration.
%
% ## References
%
% The following article discusses the grid-search method for minimizing the
% a minimum distance criterion to find good regularization weights. I use
% the same grid search to minimize the mean square error with respect to
% the true image:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
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
% See also selectWeights, selectWeightsGrid, baek2017Algorithm2,
% solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 7, 2018

    function [err, I_out] = getErr(weights)
        n = size(weights, 1);
        if n == 1
            I_out = f(...
                image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
                J_f, weights, f_args{:}...
            );
            I_diff = I_out((border + 1):(end - border), (border + 1):(end - border), :) - I_patch_gt_clipped;
            err = mean(mean(mean(I_diff.^2)));
        else
            err = zeros(n, 1);
            parfor s = 1:n
                I_s = f(...
                    image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
                    J_f, weights(s, :), f_args{:}...
                );
                I_diff = I_s((border + 1):(end - border), (border + 1):(end - border), :) - I_patch_gt_clipped;
                err(s) = mean(mean(mean(I_diff.^2)));
            end
        end
    end

nargoutchk(1, 4);
narginchk(9, 11);

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
border = options.border;
I_patch_gt_clipped = I(...
    (patch_lim(1, 1) + border):(patch_lim(2, 1) - border),...
    (patch_lim(1, 2) + border):(patch_lim(2, 2) - border), : ...
);

enabled_weights = options.enabled_weights;
n_weights = length(enabled_weights);
n_active_weights = sum(enabled_weights);
to_all_weights = find(enabled_weights);

min_weights = reshape(options.minimum_weights, 1, n_weights);
max_weights = reshape(options.maximum_weights, 1, n_weights);

% Grid search iteration
output_path = (nargout > 3);
if output_path
    search.weights = zeros(options.n_iter(1), n_weights);
    search.err = zeros(options.n_iter(1), 1);
end

grid_side_length = 4;
grid_vectors = zeros(grid_side_length, n_active_weights);
for w = 1:n_active_weights
    grid_vectors(:, w) = logspace(...
                log10(min_weights(to_all_weights(w))),...
                log10(max_weights(to_all_weights(w))),...
                grid_side_length...
    ).';
end
eval_side_length = grid_side_length - 2;
n_samples = eval_side_length ^ n_active_weights;
weights_samples = zeros(n_samples, n_weights);

converged = false;
subs = cell(n_active_weights, 1);
eval_side_length_rep = repmat(eval_side_length, 1, n_active_weights);
for iter = 1:options.n_iter(1)
    
    % Generate the grid of weights
    for w = 1:n_active_weights
        weights_samples(:, to_all_weights(w)) = repmat(...
            repelem(grid_vectors(2:(end - 1), w), eval_side_length ^ (w - 1)),...
            eval_side_length ^ (n_active_weights - w), 1 ...
        );
    end
    
    % Evaluate the response surface at the grid points
    err_samples = getErr(weights_samples);
    
    % Evaluate the minimum distance criterion
    [mse_current, min_ind] = min(err_samples);
    
    if output_path
        search.weights(iter, :) = weights_samples(min_ind, :);
        search.err(iter, :) = err_samples(min_ind, :);
    end
    
    % Check for convergence
    change = NaN;
    if iter == 1 || mse_current <= mse_prev
        if iter > 1
            change = abs(mse_current - mse_prev) ./ abs(mse_prev);
        end
        converged = (change < options.tol);
        mse_prev = mse_current;
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
        fprintf('), MSE = %g (change: %g)\n', mse_prev, change);
    end
    
    if converged && iter >= options.n_iter(2)
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

if nargout > 2
    [~, I_patch] = getErr(weights);
end

if output_path
    search.weights = search.weights(1:iter, :);
    search.err = search.err(1:iter, :);
    varargout{1} = search;
end

end
