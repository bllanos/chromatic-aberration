function [ I, weights, in, in_admm, varargout ] = weightsLowMemory(...
    admm_options, options, in, in_admm, varargin...
    )
% WEIGHTSLOWMEMORY  Select regularization weights using a grid search
%
% ## Usage
%
% This function combines the functionality of 'selectWeightsGrid()' and
% 'trainWeights()', but uses 'baek2017Algorithm2LowMemory()' for image
% estimation, so as not allocate any large arrays. Also, it operates on the
% entire image, as opposed to on a patch.
%
% ## Syntax
% [ I, weights, in, in_admm ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, I_in [, verbose]...
% )
% [ I, weights, in, in_admm, search ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, I_in [, verbose]...
% )
% [ I, weights, in, in_admm, in_penalties ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, in_penalties [, verbose]...
% )
% [ I, weights, in, in_admm, in_penalties, search ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, in_penalties [, verbose]...
% )
%
% ## Description
% [ I, weights, in, in_admm ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, I_in [, verbose]...
% )
%   Returns the regularization weights, selected using the grid search
%   method of Song et al. 2016, that minimize the true mean-squared error.
%
% [ I, weights, in, in_admm, search ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, I_in [, verbose]...
% )
%   Additionally return the search path taken to select the weights.
%
% [ I, weights, in, in_admm, in_penalties ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, in_penalties [, verbose]...
% )
%   Returns the regularization weights, selected using the grid search
%   method of Song et al. 2016, that minimize the minimum distance
%   criterion.
%
% [ I, weights, in, in_admm, in_penalties, search ] = weightsLowMemory(...
%   admm_options, options, in, in_admm, in_penalties [, verbose]...
% )
%   Additionally return the search path taken to select the weights.
%
% ## Input Arguments
%
% admm_options -- Options for latent image estimation
%   The `options` input argument of 'baek2017Algorithm2LowMemory()'. Refer
%   to the documentation of baek2017Algorithm2LowMemory.m for details.
%
% options -- Regularization weight selection options
%   There are three regularization terms which can be enabled:
%   1 - Regularization of the spatial gradient of the image, as in Equation
%       6 of Baek et al. 2017.
%   2 - Regularization of the spectral gradient of the spatial gradient of
%       the image, as in Equation 6 in Baek et al. 2017.
%   3 - A second-order gradient prior designed to penalize colour-filter
%       array artifacts, implemented in antiMosaicMatrix.m.
%
%   `options` is a structure with the following fields, controlling
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
%   as the final regularization weights. In this case, `search` will be
%   empty (`[]`).
%
% in -- Preallocated intermediate data and results
%   If `I_in` is empty, `in` is a structure with no fields, as created by
%   'initWeightsLowMemory()'. Otherwise, 'initWeightsLowMemory()' can
%   create `in` as a structure with the following fields:
%   - 'Omega_Phi': A matrix mapping the vectorized estimated latent image
%     to the reference image `I_in.I`.
%   - 'I_est': A vector containing the estimated version of `I_in.I`,
%     created by applying 'Omega_Phi' to the vectorized estimated latent
%     image.
%
% in_admm -- Preallocated intermediate data and results for ADMM
%   The `in` input/output argument of 'baek2017Algorithm2LowMemory()'.
%   Refer to the documentation of baek2017Algorithm2LowMemory.m.
%
% I_in -- True image structure
%   `I_in` is used to select regularization weights based on similarity
%   with another image (such as the true image). `I_in` is must be a
%   structure with the following fields:
%   - 'I': A 2D or 3D array containing the image, against which the
%     estimated latent image will be compared.
%   - 'spectral_weights': A 2D array, where `spectral_weights(i, j)` is the
%     sensitivity of the i-th colour channel or spectral band of 'I' to the
%     j-th colour channel or spectral band of the estimated latent image.
%     `spectral_weights` is a matrix mapping colours in the estimated
%     latent image to colours in the representation of 'I'.
%     `spectral_weights` must account for any numerical intergration that
%     is part of colour/spectral conversion.
%
% in_penalties -- Preallocated intermediate data and results for 'penalties()'
%   The `in` input/output argument of 'penalties()'. Refer to the
%   documentation of penalties.m.
%
%   When `in_penalties` is passed, regularization weights will be selected
%   using the minimum distance criterion described in Song et al. 2016.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show the progress of the
%   iterative search for `weights`.
%
% ## Output Arguments
%
% I -- Estimated latent image
%   The vectorized form of the estimated latent image.
%
% weights -- Selected regularization weights
%   The value of the `weights` input argument for
%   'baek2017Algorithm2LowMemory()', chosen by the grid search method of
%   Song et al. 2016. `weights` is a vector of the same length as
%   `options.enabled`.
%
% in_admm -- Preallocated intermediate data and results for ADMM
%   The `in` input/output argument of 'baek2017Algorithm2LowMemory()'.
%   Refer to the documentation of baek2017Algorithm2LowMemory.m.
%
% in_penalties -- Preallocated intermediate data and results for 'penalties()'
%   The `in` input/output argument of 'penalties()'. Refer to the
%   documentation of penalties.m.
%
% search -- Grid search method search path
%   A structure with the following fields:
%   - 'weights': An array of dimensions (n + 1) x length(options.enabled),
%     where 'n' is the number of iterations used by the grid search
%     algorithm to select the final regularization weights. The `weights`
%     output argument is copied to `search.weights(n + 1, :)`.
%     `search.weights(i, :)` is the estimate of the regularization weights
%     at the i-th iteration.
%   - 'criterion': The best value of the search criterion (image
%     similarity, or the minimum distance criterion) at each iteration. A
%     vector of length (n + 1), where the last element correspond to
%     `search.weights(n + 1, :)`.
%
%   If `in_penalties` is passed, `search` has the following additional
%   fields:
%   - 'err': An array of dimensions (n + 1) x (length(options.enabled) + 1)
%     containing the value of `in_penalties.err` (expanded with zeros at
%     the positions of disabled regularization terms) corresponding to the
%     smallest value of the minimum distance criterion within each
%     iteration.
%   - 'err_max': The maximum values of all objectives in the image
%     estimation problem. A vector of length (length(options.enabled) + 1),
%     where the first element corresponds to the residual, and the
%     remaining elements correspond to the regularization terms. The
%     elements corresponding to disabled regularization terms are set to
%     zero. 'err_max' and 'origin' are used to scale the response surface.
%   - 'origin': The reference point of the minimum distance criterion. A
%     vector of length (length(options.enabled) + 1), where the first
%     element corresponds to the residual, and the remaining elements
%     correspond to the regularization terms. The elements corresponding to
%     disabled regularization terms are set to zero.
%
% ## Notes
% - In contrast to 'selectWeightsGrid()', the origin of the minimum
%   distance criterion is always set using `options.minimum_weights` and
%   `options.maximum_weights`, rather than chosen semi-automatically.
% - This function imitates the behaviour of 'selectWeightsGrid()' with
%   `'normalized'` as the value of its `options.scaling` argument.
%
% ## References
%
% The following article discusses the grid-search method for minimizing the
% minimum distance function, as well as the minimum distance criterion:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
%
% The minimum distance criterion was first proposed in:
%
%   Belge, M, Kilmer, M. E., & Miller, E. L.. "Efficient determination of
%     multiple regularization parameters in a generalized L-curve
%     framework." Inverse Problems, vol. 18, pp. 1161-1183, 2002.
%     doi:10.1088/0266-5611/18/4/314
%
% See also selectWeightsGrid, penalties, baek2017Algorithm2LowMemory,
% solvePatchesADMM

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 10, 2018

    function [criterion, min_ind, I_current_best, err] = getCriterion(weights)
        n = size(weights, 1);
        criterion = Inf;
        for s = 1:n
            weights_s = zeros(1, n_weights);
            weights_s(enabled_weights) = weights(s, :);
            [in_admm, weights_s] = initBaek2017Algorithm2LowMemory(...
                in_admm, weights_s, admm_options...
            );
            in_admm = baek2017Algorithm2LowMemory(weights_s, admm_options, in_admm);
            if input_I_in
                in.I_est = in.Omega_Phi * in_admm.I;
                criterion_s = immse(I_in.I(:), in.I_est);
                if criterion_s < criterion
                    criterion = criterion_s;
                    min_ind = s;
                    I_current_best = in_admm.I;
                end
            else
                % Note: This line is needed before the check `nargout > 0`
                % because `in_penalties` is an implicit output argument.
                in_penalties = penalties(in_admm.J, in_admm.I, in_admm.M_Omega_Phi, in_admm.G, admm_options.norms, in_penalties);
                if nargout > 0
                    criterion_s = sum((...
                        (in_penalties.err(err_filter) - origin) ./ err_range...
                        ).^2);
                    if criterion_s < criterion
                        criterion = criterion_s;
                        min_ind = s;
                        I_current_best = in_admm.I;
                        err = in_penalties.err(err_filter);
                    end
                end
            end
        end
    end

narginchk(5, 6);

input_I_in = isfield(varargin{1}, 'I');
if input_I_in
    nargoutchk(1, 5);
    I_in = varargin{1};
    if size(I_in.spectral_weights, 1) ~= size(I_in.I, 3)
        error('The number of rows of `I_in.spectral_weights` must equal the size of `I_in.I` in its third dimension.');
    end
    output_path = (nargout > 4);
else
    nargoutchk(1, 6);
    in_penalties = varargin{1};
    output_path = (nargout > 5);
    varargout{1} = in_penalties;
end

verbose = false;
if length(varargin) > 1
    verbose = varargin{2};
end

enabled_weights = options.enabled;
err_filter = [true, enabled_weights];
n_weights = length(enabled_weights);
n_active_weights = sum(enabled_weights);
to_active_weights = double(enabled_weights);
to_active_weights(enabled_weights) = 1:n_active_weights;
min_weights = reshape(options.minimum_weights(enabled_weights), 1, n_active_weights);
max_weights = reshape(options.maximum_weights(enabled_weights), 1, n_active_weights);

if any(min_weights <= 0)
    error('All minimum values for enabled regularization weights must be greater than zero.');
end
if any(min_weights > max_weights)
    error('All maximum values for enabled regularization weights must be greater than the corresponding minimum values.');
end

% Handle the case where the weights are fixed
if all(min_weights == max_weights)
    weights = zeros(1, n_weights);
    weights(enabled_weights) = min_weights;
    [in_admm, weights_normalized] = initBaek2017Algorithm2LowMemory(...
        in_admm, weights, admm_options...
    );
    in_admm = baek2017Algorithm2LowMemory(weights_normalized, admm_options, in_admm);
    I = in_admm.I;
    if output_path
        if input_I_in
            varargout{1} = [];
        else
            varargout{2} = [];
        end
    end
    return;
end

low_guess = reshape(options.low_guess(enabled_weights), 1, n_active_weights);
high_guess = reshape(options.high_guess(enabled_weights), 1, n_active_weights);

if any(low_guess > high_guess)
    error('All elements for enabled terms in `options.high_guess` must be greater than the corresponding elements of `options.low_guess`.');
end
if any(low_guess < min_weights)
    error('All elements for enabled terms in `options.low_guess` must be greater than or equal to the corresponding elements of `options.minimum_weights`.');
end
if any(high_guess > max_weights)
    error('All elements for enabled terms in `options.high_guess` must be less than or equal to the corresponding elements of `options.maximum_weights`.');
end

% Select the origin of the minimum distance function
% See Section 3.4 of Belge et al. 2002 and Section IV-B of Song et al. 2016
if ~input_I_in
    
    nonneg = admm_options.nonneg;
    % Select the origin in the unconstrained version of the problem
    admm_options.nonneg = false;

    getCriterion(min_weights);
    err_filtered = in_penalties.err(err_filter);
    origin = [err_filtered(1), zeros(1, n_active_weights)];
    err_max = [0, err_filtered(2:end)];
    getCriterion(max_weights);
    err_max(1) = in_penalties.err(1);
    
    for w = 1:n_active_weights
        point = min_weights;
        point(w) = max_weights(w);
        getCriterion(point);
        err_filtered = in_penalties.err(err_filter);
        origin(w + 1) = err_filtered(w + 1);
    end
    
    err_range = err_max - origin;
    
    if verbose
        fprintf('The origin is (%d', origin(1));
        for aw = 1:n_weights
            if enabled_weights(aw)
                fprintf(', %d', origin(to_active_weights(aw) + 1));
            else
                fprintf(', _');
            end
        end
        fprintf(')\n');
    end

    admm_options.nonneg = nonneg;
end

if output_path && ~input_I_in
    search.origin = [origin(1), zeros(1, n_weights)];
    search.origin([false, enabled_weights]) = origin(2:end);
    search.err_max = [err_max(1), zeros(1, n_weights)];
    search.err_max([false, enabled_weights]) = err_max(2:end);
end

% Grid search iteration

grid_side_length = 4;
grid_vectors = zeros(grid_side_length, n_active_weights);
for w = 1:n_active_weights
    grid_vectors(:, w) = logspace(...
                log10(low_guess(w)),...
                log10(high_guess(w)),...
                grid_side_length...
    ).';
end
border = 0;
eval_side_length = grid_side_length - 2 * border;
n_samples = eval_side_length ^ n_active_weights;
weights_samples = zeros(n_samples, n_active_weights);

if output_path
    search.weights = zeros(options.n_iter(1) + 1, n_weights);
    search.criterion = zeros(options.n_iter(1) + 1, 1);
    if ~input_I_in
        n_err = n_weights + 1;
        search.err = zeros(options.n_iter(1) + 1, n_err);
        to_all_err = [1, 1 + find(enabled_weights)];
    end
end

converged = false;
weights = zeros(1, n_weights);
subs = cell(n_active_weights, 1);
eval_side_length_rep = repmat(eval_side_length, 1, n_active_weights);
for iter = 1:options.n_iter(1)
    
    % Generate the grid of weights
    for w = 1:n_active_weights
        weights_samples(:, w) = repmat(repelem(...
            grid_vectors((1 + border):(end - border), w),...
            eval_side_length ^ (w - 1)...
        ), eval_side_length ^ (n_active_weights - w), 1 );
    end
    
    % Evaluate the response surface at the grid points
    if output_path && ~input_I_in
        [criterion_samples, min_ind, I_current_best, err_samples] = getCriterion(weights_samples);
    else
        [criterion_samples, min_ind, I_current_best] = getCriterion(weights_samples);
    end
        
    if output_path
        search.weights(iter, enabled_weights) = weights_samples(min_ind, :);
        search.criterion(iter) = criterion_samples;
        if ~input_I_in
            search.err(iter, to_all_err) = err_samples;
        end        
    end
    
    % Check for convergence
    if iter > 1
        change = abs(criterion_samples - criterion_prev) ./ abs(criterion_prev);
    else
        change = NaN;
    end
    converged = (change < options.tol);
    if iter == 1 || criterion_samples <= criterion_prev
        criterion_prev = criterion_samples;
        weights(enabled_weights) = weights_samples(min_ind, :);
        I = I_current_best;
        if output_path && ~input_I_in
            err_prev = err_samples;
        end 
    end
    
    if verbose
        fprintf('%d:   active weights = (', iter);
        for w = 1:n_active_weights
            fprintf('%g', weights_samples(min_ind, w));
            if w < n_active_weights
                fprintf(', ');
            end
        end
        fprintf('), criterion = %g', criterion_samples);
        if iter > 1
            if criterion_samples <= criterion_prev
                fprintf(' (better by %g)', change);
            else
                fprintf(' (worse by %g)', change);
            end
        end
        fprintf('\n');
    end
    
    if converged && iter >= options.n_iter(2)
        break;
    end
    
    % Find the next sub-grid to sample
    [subs{:}] = ind2sub(eval_side_length_rep, min_ind);
    subs_vector = cell2mat(subs) + border;
    
    if iter == 1
        for w = 1:n_active_weights
            if subs_vector(w) == 1
                % The low guess is too high
                if verbose
                    fprintf(...
                        '%d:   Low guess of %g for the %d-th active weight is too high.\n',...
                        iter, grid_vectors(1, w), w...
                    );
                end
                grid_vectors(1, w) = min_weights(w);
                grid_vectors(3, w) = grid_vectors(2, w);
                subs_vector(w) = 2;
            elseif subs_vector(w) == grid_side_length
                % The high guess is too low
                if verbose
                    fprintf(...
                        '%d:   High guess of %g for the %d-th active weight is too low.\n',...
                        iter, grid_vectors(end, w), w...
                    );
                end
                grid_vectors(end, w) = max_weights(w);
                grid_vectors(end - 2, w) = grid_vectors(end - 1, w);
                subs_vector(w) = grid_side_length - 1;
            end
        end
    
        % Subsequent iterations can test fewer points.
        border = 1;
        eval_side_length = grid_side_length - 2 * border;
        n_samples = eval_side_length ^ n_active_weights;
        weights_samples = zeros(n_samples, n_active_weights);
        eval_side_length_rep = repmat(eval_side_length, 1, n_active_weights);
    end
    
    % Set up the grid sampling for the next iteration
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

if output_path
    output_index = iter + 1;
    search.weights(output_index, :) = weights;
    search.weights = search.weights(1:output_index, :);
    
    search.criterion(output_index, :) = criterion_prev;
    search.criterion = search.criterion(1:output_index, :);
    if input_I_in
        varargout{1} = search;
    else
        search.err(output_index, to_all_err) = err_prev;
        search.err = search.err(1:output_index, :);
        varargout{2} = search;
    end
end

end