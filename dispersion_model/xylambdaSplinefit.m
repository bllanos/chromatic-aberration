function [ splinefun, splinefun_data ] = xylambdaSplinefit(...
    X, x_field, disparity, disparity_field, smooth, varargin...
)
% XYLAMBDASPLINEFIT  Fit a spline model in two or three variables
%
% ## Syntax
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, smooth, lambda [, verbose]...
% )
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, smooth [, verbose]...
% )
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, options, lambda [, verbose]...
% )
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, options [, verbose]...
% )
%
% ## Description
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, smooth, lambda [, verbose]...
% )
%   Returns a spline model of disparity in terms of three variables X,
%   Y, and lambda.
%
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, smooth [, verbose]...
% )
%   Returns a composite model of disparity, with one spline model per
%   colour channel. Each model is in terms of the variables X and Y.
%
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, options, lambda [, verbose]...
% )
%   Returns a spline model of disparity in terms of three variables X,
%   Y, and lambda, using generalized cross validation to select the spline
%   smoothing parameter.
% 
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, options [, verbose]...
% )
%   Returns a composite model of disparity, with one spline model per
%   colour channel. Each model is in terms of the variables X and Y.
%   Generalized cross validation is used to select the spline smoothing
%   parameter.
%
% Both invocations can optionally return the data used to construct the
% function handle form of the model of disparity. (One to two output
% arguments can be requested.)
%
% ## Input Arguments
%
% X -- Data for the first two independent variables
%   A structure array of the form of the 'stats_real' or 'stats_ideal'
%   output arguments of 'doubleSphericalLensPSF()'. 'X' must have a
%   size of 1 in its third dimension. The second dimension of 'X' can
%   represent either wavelengths, or colour channels.
%
% x_field -- First two independent variables field
%   A character vector containing the fieldname in 'X' of the data to use
%   for the first two independent variables of the spline model.
%   `X(i, j, 1).(x_field)` is expected to be a two-element row vector.
%
% disparity -- Data for the dependent variables
%   A structure of the form of the 'disparity_raw' output argument of
%   'statsToDisparity()'. 'disparity.(name)' must have a size of 1 in its
%   fourth dimension.
%
% disparity_field -- Dependent variables field
%   A character vector containing the fieldname in 'disparity' of the data
%   to use for the dependent variables of the spline model.
%   `disparity.(disparity_field)(i, :, j)` is expected to be a two-element
%   row vector, containing the values of the two dependent variables.
%
% smooth -- Regularization weight
%   The smoothing weight applied during spline fitting, in the interval [0,
%   Inf).
%
% options -- Regularization weight selection options
%   A structure with the following fields, describing how to perform
%   generalized cross validation to select the spline smoothing weight:
%   - 'n_iter': The first element is the maximum number of grid refinement
%     iterations to perform in the method of Song et al. 2016. The second
%     is the minimum number of iterations to perform (which takes priority
%     over termination based on 'tol').
%   - 'grid_size': A two-element vector, where the first element is the
%     number of samples in the grid during the first iteration. The second
%     element is the number of samples in the grid on subsequent
%     iterations. The first element should be high if there is a concern
%     that the generalized cross validation objective function has multiple
%     local minima in the search region for the smoothing weight. The first
%     element should also be high if a more complete plot of the objective
%     function is desired (when `verbose` is `true`). Both elements must be
%     at least 4.
%   - 'minimum': The minimum value for the smoothing weight (greater than zero).
%   - 'maximum': The maximum value for the smoothing weight (greater than zero).
%   - 'tol': The threshold value of the relative change in the generalized
%     cross validation criterion from one iteration to the next. When the
%     change is less than this threshold, iteration terminates.
%
% lambda -- Data for the third independent variable
%   A vector of values of the third independent variable. `lambda(j)` is
%   the value corresponding to the data in `X(:, j, 1).(x_field)` and
%   `disparity(:, :, j).(disparity_field)`.
%
% verbose -- Debugging and visualization controls
%   If `verbose` is `true`, graphical output and console output will be
%   generated for debugging purposes.
%
% ## Output Arguments
%
% splinefun -- Spline model
%   A function which takes an input 2D array 'xylambda', where the three
%   columns represent the three independent variables. The output of the
%   function is a 2D array, 'disparity', where the columns represent the
%   two dependent variables. 'disparity' is an evaluation of the
%   spline model in X, Y, and lambda for the two dependent variables.
%
%   If 'xylambdaSplinefit()' is modelling colour channels instead of
%   wavelengths, the third column of 'xylambda' should contain the indices
%   of the colour channels (corresponding to the indices of the second
%   dimension of the input argument `X`).
%
% splinefun_data -- Spline model data
%   A structure which can be used to obtain a function with the same
%   behaviour as `splinefun` by calling `splinefun =
%   makeSplinefun(splinefun_data)`. `splinefun_data` has the following
%   fields:
%   - type: Type of model ('spline')
%   - T_points: The normalization transformation applied to input spatial
%     coordinates prior to evaluating the model. A 3 x 3 matrix.
%   - T_lambda: In the case of data for wavelengths, the normalization
%     transformation applied to wavelengths prior to evaluating the model.
%     In the case of colour channels, this is the identity transformation.
%     A 2 x 2 matrix.
%   - T_disparity_inv: The inverse normalization transformation applied to
%     the output of the model to obtain disparity values. A 3 x 3 matrix.
%   - xylambda_training: A 2D array, where the three columns contain
%     the normalized X, Y, and lambda values of the training datapoints.
%     When 'xylambdaSplinefit()' is modelling colour channels, this array
%     only contains X and Y values.
%   - coeff_basis: A matrix of coefficients of the Green's functions
%     corresponding to the rows of 'xylambda_training'. The columns contain
%     coefficients for the different dependent variables of the spline
%     model.
%   - coeff_affine: A matrix containing the coefficients for the affine
%     terms of the spline model. The columns contain coefficients for the
%     different dependent variables of the spline model.
%   - smooth: A copy of the `smooth` input argument
%
%   If 'xylambdaSplinefit()' is modelling wavelengths, `splinefun_data` is
%   a scalar structure. If 'xylambdaSplinefit()' is modelling colour
%   channels, `splinefun_data` contains one element per colour channel, and
%   each element has a Boolean field 'reference_channel' indicating whether
%   or not it is the reference colour channel. The element for the
%   reference channel has empty fields except for 'reference_channel',
%   which has a value of `true`.
%
% ## References
%
% Part of this code was ported from the 3D thin plate splines code in the
% Geometric Tools library, written by David Eberly, available at
% https://www.geometrictools.com/GTEngine/Include/Mathematics/GteIntpThinPlateSpline3.h
% (File Version: 3.0.0 (2016/06/19))
%
% Geometric Tools is licensed under the Boost Software License:
%
% Boost Software License - Version 1.0 - August 17th, 2003
% 
% Permission is hereby granted, free of charge, to any person or organization
% obtaining a copy of the software and accompanying documentation covered by
% this license (the "Software") to use, reproduce, display, distribute,
% execute, and transmit the Software, and to prepare derivative works of the
% Software, and to permit third-parties to whom the Software is furnished to
% do so, all subject to the following:
% 
% The copyright notices in the Software and this entire statement, including
% the above license grant, this restriction and the following disclaimer,
% must be included in all copies of the Software, in whole or in part, and
% all derivative works of the Software, unless such copies or derivative
% works are solely in the form of machine-executable object code generated by
% a source language processor.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
% SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
% FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% The grid search method for selecting a regularization parameter is
% described in:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
%
% I have adapted their method for generalized cross validation by
% substituting the generalized cross validation objective function for
% their objective function. Unfortunately, I don't think the generalized
% cross validation objective function is convex, so I have added the
% `options.grid_size` parameter to help the search avoid non-global minima.
%
% Generalized cross validation is described in Chapter 4 of:
%
%   Wahba, G.. (1990). "Spline models for observational data."
%     Philadelphia, Pa.: Society for Industrial and Applied Mathematics
%     SIAM. doi:10.1137/1.9781611970128
%
% See also statsToDisparity, normalizePointsPCA, makeSplinefun,
% xylambdaPolyfit

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 6, 2018

nargoutchk(1, 2);
narginchk(5, 7);

use_crossval = isstruct(smooth);
options = [];
if use_crossval
    options = smooth;
    smooth = [];
    if options.minimum <= 0
        error('`options.minimum` must be greater than zero.');
    elseif options.maximum <= 0
        error('`options.maximum` must be greater than zero.');
    elseif options.maximum <= options.minimum
        error('`options.maximum` must be greater than `options.minimum`.');
    end
    if any(options.grid_size < 4)
        error('All elements of `options.grid_size` must be at least four.');
    end
elseif smooth < 0
    error('`smooth` must be non-negative.');
end

sz = size(X);
n_models = sz(2);

lambda = (1:n_models).';
channel_mode = true;
verbose = false;
if ~isempty(varargin)
    if islogical(varargin{1})
        verbose = varargin{1};
        if length(varargin) > 1
            error('Expected the first optional argument to be a non-logical array, `lambda`, when two optional arguments are passed.');
        end
    else
        lambda = varargin{1};
        n_models = 1;
        channel_mode = false;
        if sz(2) ~= length(lambda)
            error('Expected as many wavelengths in `lambda` as the size of `X` in the third dimension.')
        end
        if length(varargin) > 1
            if islogical(varargin{2})
                verbose = varargin{2};
            else
                error('Expected the second optional argument to be a logical scalar, `verbose`, when two optional arguments are passed.');
            end
        end
    end
end

n_spatial_dim = 2;
if length(X(1).(x_field)) ~= n_spatial_dim
    error('Expected two independent variables in each element of `X`, "x", and "y".');
end

X_unpacked = permute(reshape([X.(x_field)], n_spatial_dim, []), [2 1]);
if size(lambda, 2) > size(lambda, 1)
    % Use a column vector instead of a row vector
    lambda = lambda.';
end
lambda_unpacked = repelem(lambda, sz(1), 1);
disparity_unpacked = reshape(permute(disparity.(disparity_field), [1 3 2]), [], 2);

% Find the reference colour channel
if channel_mode
    reference_channel = 0;
    for c = 1:n_models
        if sum(sum(disparity.(disparity_field)(:, :, c))) == 0
            reference_channel = c;
            break;
        end
    end
    if reference_channel == 0
        error('No reference colour channel detected based on zero disparity vectors.');
    end
end

dataset = [X_unpacked lambda_unpacked disparity_unpacked];

% Filter NaN values
dataset = dataset(all(isfinite(dataset), 2), :);
n_points = size(dataset, 1);

% Create one model per colour channel, or one model for all wavelengths
splinefun_data = struct(...
    'type', repmat({'spline'}, n_models, 1),...
    'T_points', cell(n_models, 1),...
    'T_lambda', cell(n_models, 1),...
    'T_disparity_inv', cell(n_models, 1),...
    'xylambda_training', cell(n_models, 1),...
    'coeff_basis', cell(n_models, 1),...
    'coeff_affine', cell(n_models, 1),...
    'smooth', repmat({smooth}, n_models, 1)...
);
if channel_mode
    [splinefun_data.reference_channel] = deal(false);
    splinefun_data(reference_channel).reference_channel = true;
end

for c = 1:n_models
    if channel_mode
        if c == reference_channel
            continue;
        end
        filter_c = (dataset(:, 3) == c);
    else
        filter_c = true(n_points, 1);
    end
    dataset_c = dataset(filter_c, :);
    n_points_c = size(dataset_c, 1);
    
    [dataset_normalized_points, splinefun_data(c).T_points] = normalizePointsPCA([dataset_c(:, 1:2), ones(n_points_c, 1)]);
    if channel_mode
        splinefun_data(c).T_lambda = eye(2);
        splinefun_data(c).xylambda_training = dataset_normalized_points(:, 1:(end-1));
    else
        [dataset_normalized_lambda, splinefun_data(c).T_lambda] = normalizePointsPCA([dataset_c(:, 3), ones(n_points_c, 1)]);
        splinefun_data(c).xylambda_training = [dataset_normalized_points(:, 1:(end-1)), dataset_normalized_lambda(:, 1:(end-1))];
    end
    [dataset_normalized_disparity, T_disparity] = normalizePointsPCA([dataset_c(:, 4:5), ones(n_points_c, 1)]);
    splinefun_data(c).T_disparity_inv = inv(T_disparity);
    
    distances = repmat(...
        permute(splinefun_data(c).xylambda_training, [1, 3, 2]), 1, n_points_c, 1 ...
    ) - repmat(...
        permute(splinefun_data(c).xylambda_training, [3, 1, 2]), n_points_c, 1, 1 ...
    );
    distances = sqrt(sum(distances .^ 2, 3));
    if channel_mode
        M = splineKernel2D(distances);
    else
        M = splineKernel3D(distances);
    end
    
    % 'B' matrix
    B = [ones(n_points_c, 1), splinefun_data(c).xylambda_training];
    
    if use_crossval
        % Spline smoothing parameter selection by generalized cross
        % validation. Use a grid search to minimize the generalized cross
        % validation objective function.
        [Q_B, ~] = qr(B);
        Q2_B = Q_B(:, (size(B, 2) + 1):end);
        
        grid_vector = logspace(...
            log10(options.minimum),...
            log10(options.maximum),...
            options.grid_size(1)...
        ).';
    
        if verbose
            smooth_samples_all = [
                grid_vector(2:(end - 1));
                zeros((options.grid_size(2) - 2) * (options.n_iter(1) - 1), 1)
            ];
            search_index_selected = zeros(options.n_iter(1), 1);
            err_samples_all = zeros(size(smooth_samples_all));
            search_index = 1;
            selected_index = 1;
            fprintf('Starting generalized cross validation for model %d.\n', c);
        end
                
        converged = false;
        for iter = 1:options.n_iter(1)
            
            % Evaluate the objective function at the grid points
            n_err = length(grid_vector) - 2;
            err = zeros(n_err, 1);
            for s = 2:(n_err + 1)
                smooth_s = grid_vector(s);
                IMinusInfluence = smooth_s * Q2_B * ((Q2_B.' * (M + smooth_s * eye(n_points_c)) * Q2_B) \ Q2_B.');
                err_s = IMinusInfluence * dataset_normalized_disparity;
                err_s = sum(sum(err_s .^ 2)) / n_points_c;
                err_s = err_s / ((trace(IMinusInfluence) / n_points_c) ^ 2);
                err(s - 1) = err_s;
                if verbose
                    smooth_samples_all(search_index) = smooth_s;
                    err_samples_all(search_index) = err_s;
                    search_index = search_index + 1;
                end
            end
            
            [err_current, min_ind] = min(err);
            
            % Check for convergence
            change = NaN;
            if iter == 1 || err_current <= err_prev
                if iter > 1
                    change = abs(err_current - err_prev) ./ abs(err_prev);
                end
                converged = (change < options.tol);
                err_prev = err_current;
                smooth_c = grid_vector(min_ind + 1, :);
                if verbose
                    search_index_selected(selected_index) = search_index - n_err - 1 + min_ind;
                    selected_index = selected_index + 1;
                end
            end
            
            if verbose
                fprintf(...
                    '%d:   smooth = %g, V(lambda) = %g (change: %g)\n',...
                    iter, smooth_c, err_prev, change...
                );
            end
            
            if converged && iter >= options.n_iter(2)
                break;
            end
            
            % Set up the grid sampling for the next iteration
            grid_vector = logspace(...
                log10(grid_vector(min_ind)),...
                log10(grid_vector(min_ind + 2)),...
                options.grid_size(2)...
            ).';
        end
        
        if verbose
            if converged
                fprintf('Convergence after %d iterations.\n', iter);
            else
                fprintf('Maximum number of iterations, %d, reached without convergence.\n', iter);
            end
            
            % Visualization
            smooth_samples_all = smooth_samples_all(1:(search_index - 1));
            search_index_selected = search_index_selected(1:(selected_index - 1));
            err_samples_all = err_samples_all(1:(search_index - 1));
            iteration_colors = jet(iter);
            log_smooth_all = log10(smooth_samples_all);
            log_smooth_selected = log_smooth_all(search_index_selected);
            err_selected = err_samples_all(search_index_selected);
            iter_selected = [
                1;
                floor(...
                    (search_index_selected(2:end) - options.grid_size(1) + 1) / ...
                    (options.grid_size(2) - 2)...
                ) + 2
            ];
            
            log_smooth_selected_diff = [diff(log_smooth_selected); 0];
            err_selected_diff = [diff(err_selected); 0];
            
            figure;
            hold on
            scatter(...
                ones(options.grid_size(1) - 2, 1),...
                log_smooth_all(1:(options.grid_size(1) - 2)),...
                'Marker', 'o', 'MarkerEdgeColor', iteration_colors(1, :),...
                'MarkerFaceColor', iteration_colors(1, :)...
            );
            for io = 2:iter
                start_io = (options.grid_size(1) - 2 + (io - 2) * (options.grid_size(2) - 2) + 1);
                scatter(...
                    repmat(io, options.grid_size(2) - 2, 1),...
                    log_smooth_all(...
                        start_io:(start_io + options.grid_size(2) - 3)...
                    ),...
                    'Marker', 'o', 'MarkerEdgeColor', iteration_colors(io, :),...
                    'MarkerFaceColor', iteration_colors(io, :)...
                );
            end
            plot(iter_selected, log_smooth_selected, 'k-');
            xlabel('Iteration number')
            ylabel('log_{10}(\lambda)')
            title('Search path for the smoothing weight')
            hold off
            
            figure;
            hold on
            scatter(...
                log_smooth_all(1:(options.grid_size(1) - 2)),...
                err_samples_all(1:(options.grid_size(1) - 2)),...
                'Marker', 'o', 'MarkerEdgeColor', iteration_colors(1, :),...
                'MarkerFaceColor', iteration_colors(1, :)...
            );
            for io = 2:iter
                start_io = (options.grid_size(1) - 2 + (io - 2) * (options.grid_size(2) - 2) + 1);
                end_io = (start_io + options.grid_size(2) - 3);
                scatter(...
                    log_smooth_all(start_io:end_io),...
                    err_samples_all(start_io:end_io),...
                    'Marker', 'o', 'MarkerEdgeColor', iteration_colors(io, :),...
                    'MarkerFaceColor', iteration_colors(io, :)...
                );
            end
            quiver(...
                log_smooth_selected, err_selected,...
                log_smooth_selected_diff, err_selected_diff,...
                'AutoScale', 'off', 'Color', zeros(1, 3)...
            );
            xlabel('log_{10}(\lambda)')
            ylabel('Generalized cross validation objective, V(\lambda)')
            title('Search path for the smoothing weight')
            hold off
        end

        splinefun_data(c).smooth = smooth_c;
    else
        smooth_c = smooth;
    end
    
    % 'A' matrix
    A = M + smooth_c * eye(n_points_c);
    
    % 'P' matrix
    P = B.' / A;
    
    % 'Q' matrix
    Q = P * B;
    
    if channel_mode
        splinefun_data(c).coeff_affine = zeros(3, n_spatial_dim);
    else
        splinefun_data(c).coeff_affine = zeros(4, n_spatial_dim);
    end
    splinefun_data(c).coeff_basis = zeros(n_points_c, n_spatial_dim);
    for dim = 1:n_spatial_dim
        % 'Pw'
        Pw = P * dataset_normalized_disparity(:, dim);

        splinefun_data(c).coeff_affine(:, dim) = Q \ Pw;
        splinefun_data(c).coeff_basis(:, dim) = A \ (dataset_normalized_disparity(:, dim) - B * splinefun_data(c).coeff_affine(:, dim));
    end
end

splinefun = makeDispersionfun(splinefun_data);
end