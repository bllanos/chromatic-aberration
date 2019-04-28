function [ polyfun, polyfun_data ] = xylambdaPolyfit(...
    X, x_field, max_degree_xy, disparity, disparity_field, varargin...
)
% XYLAMBDAPOLYFIT  Fit a polynomial model in three variables
%
% ## Syntax
% [ polyfun, polyfun_data ] = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field,...
%     lambda, max_degree_lambda [, verbose]...
% )
% [ polyfun, polyfun_data ] = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field [, verbose]...
% )
%
% ## Description
% [ polyfun, polyfun_data ] = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field,...
%     lambda, max_degree_lambda [, verbose]...
% )
%   Returns a polynomial model of disparity in terms of three variables X,
%   Y, and lambda.
%
% [ polyfun, polyfun_data ] = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field [, verbose]...
% )
%   Returns a composite model of disparity, with one polynomial model per
%   colour channel. Each model is in terms of the variables X and Y.
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
%   for the first two independent variables of the polynomial model.
%   `X(i, j, 1).(x_field)` is expected to be a two-element row vector.
%
% max_degree_xy -- Maximum degree for the first two independent variables
%   The maximum degree of each of the first two independent variables in
%   the set of permissible polynomial models.
%
% disparity -- Data for the dependent variables
%   A structure of the form of the 'disparity_raw' output argument of
%   'statsToDisparity()'. 'disparity.(name)' must have a size of 1 in its
%   fourth dimension.
%
% disparity_field -- Dependent variables field
%   A character vector containing the fieldname in 'disparity' of the data
%   to use for the dependent variables of the polynomial model.
%   `disparity.(disparity_field)(i, :, j)` is expected to be a two-element
%   row vector, containing the values of the two dependent variables.
%
% lambda -- Data for the third independent variable
%   A vector of values of the third independent variable. `lambda(j)` is
%   the value corresponding to the data in `X(:, j, 1).(x_field)` and
%   `disparity(:, :, j).(disparity_field)`.
%
% max_degree_lambda -- Maximum degree for the third independent variable
%   The maximum degree of the third independent variable in the set of
%   permissible polynomial models. It will be capped at `length(lambda) -
%   2`, but then will be set to zero, if `lambda` has length `1`, or one,
%   if `lambda` has length `2`.
%
% verbose -- Debugging and visualization controls
%   If `verbose` is `true`, graphical output will be generated for
%   debugging purposes.
%
% ## Output Arguments
%
% polyfun -- Polynomial model
%   A function which takes an input 2D array 'xylambda', where the three
%   columns represent the three independent variables. The output of the
%   function is a 2D array, 'disparity', where the columns represent the
%   two dependent variables. 'disparity' is an evaluation of the
%   polynomials in X, Y, and lambda for the two dependent variables.
%
%   If 'xylambdaPolyfit()' is modelling colour channels instead of
%   wavelengths, the third column of 'xylambda' should contain the indices
%   of the colour channels (corresponding to the indices of the second
%   dimension of the input argument `X`).
%
% polyfun_data -- Polynomial model data
%   A structure which can be used to obtain a function with the same
%   behaviour as `polyfun` by calling `polyfun =
%   makePolyfun(polyfun_data)`. `polyfun_data` has the following fields:
%   - type: Type of model ('polynomial')
%   - T_points: The normalization transformation applied to input spatial
%     coordinates prior to evaluating the model. A 3 x 3 matrix.
%   - T_lambda: In the case of data for wavelengths, the normalization
%     transformation applied to wavelengths prior to evaluating the model.
%     In the case of colour channels, this is the identity transformation.
%     A 2 x 2 matrix.
%   - T_disparity_inv: The inverse normalization transformation applied to
%     the output of the model to obtain disparity values. A 3 x 3 matrix.
%   - powers: A 3D array where the first dimension is of size one, the
%     second dimension indexes different powers in the polynomial model,
%     and the third dimension indexes independent variables (x, y,
%     wavelength). The array contains the exponents applied to each
%     independent variable in the polynomial model.
%   - n_powers: The number of terms in the polynomial model
%   - coeff_x: A vector of coefficients of the terms in the polynomial
%     model of the first component of disparity (the x-coordinate of
%     disparity)
%   - coeff_y: A vector of coefficients of the terms in the polynomial
%     model of the second component of disparity (the y-coordinate of
%     disparity)
%   - degree_xy: The degree of the polynomial model in each of the spatial
%     independent variables (x, y)
%   - degree_lambda: The degree of the polynomial model in the wavelength
%     independent variable
%
%   If 'xylambdaPolyfit()' is modelling wavelengths, `polyfun_data` is a
%   scalar structure. If 'xylambdaPolyfit()' is modelling colour
%   channels, `polyfun_data` contains one element per colour channel, and
%   each element has a Boolean field 'reference_channel' indicating whether
%   or not it is the reference colour channel. The element for the
%   reference channel has empty fields except for 'reference_channel',
%   which has a value of `true`.
%
% ## Algorithm
% This function constructs polynomials for the two dependent variables
% having degrees from zero to max_degree_xy in each of the first two
% independent variables, and having degrees from zero to max_degree_lambda
% in the third independent variable. It selects the degrees of the final
% polynomial model (one trivariate polynomial per dependent variable) by
% cross-validation. The cross-validation error estimates for the
% polynomials for the two dependent variables are pooled.
%
% The input data (both the independent and dependent variables) are
% centered and scaled using 'normalizePointsPCA()' prior to polynomial
% fitting.
%
% When modelling data for colour channels instead of wavelengths, this
% function constructs a separate polynomial model for each colour channel.
% Each model can be thought of as having a maximum degree of zero in the
% 'lambda' variable.
%
% When modelling data for wavelengths, data is partitioned during
% cross-validation first by wavelength. All data for one wavelength is set
% aside for validation. Part of the data for the other wavelengths is added
% to this validation set, while the rest is used for training. This
% procedure is necessary because there are only a limited number of
% wavelengths. The polynomials need to be well-behaved between the
% wavelengths. Omitting all data for certain wavelengths from the training
% data, and then testing on those wavelengths is the only way to estimate
% the prediction error between wavelengths.
%
% If there are two wavelengths, the polynomial model will be constrained to
% have degree one in wavelength, rather than using cross-validation to set
% the degree in wavelength.
%
% ## References
% - MATLAB Help page for the polyfit() function
% - T. Hastie and R. Tibshirani. The Elements of Statistical Learning: Data
%   Mining, Inference, and Prediction, 2nd Edition. New York: Springer,
%   2009.
% - Rudakova, V. & Monasse, P. (2014). "Precise correction of lateral
%   chromatic aberration in images" (Guanajuato). 6th Pacific-Rim Symposium
%   on Image and Video Technology, PSIVT 2013. Springer Verlag.
%   doi:10.1007/978-3-642-53842-1_2
%
% See also doubleSphericalLensPSF, statsToDisparity, normalizePointsPCA,
% makePolyfun, crossval

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 17, 2018

    function [powers, n_powers] = powersarray(deg_xy, deg_lambda)
        exponents_xy = 0:deg_xy;
        exponents_lambda = 0:deg_lambda;
        [x_powers_grid, y_powers_grid, lambda_powers_grid] = ndgrid(...
            exponents_xy, exponents_xy, exponents_lambda...
        );
        powers = cat(3,...
            reshape(x_powers_grid, 1, []),...
            reshape(y_powers_grid, 1, []),...
            reshape(lambda_powers_grid, 1, [])...
            );
        n_powers = size(powers, 2);
    end

    function y_est = predfun(vandermonde_matrix_train, y_train, vandermonde_matrix_test)
        coeff = vandermonde_matrix_train \ y_train;
        y_est = vandermonde_matrix_test * coeff;
    end

    function mse = crossvalfun(x_train, x_test)
        x_train_3d = repmat(permute(x_train(:, 1:3), [1 3 2]), 1, n_powers, 1);
        powers_3d = repmat(powers, size(x_train, 1), 1, 1);
        vandermonde_matrix_train = prod(x_train_3d .^ powers_3d, 3);
        
        x_test_3d = repmat(permute(x_test(:, 1:3), [1 3 2]), 1, n_powers, 1);
        powers_3d = repmat(powers, size(x_test, 1), 1, 1);
        vandermonde_matrix_test = prod(x_test_3d .^ powers_3d, 3);
        
        y_est_x = predfun(vandermonde_matrix_train, x_train(:, 4), vandermonde_matrix_test);
        y_est_y = predfun(vandermonde_matrix_train, x_train(:, 5), vandermonde_matrix_test);
        mse = sum([
            (x_test(:, 4) - y_est_x);
            (x_test(:, 5) - y_est_y)
            ] .^ 2 ...
            ) / (2 * length(y_est_x));
    end

nargoutchk(1, 2);
narginchk(5, 8);

n_folds = 10; % Number of cross-validation folds

sz = size(X);
n_models = sz(2);

verbose = false;
lambda = (1:n_models).';
max_degree_lambda = 0;
channel_mode = true;
if length(varargin) == 1
    verbose = varargin{1};
elseif length(varargin) > 1
    lambda = varargin{1};
    n_lambda = length(lambda);
    max_degree_lambda = max(1, min(varargin{2}, n_lambda - 2));
    if n_lambda == 1
        max_degree_lambda = 0;
    elseif isempty(lambda)
        error('There must be at least one wavelength in `lambda`.');
    end
    n_models = 1;
    channel_mode = false;
    if sz(2) ~= n_lambda
        error('Expected as many wavelengths in `lambda` as the size of `X` in the third dimension.')
    end
    if length(varargin) > 2
        verbose = varargin{3};
    end
end
use_basic_crossval = channel_mode || n_lambda < 3;

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
if ~use_basic_crossval
    lambda_filter = repelem((1:n_lambda).', sz(1), 1);
end
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

% Filter NaN values. I assume that NaN values are evenly-distributed across
% `lambda` values. Otherwise, stratified cross-validation might be
% necessary.
dataset_filter = all(isfinite(dataset), 2);
dataset = dataset(dataset_filter, :);
n_points = size(dataset, 1);
if ~use_basic_crossval
    lambda_filter = lambda_filter(dataset_filter);
    lambda_filter_ind = find([1; diff(lambda_filter); 1]);
end

% Create one model per colour channel, or one model for all wavelengths
polyfun_data = struct(...
    'type', repmat({'polynomial'}, n_models, 1),...
    'T_points', cell(n_models, 1),...
    'T_lambda', cell(n_models, 1),...
    'T_disparity_inv', cell(n_models, 1),...
    'powers', cell(n_models, 1),...
    'n_powers', cell(n_models, 1),...
    'coeff_x', cell(n_models, 1),...
    'coeff_y', cell(n_models, 1),...
    'degree_xy', cell(n_models, 1),...
    'degree_lambda', cell(n_models, 1)...
);
if channel_mode
    [polyfun_data.reference_channel] = deal(false);
    polyfun_data(reference_channel).reference_channel = true;
end

powers = [];
n_powers = [];
mean_mse = zeros(max_degree_xy + 1, max_degree_lambda + 1);
std_mse = zeros(max_degree_xy + 1, max_degree_lambda + 1);
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
    
    [dataset_normalized_points, polyfun_data(c).T_points] = normalizePointsPCA([dataset_c(:, 1:2), ones(n_points_c, 1)]);
    if channel_mode
        dataset_normalized_lambda = [dataset_c(:, 3), ones(n_points_c, 1)];
        polyfun_data(c).T_lambda = eye(2);
    else
        [dataset_normalized_lambda, polyfun_data(c).T_lambda] = normalizePointsPCA([dataset_c(:, 3), ones(n_points_c, 1)]);
    end
    [dataset_normalized_disparity, T_disparity] = normalizePointsPCA([dataset_c(:, 4:5), ones(n_points_c, 1)]);
    polyfun_data(c).T_disparity_inv = inv(T_disparity);
    dataset_normalized = [
        dataset_normalized_points(:, 1:(end-1)),...
        dataset_normalized_lambda(:, 1:(end-1)),...
        dataset_normalized_disparity(:, 1:(end-1))
        ];

    % Use cross validation to find the optimal polynomial complexity
    for deg_xy = 0:max_degree_xy
        for deg_lambda = 0:max_degree_lambda
            [powers, n_powers] = powersarray(deg_xy, deg_lambda);
            if use_basic_crossval
                mse = crossval(@crossvalfun, dataset_normalized, 'kfold', n_folds);
            else
                mse = zeros(n_folds, n_lambda);
                for b = 1:n_lambda
                    testset_b = dataset_normalized(lambda_filter_ind(b):(lambda_filter_ind(b + 1) - 1), :);
                    trainset_b = dataset_normalized([1:(lambda_filter_ind(b) - 1), lambda_filter_ind(b + 1):end], :);
                    n_train_b = size(trainset_b, 1);
                    if n_train_b <= n_folds
                        error('Insufficient data for cross-validation.');
                    else
                        fold_size = floor(n_train_b / n_folds);
                    end
                    for f = 1:n_folds
                        perm_bf = randperm(n_train_b);
                        testset_bf = [testset_b; trainset_b(perm_bf(1:fold_size), :)];
                        trainset_bf = trainset_b(perm_bf((fold_size + 1):end), :);
                        mse(f, b) = crossvalfun(trainset_bf, testset_bf);
                    end
                end
            end
            mean_mse(deg_xy + 1, deg_lambda + 1) = mean(mse(:));
            std_mse(deg_xy + 1, deg_lambda + 1) = std(mse(:));
        end
    end

    % "Often a 'one-standard error' rule is used with cross-validation, in
    % which we choose the most parsimonious model whose error is no more than
    % one standard error above the error of the best model."
    %
    % From T. Hastie and R. Tibshirani, 2009.
    [min_mse, ind] = min(mean_mse(:));
    min_mse_std = std_mse(ind);
    choice_mse = mean_mse - min_mse_std - min_mse;
    possible_choices = (choice_mse <= 0);
    [deg_xy_matrix, deg_lambda_matrix] = ndgrid(0:max_degree_xy, 0:max_degree_lambda);
    % Prioritize simplicity with respect to wavelength
    if n_lambda == 2
        polyfun_data(c).degree_lambda = 1;
    else
        polyfun_data(c).degree_lambda = min(deg_lambda_matrix(possible_choices));
    end
    possible_choices_lambda_best = possible_choices & (deg_lambda_matrix == polyfun_data(c).degree_lambda);
    polyfun_data(c).degree_xy = min(deg_xy_matrix(possible_choices_lambda_best));

    % Visualization
    if verbose
        mse_minus1std = mean_mse - min_mse_std;
        mse_plus1std = mean_mse + min_mse_std;
        figure;
        if max_degree_lambda == 0
            hold on
            plot(deg_xy_matrix, mse_minus1std, 'r-');
            plot(deg_xy_matrix, mse_plus1std, 'r-');
            plot(deg_xy_matrix(~possible_choices), mean_mse(~possible_choices), 'go');
            plot(deg_xy_matrix(possible_choices), mean_mse(possible_choices), 'bx');
            scatter(polyfun_data(c).degree_xy, mean_mse(polyfun_data(c).degree_xy + 1), 'filled')
            hold off
            xlabel('Spatial polynomial degree')
            ylabel('Mean square error')
            legend(...
                '-1 std', '+1 std',...
                'Mean square error >1 std. dev. from min',...
                'Mean square error <=1 std. dev. from min',...
                'Selected model'...
            );
        else
            hold on
            s1 = surf(deg_xy_matrix, deg_lambda_matrix, mse_minus1std);
            set(s1, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            s2 = surf(deg_xy_matrix, deg_lambda_matrix, mse_plus1std);
            set(s2, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            s3 = surf(deg_xy_matrix, deg_lambda_matrix, mean_mse);
            set(s3, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            cdata = double(cat(3, ~possible_choices, possible_choices, zeros(size(possible_choices))));
            set(s3, 'CData', cdata, 'FaceColor', 'texturemap');
            scatter3(polyfun_data(c).degree_xy, polyfun_data(c).degree_lambda, mean_mse(polyfun_data(c).degree_xy + 1, polyfun_data(c).degree_lambda + 1), 'filled')
            hold off
            xlabel('Spatial polynomial degree')
            ylabel('Wavelength polynomial degree')
            zlabel('Mean square error')
            legend('-1 std', '+1 std', 'Mean square error (<=1 std. dev. from min in green)', 'Selected model');
        end
        
        if channel_mode
            title(sprintf('(Channel %d) Cross validation error estimates', c));
        else
            title('Cross validation error estimates');
        end
    end

    % Fit the final model using all data
    [polyfun_data(c).powers, polyfun_data(c).n_powers] = powersarray(polyfun_data(c).degree_xy, polyfun_data(c).degree_lambda);
    x_train_final = repmat(permute(dataset_normalized(:, 1:3), [1 3 2]), 1, polyfun_data(c).n_powers, 1);
    powers_final_rep = repmat(polyfun_data(c).powers, n_points_c, 1, 1);
    vandermonde_matrix_final = prod(x_train_final .^ powers_final_rep, 3);
    polyfun_data(c).coeff_x = vandermonde_matrix_final \ dataset_normalized(:, 4);
    polyfun_data(c).coeff_y = vandermonde_matrix_final \ dataset_normalized(:, 5);
end

polyfun = makeDispersionfun(polyfun_data);
end