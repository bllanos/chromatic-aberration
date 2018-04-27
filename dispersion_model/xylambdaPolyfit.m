function [ polyfun ] = xylambdaPolyfit(...
    X, x_field, max_degree_xy, disparity, disparity_field, varargin...
)
% XYLAMBDAPOLYFIT  Fit a polynomial model in three variables
%
% ## Syntax
% polyfun = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field,...
%     lambda, max_degree_lambda [, verbose]...
% )
% polyfun = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field [, verbose]...
% )
%
% ## Description
% polyfun = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field,...
%     lambda, max_degree_lambda [, verbose]...
% )
%   Returns a polynomial model of disparity in terms of three variables X,
%   Y, and lambda.
%
% polyfun = xylambdaPolyfit(...
%     X, x_field, max_degree_xy, disparity, disparity_field [, verbose]...
% )
%   Returns a composite model of disparity, with one polynomial model per
%   colour channel. Each model is in terms of the variables X and Y.
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
%   permissible polynomial models.
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
% ## Algorithm
% This function constructs polynomials for the two dependent variables
% having degrees from zero to max_degree_xy in each of the first two
% independent variables, and having degrees from zero to max_degree_lambda
% in the third independent variable. It selects the degrees of the final
% polynomial model (one trivariate polynomial per dependent variable) by
% cross-validation. The cross-validation error estimates for the
% polynomials for the two dependent variables are pooled.
%
% The input data (both the independent and dependent variables) are centered
% and scaled using 'normalizePointsPCA()' prior to polynomial fitting.
%
% When modelling data for colour channels instead of wavelengths, this
% function constructs a separate polynomial model for each colour channel.
% Each model can be thought of as having a maximum degree of zero in the
% 'lambda' variable.
%
% ## References
% - MATLAB Help page for the polyfit() function
% - T. Hastie and R. Tibshirani. The Elements of Statistical Learning: Data
%   Mining, Inference, and Prediction, 2nd Edition. New York: Springer,
%   2009.
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.
%
% See also doubleSphericalLensPSF, statsToDisparity, normalizePointsPCA, crossval

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

nargoutchk(1, 1);
narginchk(5, 8);

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
    max_degree_lambda = varargin{2};
    n_models = 1;
    channel_mode = false;
    if sz(2) ~= length(lambda)
        error('Expected as many wavelengths in `lambda` as the size of `X` in the third dimension.')
    end
    if length(varargin) > 2
        verbose = varargin{3};
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
T_points = cell(n_models, 1);
T_lambda = cell(n_models, 1);
T_disparity_inv = cell(n_models, 1);
powers_final = cell(n_models, 1);
n_powers_final = cell(n_models, 1);
coeff_x = cell(n_models, 1);
coeff_y = cell(n_models, 1);

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
    
    [dataset_normalized_points, T_points{c}] = normalizePointsPCA([dataset_c(:, 1:2), ones(n_points_c, 1)]);
    if channel_mode
        dataset_normalized_lambda = [dataset_c(:, 3), ones(n_points_c, 1)];
        T_lambda{c} = eye(2);
    else
        [dataset_normalized_lambda, T_lambda{c}] = normalizePointsPCA([dataset_c(:, 3), ones(n_points_c, 1)]);
    end
    [dataset_normalized_disparity, T_disparity] = normalizePointsPCA([dataset_c(:, 4:5), ones(n_points_c, 1)]);
    T_disparity_inv{c} = inv(T_disparity);
    dataset_normalized = [
        dataset_normalized_points(:, 1:(end-1)),...
        dataset_normalized_lambda(:, 1:(end-1)),...
        dataset_normalized_disparity(:, 1:(end-1))
        ];

    % Use cross validation to find the optimal polynomial complexity
    for deg_xy = 0:max_degree_xy
        for deg_lambda = 0:max_degree_lambda
            [powers, n_powers] = powersarray(deg_xy, deg_lambda);
            mse = crossval(@crossvalfun, dataset_normalized); % Defaults to 10-fold cross validation
            mean_mse(deg_xy + 1, deg_lambda + 1)  = mean(mse);
            std_mse(deg_xy + 1, deg_lambda + 1)  = std(mse);
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
    deg_lambda_best = min(deg_lambda_matrix(possible_choices));
    possible_choices_lambda_best = possible_choices & (deg_lambda_matrix == deg_lambda_best);
    degree_xy_best = min(deg_xy_matrix(possible_choices_lambda_best));

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
            scatter(degree_xy_best, mean_mse(degree_xy_best + 1), 'filled')
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
            scatter3(degree_xy_best, deg_lambda_best, mean_mse(degree_xy_best + 1, deg_lambda_best + 1), 'filled')
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
    [powers_final{c}, n_powers_final{c}] = powersarray(degree_xy_best, deg_lambda_best);
    x_train_final = repmat(permute(dataset_normalized(:, 1:3), [1 3 2]), 1, n_powers_final{c}, 1);
    powers_final_rep = repmat(powers_final{c}, n_points_c, 1, 1);
    vandermonde_matrix_final = prod(x_train_final .^ powers_final_rep, 3);
    coeff_x{c} = vandermonde_matrix_final \ dataset_normalized(:, 4);
    coeff_y{c} = vandermonde_matrix_final \ dataset_normalized(:, 5);
end

    function disparity = modelfun(xylambda)
        n_all = size(xylambda, 1);
        disparity = zeros(n_all, n_spatial_dim);
        for d = 1:n_models
            % Apply and reverse normalization
            if channel_mode
                if d == reference_channel
                    continue;
                end
                filter_d = (xylambda(:, 3) == d);
            else
                filter_d = true(n_all, 1);
            end
            dataset_d = xylambda(filter_d, :);
            n_d = size(dataset_d, 1);

            xy_normalized = (T_points{d} * [dataset_d(:, 1:2), ones(n_d, 1)].').';
            lambda_normalized = (T_lambda{d} * [dataset_d(:, 3), ones(n_d, 1)].').';
            xylambda_normalized = [
                xy_normalized(:, 1:(end - 1)),...
                lambda_normalized(:, 1:(end - 1)),...
                ];
            xylambda_normalized_3d = repmat(permute(xylambda_normalized, [1 3 2]), 1, n_powers_final{d}, 1);
            powers_rep = repmat(powers_final{d}, n_d, 1, 1);
            vandermonde_matrix = prod(xylambda_normalized_3d .^ powers_rep, 3);

            disparity_x_normalized = vandermonde_matrix * coeff_x{d};
            disparity_y_normalized = vandermonde_matrix * coeff_y{d};
            disparity_normalized = [disparity_x_normalized disparity_y_normalized, ones(n_d, 1)];
            disparity_d = (T_disparity_inv{d} * disparity_normalized.').';
            disparity(filter_d, :) = disparity_d(:, 1:(end-1));
        end
    end

polyfun = @modelfun;
end