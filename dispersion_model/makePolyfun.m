function [ polyfun ] = makePolyfun(polyfun_data)
% MAKEPOLYFUN  Create a function to evaluate a polynomial function in three variables
%
% ## Syntax
% polyfun = makePolyfun(polyfun_data)
%
% ## Description
% polyfun = makePolyfun(polyfun_data)
%   Returns a function for evaluating polynomial model of disparity in
%   terms of three variables, X, Y, and lambda/colour channel.
%
% ## Input Arguments
%
% polyfun_data -- Polynomial model data
%   The `polyfun_data` output argument of 'xylambdaPolyfit()'.
%
% ## Output Arguments
%
% polyfun -- Polynomial model
%   A function which takes an input 2D array 'xylambda', where the three
%   columns represent the three independent variables (x, y, lambda). The
%   output of the function is a 2D array, 'disparity', where the columns
%   represent the two dependent variables (x and y-components of
%   disparity). 'disparity' is an evaluation of the polynomials in X, Y,
%   and lambda for the two dependent variables.
%
%   If 'xylambdaPolyfit()', which created `polyfun_data`, was modelling
%   colour channels instead of wavelengths, the third column of 'xylambda'
%   should contain the indices of the colour channels (corresponding to the
%   indices of the second dimension of `X`, the input argument of
%   'xylambdaPolyfit()').
%
% See also xylambdaPolyfit

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 4, 2018

nargoutchk(1, 1);
narginchk(1, 1);

n_spatial_dim = 2;
n_models = length(polyfun_data);
channel_mode = isfield(polyfun_data, 'reference_channel');
if channel_mode
    reference_channel = find([polyfun_data.reference_channel]);
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

            xy_normalized = (polyfun_data(d).T_points * [dataset_d(:, 1:2), ones(n_d, 1)].').';
            lambda_normalized = (polyfun_data(d).T_lambda * [dataset_d(:, 3), ones(n_d, 1)].').';
            xylambda_normalized = [
                xy_normalized(:, 1:(end - 1)),...
                lambda_normalized(:, 1:(end - 1)),...
                ];
            xylambda_normalized_3d = repmat(permute(xylambda_normalized, [1 3 2]), 1, polyfun_data(d).n_powers, 1);
            powers_rep = repmat(polyfun_data(d).powers, n_d, 1, 1);
            vandermonde_matrix = prod(xylambda_normalized_3d .^ powers_rep, 3);

            disparity_x_normalized = vandermonde_matrix * polyfun_data(d).coeff_x;
            disparity_y_normalized = vandermonde_matrix * polyfun_data(d).coeff_y;
            disparity_normalized = [disparity_x_normalized disparity_y_normalized, ones(n_d, 1)];
            disparity_d = (polyfun_data(d).T_disparity_inv * disparity_normalized.').';
            disparity(filter_d, :) = disparity_d(:, 1:(end-1));
        end
    end

polyfun = @modelfun;
end

