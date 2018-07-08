function [ splinefun, splinefun_data ] = xylambdaSplinefit(...
    X, x_field, disparity, disparity_field, smooth, varargin...
)
% XYLAMBDASPLINEFIT  Fit a spline model in two or three variables
%
% ## Syntax
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, smooth, lambda...
% )
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, smooth...
% )
%
% ## Description
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field, lambda...
% )
%   Returns a spline model of disparity in terms of three variables X,
%   Y, and lambda.
%
% [ splinefun, splinefun_data ] = xylambdaSplinefit(...
%     X, x_field, disparity, disparity_field...
% )
%   Returns a composite model of disparity, with one spline model per
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
%   A non-negative smoothing weight applied during spline fitting.
%
% lambda -- Data for the third independent variable
%   A vector of values of the third independent variable. `lambda(j)` is
%   the value corresponding to the data in `X(:, j, 1).(x_field)` and
%   `disparity(:, :, j).(disparity_field)`.
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
%   makeSplinefun(splinefun_data)`. `splinefun_data` has the following fields:
%   - T_points: The normalization transformation applied to input spatial
%     coordinates prior to evaluating the model. A 3 x 3 matrix.
%   - T_lambda: In the case of data for wavelengths, the normalization
%     transformation applied to wavelengths prior to evaluating the
%     polynomial. In the case of colour channels, this is the identity
%     transformation. A 2 x 2 matrix.
%   - T_disparity_inv: The inverse normalization transformation applied to
%     the output of the polynomial to obtain disparity values. A 3 x 3
%     matrix.
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
% This code was ported from the 3D thin plate splines code in the Geometric
% Tools library, written by David Eberly, available at
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
% See also statsToDisparity, normalizePointsPCA, makeSplinefun,
% xylambdaPolyfit

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 6, 2018

nargoutchk(1, 2);
narginchk(5, 6);

if smooth < 0
    error('The smoothing parameter must be greater than zero.');
end

sz = size(X);
n_models = sz(2);

lambda = (1:n_models).';
channel_mode = true;
if length(varargin) > 1
    lambda = varargin{1};
    n_models = 1;
    channel_mode = false;
    if sz(2) ~= length(lambda)
        error('Expected as many wavelengths in `lambda` as the size of `X` in the third dimension.')
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
    'T_points', cell(n_models, 1),...
    'T_lambda', cell(n_models, 1),...
    'T_disparity_inv', cell(n_models, 1),...
    'xylambda_training', cell(n_models, 1),...
    'coeff_basis', cell(n_models, 1),...
    'coeff_affine', cell(n_models, 1)...
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
    
    % 'A' matrix
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
    A = M + smooth * eye(n_points_c);
    
    % 'B' matrix
    B = [ones(n_points_c, 1), splinefun_data(c).xylambda_training];
    
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

        splinefun_data(c).coeff_affine = Q \ Pw;
        splinefun_data(c).coeff_basis = A \ (dataset_normalized_disparity(:, dim) - B * splinefun_data(c).coeff_affine);
    end
end

splinefun = makesplinefun(splinefun_data);
end