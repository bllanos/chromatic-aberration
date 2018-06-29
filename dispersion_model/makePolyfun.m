function [ polyfun ] = makePolyfun(polyfun_data, varargin)
% MAKEPOLYFUN  Create a function to evaluate a polynomial function of disparity in three variables
%
% ## Syntax
% polyfun = makePolyfun(polyfun_data [, T])
%
% ## Description
% polyfun = makePolyfun(polyfun_data [, T])
%   Returns a function for evaluating polynomial model of disparity in
%   terms of three variables, X, Y, and lambda/colour channel.
%
% ## Input Arguments
%
% polyfun_data -- Polynomial model data
%   The `polyfun_data` output argument of 'xylambdaPolyfit()'.
%
% T -- Coordinate transformation
%   A 3 x 3 transformation matrix applied to the spatial variables prior to
%   evaluating the polynomial. For instance, `T` might convert the point
%   (x, y, 1), with 'x' and 'y', in pixel coordinates, to a coordinate
%   system having its origin at the centre of the image, and with 'x' and
%   'y' measured in millimetres. `T` is applied to homogenous coordinates,
%   and is assumed to be an affine transformation.
%
%   The inverse of `T` is applied to the disparity vectors produced by the
%   polynomial model, to convert them to the coordinate frame of the input.
%   Note that disparity vectors are unaffected by the translational
%   component of the inverse of `T`.
%
% ## Output Arguments
%
% polyfun -- Polynomial model
%   A function which takes an input 2D array 'xylambda', where the three
%   columns represent two spatial coordinates and a wavelength or colour
%   channel index, (x, y, lambda). The output of the function is a 2D
%   array, 'disparity', where the columns represent the x and y-components
%   of disparity vectors. 'disparity' is an evaluation of the polynomials
%   in X, Y, and lambda for the two components of disparity vectors.
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
narginchk(1, 2);

n_spatial_dim = 2;
n_models = length(polyfun_data);
channel_mode = isfield(polyfun_data, 'reference_channel');
if channel_mode
    reference_channel = find([polyfun_data.reference_channel]);
end

if ~isempty(varargin)
    T_frame = varargin{1};
    T_frame_disparity = inv(T_frame);
    % The last column containing zeros is needed because disparity values
    % are vectors, and so cannot be translated.
    T_frame_disparity = [T_frame_disparity(:, 1:(end - 1)), zeros(size(T_frame, 1), 1)];
    for i = 1:n_models
        polyfun_data(i).T_points = polyfun_data(i).T_points * T_frame;
        % This next line relies on the assumption that both transformations
        % are affine, such that division by the homogenous coordinate is
        % not needed between the two transformations.
        polyfun_data(i).T_disparity_inv = T_frame_disparity * polyfun_data(i).T_disparity_inv;
    end
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
            xylambda_normalized_3d = permute(xylambda_normalized, [1 3 2]);
            
            disparity_normalized = ones(n_d, 3);
            for j = 1:n_d
                vandermonde_vector = prod(...
                    repmat(xylambda_normalized_3d(j, :, :), 1, polyfun_data(d).n_powers, 1)...
                    .^ polyfun_data(d).powers, 3 ...
                );
                disparity_normalized(j, 1) = dot(vandermonde_vector, polyfun_data(d).coeff_x);
                disparity_normalized(j, 2) = dot(vandermonde_vector, polyfun_data(d).coeff_y);
            end
            
            disparity_d = (polyfun_data(d).T_disparity_inv * disparity_normalized.').';
            disparity(filter_d, :) = disparity_d(:, 1:(end-1));
        end
    end

polyfun = @modelfun;
end

