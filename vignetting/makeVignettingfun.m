function [ vignettingfun ] = makeVignettingfun(data, varargin)
% MAKEVIGNETTINGFUN  Create a function to evaluate a model of vignetting
%
% ## Syntax
% vignettingfun = makeVignettingfun(polyfun_data [, T])
%
% ## Description
% vignettingfun = makeVignettingfun(polyfun_data [, T])
%   Returns a function for evaluating model of vignetting in terms of image
%   position.
%
% ## Input Arguments
%
% polyfun_data -- Model data
%   The `polyfun_data` output argument of 'vignettingPolyfit()'.
%
% T -- Coordinate transformation
%   A 3 x 3 transformation matrix applied to image positions prior to
%   evaluating the model. `T` is applied to homogenous coordinates, and is
%   assumed to be an affine transformation, as the homogenous coordinate of
%   image positions will be dropped after applying `T`. `T` might
%   compensate for cropping the image, for example.
%
% ## Output Arguments
%
% vignettingfun -- Vignetting model
%   A function which takes an input 2D array 'xy', where the two columns
%   represent image x and y-coordinates. The output of the function is a
%   column vector, 'factor', containing multiplicative factors representing
%   vignetting at the given image positions. 'factor' is an evaluation of
%   the polynomial in 'x' and 'y' for the vignetting factor. An image can
%   be corrected for vignetting by dividing pixels by the corresponding
%   elements of 'factor'.
%
% See also vignettingPolyfit, makeDispersionfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 21, 2019

nargoutchk(1, 1);
narginchk(1, 2);

if ~isempty(varargin)
    T_frame = varargin{1};
    data.T_points = data.T_points * T_frame;
end

    function factors = modelfun(xy)
        n = size(xy, 1);
        xy_normalized = (data.T_points * [xy(:, 1:2), ones(n, 1)].').';
        xy_normalized_3d = permute(xy_normalized(:, 1:(end - 1)), [1 3 2]);
        
        factors = ones(n, 2);
        for j = 1:n
            vandermonde_vector = prod(...
                repmat(xy_normalized_3d(j, :, :), 1, data.n_powers, 1)...
                .^ data.powers, 3 ...
                );
            factors(j, 1) = dot(vandermonde_vector, data.coeff);
        end
        factors = (data.T_factor_inv * factors.').';
        factors = factors(:, 1:(end-1));
    end

vignettingfun = @modelfun;
end

