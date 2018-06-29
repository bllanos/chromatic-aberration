function [w] = integrationWeights(x, method)
% INTEGRATIONWEIGHTS Calculate weights for numerical integration
%
% ## Syntax
% w = integrationWeights(x, method)
%
% ## Description
% w = integrationWeights(x, method)
%   Returns weights to use in a weighted sum for numerical integration
%   of a function, given the sampling locations in its 1D domain
%
% ## Input Arguments
%
% x -- Sampling locations
%   A vector of values at which the function to be integrated is sampled.
%
% method -- Integration method
%   A character vector or a string scalar indicating which numerical
%   integration method to use. If `method` is:
%   - 'rect', then weights will be produced according to a rectangle rule
%     with the value of the function within each subinterval taken as its
%     value on the left end of the subinterval.
%   - 'trap', then weights will be produced according to the trapezoid rule
%
% ## Output Arguments
%
% w -- Integration weights
%   A column vector with the same length as `x`, containing weights for
%   numerical integration. A function, represented as a vector, `f`, of
%   values sampled at the locations in `x`, has an integral in the interval
%   [x(1), x(end)] approximated by `dot(f, w)`.
%
% See also trapz

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 28, 2018

nargoutchk(1, 1);
narginchk(2, 2);

if ~isvector(x)
    error('The domain of integration must be one-dimensional.');
end

if size(x, 1) < size(x, 2)
    x = x.';
end

if isStringScalar(method) || ischar(method)
    dx = diff(x);
    if strcmp(method, 'rect')
        w = [dx; 0];
    elseif strcmp(method, 'trap')
        w = ([dx; 0] + [0; dx] ) / 2;
    else
        error('Unexpected value of `method`.');
    end
else
    error('`method` must be a character vector or a string scalar.');
end

end