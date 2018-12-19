function w = triangle(x)
% TRIANGLE  Triangle interpolation kernel (linear interpolation)
%
% ## Syntax
% w = triangle(x)
%
% ## Description
% w = triangle(x)
%   Returns the interpolation weight for a sample at location `x`.
%
% ## Input Arguments
%
% x -- Sample location
%   An array containing the locations of samples being interpolated. The
%   location of the interpolation point is the origin, `x = 0`.
%
% ## Output Arguments
%
% w -- Interpolation weights
%   An array with the same dimensions as `x` containing interpolation
%   weights for the samples at the locations in `x`. The weights are such
%   that, when convolved with samples spaced by increments of one unit,
%   this function performs linear interpolation of the samples.
%
% See also resamplingWeights

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 9, 2018

narginchk(1, 1);
nargoutchk(1, 1);

w = max(1 - abs(x), 0);
    
end
