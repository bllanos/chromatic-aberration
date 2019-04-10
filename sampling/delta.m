function w = delta(x)
% DELTA  Delta distribution interpolation kernel (i.e. no interpolation)
%
% ## Syntax
% w = delta(x)
%
% ## Description
% w = delta(x)
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
%   that pointwise sampling occurs instead of interpolation.
%
% See also resamplingWeights

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 9, 2019

narginchk(1, 1);
nargoutchk(1, 1);

w = zeros(size(x));
w(x == 0) = 1;
    
end
