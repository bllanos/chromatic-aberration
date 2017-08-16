function [ div ] = censusDivergence( I_1, I_2, weights )
% CENSUSDIVERGENCE  Compare image windows using the Census transform
%
% ## Syntax
% div = censusDivergence( I_1, I_2, weights )
%
% ## Description
% div = censusDivergence( I_1, I_2, weights )
%   Returns the difference between the two image regions, quantified using
%   the Census transform.
%
% ## Input Arguments
%
% I_1 -- First image window
%   An h x w array, representing a greyscale image window.
%
% I_2 -- Second image window
%   An h x w array, representing a greyscale image window.
%
% weights -- Window weights
%   An h x w array of weights to apply to the per-pixel divergences within
%   the window.
%
% ## Output Arguments
%
% div -- Divergence value
%   A scalar equal to the Hamming distance between the Census transforms of
%   the two image windows.
%
% ## Notes
% - `w` and `h`, the window dimensions, must be odd integers
%
% ## References
% - J. Mustaniemi, J. Kannala, and J. Heikkila, Parallax correction via
%   disparity estimation in a multi-aperture camera," Machine Vision and
%   Applications, vol. 27, pp. 1313–1323, 2016. doi: 10.1007/s00138-016-0773-7.
%
% See also divergenceMap

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2017

x = (size(I_1, 2) + 1) / 2;
y = (size(I_1, 1) + 1) / 2;
center_1 = I_1(y, x);
center_2 = I_2(y, x);

I_1 = I_1 < center_1;
I_2 = I_2 < center_2;

div = (I_1 ~= I_2);
% Note that the centre pixels will never contribute to the Hamming distance
div = sum(sum(div .* weights));

end

