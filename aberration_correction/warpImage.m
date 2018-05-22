function [ J ] = warpImage(I, W, varargin)
% WARPIMAGE  Warp an image using a warp matrix
%
% ## Syntax
% J = warpImage(I, W [, sz])
%
% ## Description
% J = warpImage(I, W [, sz])
%   Returns the image warped with the given warp matrix
%
% ## Input Arguments
%
% I -- Input image
%   A 2D or 3D array containing an image.
%
% W -- Warp matrix
%   A m x n array which warps the input image to the output image,
%   according to the equation:
%     `J_vector = W * I_vector`
%   `I_vector` is a vectorized form of `I` where all pixels have been
%   rearranged from columnwise order into a column vector. If `W`
%   represents a warp, `I` and `J` would have the same number of channels,
%   but `W` can be any matrix with compatible dimensions for multiplication
%   with `I`. `W` may be the first output argument of 'polyfunToMatrix()',
%   for example.
%
% sz -- Image dimensions
%   A two-element or three element vector containing the image height, and
%   width, or the image height, image width, and number of channels,
%   respectively, outupt image `J`. If `sz` has two elements, it will be
%   extended with a third element equal to the size of `I` in the third
%   dimension. If `sz` is not passed, it will be set equal to the
%   dimensions of `I`.
%
% ## Output Arguments
%
% J -- Warped image
%   The image, of dimensions `sz`, created by applying the matrix
%   transformation `W` to the input image `I`.
%
% See also polyfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 22, 2018

nargoutchk(1, 1);
narginchk(2, 3);

if numel(I) ~= size(W, 2)
    error('Input image and warp matrix have mismatched dimensions.')
end

if ~isempty(varargin)
    sz = varargin{1};
    if length(sz) == 2
        sz(3) = size(I, 3);
    end
else
    sz = size(I);
end

I_vector = I(:);
J_vector = W * I_vector;
J = full(reshape(J_vector, sz));

end

