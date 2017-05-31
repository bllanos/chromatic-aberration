function [ I ] = uniformImage( sz, c )
% UNIFORMIMAGE  Generate a uniformly-coloured image
%
% ## Syntax
% [ I ] = uniformImage( sz, c )
%
% ## Description
% [ I ] = uniformImage( sz, c )
%   Returns an image with the given size and colour.
%
% ## Input Arguments
%
% sz -- Image dimensions
%   A two-element vector containing the image height and width.
%
% c -- Image colour
%   If `c` is a scalar value between 0 and 1, it is interpreted as a
%   greyscale value, and will be assigned to each of the Red, Green, and
%   Blue channels of the output image.
%
%   If `c` is an RGB colour vector (3 elements) with elements in the range
%   between 0 and 1, the output image will be filled with this colour.
%
% ## Output Arguments
%
% I -- Output Image
%   An image of dimensions `sz` with colour `c`.
%
% ## References
% - Mosleh, A., Green, P., Onzon, E., Begin, I., & Langlois, J. M. P. (2015).
%   "Camera intrinsic blur kernel estimation: A reliable framework." In
%   IEEE Conference on Computer Vision and Pattern Recognition, CVPR
%   2015 (Vol. 07-12-June-2015, pp. 4961-4968). IEEE Computer
%   Society. doi:10.1109/CVPR.2015.7299130

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 31, 2017

nargoutchk(1, 1);
narginchk(2, 2);

if length(c) == 1
    c = repmat(c, 1, 3);
end

I = ones(sz);
I = cat(3, I * c(1), I * c(2), I * c(3));

end