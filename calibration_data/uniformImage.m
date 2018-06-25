function [ I, I_rgb, varargout ] = uniformImage( sz, c )
% UNIFORMIMAGE  Generate a uniformly-coloured image
%
% ## Syntax
% [ I ] = uniformImage( sz, c )
% [ I, I_rgb ] = uniformImage( sz, c )
% [ I, I_rgb, fg ] = uniformImage( sz, c )
%
% ## Description
% [ I ] = uniformImage( sz, c )
%   Returns an image with the given size and colour.
% [ I, I_rgb ] = uniformImage( sz, c )
%   Additionally returns individual colour channels as full images.
% [ I, I_rgb, fg ] = uniformImage( sz, c )
%   Additionally returns a handle to a figure for a print version of the
%   image.
%
% ## Input Arguments
%
% sz -- Image dimensions
%   A two-element vector containing the image height and width in pixels.
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
% I -- Output image
%   An image of dimensions `sz` with colour `c`.
%
% I_rgb -- Output image channels
%   A 3-element cell vector where the elements are versions of `I` with the
%   same Red, Green, and Blue channels as `I`, respectively, but with the
%   two other colour channels set to zero.
%
% fg -- Output figure for printing
%   A handle to a figure containing a vector graphics version of `I` for
%   printing.
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

narginchk(2, 3);
narginchk(2, 2);
print_requested = nargout > 2;

if length(c) == 1
    c = repmat(c, 1, 3);
end

I = ones(sz);
I = cat(3, I * c(1), I * c(2), I * c(3));

if nargout > 1
    I_rgb = {I, I, I};
    I_rgb{1}(:, :, 2:3) = 0;
    I_rgb{2}(:, :, [1 3]) = 0;
    I_rgb{3}(:, :, 1:2) = 0;
end

if print_requested
    fg = figure;
    ax = axes(fg, 'Position', [0, 0, 1, 1], 'OuterPosition', [0, 0, 1, 1]);
    ax.XLim = [0, 1];
    ax.YLim = [0, 1];
    ax.Visible = 'off';
    rectangle(...
        'Position', [0, 0, 1, 1], 'FaceColor', c,...
        'LineStyle', 'none'...
    );
    varargout{1} = fg;
end

end