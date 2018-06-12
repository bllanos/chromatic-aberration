function [T] = pixelsToWorldTransform(image_size, pixel_size)
% PIXELSTOWORLDTRANSFORM  Create a transformation from pixel to world coordinates
%
% ## Syntax
% T = pixelsToWorldTransform(image_size, pixel_size)
%
% ## Description
% T = pixelsToWorldTransform(image_size, pixel_size)
%   Returns a transformation matrix for converting image pixel coordinates
%   to world coordinates on the image plane.
%
% ## Input Arguments
%
% image_size -- Image pixel dimensions
%   A two-element vector containing the image height and width in pixels.
%
% pixel_size -- Pixel size
%   The side length of a pixel in world units. (Pixels are assumed to be
%   square.)
%
% ## Output Arguments
%
% T -- Coordinate transformation
%   A 3 x 3 transformation matrix for converting pixel coordinates to world
%   coordinates. The world origin is at the centre of the image, and the
%   world y-coordinate increases in the opposite direction from the image
%   y-coordinate. The homogenous vector `[x; y; 1]`, with 'x' and 'y', in
%   pixel coordinates, maps to the homogenous vector `T * [x; y; 1]` in
%   world coordinates.
%
% See also makePolyfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 12, 2018

T = [
    pixel_size, 0, 0;
    0, pixel_size, 0;
    0, 0, 1
    ] * [
    1, 0, -image_size(2) / 2;
    0, -1, image_size(1) / 2;
    0, 0, 1
    ];

end

