function [ I ] = diskImage( sz, res, r, sep, c_bg, c_fg )
% DISKIMAGE  Generate a image of a grid of disks
%
% ## Syntax
% [ I ] = diskImage( sz, res, r, sep, c_bg, c_fg )
%
% ## Description
% [ I ] = diskImage( sz, res, r, sep, c_bg, c_fg )
%   Returns an image containing a grid of sharp-edged disks
%
% ## Input Arguments
%
% sz -- Image dimensions
%   A two-element vector containing the image height and width.
%
% res -- Image resolution
%   The resolution of the image, in pixels per millimetre. (Pixels are
%   assumed to be square.)
%
% r -- Disk radius
%   The radius of a disk in millimetres.
%
% sep -- Disk separation
%   The separation between the centres of adjacent disks, in millimetres.
%   Disks are separated by the same amount in both the horizontal and
%   vertical directions.
%
% c_bg -- Image background colour
%   If `c_bg` is a scalar value between 0 and 1, it is interpreted as a
%   greyscale value, and will be assigned to each of the Red, Green, and
%   Blue channels of the output image.
%
%   If `c_bg` is an RGB colour vector (3 elements) with elements in the
%   range between 0 and 1, the output image will be filled with this
%   colour.
%
%  `c_bg` is the background colour surrounding the disks.
%
% c_fg -- Image foreground colour
%   Similar to `c_bg`, but corresponds to the colour of the disks.
%
% ## Output Arguments
%
% I -- Output Image
%   An image of dimensions `sz` with a centered grid of disks meeting the
%   requirements of the input arguments.
%
% ## References
% - Mannan, F. & Langer, M. S. (2016a). "Blur calibration for depth from
%   defocus." In J. Guerrero (Ed.), 13th Conference on Computer and Robot
%   Vision, CRV 2016 (pp. 281-288). Institute of Electrical and Electronics
%   Engineers Inc. doi:10.1109/CRV.2016.62
% - Rudakova, V. & Monasse, P. (2014). "Precise correction of lateral
%   chromatic aberration in images" (Guanajuato). 6th Pacific-Rim Symposium
%   on Image and Video Technology, PSIVT 2013. Springer Verlag.
%   doi:10.1007/978-3-642-53842-1_2

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 31, 2017

nargoutchk(1, 1);
narginchk(6, 6);

if length(c_bg) == 1
    c_bg = repmat(c_bg, 1, 3);
end

if length(c_fg) == 1
    c_fg = repmat(c_fg, 1, 3);
end

I = ones(sz);
sz = flip(sz);
I = cat(3, I * c_bg(1), I * c_bg(2), I * c_bg(3));

% Disk dimensions in pixels
r = r * res;
sep = sep * res;

% Number of disks in each direction
n = floor(((sz - r * 2) / sep) + 1);

% Disk center positions
image_center = sz / 2;
disk_origin = image_center - (((n-1) / 2) * sep);
disk_centers_x = disk_origin(1) + ((1:n(1)) - 1) * sep;
disk_centers_y = disk_origin(2) + ((1:n(2)) - 1) * sep;

[disk_centers_x, disk_centers_y] = meshgrid(disk_centers_x, disk_centers_y);

% Draw circles without anti-aliasing
I = insertShape(...
    I, 'FilledCircle',...
    [disk_centers_x(:) disk_centers_y(:) repmat(r, prod(n), 1)],...
    'Color', c_fg, 'Opacity', 1, 'SmoothEdges', false...
);

end