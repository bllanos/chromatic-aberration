function [ fg ] = diskImagePrint( sz, r, sep, c_bg, c_fg )
% DISKIMAGEPRINT  Generate an image of a grid of disks for printing
%
% ## Syntax
% fg = diskImagePrint( sz, r, sep, c_bg, c_fg )
%
% ## Description
% fg = diskImagePrint( sz, r, sep, c_bg, c_fg )
%   Returns a handle to a figure of an image an image containing a grid of
%   sharp-edged disks
%
% ## Input Arguments
%
% sz -- Image dimensions
%   A two-element vector containing the image width and height,
%   respectively, in inches. Note that these dimensions are flipped
%   relative to the usual ordering (height, width) of image dimensions,
%   reflecting the difference in conventions between array dimensions and
%   paper sizes.
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
% fg -- Output figure for printing
%   A handle to a figure containing the image for printing. The image is a
%   centered grid of disks meeting the requirements of the input arguments.
%
% ## Notes
%
% MATLAB is not designed for vector image generation. This function uses
% the 'scatter()' plotting function to produce disks, as it is the only
% function I could find which produces filled vector graphics circles.
% Unfortunately, the radius of a plotted circle is somewhat smaller than
% the radius value given as input to 'diskImagePrint()'.
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
% File created June 25, 2018

nargoutchk(1, 1);
narginchk(5, 5);

if length(c_bg) == 1
    c_bg = repmat(c_bg, 1, 3);
end

if length(c_fg) == 1
    c_fg = repmat(c_fg, 1, 3);
end

% Background
max_size_px = 2048;
max_size_inches = max(sz);
sz_pixels = ceil(flip(sz) * max_size_px / max_size_inches);

I = ones(sz_pixels);
I = cat(3, I * c_bg(1), I * c_bg(2), I * c_bg(3));

fg = figure;
imshow(I, 'Border', 'tight');

% Disk quantities in pixels
pixelsPerMM = (max_size_px / max_size_inches) * unitsratio('inch', 'mm');
sep = sep * pixelsPerMM;

% According to MATLAB's documentation, scatter markers have sizes in units
% of points. A typography point equals 1/72 inch.
pointsPerMM = 72 * unitsratio('inch', 'mm');
r = r * pointsPerMM;

% Number of disks in each direction
sz_pixels = flip(sz_pixels);
n = floor((sz_pixels - r * 2) / sep);

% Disk center positions
image_center = sz_pixels / 2;
disk_origin = image_center - (((n-1) / 2) * sep);
disk_centers_x = disk_origin(1) + ((1:n(1)) - 1) * sep;
disk_centers_y = disk_origin(2) + ((1:n(2)) - 1) * sep;

[disk_centers_x, disk_centers_y] = meshgrid(disk_centers_x, disk_centers_y);

% Draw circles without anti-aliasing
hold on
scatter(...
    disk_centers_x(:), disk_centers_y(:),...
    pi * (r ^ 2),...
    c_fg,...
    'filled'...
    );
hold off

end