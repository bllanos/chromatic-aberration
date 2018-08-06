function align_new = offsetBayerPattern(corner, align)
% OFFSETBAYERPATTERN  Determine the Bayer pattern for an offset image
%
% ## Syntax
% align = offsetBayerPattern(corner, align)
%
% ## Description
% align = offsetBayerPattern(corner, align)
%   Returns a Bayer pattern string for an image with the given origin
%
% ## Input Arguments
%
% corner -- Offset origin location
%   The `[1, 1]` pixel location of an offset image in the pixel index space
%   of the current image. The first element is the row offset, whereas the
%   second element is the column offset.
%
% align -- Bayer pattern format
%   A four-character character vector, specifying the Bayer tile pattern of
%   the current image. For example, 'gbrg'.
%
%   This argument has the same form as the `sensorAlignment` input argument
%   of `demosaic()`.
%
% ## Output Arguments
%
% align -- Bayer pattern format
%   The Bayer tile pattern of an image taken from the current Bayer pattern
%   by choosing, as its top left pixel, the (row, column) indices of
%   `corner`.
%
% See also bayerMask, demosaic

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 17, 2018

if mod(corner(1), 2) == 0
    if mod(corner(2), 2) == 0
        align_new = align([4, 3, 2, 1]);
    else
        align_new = align([3, 4, 1, 2]);
    end
else
    if mod(corner(2), 2) == 0
        align_new = align([2, 1, 4, 3]);
    else
        align_new = align;
    end
end
end