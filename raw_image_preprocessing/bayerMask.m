function [ mask ] = bayerMask( image_height, image_width, align )
% BAYERMASK  Create a logical array to index colour channels in a raw image
%
% ## Syntax
% mask = bayerMask( image_height, image_width, align )
%
% ## Description
% mask = bayerMask( image_height, image_width, align )
%   Returns logical arrays indexing colour channels in a Bayer pattern
%
% ## Input Arguments
%
% image_height -- Image height
%   The number of rows in the raw image
%
% image_width -- Image width
%   The number of columns in the raw image
%
% align -- Bayer pattern format
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'.
%
%   This argument has the same form as the `sensorAlignment` input argument
%   of `demosaic()`.
%
% ## Output Arguments
%
% mask -- Colour channel indices
%   An image_height x image_width x 3 array, where `mask(i, j, k)` is
%   `true` if the pixel at subscripts `i` and `j` belongs to the k-th
%   colour channel. (`k` is 1, 2, or 3, for Red, Green, or Blue,
%   respectively.)
%
% ## References
% - "Processing RAW Images in MATLAB", by Rob Sumner:
%   http://rcsumner.net/raw_guide/RAWguide.pdf
%
% See also demosaic

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 4, 2017

nargoutchk(1, 1);
narginchk(3, 3);

mask = false(image_height, image_width, 3);

switch align
    case 'rggb'
        mask(1:2:end,1:2:end, 1) = true; %r
        mask(2:2:end,2:2:end, 3) = true; %b
    case 'bggr'
        mask(2:2:end,2:2:end, 1) = true; %r
        mask(1:2:end,1:2:end, 3) = true; %b
    case 'grbg'
        mask(1:2:end,2:2:end, 1) = true; %r
        mask(2:2:end,1:2:end, 3) = true; %b
    case 'gbrg'
        mask(2:2:end,1:2:end, 1) = true; %r
        mask(1:2:end,2:2:end, 3) = true; %b
    otherwise
        error('Unrecognized Bayer pattern format string.');
end

mask(:, :, 2) = not(or(mask(:, :, 1), mask(:, :, 3)));

end

