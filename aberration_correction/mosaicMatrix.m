function [ M ] = mosaicMatrix(image_sampling, align)
% MOSAICMATRIX  Create a sparse matrix to mosaic an image
%
% ## Syntax
% M = mosaicMatrix(image_sampling, align)
%
% ## Description
% M = mosaicMatrix(image_sampling, align)
%   Returns a matrix for mosaicing an image.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the image height and width, respectively.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'. `align` has the same form as the `sensorAlignment`
%   input argument of `demosaic()`.
%
% ## Output Arguments
%
% M -- Mosaicing matrix
%   A n_px x (n_px * 3) array, where `n_px = prod(image_sampling)`. `M`
%   converts 3-channel (RGB) images to RAW (colour-filter array) images:
%     `J = M * I`
%   `I` is a vectorized form of an image where all pixels have been
%   rearranged from columnwise order into a column vector. Specifically, if
%   the original image had a height of `image_sampling(1)`, a width of
%   `image_sampling(2)`, and 3 colour channels, then `I` contains the data
%   from the image in order first by row, then by column, then by colour
%   channel. `J` is a vectorized form of the 2D RAW image.
%
% See also mosaic, bayerMask, channelConversionMatrix, bilinearDemosaic,
% demosaic

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 27, 2018

nargoutchk(1, 1);
narginchk(2, 2);

if length(image_sampling) ~= 2
    error('The `image_sampling` input argument must contain an image height and width only.');
end

mask = bayerMask(image_sampling(1), image_sampling(2), align);

n_px = prod(image_sampling);
n_channels = size(mask, 3);
n_px_c = n_px * n_channels;

% Column indices
columns = find(mask);

% Row indices
rows = columns;
for c = 1:n_channels
    filter = (rows > (n_px * (c - 1))) & (rows <= (n_px * c));
    rows(filter) = rows(filter) - (n_px * (c - 1));
end

% Matrix values
elements = ones(n_px, 1);

% Assemble the sparse matrix
M = sparse(...
    rows,...
    columns,...
    elements,...
    n_px, n_px_c...
);
end

