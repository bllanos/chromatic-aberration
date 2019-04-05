function [ Omega ] = channelConversionMatrix(image_sampling, sensitivity)
% CHANNELCONVERSIONMATRIX  Create a sparse matrix to convert images between colour spaces
%
% ## Syntax
% Omega = channelConversionMatrix(image_sampling, sensitivity)
%
% ## Description
% Omega = channelConversionMatrix(image_sampling, sensitivity)
%   Returns a matrix for changing the spectral representation of an image.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the image height and width, respectively.
%
% sensitivity -- Colour space conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   output colour channel or spectral band to the j-th input colour channel
%   or spectral band. `sensitivity` is a matrix mapping colours in the
%   input colour space or hyperspectral representation to colours in the
%   output colour space or hyperspectral representation. `sensitivity` must
%   account for any numerical intergration that is part of colour
%   conversion.
%
% ## Output Arguments
%
% Omega -- Colour channel conversion matrix
%   A (n_px x c2)-by-(n_px x c1) array, where `n_px = prod(image_sampling)`,
%   `c1 = size(sensitivity, 2)`, and `c2 = size(sensitivity, 1)`. `Omega`
%   converts images from one spectral representation to another:
%     `J = Omega * I`
%   `I` is a vectorized form of an image where all pixels have been rearranged
%   from columnwise order into a column vector. Specifically, if the original
%   image had a height of `image_sampling(1)`, a width of `image_sampling(2)`,
%   and `c1` colour channels or wavelength bands, then `I` contains the data
%   from the image in order first by row, then by column, then by colour channel.
%   `J` is a vectorized form of the image in the output colour space, with `c2`
%   colour channels.
%
% See also findSampling, integrationWeights, sonyQuantumEfficiency

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 24, 2018

nargoutchk(1, 1);
narginchk(2, 2);

if length(image_sampling) ~= 2
    error('The `image_sampling` input argument must contain an image height and width only.');
end

n_px = prod(image_sampling);
c1 = size(sensitivity, 2);
c2 = size(sensitivity, 1);
n_px_c1 = n_px * c1;
n_px_c2 = n_px * c2;

% Row indices
% For each output pixel in each output colour channel, iterate over the colour
% channels of the input colour space
rows = repelem((1:n_px_c2).', c1, 1);

% Column indices
% Iterate over the input colour channels for each pixel, and repeat for each
% output colour channel
columns = repmat(repelem((1:n_px).', c1, 1), c2, 1) +...
    repmat((0:(c1 - 1)).' * n_px, n_px_c2, 1);

% Matrix values
elements = reshape(repmat(sensitivity.', n_px, 1), [], 1);

% Assemble the sparse matrix
Omega = sparse(...
    rows,...
    columns,...
    elements,...
    n_px_c2, n_px_c1...
);
end
