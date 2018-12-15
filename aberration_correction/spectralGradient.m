function [ G_lambda ] = spectralGradient(image_sampling, replicate)
% SPECTRALGRADIENT  Create a sparse matrix acting as an image spectral gradient operator
%
% ## Syntax
% G_lambda = spectralGradient(image_sampling, replicate)
%
% ## Description
% G_lambda = spectralGradient(image_sampling, replicate)
%   Returns a matrix representing the spectral gradient operator.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A three-element vector containing the height, width, and number of
%   colour channels or wavelength bands, respectively, of the image.
%
% replicate -- Boundary conditions
%   A Boolean value indicating how to treat the last colour channel (or
%   wavelength band) of the image. If `true`, the last colour channel will
%   be replicated prior to taking the gradient. Otherwise, the gradient
%   will not be calculated for the last colour channel.
%
% ## Output Arguments
%
% G_lambda -- Spectral gradient matrix
%   A (n_px x c_out)-by-(n_px x c) array, where `n_px =
%   prod(image_sampling(1:2))`, and  `c = image_sampling(3)`. `G_lambda`
%   produces the spectral gradient (the image gradient taken along the
%   spectral dimension - the third dimension) as follows:
%     `gradient_lambda = G_lambda * I`
%   `I` is a vectorized form of an image where all pixels have been
%   rearranged from columnwise order into a column vector. Specifically, if
%   the original image had a height of `image_sampling(1)`, a width of
%   `image_sampling(2)`, and `c = image_sampling(3)` colour channels or
%   wavelength bands, then `I` contains the data from the image in order
%   first by row, then by column, then by colour channel. `c_out` is equal
%   to `c` if `replicate` is `true`. Otherwise, `c_out` is equal to `c -
%   1`.
%
% ## Algorithm
%
% The forward difference formula is used to produce the gradient. Boundary
% colour channels are replicated to allow the formula to be applied to
% every colour channel, if `replicate` is `true`. Otherwise, the gradient
% is only calculated for colour channels with a following colour channel.
%
% See also spatialGradient

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 24, 2018

nargoutchk(1, 1);
narginchk(2, 2);

if length(image_sampling) ~= 3
    error('The `image_sampling` input argument must contain the image height, width, and number of channels/bands');
end

c = image_sampling(3);
c_nonzero = c - 1;
n_px = prod(image_sampling(1:2));
n_px_c = prod(image_sampling);
n_px_c_nonzero = n_px * c_nonzero;

% Row indices
% - Go through all pixels once, and create two matrix elements per pixel
offsets = [0; 1];
n_offsets = size(offsets, 1);
rows = repelem((1:n_px_c_nonzero).', n_offsets);

% Column indices
% Each element of the forward difference gradient is calculated from two
% colour channel values
channel_indices = repmat(offsets, n_px_c_nonzero, 1) +...
    reshape(repelem((1:c_nonzero).', n_px * n_offsets), [], 1);

% Convert to linear indices
columns = repmat(repelem((1:n_px).', n_offsets), c_nonzero, 1) +...
    (channel_indices - 1) * n_px;

% Matrix values
elements = repmat([-1; 1], n_px_c_nonzero, 1);

% Assemble the sparse matrix
if replicate
    sz_1 = n_px_c;
else
    sz_1 = n_px_c_nonzero;
end
G_lambda = sparse(...
    rows,...
    columns,...
    elements,...
    sz_1, n_px_c...
);
end
