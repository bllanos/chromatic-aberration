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
%   A Boolean value indicating how to treat the first and last colour
%   channels (or wavelength bands) of the image. If `true`, the first and
%   last colour channels will be replicated prior to taking the gradient.
%   Otherwise, the gradient will not be calculated for the first and last
%   colour channels.
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
%   2`.
%
% ## Algorithm
%
% The central difference formula is used to produce the gradient. Boundary
% colour channels are replicated to allow the formula to be applied to
% every colour channel, if `replicate` is `true`. Otherwise, the gradient
% is only calculated for colour channels with two adjacent colour channels.
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
if replicate
    c_out = c;
else
    c_out = c - 2;
end
if c_out < 1
    error('Cannot output less than one gradient channel.');
end
n_px = prod(image_sampling(1:2));
n_px_c = prod(image_sampling);
n_px_cout = n_px * c_out;

% Row indices
% - Go through all pixels once, and create two matrix elements per pixel
offsets = [1; 3];
n_offsets = size(offsets, 1);
rows = repelem((1:n_px_cout).', n_offsets);

% Column indices
% Each element of the central difference gradient is calculated from two colour
% channel values
channel_indices = repmat(offsets, n_px * (c - 2), 1) +...
    repelem((0:(c - 3)).', n_px * n_offsets);
if replicate
    channel_indices = [
        ones(n_px * n_offsets, 1);
        channel_indices;
        repmat(c, n_px * n_offsets, 1)
    ];
end

% Convert to linear indices
columns = repmat(repelem((1:n_px).', 2), c_out, 1) + (channel_indices - 1) * n_px;

% Matrix values
elements = repmat([-0.5; 0.5], n_px_cout);

% Assemble the sparse matrix
G_lambda = sparse(...
    rows,...
    columns,...
    elements,...
    n_px_cout, n_px_c...
);
end
