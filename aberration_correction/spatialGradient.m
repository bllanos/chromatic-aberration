function [ G_xy ] = spatialGradient(image_sampling)
% SPATIALGRADIENT  Create a sparse matrix acting as an image spatial gradient operator
%
% ## Syntax
% G_xy = spatialGradient(image_sampling)
%
% ## Description
% G_xy = spatialGradient(image_sampling)
%   Returns a matrix representing the spatial gradient operator.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A three-element vector containing the height, width, and number of colour
%   channels or wavelength bands, respectively, of the image.
%
% ## Output Arguments
%
% G_xy -- Spatial gradient matrix
%   A (n_px x c x 2)-by-(n_px x c) array, where `n_px = prod(image_sampling(1:2))`,
%   and  `c = image_sampling(3)`. `G_xy` produces the image gradient as follows:
%     `gradient_xy = G_xy * I`
%   `I` is a vectorized form of an image where all pixels have been rearranged
%   from columnwise order into a column vector. Specifically, if the original
%   image had a height of `image_sampling(1)`, a width of `image_sampling(2)`,
%   and `c = image_sampling(3)` colour channels or wavelength bands, then `I`
%   contains the data from the image in order first by row, then by column, then
%   by colour channel.
%   `gradient_xy` has twice the number of elements as `I` because it contains
%   the x-gradient at each pixel, followed by the y-gradient at each pixel.
%
% ## Algorithm
%
% The intermediate difference formula is used to produce the image
% gradient. Boundary pixels are replicated to allow the formula to be
% applied at the image boundaries.
%
% ## References
% - The intermediate difference formula for the image gradient is given in
%   the MATLAB help page for the 'imgradientxy()' function.
%   (https://www.mathworks.com/help/images/ref/imgradientxy.html)
%
% See also spectralGradient

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 24, 2018

nargoutchk(1, 1);
narginchk(1, 1);

if length(image_sampling) ~= 3
    error('The `image_sampling` input argument must contain the image height, width, and number of channels/bands');
end

image_height = image_sampling(1);
image_width = image_sampling(2);
c = image_sampling(3);
n_px_c = prod(image_sampling);
% Each element of the intermediate difference gradient is calculated from
% two pixels
n_offsets = 2;
% x-gradients, and y-gradients
n_gradients = 2;

% Create filters to enforce the boundary conditions
% - Omit the last column for the x-gradient, and the last row for the
%   y-gradient, as these pixels get zero values from the intermediate
%   difference formula, according to the boundary conditions (pixel
%   replication).
nonzero_gx = sub2ind(...
    image_sampling,...
    repmat((1:image_height).', (image_width - 1) * c, 1),...
    repmat(repelem((1:(image_width - 1)).', image_height), c, 1),...
    repelem((1:c).', image_height * (image_width - 1))...
    );
nonzero_gy = sub2ind(...
    image_sampling,...
    repmat((1:(image_height - 1)).', image_width * c, 1),...
    repmat(repelem((1:image_width).', (image_height - 1)), c, 1),...
    repelem((1:c).', (image_height - 1) * image_width)...
    );

% Row indices
rows = repelem([
    nonzero_gx;
    nonzero_gy + n_px_c
    ], n_offsets);

% Column indices
col_neighbour_gx = sub2ind(...
    image_sampling,...
    repmat((1:image_height).', (image_width - 1) * c, 1),...
    repmat(repelem((2:image_width).', image_height), c, 1),...
    repelem((1:c).', image_height * (image_width - 1))...
    );
col_neighbour_gy = sub2ind(...
    image_sampling,...
    repmat((2:image_height).', image_width * c, 1),...
    repmat(repelem((1:image_width).', (image_height - 1)), c, 1),...
    repelem((1:c).', (image_height - 1) * image_width)...
    );

n_px_x = length(col_neighbour_gx);
n_px_y = length(col_neighbour_gy);
cols = zeros((n_px_x + n_px_y) * n_offsets, 1);
cols(1:2:(n_px_x * 2)) = nonzero_gx;
cols(2:2:(n_px_x * 2)) = col_neighbour_gx;
cols(((n_px_x * 2) + 1):2:end) = nonzero_gy;
cols(((n_px_x * 2) + 2):2:end) = col_neighbour_gy;

% Matrix values
elements = repmat([-1; 1], n_px_x + n_px_y, 1);

% Assemble the sparse matrix
G_xy = sparse(...
    rows,...
    cols,...
    elements,...
    n_px_c * n_gradients, n_px_c...
);
end
