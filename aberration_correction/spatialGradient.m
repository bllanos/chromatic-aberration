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
% The central difference formula is used to produce the image gradient. Boundary
% pixels are replicated to allow the formula to be applied at the image boundaries.
%
% ## References
% - The central difference formula for the image gradient is given in the MATLAB
%   help page for the 'imgradientxy()' function.
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
n_px = prod(image_sampling(1:2));
n_px_c = prod(image_sampling);
% Each element of the central difference gradient is calculated from two pixels
n_offsets = 2;
% x-gradients, and y-gradients
n_gradients = 2;

% Row indices
% - Go through all pixels once per gradient direction
% - Within each gradient direction, create two matrix elements per pixel
rows = repmat(repelem((1:n_px_c).', n_offsets), n_gradients, 1);

% Column indices

% Enumerate the positions of all pixels in the image
[X, Y] = meshgrid(1:image_width, 1:image_height);
x = repmat(X(:), c, 1);
y = repmat(Y(:), c, 1);

neighbour_x = zeros(n_px_c, n_offsets, n_gradients);
neighbour_y = zeros(n_px_c, n_offsets, n_gradients);
for g = 1:n_gradients
    if g == 1
        offsets = [-1 0; 1 0];
    else
        offsets = [0 -1; 0 1];
    end
    for i = 1:n_offsets
        neighbour_x(:, i, g) = x + offsets(i, 1);
        neighbour_y(:, i, g) = y + offsets(i, 2);
    end
end
neighbour_x = neighbour_x(:);
neighbour_y = neighbour_y(:);
neighbour_c = repmat(repelem((1:c).', n_px), n_gradients, 1);

% Replicate pixels outside the image boundaries
neighbour_x(neighbour_x < 1) = 1;
neighbour_x(neighbour_x > image_width) = image_width;
neighbour_y(neighbour_y < 1) = 1;
neighbour_y(neighbour_y > image_height) = image_height;

% Convert to linear indices
neighbour_index = sub2ind(...
    image_sampling,...
    neighbour_y,...
    neighbour_x,...
    neighbour_c...
);

% Matrix values
elements = repmat([-0.5; 0.5], n_px_c * n_gradients);

% Assemble the sparse matrix
G_xy = sparse(...
    rows,...
    neighbour_index,...
    elements,...
    n_px_c * n_gradients, n_px_c...
);
end
