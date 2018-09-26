function [ G_xy, G_diag ] = spatialGradient(image_sampling)
% SPATIALGRADIENT  Create sparse matrices acting as image spatial gradient operators
%
% ## Syntax
% G_xy = spatialGradient(image_sampling)
% [G_xy, G_diag] = spatialGradient(image_sampling)
%
% ## Description
% G_xy = spatialGradient(image_sampling)
%   Returns a matrix representing the spatial gradient operators.
% [G_xy, G_diag] = spatialGradient(image_sampling)
%   Additionally returns a matrix representing the diagonal spatial
%   gradient operators
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
% G_diag -- Spatial gradient matrix for diagonal directions
%   A matrix similar to `G_xy`, except that it produces gradients taken
%   along the diagonal directions in the image. The result of multiplying a
%   vectorized image with this matrix is a vector containing the gradient
%   in the top right direction for every pixel, followed by the gradient in
%   the bottom right direction for every pixel.
%
%   Note that the derivatives in the diagonal directions are computed from
%   pixel differences, rather than derived from the image x- and
%   y-gradients.
%
% ## Algorithm
%
% The intermediate difference formula is used to produce the image
% gradients. Boundary pixels are replicated to allow the formula to be
% applied at the image boundaries.
%
% ## References
% - The intermediate difference formula for the image gradient is given in
%   the MATLAB help page for the 'imgradientxy()' function.
%   (https://www.mathworks.com/help/images/ref/imgradientxy.html)
%
% See also spatialGradient2, spectralGradient

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 24, 2018

nargoutchk(1, 2);
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
% x-gradients, and y-gradients, or two diagonal gradients
n_gradients = 2;

% Create filters to enforce the boundary conditions
% - Omit the last column for the x-gradient, and the last row for the
%   y-gradient, as these pixels get zero values from the intermediate
%   difference formula, according to the boundary conditions (pixel
%   replication).
nonzero_gx = sub2ind(...
    image_sampling,...
    repmat((1:image_height).', (image_width - 1) * c, 1),...
    repmat(repelem((1:(image_width - 1)), image_height).', c, 1),...
    repelem((1:c), image_height * (image_width - 1)).'...
    );
nonzero_gy = sub2ind(...
    image_sampling,...
    repmat((1:(image_height - 1)).', image_width * c, 1),...
    repmat(repelem((1:image_width), (image_height - 1)).', c, 1),...
    repelem((1:c), (image_height - 1) * image_width).'...
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
    repmat(repelem((2:image_width), image_height).', c, 1),...
    repelem((1:c), image_height * (image_width - 1)).'...
    );
col_neighbour_gy = sub2ind(...
    image_sampling,...
    repmat((2:image_height).', image_width * c, 1),...
    repmat(repelem((1:image_width), (image_height - 1)).', c, 1),...
    repelem((1:c), (image_height - 1) * image_width).'...
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

% Compute diagonal gradients
if nargout > 1
    % Create filters to enforce the boundary conditions
    % - Omit the last column and first row for the up-right-gradient, and
    %   the last column and last row for the down-right-gradient, as these
    %   pixels get zero values from the intermediate difference formula,
    %   according to the boundary conditions (pixel replication).
    nonzero_g1 = sub2ind(...
        image_sampling,...
        repmat((2:image_height).', (image_width - 1) * c, 1),...
        repmat(repelem((1:(image_width - 1)).', (image_height - 1)), c, 1),...
        repelem((1:c).', (image_height - 1) * (image_width - 1))...
        );
    nonzero_g2 = sub2ind(...
        image_sampling,...
        repmat((1:(image_height - 1)).', (image_width - 1) * c, 1),...
        repmat(repelem((1:(image_width - 1)).', (image_height - 1)), c, 1),...
        repelem((1:c).', (image_height - 1) * (image_width - 1))...
        );

    % Row indices
    rows = repelem([
        nonzero_g1;
        nonzero_g2 + n_px_c
        ], n_offsets);

    % Column indices
    col_neighbour_g1 = sub2ind(...
        image_sampling,...
        repmat((1:(image_height-1)).', (image_width - 1) * c, 1),...
        repmat(repelem((2:image_width).', (image_height - 1)), c, 1),...
        repelem((1:c).', (image_height - 1) * (image_width - 1))...
        );
    col_neighbour_g2 = sub2ind(...
        image_sampling,...
        repmat((2:image_height).', (image_width - 1) * c, 1),...
        repmat(repelem((2:image_width).', (image_height - 1)), c, 1),...
        repelem((1:c).', (image_height - 1) * (image_width - 1))...
        );

    n_px_1 = length(col_neighbour_g1);
    n_px_2 = length(col_neighbour_g2);
    cols = zeros((n_px_1 + n_px_2) * n_offsets, 1);
    cols(1:2:(n_px_1 * 2)) = nonzero_g1;
    cols(2:2:(n_px_1 * 2)) = col_neighbour_g1;
    cols(((n_px_1 * 2) + 1):2:end) = nonzero_g2;
    cols(((n_px_1 * 2) + 2):2:end) = col_neighbour_g2;

    % Matrix values
    % Divide by `sqrt(2)` to account for the larger diagonal step size.
    elements = repmat([-1; 1] ./ sqrt(2), n_px_1 + n_px_2, 1);

    % Assemble the sparse matrix
    G_diag = sparse(...
        rows,...
        cols,...
        elements,...
        n_px_c * n_gradients, n_px_c...
    );
end
end
