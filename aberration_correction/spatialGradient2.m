function [ G_xy2, G_diag2 ] = spatialGradient2(image_sampling)
% SPATIALGRADIENT2  Create sparse matrices acting as image second-order spatial gradient operators
%
% ## Syntax
% G_xy2 = spatialGradient2(image_sampling)
% [G_xy2, G_diag2] = spatialGradient2(image_sampling)
%
% ## Description
% G_xy2 = spatialGradient2(image_sampling)
%   Returns a matrix representing the second-order spatial gradient
%   operators.
% [G_xy2, G_diag2] = spatialGradient2(image_sampling)
%   Additionally returns a matrix representing the second-order diagonal
%   spatial gradient operators
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A three-element vector containing the height, width, and number of colour
%   channels or wavelength bands, respectively, of the image.
%
% ## Output Arguments
%
% G_xy2 -- Second-order spatial gradients matrix
%   A (n_px x c x 2)-by-(n_px x c) array, where `n_px =
%   prod(image_sampling(1:2))`, and  `c = image_sampling(3)`. `G_xy2`
%   produces the image second-order x- and y-gradients as follows:
%     `gradient_xy2 = G_xy2 * I`
%   `I` is a vectorized form of an image where all pixels have been
%   rearranged from columnwise order into a column vector. Specifically, if
%   the original image had a height of `h`, a width of
%   `w`, and `c = image_sampling(3)` colour channels or
%   wavelength bands, then `I` contains the data from the image in order
%   first by row, then by column, then by colour channel. `gradient_xy2`
%   has twice the number of elements as `I` because it contains the
%   centered difference second-order x-gradient at each pixel, followed by
%   the centered difference second-order y-gradient at each pixel.
%
% G_diag2 -- Second order spatial gradients matrix for diagonal directions
%   A matrix similar to `G_xy2`, except that it produces second-order
%   gradients taken along the diagonal directions in the image. The result
%   of multiplying a vectorized image with this matrix is a vector
%   containing the gradient in the top right direction for every pixel,
%   followed by the gradient in the bottom right direction for every pixel.
%
%   Note that the derivatives in the diagonal directions are computed from
%   pixel differences, rather than derived from the image x- and
%   y-gradients.
%
% ## Algorithm
%
% The centered difference formula is used to produce the image gradients.
% At the image boundaries, the center pixel is substituted for the values
% of pixels outside the image boundaries, resulting in first order x- and
% y-gradients, and scaled first order diagonal gradients, as opposed to
% second order gradients.
%
% The second-order x-gradient kernel is: [-1 2 -1]
% The second-order y-gradient kernel is: [-1; 2; -1]
% The second-order top-right diagonal gradient kernel is:
%   [0 0 -0.5; 0 1 0; -0.5 0 0]
% The second-order bottom-right diagonal gradient kernel is:
%   [-0.5 0 0; 0 1 0; 0 0 -0.5]
%
% See also spatialGradient, spectralGradient

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 12, 2018

nargoutchk(1, 2);
narginchk(1, 1);

if length(image_sampling) ~= 3
    error('The `image_sampling` input argument must contain the image height, width, and number of channels/bands');
end

h = image_sampling(1);
w = image_sampling(2);
c = image_sampling(3);
n_px = h * w;
n_px_c = prod(image_sampling);
n_px_cNeg1 = n_px_c - n_px;
% x-gradients, and y-gradients, or two diagonal gradients
n_gradients = 2;
% Spatial indices of pixels in the image
indices = reshape(1:n_px, h, w);

% Second-order x-gradient
rows_gx2 = [
    repelem(indices(:, 1), 2);
    reshape(repelem(indices(:, 2:(end - 1)), 3, 1), [], 1);
    repelem(indices(:, end), 2)
];
col_gx2 = [
    repmat([0; h], h, 1);
    repmat([-h; 0; h], h * (w - 2), 1);
    repmat([-h; 0], h, 1)
] + rows_gx2;
channel_offsets = repelem((0:n_px:n_px_cNeg1).', length(col_gx2), 1);
col_gx2 = repmat(col_gx2, c, 1) + channel_offsets;
rows_gx2 = repmat(rows_gx2, c, 1) + channel_offsets;

elements_gx2 = repmat([
    repmat([1; -1], h, 1);
    repmat([-1; 2; -1], h * (w - 2), 1);
    repmat([-1; 1], h, 1)
    ], c, 1);

% Second-order y-gradient
rows_gy2 = reshape([
    repmat(indices(1, :), 2, 1);
    repelem(indices(2:(end-1), :), 3, 1);
    repmat(indices(end, :), 2, 1)
    ], [], 1 ...
);
col_gy2 = repmat([0; 1; repmat([-1; 0; 1], h - 2, 1); -1; 0], w, 1) + rows_gy2;
channel_offsets = repelem((0:n_px:n_px_cNeg1).', length(col_gy2), 1);
col_gy2 = repmat(col_gy2, c, 1) + channel_offsets;
rows_gy2 = repmat(rows_gy2, c, 1) + channel_offsets;

elements_gy2 = repmat([1; -1; repmat([-1; 2; -1], h - 2, 1); -1; 1], w * c, 1);

% Assemble the sparse matrix
G_xy2 = sparse(...
    [rows_gx2; (rows_gy2 + n_px_c)],...
    [col_gx2; col_gy2],...
    [elements_gx2; elements_gy2],...
    n_px_c * n_gradients, n_px_c...
);

% Compute diagonal gradients
if nargout > 1
    
    % Up-right second-order diagonal gradient
    rows_g12 = repmat(indices, 1, 1, 3);
    col_g12 = rows_g12 + repmat(reshape([-h + 1; 0; h - 1], 1, 1, 3), h, w, 1);
    col_g12(:, 1, 1) = 0;
    col_g12(1, :, 3) = 0;
    col_g12(end, :, 1) = 0;
    col_g12(:, end, 3) = 0;
    
    col_g12 = reshape(permute(col_g12, [3, 1, 2]), [], 1);
    col_g12_filter = (col_g12 ~= 0);
    col_g12 = col_g12(col_g12_filter);
    channel_offsets = repelem((0:n_px:n_px_cNeg1).', length(col_g12), 1);
    col_g12 = repmat(col_g12, c, 1) + channel_offsets;
    
    rows_g12 = reshape(permute(rows_g12, [3, 1, 2]), [], 1);
    rows_g12 = repmat(rows_g12(col_g12_filter), c, 1) + channel_offsets;
    
    elements_g12 = [
        0; repmat([0.5; -0.5], h - 1, 1); % First column
        repmat([-0.5; 0.5; repmat([-0.5; 1; -0.5], h - 2, 1); 0.5; -0.5], w - 2, 1); % Middle columns
        repmat([-0.5; 0.5], h - 1, 1); 0 % Last column
    ];
    elements_g12 = repmat(elements_g12, c, 1);
    
    % Down-right second-order diagonal gradient
    rows_g22 = repmat(indices, 1, 1, 3);
    col_g22 = rows_g22 + repmat(reshape([-h - 1; 0; h + 1], 1, 1, 3), h, w, 1);
    col_g22(:, 1, 1) = 0;
    col_g22(1, :, 1) = 0;
    col_g22(end, :, 3) = 0;
    col_g22(:, end, 3) = 0;
    
    col_g22 = reshape(permute(col_g22, [3, 1, 2]), [], 1);
    col_g22_filter = (col_g22 ~= 0);
    col_g22 = col_g22(col_g22_filter);
    channel_offsets = repelem((0:n_px:n_px_cNeg1).', length(col_g22), 1);
    col_g22 = repmat(col_g22, c, 1) + channel_offsets;
    
    rows_g22 = reshape(permute(rows_g22, [3, 1, 2]), [], 1);
    rows_g22 = repmat(rows_g22(col_g22_filter), c, 1) + channel_offsets;
    
    elements_g22 = [
        repmat([0.5; -0.5], h - 1, 1); 0; % First column
        repmat([0.5; -0.5; repmat([-0.5; 1; -0.5], h - 2, 1); -0.5; 0.5], w - 2, 1); % Middle columns
        0; repmat([-0.5; 0.5], h - 1, 1) % Last column
    ];
    elements_g22 = repmat(elements_g22, c, 1);
    
    % Assemble the sparse matrix
    G_diag2 = sparse(...
        [rows_g12; (rows_g22 + n_px_c)],...
        [col_g12; col_g22],...
        [elements_g12; elements_g22],...
        n_px_c * n_gradients, n_px_c...
    );
    
end


end