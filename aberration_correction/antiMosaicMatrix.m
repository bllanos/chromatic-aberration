function [ M ] = antiMosaicMatrix(image_sampling, align)
% ANTIMOSAICMATRIX  Create a sparse matrix to penalize colour-filter array-like images
%
% ## Syntax
% M = antiMosaicMatrix(image_sampling, align)
%
% ## Description
% M = antiMosaicMatrix(image_sampling, align)
%   Returns a matrix for quantifying how similar an image is to a raw
%   colour-filter array pattern.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the image height and width,
%   respectively. The image dimensions must be even integers.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'. `align` has the same form as the `sensorAlignment`
%   input argument of `demosaic()`.
%
% ## Output Arguments
%
% M -- Penalty matrix
%   A (n_px x 3) x (n_px * 3) array, where `n_px = prod(image_sampling)`.
%   `M` produces penalty values for each pixel in a 3-channel (RGB) image
%   that is not measured in the corresponding RAW (color-filter array)
%   image:
%     `P = M * I`
%   `I` is a vectorized form of an image where all pixels have been
%   rearranged from columnwise order into a column vector. Specifically, if
%   the original image had a height of `image_sampling(1)`, a width of
%   `image_sampling(2)`, and 3 colour channels, then `I` contains the data
%   from the image in order first by row, then by column, then by colour
%   channel. `P` is a vectorized form of the penalties on the image, with
%   the penalties on specific colour channels in the same order as the
%   colour channels are given in the image `I`.
%
% ## Algorithm
% In a region of constant colour, a colour-filter array image exhibits
% periodic variations according to its colour filter pattern. Such
% variations can be distinguished from edges using second-order image
% gradients: Edges correspond to high first-order image gradients, but not
% necessarily to high second-order image gradients, whereas colour filter
% patterns produce very high second-order image gradients.
%
% This function produces a matrix which measures second-order image
% gradients, but which only measures gradients between the appropriate
% colour channel values of the appropriate pixels. The following twelve
% gradients are measured, and are output as `P` in the the order given:
% - Red channel gradients:
%   - At pixels corresponding to Green filters in the colour-filter array:
%     - Second order derivative in the horizontal or vertical direction,
%       depending on the location of the pixel. (Given a double weight.)
%   - At pixels corresponding to Blue filters in the colour-filter array:
%     - Second order derivative in one diagonal direction (down-left to
%       up-right).
%     - Second order derivative in the other diagonal direction (up-left to
%       down-right).
% - Green channel gradients:
%   - At pixels corresponding to Red filters in the colour-filter array:
%     - Second order derivative in the horizontal direction
%     - Second order derivative in the vertical direction
%   - At pixels corresponding to Blue filters in the colour-filter array:
%     - Second order derivative in the horizontal direction
%     - Second order derivative in the vertical direction
% - Blue channel gradients:
%   - At pixels corresponding to Red filters in the colour-filter array:
%     - Second order derivative in one diagonal direction (down-left to
%       up-right).
%     - Second order derivative in the other diagonal direction (up-left to
%       down-right).
%   - At pixels corresponding to Green filters in the colour-filter array:
%     - Second order derivative in the horizontal or vertical direction,
%       depending on the location of the pixel. (Given a double weight.)
%
% See also mosaic, mosaicMatrix, bayerMask, spatialGradient2

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 10, 2018

nargoutchk(1, 1);
narginchk(2, 2);

if length(image_sampling) ~= 2
    error('The `image_sampling` input argument must contain an image height and width only.');
elseif any(mod(image_sampling, 2) ~= 0)
    error('The image dimensions must be even integers in order for the image to be a valid color filter array.');
end

% Does the first Red value at the first Green pixel use a horizontal or
% vertical gradient? (Note that pixels are stored in column-major order in
% MATLAB.)
horizontal_red = false;

% Enumerate the row and column indices of pixels in the colour filter array
% pattern. The order is Red, first Green, second Green, Blue:
rows = {
    1:2:image_sampling(1);
    2:2:image_sampling(1);
    1:2:image_sampling(1);
    2:2:image_sampling(1)
};
cols = {
    1:2:image_sampling(2);
    1:2:image_sampling(2);
    2:2:image_sampling(2);
    2:2:image_sampling(2)
};
switch align
    case 'rggb'
    case 'bggr'
        horizontal_red = true;
        rows = rows([4, 2, 3, 1]);
        cols = cols([4, 2, 3, 1]);
    case 'grbg'
        horizontal_red = true;
        rows = rows([3, 1, 4, 2]);
        cols = cols([3, 1, 4, 2]);
    case 'gbrg'
        rows = rows([2, 1, 4, 3]);
        cols = cols([2, 1, 4, 3]);
    otherwise
        error('Unrecognized Bayer pattern format string.');
end

% Convert subscripts to indices
rows_full = cell(length(rows), 1);
cols_full = cell(length(cols), 1);
for i = 1:length(rows)
    [cols_full{i}, rows_full{i}] = meshgrid(cols{i}, rows{i});
    rows_full{i} = reshape(rows_full{i}, [], 1);
    cols_full{i} = reshape(cols_full{i}, [], 1);
end
n_channels = 3;
image_sampling3 = [image_sampling n_channels];
n_px = prod(image_sampling);
n_px_c = n_px * n_channels;
n_px_div4 = n_px / 4;
red_at_first_green = sub2ind(...
    image_sampling3, rows_full{2}, cols_full{2}, repelem(1, n_px_div4, 1)...
);
red_at_second_green = sub2ind(...
    image_sampling3, rows_full{3}, cols_full{3}, repelem(1, n_px_div4, 1)...
);
red_at_blue = sub2ind(...
    image_sampling3, rows_full{4}, cols_full{4}, repelem(1, n_px_div4, 1)...
);
green_at_red = sub2ind(...
    image_sampling3, rows_full{1}, cols_full{1}, repelem(2, n_px_div4, 1)...
);
green_at_blue = sub2ind(...
    image_sampling3, rows_full{4}, cols_full{4}, repelem(2, n_px_div4, 1)...
);
blue_at_red = sub2ind(...
    image_sampling3, rows_full{1}, cols_full{1}, repelem(3, n_px_div4, 1)...
);
blue_at_first_green = sub2ind(...
    image_sampling3, rows_full{2}, cols_full{2}, repelem(3, n_px_div4, 1)...
);
blue_at_second_green = sub2ind(...
    image_sampling3, rows_full{3}, cols_full{3}, repelem(3, n_px_div4, 1)...
);

% Compute the image second order derivatives
[ G_xy2, G_diag2 ] = spatialGradient2(image_sampling3);
G_x2 = G_xy2(1:n_px_c, :);
G_y2 = G_xy2((n_px_c + 1):end, :);
G_12 = G_diag2(1:n_px_c, :);
G_22 = G_diag2((n_px_c + 1):end, :);

% Select the appropriate second order derivatives for the penalty
% generation matrix
RB_weight = 2;
if horizontal_red
    M_RG = [
        G_x2(red_at_first_green, :);
        G_y2(red_at_second_green, :)
    ];
    M_BG = [
        G_y2(blue_at_first_green, :);
        G_x2(blue_at_second_green, :)
    ];
else
    M_RG = [
        G_y2(red_at_first_green, :);
        G_x2(red_at_second_green, :)
    ];
    M_BG = [
        G_x2(blue_at_first_green, :);
        G_y2(blue_at_second_green, :)
    ];
end
M = [
    RB_weight * M_RG;
    G_12(red_at_blue, :);
    G_22(red_at_blue, :);
    G_x2(green_at_red, :);
    G_y2(green_at_red, :);
    G_x2(green_at_blue, :);
    G_y2(green_at_blue, :);
    G_12(blue_at_red, :);
    G_22(blue_at_red, :);
    RB_weight * M_BG
];

end

