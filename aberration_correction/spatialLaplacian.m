function [ L ] = spatialLaplacian(image_sampling)
% SPATIALLAPLACIAN  Create a sparse matrix acting as an image spatial Laplacian operator
%
% ## Syntax
% L = spatialLaplacian(image_sampling)
%
% ## Description
% L = spatialLaplacian(image_sampling)
%   Returns a matrix representing the image spatial-domain Laplacian operator.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A three-element vector containing the height, width, and number of colour
%   channels or wavelength bands, respectively, of the image.
%
% ## Output Arguments
%
% L -- Spatial Laplacian matrix
%   A (n_px x c)-by-(n_px x c) array, where `n_px = prod(image_sampling(1:2))`,
%   and  `c = image_sampling(3)`. `L` produces the image spatial Laplacian as
%   follows:
%     `laplacian_xy = L * I`
%   `I` is a vectorized form of an image where all pixels have been rearranged
%   from columnwise order into a column vector. Specifically, if the original
%   image had a height of `image_sampling(1)`, a width of `image_sampling(2)`,
%   and `c = image_sampling(3)` colour channels or wavelength bands, then `I`
%   contains the data from the image in order first by row, then by column, then
%   by colour channel. `laplacian_xy(i)` is the spatial Laplacian of the pixel
%   corresponding to `I(i)`.
%
% ## Algorithm
%
% The image spatial Laplacian is given by the following convolution kernel:
% ```
% [
%    0 -1  0
%   -1  4 -1
%    0 -1  0
% ]
% ```
%
% Boundary pixels are replicated to allow the kernel to be applied to the image
% boundaries.
%
% ## References
%
% The Laplacian energy of an image is a measure of roughness. Bilinear
% interpolation minimizes Laplacian energy:
%
%   Kiku, D., Monno, Y., Tanaka, M, & Okutomi, M. (2016). "Beyond Color
%     Difference: Residual Interpolation for Color Image Demosaicking." IEEE
%     Transactions on Image Processing, 25(3), 1288-1300.
%     doi:10.1109/TIP.2016.2518082
%
% The Laplacian energy is one of the penalty terms used by:
%
%   Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
%
% See also spatialGradient, spectralGradient

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 4, 2019

nargoutchk(1, 1);
narginchk(1, 1);

if length(image_sampling) ~= 3
    error('The `image_sampling` input argument must contain the image height, width, and number of channels/bands');
end

image_height = image_sampling(1);
image_width = image_sampling(2);
c = image_sampling(3);
n_px = image_height * image_width;
n_px_c = prod(image_sampling);
% Each element of the Laplacian is calculated from
% five pixels
n_offsets = 5;

% Handle special cases
if image_height < 2 && image_width < 2
    L = spalloc(n_px_c, n_px_c, 0);
    return;
end
has_bottom = 1;
has_right = 1;
if image_height == 1
    has_bottom = 0;
end
if image_width == 1
    has_right = 0;
end
h1 = image_height - 1;
h2 = max(image_height - 2, 0);
w2 = max(image_width - 2, 0);
has_center = image_width > 2;

n_values_per_channel = h2 * w2 * n_offsets +... % Centre, non-edge case
    1 + has_bottom + has_right +... % Top left
    has_bottom * (2 + has_right) +... % Bottom left
    has_right * (2 + has_bottom) +... % Top right
    has_bottom * has_right * 3 +... % Bottom right
    h2 * (3 + has_right) * (1 + has_right) +... % Left and right
    w2 * (3 + has_bottom) * (1 + has_bottom); % Top and bottom
    
rows = repmat([
    ones(2 + has_right * has_bottom, 1); % Top left pixel
    repelem((2:h1).', 3 + has_right, 1); % First column
    repelem(image_height, has_bottom * 2 + has_right * has_bottom * 1, 1); % Bottom left pixel
    repmat(...
        [
            ones(3 + has_bottom, 1); % Top pixel
            repelem((2:h1).', 5, 1); % Middle pixels
            repelem(image_height, has_bottom * 4, 1); % Bottom pixel
        ], has_center * w2, 1 ...
    ) + repelem(...
        (1:w2).' * image_height,...
        has_center * (3 + has_bottom + h2 * 5 + has_bottom * 4), 1 ...
    );
    repelem(n_px - h1, has_right * 2 + has_right * has_bottom, 1); % Top right pixel
    repmat(repelem(((n_px - h2):(n_px - 1)).', 4, 1), has_right, 1); % Last column
    repelem(n_px, has_bottom * has_right * 3, 1) % Bottom right pixel
], c, 1) + repelem((0:n_px:(n_px * (c - 1))).', n_values_per_channel, 1); % Replicate for all channels

% For the filter kernel
% ```
% [
%    0 -1  0
%   -1  4 -1
%    0 -1  0
% ]
% ```
% let the offsets be numbered as follows:
% ```
% [
%    0  2  0
%    1  3  5
%    0  4  0
% ]
% ```

columns = repmat([
    [1; repmat(2, has_bottom, 1); repmat(image_height + 1, has_right, 1)]; % Top left pixel
    repmat([-1; 0; 1; repmat(image_height, has_right, 1)], h2, 1) + repelem((2:h1).', 3 + has_right, 1); % First column
    repmat([h1; image_height; repmat(image_height * 2, has_right, 1)], has_bottom, 1); % Bottom left pixel
    repmat(...
        [
            [-image_height; 0; ones(has_bottom, 1); image_height]; % Top pixel
            repmat([-image_height; -1; 0; 1; image_height], h2, 1) + repelem((1:h2).', 5, 1); % Middle pixels
            repmat([-image_height; -1; 0; image_height] + h1, has_bottom, 1) % Bottom pixel
        ], has_center * w2, 1 ...
    ) + repelem((1:w2).' * image_height + 1, has_center * (3 + has_bottom + h2 * 5 + has_bottom * 4), 1);
    repmat([n_px - image_height * 2 + 1; n_px - h1; repmat(n_px - h2, has_bottom, 1)], has_right, 1); % Top right pixel
    repmat([-image_height; -1; 0; 1], has_right * h2, 1) + repelem(((n_px - h2):(n_px - 1)).', has_right * 4, 1); % Last column
    repmat([n_px - image_height; n_px - 1; n_px], has_bottom * has_right, 1) % Bottom right pixel
], c, 1) + repelem((0:n_px:(n_px * (c - 1))).', n_values_per_channel, 1); % Replicate for all channels

values = repmat([
    [1 + has_right * has_bottom; repmat(-1, has_right, 1); repmat(-1, has_bottom, 1)]; % Top left pixel
    repmat([-1; 2 + has_right; -1; repmat(-1, has_right, 1)], h2, 1); % First column
    repmat([-1; 1 + has_right; repmat(-1, has_right, 1)], has_bottom, 1); % Bottom left pixel
    repmat(...
        [
            [-1; 2 + has_bottom; repmat(-1, has_bottom, 1); -1]; % Top pixel
            repmat([-1; -1; 4; -1; -1], h2, 1); % Middle pixels
            repmat([-1; -1; 3; -1], has_bottom, 1); % Bottom pixel
        ], has_center * w2, 1 ...
    );
    repmat([-1; 1 + has_bottom; repmat(-1, has_bottom, 1)], has_right, 1); % Top right pixel
    repmat([-1; -1; 3; -1], has_right * h2, 1); % Last column
    repmat([-1; -1; 2], has_bottom * has_right, 1) % Bottom right pixel
], c, 1); % Replicate for all channels

% Assemble the sparse matrix
L = sparse(...
    rows,...
    columns,...
    values,...
    n_px_c, n_px_c...
);

end