function [ I_rgb ] = bilinearDemosaic(I_raw, align, varargin)
% BILINEARDEMOSAIC  Demosaic an image by bilinear interpolation
%
% ## Syntax
% I_rgb = bilinearDemosaic(I_raw, align [, channels])
%
% ## Description
% I_rgb = bilinearDemosaic(I_raw, align [, channels])
%   Converts a colour-filter array image to a colour image by bilinear
%   interpolation.
%
% ## Input Arguments
%
% I_raw -- RAW image
%   A image_height x image_width array storing the raw colour filter array
%   data of an image.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'. `align` has the same form as the `sensorAlignment`
%   input argument of `demosaic()`.
%
% channels -- Requested channels
%   A 3-element logical vector, where `channels(i)` is a flag indicating
%   whether or not to produce the i-th demosaiced colour channel in
%   `I_rgb`. If empty or not passed, `channels` defaults to `true(3, 1)`.
%
% ## Output Arguments
%
% I_rgb -- Full-colour image
%   An image_height x image_width x sum(channels) array, containing the
%   demosaicked image, produced by per-channel bilinear interpolation of
%   `I_raw`. Only the colour channels corresponding to `true` values in
%   `channels` will be output. Consequently, while the order of the colour
%   channels in `I_rgb` is always Red, Green, Blue, the indices of the
%   channels depend on the contents of `channels`. For instance, `Blue`
%   will be the first channel if it is the only channel requested.
%
% ## Notes
% - Pixels at the edges of the image will be given the values of their
%   neighbours for colour channels whose pixels do not extend to the
%   corresponding edges, or which omit image corners.
% - The image dimensions are assumed to be even numbers.
%
% ## References
% - Bilinear interpolation is mentioned as a rudimentary demosaicking method
%   in D. Menon and G. Calvagno. "Color image demosaicking: An overview,"
%   Signal Processing: Image Communication, vol. 26, pp. 518-533, 2011.
% - Bilinear interpolation formulae:
%   https://en.wikipedia.org/wiki/Bilinear_interpolation
%
% See also demosaic, interp2, interp1

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 20, 2018

nargoutchk(1, 1);
narginchk(2, 3);

n_channels_max = 3;
if ~isempty(varargin)
    channels = varargin{1};
    if length(channels) ~= n_channels_max
        error('`channels` must have three elements, for Red, Green, and Blue.');
    end
else
    channels = true(n_channels_max, 1);
end

I_raw = im2double(I_raw);

image_height = size(I_raw, 1);
image_width = size(I_raw, 2);
mask = bayerMask( image_height, image_width, align );

n_channels = sum(channels);
I_rgb = nan(image_height, image_width, n_channels);

% Red and Blue
[X, Y] = meshgrid(1:image_width, 1:image_height);
for c = [1 n_channels_max]
    if channels(c)
        mask_c = mask(:, :, c);
        I_raw_c = reshape(I_raw(mask_c), image_height / 2, image_width / 2);
        if image_height <= 2 || image_width <= 2
            I_c = repelem(I_raw_c, 2, 2);
        else
            x_c = reshape(X(mask_c), image_height / 2, image_width / 2);
            y_c = reshape(Y(mask_c), image_height / 2, image_width / 2);
            I_c = interp2(x_c, y_c, I_raw_c, X, Y, 'linear');
        end
        if c == 1
            I_rgb(:, :, c) = I_c;
        else
            I_rgb(:, :, n_channels) = I_c;
        end
    end
end

% Green
green_index = 2;
if channels(green_index)
    mask_g = mask(:, :, green_index);
    I_g = nan(image_height, image_width);
    I_g(mask_g) = I_raw(mask_g);
    center_mask_g = ~mask_g;
    center_mask_g(1, :) = false;
    center_mask_g(end, :) = false;
    center_mask_g(:, 1) = false;
    center_mask_g(:, end) = false;
    missing_ind_g = sub2ind([image_height, image_width], Y(center_mask_g), X(center_mask_g));
    I_g(missing_ind_g) = (...
        I_g(missing_ind_g + 1) +...
        I_g(missing_ind_g - 1) +...
        I_g(missing_ind_g + image_height) +...
        I_g(missing_ind_g - image_height)...
        ) / 4; % Central region of image
    
    if image_height <= 2
        I_g(:, 1) = I_g(mask_g(:, 1), 1);
        I_g(:, end) = I_g(mask_g(:, end), end);
    else
        I_g(:, 1) = interp1(Y(mask_g(:, 1), 1), I_g(mask_g(:, 1), 1), Y(:, 1), 'linear'); % Left
        I_g(:, end) = interp1(Y(mask_g(:, end), end), I_g(mask_g(:, end), end), Y(:, end), 'linear'); % Right
    end
    if image_width <= 2
        I_g(1, :) = I_g(1, mask_g(1, :));
        I_g(end, :) = I_g(end, mask_g(end, :));
    else
        I_g(1, :) = interp1(X(1, mask_g(1, :)), I_g(1, mask_g(1, :)), X(1, :), 'linear'); % Top
        I_g(end, :) = interp1(X(end, mask_g(end, :)), I_g(end, mask_g(end, :)), X(end, :), 'linear'); % Bottom
    end

    if channels(green_index - 1)
        I_rgb(:, :, green_index) = I_g;
    else
        I_rgb(:, :, green_index - 1) = I_g;
    end
end

% Fill in boundary values by copying adjacent pixels
ind_unknown = find(~isfinite(I_rgb));
[row,col,channel] = ind2sub(size(I_rgb), ind_unknown);
row_adj = row;
col_adj = col;
% Left
col_adj(col == 1) = col(col == 1) + 1;
% Top
row_adj(row == 1) = row(row == 1) + 1;
% Bottom
row_adj(row == image_height) = row(row == image_height) - 1;
% Right
col_adj(col == image_width) = col(col == image_width) - 1;

ind_adj = sub2ind(size(I_rgb), row_adj, col_adj, channel);
I_rgb(ind_unknown) = I_rgb(ind_adj);

end