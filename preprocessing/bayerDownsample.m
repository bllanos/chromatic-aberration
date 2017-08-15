function [...
    I_channels, downsample_map, upsample_map, upsample_map_centers...
] = bayerDownsample( I_raw, align )
% BAYERDOWNSAMPLE  Create subresolution image from colour-filter array data
%
% ## Syntax
% I_channels = bayerDownsample( I_raw, align )
% [ I_channels, downsample_map ] = bayerDownsample( I_raw, align )
% [ I_channels, downsample_map, upsample_map ] = bayerDownsample( I_raw, align )
% [...
%   I_channels, downsample_map, upsample_map, upsample_map_centers...
% ] = bayerDownsample( I_raw, align )
%
% ## Description
% I_channels = bayerDownsample( I_raw, align )
%   Returns a half-size RGB image created by downsampling the input raw
%   image.
% [ I_channels, downsample_map ] = bayerDownsample( I_raw, align )
%   Additionally returns a map from pixels in the raw image to pixels in
%   the output image.
% [ I_channels, downsample_map, upsample_map ] = bayerDownsample( I_raw, align )
%   Additionally returns a map from pixels in the output image to pixels in
%   the raw image.
% [...
%   I_channels, downsample_map, upsample_map, upsample_map_centers...
% ] = bayerDownsample( I_raw, align )
%   Additionally returns a map from pixel centers in the output image to
%   coordinates in the raw image.
%
% ## Input Arguments
%
% I_raw -- Colour-filter array data
%   An image_height x image_width array, containing colour filter array
%   data (Bayer pattern).
%
% align -- Bayer pattern format
%   A four-character character vector, specifying the Bayer tile pattern
%   corresponding to `I_raw`. For example, 'gbrg'. The function assumes
%   that there are two Green pixels per tile, and that they do not share
%   the same column.
%
%   This argument has the same form as the `sensorAlignment` input argument
%   of `demosaic()`.
%
% ## Output Arguments
%
% I_channels -- Full-colour downsampled image
%   An (image_height / 2) x (image_width / 2) x 3 array, containing a
%   downsampled version of `I_raw`, separated into colour channels. Each
%   pixel in `I_channels` corresponds to a Bayer tile in `I_raw`.
%   `I_channels(i,j,1)` is the Red value from the (i,j)-th Bayer tile, and
%   `I_channels(i,j,3)` is the Blue value from the (i,j)-th Bayer tile.
%   `I_channels(i,j,2)` is the average of the two Green values from the
%   (i,j)-th Bayer tile.
%
% downsample_map -- Downsampling pixel coordinate conversion
%   An image_height x image_width array. `downsample_map(i,j)` is the
%   linear index of a colour channel value in `I_channels`, corresponding
%   to the pixel `I_raw(i,j)`. For example, if `I_raw(i,j)` represents a
%   Red pixel, then this pixel maps to `I_channels(downsample_map(i,j))`,
%   where `downsample_map(i,j)` will correspond to an index into
%   `I_channels(:,:,1)`.
%
% upsample_map -- Upsampling pixel coordinate conversion
%   An image_height x image_width x 4 array. `upsample_map(i,j,:)`
%   contains the linear indices of the pixels in `I_raw` corresponding to
%   the pixel `I_channels(i,j)`. Each pixel in `I_channels` maps to a Bayer
%   tile (4 pixels) in `I_raw`. `upsample_map(i,j,k)` maps to a Red, Green,
%   Green, or Blue pixel in `I_raw`, for `k` from one to four,
%   respectively.
%
% upsample_map_centers -- Coordinate frame conversion
%   An image_height x image_width array x 2 array, where
%   `upsample_map_centers(i, j, :)` is the `(x, y)` position in `I_raw`
%   corresponding to the centre of the pixel `I_channels(i,j,:)`.
%
% ## Notes
% - `upsample_map_centers` is technically valid only for the Green channel
%   of the output image. In the raw image, Red and Blue channels are
%   sampled at locations offset by a pixel in both the horizontal and
%   vertical directions (the signs of the displacements being opposite for
%   Red vs. Blue) from the average locations of the two Green pixels in the
%   Bayer tiles. This is a design decision: A uniform bias in chromatic
%   aberration was thought to be preferable to interpolation artifacts
%   introduced by resampling the Red and Blue channels at the averaged
%   positions of the pixels in the Green channel.
%
% See also bayerMask, demosaic

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 10, 2017

nargoutchk(1, 4);
narginchk(2, 2);

image_height = size(I_raw, 1);
image_width = size(I_raw, 2);

if mod(image_height, 2) || mod(image_width, 2)
    error('A valid Bayer pattern image has dimensions which are even integers.')
end

mask = bayerMask( image_height, image_width, align );

image_height2 = image_height / 2;
image_width2 = image_width / 2;
n_px = image_height * image_width;
n_px2 = image_height2 * image_width2;
n_channels = size(mask, 3);

if nargout > 1
    downsample_map = reshape(1:n_px2, image_height2, image_width2);
    downsample_map = repelem(downsample_map, 2, 2);
    for i = 2:n_channels
        downsample_map(mask(:, :, i)) = downsample_map(mask(:, :, i)) + (i - 1) * n_px2;
    end
end

upsample_indices = reshape(1:n_px, image_height, image_width);
n_px_in_bayer_tile = 4;
upsample_map = zeros(image_height2, image_width2, n_px_in_bayer_tile);
for i = 1:n_px_in_bayer_tile
    if i < 3
        k = i;
    else
        k = i - 1;
    end
    upsample_indices_k = upsample_indices(mask(:, :, k));
    if i == 2 || i == 3
        % Assuming the Bayer pattern will never have two Green pixels in
        % the same column
        upsample_indices_k_filter = repmat(...
            [true(image_height2, 1); false(image_height2, 1)], image_width2, 1 ...
            );
    end
    if i == 2
        upsample_indices_k = upsample_indices_k(upsample_indices_k_filter);
    elseif i == 3
        upsample_indices_k = upsample_indices_k(~upsample_indices_k_filter);
    end
    upsample_map(:, :, i) = reshape(upsample_indices_k, image_height2, image_width2);
end

if nargout > 3
    [upsample_map_centersX, upsample_map_centersY] = meshgrid(...
        1:2:image_width, 1:2:image_height...
    );
    upsample_map_centers = cat(3, upsample_map_centersX, upsample_map_centersY);
    upsample_map_centers = upsample_map_centers + 0.5;
end

I_channels = zeros(image_height2, image_width2, n_channels);
I_channels(:, :, 1) = reshape(...
    I_raw(upsample_map(:, :, 1)), image_height2, image_width2...
    );
I_channels(:, :, 2) = reshape(...
    0.5 * (I_raw(upsample_map(:, :, 2)) + I_raw(upsample_map(:, :, 3))),...
    image_height2, image_width2...
    );
I_channels(:, :, 3) = reshape(...
    I_raw(upsample_map(:, :, 4)), image_height2, image_width2...
    );

end