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
%   corresponding to `I_raw`. For example, 'gbrg'.
%
%   This argument has the same form as the `sensorAlignment` input argument
%   of `demosaic()`.
%
% ## Output Arguments
%
% I_channels -- Full-colour downsampled image
%   An image_height / 2 x image_width / 2 x 3 array, containing a
%   downsampled version of `I_raw`, separated into colour channels. Each
%   pixel in `I_channels` corresponds to a Bayer tile in `I_raw`.
%   `I_channels(i,j,1)` is the Red value from the (i,j)-th Bayer tile, and
%   `I_channels(i,j,3)` is the Blue value from the (i,j)-th Bayer tile.
%   `I_channels(i,j,2)` is the average of the two Green values from the
%   (i,j)-th Bayer tile.
%
% downsample_map -- Downsampling pixel coordinate conversion
%   An image_height x image_width array. `downsample_map(i,j)` is the
%   linear index of a pixel in `I_channels` corresponding to the pixel
%   `I_raw(i,j)`. For example, if `I_raw(i,j)` represents a Red pixel, then
%   this pixel maps to `I_channels(downsample_map(i,j))`, where
%   `downsample_map(i,j)` will be the index of a value in `I_channels(:, :, 1)`.
%
% upsample_map -- Upsampling pixel coordinate conversion
%   An image_height x image_width x 4 array. `upsample_map(i,j, :)`
%   contains the linear indices of the pixels in `I_raw` corresponding to
%   the pixel `I_channels(i,j)`. Each pixel in `I_channels` maps to a Bayer
%   tile (4 pixels) in `I_raw`.
%
% upsample_map_centers -- Coordinate frame conversion
%   An image_height x image_width array x 2 array, where
%   `upsample_map_centers(i, j, :)` is the `(x, y)` position in `I_raw`
%   corresponding to the centre of the pixel `I_channels(i,j,:)`.
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
downsample_map = reshape(1:n_px2, image_height2, image_width2);
downsample_map = repelem(downsample_map, 2, 2);
for i = 2:n_channels
    downsample_map(mask(:, :, i)) = downsample_map(mask(:, :, i)) + (i - 1) * n_px2;
end

upsample_map = reshape(1:n_px, image_height, image_width);
upsample_map = cat(3,...
    upsample_map(1:2:end, 1:2:end),...
    upsample_map(2:2:end, 1:2:end),...
    upsample_map(2:2:end, 1:2:end),...
    upsample_map(2:2:end, 2:2:end)...
);

[upsample_map_centersY, upsample_map_centersX] = ind2sub(...
    [image_height, image_width], upsample_map(:, :, 1)...
    );
upsample_map_centers = cat(3, upsample_map_centersX, upsample_map_centersY);
upsample_map_centers = upsample_map_centers + 0.5;

I_channels = zeros(image_height, image_width, n_channels);
for i = 1:n_channels
    if i ~= 2
        I_channels(:, :, i) = I_raw(mask(:, :, i));
    else
        % Green channel contains 2/4 pixels
        I_raw_green = I_raw .* mask(:, :, i) / 2;
        for j = 1:size(upsample_map, 3)
            I_channels(:, :, i) = I_channels(:, :, i) +...
                reshape(I_raw_green(upsample_map(:, :, i)), image_height2, image_width2);
        end
    end
end