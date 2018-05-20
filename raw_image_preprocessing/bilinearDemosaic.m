function [ I_rgb ] = bilinearDemosaic(I_raw, align)
% BILINEARDEMOSAIC  Demosaic an image by bilinear interpolation
%
% ## Syntax
% I_rgb = bilinearDemosaic(I_raw, align)
%
% ## Description
% I_rgb = bilinearDemosaic(I_raw, align)
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
% ## Output Arguments
%
% I_rgb -- Full-colour image
%   An image_height x image_width x 3 array, containing the demosaicked
%   image, produced by per-channel bilinear interpolation of `I_raw`.
%
% ## Notes
% - Pixels at the edges of the image will be given values of 'NaN' for
%   colour channels whose pixels do not extend to the corresponding edges.
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
narginchk(2, 2);

image_height = size(I_raw, 1);
image_width = size(I_raw, 2);
mask = bayerMask( image_height, image_width, align );

I_rgb = nan(size(mask));

% Red and Blue
[X, Y] = meshgrid(1:image_width, 1:image_height);
for c = [1 3]
    mask_c = mask(:, :, c);
    x_c = reshape(X(mask_c), image_height / 2, image_width / 2);
    y_c = reshape(Y(mask_c), image_height / 2, image_width / 2);
    I_raw_c = reshape(I_raw(mask_c), image_height / 2, image_width / 2);
    I_rgb(:, :, c) = interp2(x_c, y_c, I_raw_c, X, Y, 'linear');
end

% Green
green_index = 2;
mask_g = mask(:, :, green_index);
I_g = nan(image_height, image_width);
I_g(mask_g) = I_raw(mask_g);
center_mask_g = ~mask_g(2:(end - 1), 2:(end - 1));
missing_ind_g = sub2ind([image_height, image_width], Y(center_mask_g), X(center_mask_g));
I_g(missing_ind_g) = (...
    I_g(missing_ind_g + 1) +...
    I_g(missing_ind_g - 1) +...
    I_g(missing_ind_g + image_height) +...
    I_g(missing_ind_g - image_height)...
    ) / 4; % Central region of image
I_g(1, :) = interp1(X(1, mask_g(1, :)), I_g(1, mask_g(1, :)), X(1, :), 'linear'); % Top
I_g(end, :) = interp1(X(end, mask_g(end, :)), I_g(end, mask_g(end, :)), X(end, :), 'linear'); % Bottom
I_g(:, 1) = interp1(X(mask_g(:, 1), 1), I_g(mask_g(:, 1), 1), X(:, 1), 'linear'); % Left
I_g(:, end) = interp1(X(mask_g(:, end), end), I_g(mask_g(:, end), end), X(:, end), 'linear'); % Right

I_rgb(:, :, green_index) = I_g;

end

