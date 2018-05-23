function [ I_raw ] = mosaic(I_color, align)
% MOSAIC  Mosaic colour images to produce colour-filter array images
%
% ## Syntax
% I_raw = mosaic(I_color, align)
%
% ## Description
% I_raw = mosaic(I_color, align)
%   Converts colour images to colour-filter array images
%
% ## Input Arguments
%
% I_color -- Full-colour images
%   Either an image_height x image_width x 3 x n array, where the last
%   dimension indexes colour images, or a cell vector of length 'n', where
%   each cell contains an image_height x image_width x 3 colour image.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'. `align` has the same form as the `sensorAlignment`
%   input argument of `demosaic()`.
%
% ## Output Arguments
%
% I_raw -- RAW images
%   Either an image_height x image_width x n array where the last dimension
%   indexes the raw colour-filter array versions of the input images, or a
%   cell vector of length 'n', where each cell contains an image_height x
%   image_width colour-filter array image. `I_raw` takes the same form
%   (array or cell vector) as `I_color`.
%
% See also bilinearDemosaic, demosaic, bayerMask

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 23, 2018

nargoutchk(1, 1);
narginchk(2, 2);

cell_input = iscell(I_color);
if cell_input
    I_color_i = I_color{1};
else
    I_color_i = I_color(:, :, :, 1);
end

sz = size(I_color_i);
mask = bayerMask(sz(1), sz(2), align);

if cell_input
    n_images = length(I_color);
    I_raw = cell(n_images, 1);
else
    n_images = size(I_color, 4);
    I_raw = zeros(sz(1), sz(2), n_images);
end

for i = 1:n_images
    if cell_input
        I_color_i = I_color{i};
    else
        I_color_i = I_color(:, :, :, i);
    end
    if any(sz ~= size(I_color_i))
        error('Not all images have the same dimensions.');
    end
    I_raw_i = zeros(sz(1), sz(2));
    for c = 1:sz(3)
        mask_c = mask(:, :, c);
        I_color_i_c = I_color_i(:, :, c);
        I_raw_i(mask_c) = I_color_i_c(mask_c);
    end
    if cell_input
        I_raw{i} = I_raw_i;
    else
        I_raw(:, :, i) = I_raw_i;
    end
end

end

