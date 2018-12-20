function [ sensitivity ] = relativeSensitivity(...
    var_name, paths, range, align...
)
% RELATIVESENSITIVITY  Fit scaling factors between colour channels
%
% ## Syntax
% sensitivity = relativeSensitivity(var_name, paths, range, align)
%
% ## Description
% sensitivity = relativeSensitivity(var_name, paths, range, align)
%   Load each group of images, and find scaling factors between colour channels
%   in those images.
%
% ## Input Arguments
%
% var_name -- Image variable name
%   A character vector containing the variable name of the variable used to
%   store images in the input '.mat' files. (Images are stored as '.mat' files,
%   not as image files.)
%
% paths -- Image filenames and paths
%   A cell vector of cell vectors of character vectors containing the names and
%   paths of images used to calibrate scaling factors between colour channels.
%   Separate scaling factors are to be found for groups of images corresponding
%   to different top-level cells in `paths`.
%
% range -- Clipping range
%   A two-element vector containing the values at or below which a pixel is
%   considered to be zero, and at or above which a pixel is considered to
%   be saturated, respectively. Only pixels between these two values will
%   be used to calibrate scaling factors between channels.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the images. For example, 'gbrg'.
%
%   `align` has the same form as the `sensorAlignment` input argument of
%   `demosaic()`
%
% ## Output Arguments
%
% sensitivity -- Colour channel scaling factors
%   A 3 x length(paths) matrix, containing scaling factors for the three colour
%   channels (corresponding to the rows) relative to Green (the second row).
%   `sensitivity(2, :)` is a vector of ones, therefore. Each column of
%   `sensitivity` corresponds to a cell of `paths`. In other words, a set of
%   scaling factors is computed for each group of input images.
%
% ## Notes
% - All images are expected to be raw images (colour-filter array images
%   stored as a single channel). The input images should be linearized and dark
%   corrected (dark frame subtracted).
% - The approach for computing relative sensitivities works well only in the
%   event that all colour channels have non-negligible sensitivities under the
%   imaging conditions in each group of input images.
%
% ## References
%
% The approach taken for computing relative spectral sensitivities is inspired
% by the approach taken to merge exposures in:
%
%   Darrodi, M. M., Finlayson, G., Goodman, T., & Mackiewicz, M. (2015).
%   Reference data set for camera spectral sensitivity estimation. Journal
%   of the Optical Society of America A: Optics and Image Science, and
%   Vision, 32(3), 381-391. doi:10.1364/JOSAA.32.000381
%
% See also blendExposures, darkSubtract

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 19, 2018

nargoutchk(1, 1)
narginchk(4, 4)

n_channels = 3; % RGB channels
channel_indices = 1:n_channels;
reference_channel = 2;
reference_filter = (channel_indices == reference_channel);
other_channels = channel_indices(~reference_filter);
n_other_channels = length(other_channels);
n_groups = length(paths);
sensitivity = ones(n_channels, n_groups);

offsets = [
    0, 1;
    -1, 0;
    0, -1;
    1, 0
];
n_offsets = size(offsets, 1);

for g = 1:n_groups
    paths_g = paths{g};
    n_images = length(paths_g);
    pixels = cell(n_images, n_other_channels);
    for i = 1:n_images
        I = loadImage(paths_g{i}, var_name);
        if size(I, 3) ~= 1
            error('The image "%s" is not a RAW image.', paths_g{i});
        end
        image_sampling = size(I);
        mask = reshape(bayerMask(image_sampling(1), image_sampling(2), align), [], n_channels);
        
        % Filter to well-exposed pixels in the reference channel
        mask = mask & repmat(...
            reshape((I > range(1)) & (I < range(2)), [], 1),...
            1, n_channels...
        );
    
        % Filter to pixels with low gradients
        G = bilinearDemosaic(I, align, reference_filter);
        G = imgradient(G);
        level = graythresh(G);
        gradient_mask = ~imbinarize(G, level);
        mask = mask & repmat(reshape(gradient_mask, [], 1), 1, n_channels);
        
        for c_ind = 1:n_other_channels
            c = other_channels(c_ind);
            % Find neighbours in the reference channel
            other_ind = find(mask(:, c));
            n_other = length(other_ind);
            [other_row, other_col] = ind2sub(image_sampling, other_ind);
            pixels_ref = zeros(n_other, n_offsets);
            filter_c = false(n_other, n_offsets);
            for f = 1:n_offsets
                % Select pixels with good neighbours at all offsets
                row = other_row + offsets(f, 1);
                col = other_col + offsets(f, 2);
                filter_c(:, f) = (row > 0) & (row <= image_sampling(1)) &...
                                 (col > 0) & (col <= image_sampling(2));
                row = row(filter_c(:, f));
                col = col(filter_c(:, f));
                ref_ind = sub2ind(image_sampling, row, col);
                ref_mask = mask(ref_ind, reference_channel);
                ref_ind = ref_ind(ref_mask);
                filter_c(filter_c(:, f), f) = ref_mask;
                pixels_ref(filter_c(:, f), f) = I(ref_ind);
            end
            filter_c = all(filter_c, 2);
            pixels_other = I(other_ind(filter_c));
            pixels{i, c_ind} = [
                repmat(pixels_other, n_offsets, 1),...
                reshape(pixels_ref(filter_c, :), [], 1)
            ];
        end
    end
    % Find per-channel conversion factors, from all pairs of good other
    % channel-reference channel pixels
    for c_ind = 1:n_other_channels
        pixels_mat = cell2mat(pixels(:, c_ind));
        component = pca(pixels_mat, 'NumComponents', 1);
        sensitivity(other_channels(c_ind), g) = component(1) ./ component(end);
    end
end

end