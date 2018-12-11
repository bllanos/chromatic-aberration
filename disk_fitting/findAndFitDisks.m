function [ centers ] = findAndFitDisks(...
    I, mask, align, image_bounds, radius, k0, options, varargin...
)
% FINDANDFITDISKS  Fit ellipses to blobs in an image
%
% ## Syntax
% centers = findAndFitDisks(...
%   I, mask, align, image_bounds, radius, k0, options [, verbose]...
% )
%
% ## Description
% centers = findAndFitDisks(...
%   I, mask, align, image_bounds, radius, k0, options [, verbose]...
% )
%   Returns the centres of ellipses fit to blobs in the image
%
% ## Input Arguments
%
% I -- Image
%   A 2D array representing either a RAW image, or an image from a
%   monochromatic sensor (including a non-mosaicked image for a narrow
%   wavelength band).
%
% mask -- Region of interest
%   A 2D logical array the same size as `I`, indicating the region in which
%   the function will look for blobs.
%
%   If `mask` is empty, the entire image will be processed.
%
% align -- Bayer pattern format
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'. If the image is not mosaicked, an empty array
%   (`[]`) should be passed.
%
%   This argument has the same form as the `sensorAlignment` input argument
%   of `demosaic()`.
%
% image_bounds -- Image domain
%   The rectangular domain of the image in world coordinates.
%   `image_bounds` is a vector containing the following elements:
%   1 - The x-coordinate of the bottom left corner of the image
%   2 - The y-coordinate of the bottom left corner of the image
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
%   `image_bounds` is used to convert the coordinates in `centres` from
%   image space to world space. `image_bounds` is useful if this function
%   is being called on different sub-images, as the results of the
%   different calls will be in a common frame of reference if
%   `image_bounds` is updated appropriately between calls. If
%   `image_bounds` is empty, no coordinate conversion will be performed.
%
% radius -- Binary image cleanup radius
%   The radius of the disk structuring element that will be used during
%   morphological operations to clean up an intermediate binary image prior
%   to identifying blobs.
%
% k0 -- Initial guess for disk boundary width
%   The approximate width in pixels of the transition region between the
%   intensity of blobs in the image with the intensity of their
%   surroundings. A higher value should be passed if the sensor resolution
%   is higher, or if the blobs are more defocused.
%
% options -- Data processing options
%   A structure with the following fields:
%   - mask_as_threshold: If `true`, `mask` will be used as the initial
%     binary image containing the blobs themselves, rather than indicating
%     the region in which to search for blobs.
%   - bright_disks: If `true`, the function will look for bright blobs on a
%     dark background. Otherwise, the function will look for dark blobs on
%     a bright background.
%   - area_outlier_threshold: The number of sample standard deviations away
%     from the mean blob area for blobs to be rejected as outliers. If
%     zero, no filtering based on area will be performed.
%   - group_channels: If `true`, all colour channels in a mosaicked image
%     will be processed together. This situation corresponds to a colour
%     camera taking a picture through a narrowband colour filter, for
%     example. If `false`, separate ellipses will be fit for each colour
%     channel, to allow for studying displacements between colour channels.
%     This field is optional if `align` is empty.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% centers -- Ellipse centres
%   A structure array, where each element has a field, 'center', storing a
%   two-element vector of the x and y-coordinates of an ellipse. The
%   ellipses have been fit to blobs in the image.
%
%   The coordinates in `centers` are pixel coordinates, if `image_bounds`
%   is empty, or world coordinates, if `image_bounds` is not empty.
%
%   `centers` has dimensions n x m, where 'n' is the number of blobs. 'm'
%   is `1`, if `align` is empty, or if `align` is not empty, but
%   `options.group_channels` is `true`. Otherwise, 'm' is three (for the
%   three colour channels of an RGB image). Note that the elements in a
%   given row of `centers` correspond, by design, to the same blob,
%   because, prior to the refinement stage, the initial blobs are detected
%   in binary images merged across colour channels.
%
% ## Algorithm
%
% - Blobs are detected by analyzing the binary image produced by Otsu
%   thresholding of the region of interest (`mask`) in the image. A
%   different threshold is calculated for each colour channel, in the case
%   of a mosaicked image, but the binary images from all colour channels
%   are merged prior to blob detection.
% - Initial ellipses are fit using the MATLAB 'regionprops()' function.
% - Ellipses are refined using 'refineDisk()'
%
% See also refineDisk, ellipseModel, plotEllipse, regionprops, otsuthresh, imopen, imclose

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 20, 2018

nargoutchk(1,1);
narginchk(7,8);

% Parse input arguments
if ~isempty(varargin)
    verbose = varargin{1};
    verbose_disk_search = verbose.verbose_disk_search;
    verbose_disk_refinement = verbose.verbose_disk_refinement;
    display_final_centers = verbose.display_final_centers;
else
    verbose_disk_search = false;
    verbose_disk_refinement = false;
    display_final_centers = false;
end

I = im2double(I);

image_height = size(I, 1);
image_width = size(I, 2);
if size(I, 3) ~= 1
    error('This function processes RAW RGB images or greyscale images only, not demosaicked images.');
end

% Binarize the image
single_channel = isempty(align);
if single_channel
    n_channels = 1;
    channel_mask = true(image_height, image_width);
else
    n_channels = 3;
    channel_mask = bayerMask(image_height, image_width, align);
end

bw = false(image_height, image_width, n_channels);
for c = 1:n_channels
    if isempty(mask)
        mask_c = channel_mask(:, :, c);
    else
        mask_c = mask & channel_mask(:, :, c);
    end
    if ~isempty(mask) && options.mask_as_threshold
        bw(:, :, c) = mask_c;
    elseif isempty(mask) && options.mask_as_threshold
        error('`mask` is empty, but `options.mask_as_threshold` is true.');
    else
        [counts_c, edges_c] = histcounts(I(mask_c));
        threshold_c = otsuthresh(counts_c);
        threshold_c = threshold_c * (edges_c(end) - edges_c(1)) + edges_c(1);
        if options.bright_disks
            bw(:, :, c) = imbinarize(I,threshold_c) & mask_c;
        else
            bw(:, :, c) = ~imbinarize(I,threshold_c) & mask_c;
        end
    end
    if verbose_disk_search
        figure;
        imshow(bw(:, :, c));
        title(sprintf('Binary image for colour channel %d', c));
    end
end

% Extract binary regions across all colour channels
bw_fused = any(bw, 3);
% Morphologial operations to remove small regions and fill holes
if radius > 0
    se = strel('disk',radius);
    bw_fused_cleaned = imopen(bw_fused, se);
    bw_fused_cleaned = imclose(bw_fused_cleaned, se);
    if verbose_disk_search
        figure;
        imshowpair(bw_fused, bw_fused_cleaned, 'montage');
        title('Binary image for all channels, before (left) and after (right) morphological cleanup');
    end
else
    bw_fused_cleaned = bw_fused;
end

ellipse_stats = regionprops(...
    bw_fused_cleaned,...
    'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Area'...
    );
n_ellipses = length(ellipse_stats);
if verbose_disk_search
    fg = figure;
    imshow(I);
    for i = 1:n_ellipses
        [~, ~, ellipse_to_world_i] = ellipseModel(...
            [ellipse_stats(i).MajorAxisLength ellipse_stats(i).MinorAxisLength] / 2,...
            deg2rad(ellipse_stats(i).Orientation),...
            ellipse_stats(i).Centroid,...
            0,...
            [0 1],...
            true...
        );
        plotEllipse(ellipse_to_world_i, fg);
    end
    title('Initial detected ellipses')
end
centers_matrix_initial = vertcat(ellipse_stats.Centroid);

% Filter out outlier ellipses based on size
% Reference: MATLAB documentation page on "Inconsistent Data"
if n_ellipses > 2 && options.area_outlier_threshold > 0
    areas = vertcat(ellipse_stats.Area);
    sigma_areas = std(areas);
    if sigma_areas > 0
        mu_areas = mean(areas);
        mu_areas_rep = repmat(mu_areas, n_ellipses, 1);
        sigma_areas_rep = repmat(sigma_areas, n_ellipses, 1);
        inliers_filter = (abs(areas - mu_areas_rep) < options.area_outlier_threshold * sigma_areas_rep);
        ellipse_stats = ellipse_stats(inliers_filter);
        centers_matrix_initial = centers_matrix_initial(inliers_filter, :);
        n_ellipses = length(ellipse_stats);
    end
    
    if verbose_disk_search
        figure(fg);
        hold on
        scatter(centers_matrix_initial(:, 1), centers_matrix_initial(:, 2), 'go');
        hold off
        title('Initial detected ellipses (blue), filtered by area (green)')
    end
end

% For later conversion from pixel to world coordinates
convert_coordinates = ~isempty(image_bounds);
if convert_coordinates
    pixel_width = image_bounds(3) / image_width;
    pixel_height = image_bounds(4) / image_height;
end

%% Refine disk fitting
split_channels = (~single_channel && ~options.group_channels);
if split_channels
    n_channels_out = n_channels;
    channel_masks_refinement = cell(n_channels_out, 1);
    for c = 1:n_channels_out
        channel_masks_refinement{c} = channel_mask(:, :, c);
    end
else
    n_channels_out = 1;
    channel_masks_refinement = {channel_mask};
end
centers = struct('center', cell(n_ellipses, n_channels_out));
centers_matrix = zeros(n_ellipses, 2, n_channels_out);

% Initialize blob and non-blob lightnesses to be the modes of the
% histograms of the corresponding binary regions. Perform the calculations
% for each colour channel.
lightness0 = zeros(n_channels, 2);
for c = 1:n_channels
    for b = 1:2
        if b == 1
            channel_values = I(bw(:, :, c));
        else
            if isempty(mask)
                mask_c = channel_mask(:, :, c);
            else
                mask_c = mask & channel_mask(:, :, c);
            end
            channel_values = I(~bw(:, :, c) & mask_c);
        end
        [counts_cb,edges_cb] = histcounts(channel_values);
        [~, ind_cb] = max(counts_cb);
        lightness0(c, b) = mean(edges_cb(ind_cb:(ind_cb + 1)));
    end
end

if split_channels
    lightness0_refinement = cell(n_channels_out, 1);
    for c = 1:n_channels_out
        lightness0_refinement{c} = lightness0(c, :);
    end
else
    lightness0_refinement = {lightness0};
end

% Find the average spacing between disks
if n_ellipses > 1
    sep_mean = 0;
    for i = 1:n_ellipses
        sep = centers_matrix_initial([1:(i-1), (i+1):end], :)...
            - repmat(centers_matrix_initial(i, :), n_ellipses - 1, 1);
        sep = sqrt(dot(sep, sep, 2));
        sep_mean = sep_mean + min(sep);
    end
    r_max = (sep_mean / n_ellipses) / 2;
else
    r_max = Inf;
end

% Refine the disks
refine_filter = true(n_ellipses, 1);
for i = 1:n_ellipses
    % Assume that no colour channel has an entirely zero response in bright
    % areas of other colour channels
    for c = 1:n_channels_out
        [...
            ~, ~, center_ic, ~, ~, result...
        ] = refineDisk(...
            I,...
            channel_masks_refinement{c},...
            [ellipse_stats(i).MajorAxisLength ellipse_stats(i).MinorAxisLength] / 2,...
            deg2rad(ellipse_stats(i).Orientation),...
            ellipse_stats(i).Centroid,...
            k0,...
            lightness0_refinement{c},...
            r_max,...
            verbose_disk_refinement...
        );
        refine_filter(i) = result && refine_filter(i);
        centers_matrix(i, :, c) = center_ic;
        if convert_coordinates
            center_ic = image_bounds(1:2) + [
                pixel_width * center_ic(1),...
                pixel_height * (image_height - center_ic(2))
                ];
        end
        centers(i, c).center = center_ic;
    end
end

if display_final_centers
    for c = 1:n_channels_out
        figure;
        imshow(I);
        hold on
        scatter(centers_matrix(:, 1, c), centers_matrix(:, 2, c), 'bo');
        scatter(...
            centers_matrix(refine_filter, 1, c),...
            centers_matrix(refine_filter, 2, c),...
            'go'...
        );
        hold off
        if split_channels
            title(sprintf('Refined disk centres (blue) and final centres (green) for channel %d', c));
        else
            title('Refined disk centres (blue) and final centres (green) ');
        end
    end
end

centers = centers(refine_filter, :);

end

