function [centers, I_segmented, I_soft] = segmentColors(I, K, varargin)
% SEGMENTCOLORS  Soft segmentation of image chromaticity values
%
% ## Syntax
% centers = segmentColors(I, K [, verbose])
% [centers, I_segmented] = segmentColors(I, K [, verbose])
% [centers, I_segmented, I_soft] = segmentColors(I, K [, verbose])
% [____] = segmentColors(I, 'classes' [, verbose])
%
% ## Description
% centers = segmentColors(I, K [, verbose])
%   Returns the centroids of K chromaticity clusters in the image.
%
% [centers, I_segmented] = segmentColors(I, K [, verbose])
%   Additionally returns a hard segmentation of the image into K
%   chromaticity clusters.
%
% [centers, I_segmented, I_soft] = segmentColors(I, K [, verbose])
%   Additionally returns a membership probabilities for K chromaticity
%   clusters.
%
% [____] = segmentColors(I, 'classes' [, verbose])
%   Treats the input image as containing class indices instead of colours,
%   allowing for direct control over clustering (but only for hard
%   clustering). An error will be thrown if `I` has more than one channel.
%
% ## Input Arguments
%
% I -- Input image
%   A greyscale or sRGB image.
%
% K -- Number of clusters
%   The desired number of clusters into which to divide the chromaticity
%   values of the image.
%
% verbose -- Verbosity flag
%   If `true`, graphical output will be produced for debugging and
%   visualization. Defaults to `false` if not passed.
%
% ## Output Arguments
%
% centers -- Cluster centroids
%   The centroids of each of the clusters of pixels in the segmentation of
%   the image. While the pixels are clustered based on their chromaticity
%   values, the centroids are computed as the mean colours in the original
%   colour space of the image. The i-th row of `centers` is the centroid of
%   the i-th cluster.
%
% I_segmented -- Cluster labels
%   A greyscale image, where pixels store the indices of their
%   corresponding chromaticity clusters.
%
% I_soft -- Cluster probabilities
%   An image_height x image_width x K array, where `I_soft(:, :, i)` stores
%   the probabilities that the pixels in the input image are members of the
%   i-th chromaticity cluster. For the syntax where `'classes'` is passed
%   instead of `K`, `I_soft` will be a binary image, as the input image `I`
%   represents a hard clustering.
%
% ## Notes
% - If `I` is a greyscale image, its values are used directly for
%   clustering. If `I` is a colour image, the a* and b* channels of its
%   representation in the CIE 1976 L*a*b* colour space are used for
%   clustering. The clustering is performed in a one-dimensional space
%   corresponding to the first principal component of the a* and b*
%   channels.
% - Clustering is performed by multilevel thresholding using Otsu’s method,
%   as implemented in MATLAB's 'multithresh()' function.
% - Cluster membership probabilities are computed from a Gaussian
%   Naive Bayes model fitted to the cluster labels of the pixels.
%
% See also multithresh, rgb2lab, fitcnb

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 13, 2018

nargoutchk(1, 3);
narginchk(2, 3);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

is_classes = false;
if isStringScalar(K) || ischar(K)
    if strcmp(K, 'classes')
        is_classes = true;
    else
        error('The second argument must be either an integer, or `''classes''`.');
    end
end

I = im2double(I);

% Convert to a single-channel image
image_height = size(I, 1);
image_width = size(I, 2);
n_channels = size(I, 3);
if n_channels == 1
    I_mono = I;
    I_mono_linear = I_mono(:);
elseif n_channels == 3
    if is_classes
        error('The second argument cannot be `''classes''` when `I` has three channels.');
    end
    
    I_lab = rgb2lab(I);
    I_a = I_lab(:, :, 2);
    I_b = I_lab(:, :, 3);
    [coeff, scores, latent] = pca([I_a(:), I_b(:)]);
    
    if verbose
        % Practical ranges for the CIE L*a*b* coordinates:
        % https://stackoverflow.com/questions/19099063/what-are-the-ranges-of-coordinates-in-the-cielab-color-space
        a = [-86.185 98.254];
        b = [-107.863 94.482];
        
        [A, B] = meshgrid(a(1):a(2), b(1):b(2));
        lab = cat(3, 100 * ones(size(A)), A, B);
        rgb = lab2rgb(lab);
        figure;
        image('XData', a, 'YData', b, 'CData', rgb)
        coeff_scaled = coeff .* repmat(latent.', size(coeff, 1), 1);
        hold on
        quiver(...
            zeros(size(coeff, 2), 1), zeros(size(coeff, 2), 1),...
            coeff_scaled(1, :).', coeff_scaled(2, :).' ...
            );
        hold off
        title('Principal components of the image a* and b* channels')
        xlabel('a*')
        ylabel('b*')
    end
    
    I_mono = reshape(scores(:, 1), image_height, image_width);
    % Rescale to the range [0, 1]
    I_mono_linear = I_mono(:);
    I_mono_min = min(I_mono_linear);
    I_mono = (I_mono - I_mono_min) ./ (max(I_mono_linear) - I_mono_min);
    
    if verbose
        figure;
        imshow(I_mono);
        title('Scores along the first principal component of the image a* and b* channels');
    end
else
    error('The input image `I` must have one or three colour channels.')
end

if is_classes
    % Remap class indices to consecutive integers
    [ids_old, ~, map] = unique(I_mono_linear);
    K = length(ids_old);
    ids_new = (1:K).';
    I_segmented_linear = ids_new(map);
    I_segmented = reshape(I_segmented_linear, image_height, image_width);
else
    % Threshold into classes
    t = multithresh(I, K - 1);
    I_segmented = imquantize(I_mono, t);
    I_segmented_linear = I_segmented(:);
end

class_instances_color = cell(K, 1);
centers = zeros(K, n_channels);
for k = 1:K
    filter = (I_segmented == k);
    class_instances_color_k = zeros(sum(sum(filter)), n_channels);
    for c = 1:n_channels
        I_c = I(:, :, c);
        class_instances_color_k(:, c) = I_c(filter);
    end
    class_instances_color{k} = class_instances_color_k;
    centers(k, :) = mean(class_instances_color_k, 1);
end

if verbose
    rgb = zeros(image_height * image_width, n_channels);
    for c = 1:n_channels
        centers_c = centers(:, c);
        rgb(:, c) = centers_c(I_segmented_linear);
    end
    rgb = reshape(rgb, size(I));
    figure;
    imshow(rgb)
    title('Hard-thresholded image with classes labelled with mean colours')
end

% Soft threshold by fitting a multinomial regression model
% B = mnrfit(I_mono_linear, I_segmented_linear);
% I_soft = mnrval(B, I_mono_linear);

if is_classes
    I_soft = zeros(image_height, image_width, K);
    n_px = image_height * image_width;
    I_soft((I_segmented_linear - 1) * n_px + (1:n_px).') = 1;
else
    % Soft threshold by fitting a Naive Bayes model, replacing class labels by
    % class probabilities
    B = fitcnb(...
        I_mono_linear, I_segmented_linear,...
        'DistributionNames', 'normal', 'CrossVal', 'off');
    [ ~, I_soft ] = predict(B, I_mono_linear);
    I_soft = reshape(I_soft, image_height, image_width, K);
end

if verbose
    for k = 1:K
        figure;
        imshow(I_soft(:, :, k))
        title(sprintf('Probability of the colour class %d', k));
    end
end

end

