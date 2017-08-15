function [ div_map ] = divergenceMap( I_1, I_2, f_div, std_weights, p )
% DIVERGENCEMAP  Generalized window-based comparison of images
%
% ## Syntax
% div_map = divergenceMap( I_1, I_2, f_div, std_weights, p )
% div_map = divergenceMap( I_1, I_2, f_div, [w, h], 0 )
%
% ## Description
% div_map = divergenceMap( I_1, I_2, f_div, std_weights, p )
%   Returns a difference map between the two images, quantified using the
%   given divergence operator. A Gaussian window is constructed around each
%   pixel.
% div_map = divergenceMap( I_1, I_2, f_div, [w, h], 0 )
%
% ## Input Arguments
%
% I_1 -- First image
%   An image_height x image_width array, representing a greyscale image.
%
% I_2 -- Second image
%   An image_height x image_width array, representing a greyscale image.
%
% f_div -- Divergence operator
%   A function handle to a function with a signature of the form
%   `div = f_div( I_1, I_2, weights )`. The arguments of `f_div` are as
%   follows:
%   - I_1: First image window. An h x w array, representing a greyscale
%     image window.
%   - I_2: Second image window. An h x w array, representing a greyscale
%     image window.
%   - weights: Window weights. An h x w array of weights to apply to the
%     positions within the window. `weights` allows for using "soft"
%     windows.
%
%   `h` and `w` will be odd integers. `div` is a divergence measurement.
%
% std_weights -- Gaussian window standard deviation
%   For a window defined by a Gaussian weighting function, the standard
%   deviation of the Gaussian with respect to the distance from the window
%   centre.
%
% p -- Gaussian window cutoff value
%   For a window defined by a Gaussian weighting function, the window will
%   be of width (and height) such that a vertical or horizontal cross
%   section of the window, through the centre of the window, will contain
%   `1 - 2 * p` of a Gaussian distribution with standard deviation
%   `std_weights`.
%
%   If `p` is zero, a box distribution is used instead of a Gaussian
%   distribution for defining the window.
%
% [w, h] -- Window width and height
%   The width (columns) and height (rows) of a window defined by a box
%   distribution. In other words, when `p` is zero, the soft Gaussian
%   window is replaced by a conventional window with equal weights for all
%   neighbouring pixel locations.
%
%   `w` and `h` must be odd integers.
%
% ## Output Arguments
%
% div_map -- Divergence map
%   An image_height x image_width array, where `div_map(i,j)` is the value
%   of the divergence operator `f_div` evaluated on windows surrounding the
%   pixels at indices `(i,j)` in the input images. The windows are either
%   Gaussian weightings of neighbouring pixels, or simple rectangular
%   windows, as defined by the `std_weights`, `p`, and `[w, h]` input
%   arguments.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2017

    function [ weights, w, h ] = gaussianKernel(s_distance, p)
        % Returns a filtering kernel, where filter weights are obtained
        % from a Gaussian function of distance from the center of the
        % kernel.
        
        % Define the size of the filter so that it encompasses `1 - 2 * p`
        % of the normal distribution
        radial_limit = norminv([p, 1 - p], 0, s_distance);
        radial_limit = ceil(radial_limit(2));
        w = radial_limit * 2 + 1;
        h = w;
        weights = zeros(h, w);
        center = [radial_limit, radial_limit] + 1;
        for i_inner = 1:h
            for j_inner = 1:w
                dx = j_inner - center(2);
                dy = i_inner - center(1);
                r = sqrt(dx ^ 2 + dy ^ 2);
                weights(i_inner, j_inner) = normpdf(r, 0, s_distance);
            end
        end
        % Normalize weights
        weights = weights ./ sum(sum(weights));
    end

% Create the window
if p == 0
    w = std_weights(1);
    h = std_weights(2);
    weights = ones(h, w);
else
    if length(std_weights) > 1
        error('`std_weights` must be a scalar if `p > 0`')
    end
    [ weights, w, h ] = gaussianKernel(std_weights, p);
end

% Fill in the divergence map
div_map = zeros(size(I_1));
end_i = size(I_1, 2) - w;
end_j = size(I_1, 1) - h;
for i = 1:end_i
    iw = i + w - 1;
    for j = 1:end_j
        jh = j + h - 1;
        div_map(j, i) = f_div( I_1(j:jh, i:iw), I_2(j:jh, i:iw), weights );
    end
end
end

