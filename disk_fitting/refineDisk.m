function [...
    lambda, theta, t, k, lightness...
] = refineDisk(...
    I, channel_mask, lambda0, theta0, t0, k0, lightness0, r_max, varargin...
)
% REFINEDISK  Improve an ellipse fit to an image blob
%
% ## Syntax
% [...
%    lambda, theta, t, k, lightness...
% ] = refineDisk(...
%    I, channel_mask, lambda0, theta0, t0, k0, lightness0, r_max [, verbose]...
% )
%
% ## Description
% [...
%    lambda, theta, t, k, lightness...
% ] = refineDisk(...
%    I, channel_mask, lambda0, theta0, t0, k0, lightness0, r_max [, verbose]...
% )
%   Returns one to all parameters of a refined parametric ellipse
%
% ## Input Arguments
%
% I -- Image
%   A 2D array representing either a RAW image, or an image from a
%   monochromatic sensor (including a non-mosaicked image for a narrow
%   wavelength band).
%
% channel_mask -- Image colour channel indices
%   For a monochromatic image, `channel_mask` should be `true(size(I))`.
%   For a RAW image, it should be a 3D array of logical values, such as
%   produced by 'bayerMask()'.
%
% lambda0 -- Initial guess for ellipse dimensions
%   A two-element vector containing the major and minor semi-axis lengths
%   of the ellipse, respectively.
%
% theta0 -- Initial guess for ellipse orientation
%   The angle in radians from the positive x-axis to the major axis of the
%   ellipse.
%
% t0 -- Initial guess for ellipse centre
%   A two-element vector containing the coordinates of the centre of the
%   ellipse.
%
% k0 -- Initial guess for ellipse edge width
%   The linear lightness transition between the ellipse and its
%   surroundings extends 1/k units both inside and outside of the geometric
%   boundary of the ellipse.
%
% lightness0 -- Initial guess for ellipse values
%   A c x 2 matrix where each row contains the lightness within, and
%   outside the ellipse, respectively. The different rows correspond to the
%   colour channels in `channel_mask`.
%
% r_max -- Bounding radius
%   An upper bound on the major semi-axis length of the ellipse. If `r_max`
%   is infinite, the entire image `I` will be the domain used to evaluate
%   the fitting error of the ellipse.
%
% verbose -- Debugging and visualization controls
%   If `verbose` is `true`, graphical and console output will be generated
%   for debugging purposes.
%
% ## Output Arguments
%
% lambda -- Refined ellipse dimensions
%   A two-element vector containing the major and minor semi-axis lengths
%   of the ellipse, respectively.
%
% theta -- Refined ellipse orientation
%   The angle in radians from the positive x-axis to the major axis of the
%   ellipse.
%
% t -- Refined ellipse centre
%   A two-element vector containing the coordinates of the centre of the
%   ellipse.
%
% k -- Refined ellipse edge width
%   The linear lightness transition between the ellipse and its
%   surroundings extends 1/k units both inside and outside of the geometric
%   boundary of the ellipse.
%
% lightness -- Refined ellipse values
%   A c x 2 matrix where each row contains the lightness within, and
%   outside the ellipse, respectively. The different rows correspond to the
%   colour channels in `channel_mask`.
%
% ## Algorithm
%
% The refined ellipse is obtained using the Levenberg-Marquardt algorithm
% to minimize the sum of squared differences between the parametric model
% of the ellipse and the image values. A separate set of interior and
% exterior lightnesses is used in the parametric model of the ellipse for
% each colour channel.
%
% ## References
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.
%
% See also findAndFitDisks, bayerMask, ellipseModel, plotEllipse

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 24, 2018

nargoutchk(1, 5);
narginchk(8, 9);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

n_channels = size(channel_mask, 3);

% Initial Guess
p0 = zeros(1, 8);
if verbose
    disp('Initial guess:')
    lambda0 %#ok<NOPRT>
    theta0 %#ok<NOPRT>
    t0 %#ok<NOPRT>
    k0 %#ok<NOPRT>
    lightness0 %#ok<NOPRT>
    disp('Perturbation (the actual initial guess for the optimizer):')
    disp(p0)
end

% Extract the region from which to evaluate the fitting cost function
% Use a circle of radius `r_max`
px_values = cell(n_channels, 1);
px_coords = cell(n_channels, 1);
image_height = size(I, 1);
image_width = size(I, 2);
n_px = image_height * image_width;
if isfinite(r_max)
    bounding_box_x = floor(t0(1) - r_max):ceil(t0(1) + r_max);
    bounding_box_x = bounding_box_x(bounding_box_x >= 1 & bounding_box_x <= image_width);
    bounding_box_y = floor(t0(2) - r_max):ceil(t0(2) + r_max);
    bounding_box_y = bounding_box_y(bounding_box_y >= 1 & bounding_box_y <= image_height);
    [bounding_box_x_grid, bounding_box_y_grid] = meshgrid(bounding_box_x,bounding_box_y);
    bounding_box = [bounding_box_x_grid(:), bounding_box_y_grid(:), ones(numel(bounding_box_x_grid), 1)];
    [ ~, boundaryFun0 ] = ellipseModel(...
        lambda0, theta0, t0, k0, lightness0(1:2), true...
    );
    bounding_flags = boundaryFun0(bounding_box);
    bounding_coords = bounding_box(bounding_flags < 1, :);
    bounding_ind = sub2ind(...
        [image_height, image_width],...
        bounding_coords(:, 2),...
        bounding_coords(:, 1)...
    );
    for i = 1:n_channels
        bounding_ind_i = bounding_ind + (i - 1) * n_px;
        filter_i = channel_mask(bounding_ind_i);
        px_values{i} = I(bounding_ind_i(filter_i));
        px_coords{i} = bounding_coords(filter_i, :);
    end
else
    [bounding_coords_y, bounding_coords_x] = ind2sub([image_height, image_width], 1:n_px);
    bounding_coords = [bounding_coords_x, bounding_coords_y, ones(n_px, 1)];
    for i = 1:n_channels
        filter_i = channel_mask(:, :, i);
        px_values{i} = I(filter_i);
        px_coords{i} = bounding_coords(filter_i(:), :);
    end
end

n_points_per_channel = zeros(n_channels, 1);
for i = 1:n_channels
    n_points_per_channel(i) = size(px_values{i}, 1);
end
n_points = sum(n_points_per_channel);

    function parseSolution( p )
        lambda = lambda0 + p(1:2);
        theta = theta0 + p(3);
        t = t0 + p(4);
        k = k0 + p(5);
        lightness = lightness0 + reshape(p(6:end), n_channels, []);
    end

    function [ err ] = costFunction( p )
        parseSolution(p);
        err = zeros(n_points, 1);
        ind = 1;
        for j = 1:n_channels
            valueFun_j = ellipseModel(...
                lambda, theta, t, k, lightness(j, :), true...
            );
            ellipse_values = valueFun_j(px_coords{j});
            err(ind:(ind + n_points_per_channel(j))) = px_values{j} - ellipse_values;
            ind = ind + n_points_per_channel(j);
        end
    end

% Invoke optimizer
options = optimoptions('lsqnonlin',...
    'FunValCheck', 'on', 'Algorithm', 'levenberg-marquardt'...
);
if verbose
    options.Display = 'iter-detailed';
else
    options.Display = 'none';
end
p = lsqnonlin(@costFunction, p0, [], [], options);
parseSolution(p);

if any(lambda > r_max)
    error('Ellipse exceeded bounding radius.')
end

if verbose
    disp('Final guess from optimizer:')
    disp(p)
    lambda %#ok<NOPRT>
    theta %#ok<NOPRT>
    t %#ok<NOPRT>
    k %#ok<NOPRT>
    lightness %#ok<NOPRT>
    
    fg = figure;
    imshow(I);
    [~, ~, ellipse_to_world0] = ellipseModel(...
        lambda0, theta0, t0, k0, lightness0(1:2), true...
    );
    plotEllipse(ellipse_to_world0, fg);
    title('Initial ellipse')
    
    fg = figure;
    imshow(I);
    [~, ~, ellipse_to_world] = ellipseModel(...
        lambda, theta, t, k, lightness(1:2), true...
    );
    plotEllipse(ellipse_to_world, fg);
    title('Final ellipse')
end

end

