function [ stats ] = analyzePSF( varargin )
% DENSIFYRAYS  Find image intensities from discrete samples of ray irradiance
%
% ## Syntax
% empty_stats = analyzePSF( sz )
% stats = analyzePSF( psf_spline, image_position, v_adj [, verbose] )
%
% ## Description
% empty_stats = analyzePSF( sz )
%   Returns a structure array, for preallocation.
% stats = analyzePSF( psf_spline, image_position, v_adj [, verbose] )
%   Returns statistics describing the density function `psf_spline`.
%
% ## Input Arguments
%
% sz -- Structure array dimensions
%   A vector containing the dimensions of the `empty_stats` output array,
%   such that `all(size(empty_stats) == sz)` is true.
%
% psf_spline -- Thin-plate spline model of image intensities
%   A thin-plate smoothing splines modeling image intensity (irradiance) as
%   a function of 2D position on the image plane. Presently, the
%   interpretation of the spline model is unimportant, but this function
%   was written in a context where `psf_spline` represented a point spread
%   function.
%
%   `psf_spline` can be the `image_spline` output argument of
%   'densifyRays()', for example.
%
% image_position -- Image sampling coordinates
%   The 2D coordinates on the image plane corresponding to the raw data
%   used to create `psf_spline`. `image_position` is a two-column array,
%   with the columns containing image x, and y coordinates, respectively.
%
% v_adj -- Sample adjacency lists
%   The indices of the rows, in `image_position`, which are connected to
%   the given row, in `image_position`, in the Delaunay triangulation of
%   `image_position`. `v_adj` is a cell vector of length
%   `size(image_position, 1)`. `v_adj{i}` is a column vector, containing
%   the indices of the image points (rows in `image_position`) adjacent to
%   the `image_position(i, :)`, in the triangulation.
%
%   `v_adj` can be the `v_adj` output argument of 'densifyRays()', for
%   example.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% empty_stats -- Empty statistics structure array
%   A structure array with empty field values, but with the same fields as
%   `stats`. `empty_stats` is useful for preallocating a structure array,
%   to be filled with subsequent calls to 'analyzePSF()'.
%
% stats -- Distribution statistics
%   A structure describing `psf_spline` more concisely, using the following
%   fields:
%   - mean_position: The centroid of `psf_spline`, computed by weighting
%     each position in `image_position` by the evaluation of `psf_spline`
%     at that position.
%   - mean_value: The evaluation of `psf_spline` at `mean_position`; The
%     expected value of `psf_spline` over the domain of samples in
%     `image_position`.
%   - max_position: A two-element row vector containing the x and
%     y-coordinates of the peak value in `psf_spline`.
%
%     If there are multiple local maxima in `psf_spline`, and their average
%     location has a value which is lower than more than one of them, then
%     `max_position` is an NaN 1 x 2 array. `max_position` is also NaN if
%     the location of the peak value is further than `radius` from
%     `mean_position`, as this is assumed to be a location in the
%     extrapolation region of the spline.
%   - max_value: A scalar containing the peak value corresponding to
%     `max_position`. `max_value` is the evaluation of `psf_spline` at
%     `max_position`, and is NaN if `max_position` is NaN.
%   - radius: The weighted mean of the distances of the points in
%     `image_position` from `mean_position`, where the weights are the
%     evaluations of `psf_spline` at the points. `radius` is related to the
%     second moments of `psf_spline`.
%
% See also densifyRays, neighborVertices, tpaps, fmincon

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 28, 2017

    function [f, g] = splineValueAndGradient(x)
        f = -fnval(psf_spline,[x(1);x(2)]);
        if nargout > 1 % gradient required
            g = -[
                fnval(psf_spline_dx,[x(1);x(2)]);
                fnval(psf_spline_dy,[x(1);x(2)])
                ];
        end
    end

nargoutchk(1, 1);

if nargin == 1
    sz = varargin{1};
    stats = struct(...
        'mean_position', cell(sz),...
        'mean_value', cell(sz),...
        'max_position', cell(sz),...
        'max_value', cell(sz),...
        'radius', cell(sz)...
    );
    return
else
    narginchk(3, 4);
    psf_spline =  varargin{1};
    image_position =  varargin{2};
    v_adj =  varargin{3};
end

if nargin > 3
    verbose = varargin{4};
else
    verbose = false;
end

n_points = size(image_position, 1);

% Irradiance-weighted centroid
image_irradiance_spline = fnval(psf_spline,image_position.');
image_irradiance_spline = image_irradiance_spline.';
stats.mean_position = sum(...
    repmat(image_irradiance_spline, 1, 2) .* image_position, 1 ...
) ./ sum(image_irradiance_spline);
stats.mean_value = fnval(psf_spline, stats.mean_position.');

% Irradiance-weighted uncorrected standard deviation of the distance from
% the centroid
mean_position_rep = repmat(stats.mean_position, n_points, 1);
radius = image_position - mean_position_rep;
radius = sqrt(dot(radius, radius, 2));
radius = sum(...
    image_irradiance_spline .* radius ...
) / sum(image_irradiance_spline);
stats.radius = radius;

% Quickly search for local maxima
% `fmincon` would be more accurate, but I think it is overkill
local_maxima_filter = false(n_points, 1);
for i = 1:n_points
    v_adj_i = v_adj{i};
    differences = image_irradiance_spline(i) - image_irradiance_spline(v_adj_i);
    local_maxima_filter(i) = all(differences >= 0);
end
local_maxima_irradiance = image_irradiance_spline(local_maxima_filter);
local_maxima_position = image_position(local_maxima_filter, :);

% Estimate whether there is a true peak irradiance
has_peak = true;
if sum(local_maxima_filter) > 1
    local_maxima_position_mean = sum(...
            repmat(local_maxima_irradiance, 1, 2) .* local_maxima_position, 1 ...
        ) ./ sum(local_maxima_irradiance);
    local_maxima_irradiance_mean = fnval(psf_spline, local_maxima_position_mean.');
    if sum(local_maxima_irradiance >= local_maxima_irradiance_mean) > 1
        % Probably no peak irradiance
        has_peak = false;
    end
end

if has_peak
    % Find the peak irradiance
    psf_spline_dx = fnder(psf_spline,[1, 0]);
    psf_spline_dy = fnder(psf_spline,[0, 1]);
    min_x = min(image_position(:, 1));
    max_x = max(image_position(:, 1));
    min_y = min(image_position(:, 2));
    max_y = max(image_position(:, 2));
    lb = [min_x, min_y];
    ub = [max_x, max_y];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    [~, max_ind] = max(image_irradiance_spline);
    x0 = image_position(max_ind, :);
    nonlcon = [];
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', true);
    [max_position, max_value] = fmincon(...
        @splineValueAndGradient, x0, A, b, Aeq, beq, lb, ub, nonlcon, options...
        );
    max_value = -max_value;
    
    % Validate the solution:
    % If the estimated peak is further from the centroid than the average
    % distance from the centroid, the solution is probably spurious
    r_peak = norm(max_position - stats.mean_position);
    if r_peak > stats.radius
        has_peak = false;
    else
        stats.max_position = max_position;
        stats.max_value = max_value;
    end
end

if ~has_peak
    stats.max_position = nan(1, 2);
    stats.max_value = nan;
end

if verbose
    figure
    pts = fnplt(psf_spline); % I don't like the look of the plot, so I will plot manually below
    surf(pts{1}, pts{2}, pts{3}, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Irradiance');
    colormap summer
    c = colorbar;
    c.Label.String = 'Irradiance';
    
    hold on
    legend_str = {'Thin plate spline', 'Local maxima', 'Mean value', 'Radius'};
    plot3(...
        local_maxima_position(:, 1), local_maxima_position(:, 2),...
        local_maxima_irradiance,...
        'ko','markerfacecolor','c'...
    )
    plot3(...
        stats.mean_position(1), stats.mean_position(2), stats.mean_value,...
        'ko','markerfacecolor','b'...
    )
    radius_points = linspace(0, 2 * pi * (1 - 1 / n_points), n_points);
    radius_points = stats.radius .* [cos(radius_points); sin(radius_points)];
    radius_points = radius_points + mean_position_rep.';
    radius_points(3, :) = fnval(psf_spline, radius_points);
    plot3(...
        radius_points(1, :), radius_points(2, :),...
        radius_points(3, :),...
        '-m'...
    )
    if has_peak
        title('PSF Statistics')
        plot3(...
            stats.max_position(1), stats.max_position(2), stats.max_value,...
            'ko','markerfacecolor','r'...
        )
        legend_str{end + 1} = 'Peak';
    else
        title('PSF Statistics [No peak found]')
    end
    legend(legend_str);
    hold off
end

end

