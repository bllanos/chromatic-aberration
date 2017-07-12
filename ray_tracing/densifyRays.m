function [ image_irradiance, v_adj, image_spline, I ] = densifyRays(...
    incident_position_cartesian, r_front, image_position, ray_irradiance,...
    varargin...
)
% DENSIFYRAYS  Model image intensities from discrete samples of ray irradiance
%
% ## Syntax
% image_irradiance = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance [, verbose]...
% )
% [ image_irradiance, v_adj ] = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance [, verbose]...
% )
% [ image_irradiance, v_adj, image_spline ] = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance [, verbose]...
% )
% [ image_irradiance, v_adj, image_spline, I ] = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance, image_bounds, image_sampling [, verbose]...
% )
%
% ## Description
% image_irradiance = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance [, verbose]...
% )
%   Estimate the image intensities at the ray intersections with the image
%   plane.
% [ image_irradiance, v_adj ] = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance [, verbose]...
% )
%   Additionally returns adjacency information about the rays, for further
%   analysis.
% [ image_irradiance, v_adj, image_spline ] = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance [, verbose]...
% )
%   Additionally returns a spline model of the image intensities,
%   approximating
% [ image_irradiance, v_adj, image_spline, I ] = densifyRays(...
%    incident_position_cartesian, r_front, image_position,...
%    ray_irradiance, image_bounds, image_sampling [, verbose]...
% )
%   Additionally returns the estimated image formed by the incident rays.
%
% ## Input Arguments
%
% incident_position_cartesian -- Cartesian sampling positions on the front aperture
%   The positions, expressed in cartesian coordinates (x, y, z), of the
%   incident rays on the front aperture of the lens system. Each row of
%   `incident_position_cartesian` corresponds to a row of `image_position`.
%
%   The front aperture is assumed to be a section of a spherical surface or
%   radius `r_front`, where the centre of the sphere is located at the
%   origin.
%
%   For example, `incident_position_cartesian` is an output argument of
%   'doubleSphericalLens.m'
%
% r_front -- Front lens radius
%   The radius of curvature of the front surface (aperture) of the lens
%   system.
%
% image_position -- Image coordinates of rays
%   The points of intersection of the light paths traced through the lens
%   system with the image plane. `image_position` is a two-column array,
%   with the columns containing image x, and y coordinates, respectively.
%
% ray_irradiance -- Ray irradiance
%   The incident irradiance produced by individual rays, at the points of
%   intersection of the light paths with the image plane.
%
%   `ray_irradiance` is a vector, where the i-th element corresponds to the
%   i-th row of `image_position`.
%
%   Refer to the documentation of 'doubleSphericalLens.m' for more details.
%
% image_bounds -- Image domain
%   The rectangular domain of the image to be produced. `image_bounds` is a
%   vector containing the following elements:
%   1 - The x-coordinate of the bottom left corner of the image
%   2 - The y-coordinate of the bottom left corner of the image
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
% image_sampling -- Image resolution
%   The number of pixels at which to sample in the domain specified by
%   `image_bounds`. `image_sampling` is a two-element vector containing the
%   image height, and width, respectively, measured in pixels.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% image_irradiance -- Discrete image intensities
%   The incident irradiance at the points of intersection of the light
%   paths with the image plane. `image_irradiance` is a transformation of
%   `ray_irradiance` which takes into account changes in the density of
%   rays between `incident_position_cartesian` and `image_position`.
%
%   `image_irradiance` is a vector, where the i-th element corresponds to
%   the i-th row of `image_position`.
%
% v_adj -- Ray sample adjacency lists
%   The indices of the rows, in `image_position`, which are connected to
%   the given row, in `image_position`, in the Delaunay triangulation of
%   `image_position`. `v_adj` is a cell vector of length
%   `size(image_position, 1)`. `v_adj{i}` is a column vector, containing
%   the indices of the image points (rows in `image_position`) adjacent to
%   the `image_position(i, :)`, in the triangulation.
%
%   The Delaunay triangulation of `image_position` is intermediate data
%   used to create `image_spline`, but is useful for other applications,
%   such as in 'analyzePSF()'.
%
% image_spline -- Thin-plate spline model of image intensities
%   A thin-plate smoothing spline modeling image intensity (irradiance) as
%   a function of 2D position on the image plane.
%
%   `irradiance = fnval(image_spline,[x; y])` evaluates the spline model at
%   the given image x, and y positions in the row vectors `x`, and `y`,
%   respectively.
%
% I -- Image
%   The irradiance pattern produced on the image plane by the rays traced
%   through the lens system. The image boundaries, and pixel resolution
%   are defined by the `image_bounds`, and `image_sampling` input
%   arguments, respectively.
%
% ## Notes
% - `image_spline` is expensive to generate for large numbers of rays.
% - `I` is generated from `image_spline`. In the future, it could also be
%   generated (but with more noise) from scattered data interpolation of
%   `image_irradiance` (see also 'scatteredInterpolant'). Note that
%   'analyzePSF()' already uses `scatteredInterpolant()` to interpolate
%   sparse image irradiances. Therefore, if `scatteredInterpolant()` was
%   also used in this function, the interpolant should be passed to
%   'analyzePSF()' to save computation time.
%
% See also doubleSphericalLens, sphereSection, refract, tpaps, delaunayTriangulation, neighborVertices2, analyzePSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 8, 2017

nargoutchk(1, 4);
narginchk(4, 7);

verbose = false;

output_image = false;

if ~isempty(varargin)
    n_varargs = length(varargin);
    if n_varargs == 1
        verbose = varargin{1};
    elseif n_varargs == 2
        image_bounds = varargin{1};
        image_sampling = varargin{2};
    elseif n_varargs == 3
        image_bounds = varargin{1};
        image_sampling = varargin{2};
        verbose = varargin{3};
    end
    
    image_output_requested = (nargout >= 4);
    if image_output_requested && n_varargs < 2
        error('The output image, `I`, cannot be calculated without the input arguments `image_bounds`, and `image_sampling`.');
    elseif ~image_output_requested && n_varargs > 1
        error('More input arguments were passed than are required to produce the first two output arguments.')
    elseif image_output_requested
        output_image = true;
    end    
end

dt_in = delaunayTriangulation(incident_position_cartesian(:, 1:2));
dt_out = delaunayTriangulation(image_position);
n_rays = length(ray_irradiance);

% Plot the input data
if verbose
    figure
    triplot(dt_in);
    hold on
    scatter(...
        incident_position_cartesian(:, 1),...
        incident_position_cartesian(:, 2),...
        [], ray_irradiance, 'filled'...
    )
    xlabel('X');
    ylabel('Y');
    c = colorbar;
    c.Label.String = 'Exit ray irradiance';
    title('Sampling points on the first aperture')
    hold off
    
    figure
    triplot(dt_out);
    hold on
    scatter(...
        image_position(:, 1),...
        image_position(:, 2),...
        [], ray_irradiance, 'filled'...
    )
    xlabel('X');
    ylabel('Y');
    c = colorbar;
    c.Label.String = 'Exit ray irradiance';
    title('Intersection points with the image')
    hold off
end

% Find the neighbouring points of each vertex in each triangulation
[~, v_adj_in] = neighborVertices2(dt_in);
[~, v_adj_out] = neighborVertices2(dt_out);
v_adj = v_adj_out;

% Find the average areas of circles with radii defined between neighbouring
% sample points on the front aperture
mean_areas_in = zeros(n_rays, 1);
r_front_sq = r_front ^ 2;
for i = 1:n_rays
    v_adj_i = v_adj_in{i};
    n_adj = length(v_adj_i);
    v_rep = repmat(incident_position_cartesian(i, :), n_adj, 1);
    cosines = dot(v_rep, incident_position_cartesian(v_adj_i, :), 2);
    areas = 2 * pi * ( r_front_sq - cosines );
    mean_areas_in(i) = mean(areas);
end

% Find the average areas of circles with radii defined between neighbouring
% sample points on the image
mean_areas_out = zeros(n_rays, 1);
for i = 1:n_rays
    v_adj_i = v_adj_out{i};
    n_adj = length(v_adj_i);
    v_rep = repmat(image_position(i, :), n_adj, 1);
    distances_sq = v_rep - image_position(v_adj_i, :);
    distances_sq = dot(distances_sq, distances_sq, 2);
    areas = pi * distances_sq;
    mean_areas_out(i) = mean(areas);
end

% Find change in areas
density_ratios = mean_areas_in ./ mean_areas_out;
if verbose
    figure
    triplot(dt_out);
    hold on
    scatter(...
        image_position(:, 1),...
        image_position(:, 2),...
        [], density_ratios, 'filled'...
    )
    xlabel('X');
    ylabel('Y');
    c = colorbar;
    c.Label.String = 'Ray density change';
    title('Image ray density relative to incident ray density')
    hold off
end

% Adjust ray irradiance values based on density change
image_irradiance = ray_irradiance .* density_ratios;
if verbose
    figure
    triplot(dt_out);
    hold on
    scatter(...
        image_position(:, 1),...
        image_position(:, 2),...
        [], image_irradiance, 'filled'...
    )
    xlabel('X');
    ylabel('Y');
    c = colorbar;
    c.Label.String = 'Image irradiance';
    title('Image irradiance calculated from ray irradiance and density change')
    hold off
end

% Interpolate image irradiance values
if nargout > 2
    image_spline = tpaps(image_position.',image_irradiance.');
    if verbose
        figure
        pts = fnplt(image_spline); % I don't like the look of the plot, so I will plot manually below
        surf(pts{1}, pts{2}, pts{3}, 'EdgeColor', 'none');
        colorbar
        xlabel('X');
        ylabel('Y');
        zlabel('Irradiance');
        colormap summer
        c = colorbar;
        c.Label.String = 'Irradiance';
        title('Interpolation of image irradiance values')
    end

    % Sample on a grid to produce an image
    if output_image
        I = zeros(image_sampling);
        % Avoid extrapolation by only estimating values within the convex hull
        % of the rays on the image plane
        convex_hull_indices = convhull(image_position);
        image_position_convhull = image_position(convex_hull_indices, :);
        % Convert world coordinates to pixel indices
        image_position_convhull_px =[
            image_sampling(2) * ...
                (image_position_convhull(:, 1) - image_bounds(1)) / image_bounds(3),...
            image_sampling(1) * ...
                (image_position_convhull(:, 2) - image_bounds(2)) / image_bounds(4)
        ];
        mask = roipoly(...
            I, image_position_convhull_px(:, 1), image_position_convhull_px(:, 2)...
        );
        x = linspace(image_bounds(1), image_bounds(1) + image_bounds(3), image_sampling(2));
        y = linspace(image_bounds(2), image_bounds(2) + image_bounds(4), image_sampling(1));
        [X,Y] = meshgrid(x,y);
        xy = [X(mask).'; Y(mask).'];
        I(mask) = fnval(image_spline,xy);

        if verbose
            figure
            surf(X, Y, I, 'EdgeColor', 'none');
            colorbar
            xlabel('X');
            ylabel('Y');
            zlabel('Irradiance');
            c = colorbar;
            c.Label.String = 'Irradiance';
            title('Estimated output image pixels') 
        end
    end
end

end

