function [ max_position, max_irradiance, I ] = densifyRays(...
    incident_position_cartesian, r_front, image_position, ray_irradiance,...
    varargin...
)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Parameters: 
% Varargin = image_bounds crop area (rectangle of world coordinates), image_sampling, verbose
% Assume `incident_position_cartesian` is on a sphere with radius `r_front`, centered at the
% origin

% See also doubleSphericalLens, tpaps, fmincon

nargoutchk(1, 3);
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
    
    if nargout == 3 && n_varargs < 2
        error('The output image, `I`, cannot be calculated without the input arguments `image_bounds`, and `image_sampling`.');
    elseif nargout == 3
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
    function v_adj = neighborVertices(TR, vi)
        ti = vertexAttachments(TR,vi);
        n_vi = length(vi);
        v_adj = cell(n_vi, 1);
        tri = TR.ConnectivityList;
        for k = 1:n_vi
            v_adj_k = tri(ti{k}, :);
            v_adj_k = unique(v_adj_k(:));
            v_adj{k} = v_adj_k(v_adj_k ~= vi(k));
        end
    end

vi = (1:n_rays).';
v_adj_in = neighborVertices(dt_in, vi);
v_adj_out = neighborVertices(dt_out, vi);

% Find the average areas of circles with radii defined between neighbouring
% sample points on the front aperture
mean_areas_in = zeros(n_rays, 1);
r_front_sq = r_front ^ 2;
for i = 1:n_rays
    v_adj_i = v_adj_in{i};
    n_adj = length(v_adj_i);
    v_rep = repmat(incident_position_cartesian(i, :), n_adj, 1);
    cosines = dot(v_rep, incident_position_cartesian(v_adj_i, :), 2) / r_front_sq;
    areas = 2 * pi * r_front_sq * (1 - cosines);
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
thin_plate_spline = tpaps(image_position.',image_irradiance.');
if verbose
    figure
    fnplt(thin_plate_spline)
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Irradiance');
    c = colorbar;
    c.Label.String = 'Irradiance';
    title('Interpolation of image irradiance values')
end

% Find the peak intensity
thin_plate_spline_dx = fnder(thin_plate_spline,[1, 0]);
thin_plate_spline_dy = fnder(thin_plate_spline,[0, 1]);
% if verbose
%     figure
%     fnplt(thin_plate_spline_dx)
%     colorbar
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Irradiance X-gradient');
%     title('Image irradiance interpolant X-gradient')
%     
%     figure
%     fnplt(thin_plate_spline_dy)
%     colorbar
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Irradiance Y-gradient');
%     title('Image irradiance interpolant Y-gradient')
% end
min_x = min(image_position(:, 1));
max_x = max(image_position(:, 1));
min_y = min(image_position(:, 2));
max_y = max(image_position(:, 2));
    function [f, g] = splineValueAndGradient(x)
        f = -fnval(thin_plate_spline,[x(1);x(2)]);
        if nargout > 1 % gradient required
            g = [
                fnval(thin_plate_spline_dx,[x(1);x(2)]);
                fnval(thin_plate_spline_dy,[x(1);x(2)])
                ];
        end
    end
lb = [min_x, min_y];
ub = [max_x, max_y];
A = [];
b = [];
Aeq = [];
beq = [];
% Initial guess for maximum is the maximum of `ray_irradiance`, not
% `image_irradiance`, as `image_irradiance` is noisy.
[~, max_ind] = max(ray_irradiance);
x0 = image_position(max_ind, :);
nonlcon = [];
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
[max_position, max_irradiance] = fmincon(...
    @splineValueAndGradient,x0,A,b,Aeq,beq,lb,ub, nonlcon, options...
);
max_irradiance = -max_irradiance;

% Sample on a grid to produce an image
if output_image
    x = linspace(image_bounds(1), image_bounds(1) + image_bounds(3), image_sampling(2));
    y = linspace(image_bounds(2), image_bounds(2) + image_bounds(4), image_sampling(1));
    [X,Y] = meshgrid(x,y);
    xy = [X(:).'; Y(:).'];
    I = fnval(thin_plate_spline,xy);
    I = reshape(I, image_sampling);
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

