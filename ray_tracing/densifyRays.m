function [ I ] = densifyRays( incident_position_cartesian, r_front, image_position, ray_irradiance, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Parameters: 
image_sampling = [100, 200]; % Y, X sample counts
% Image crop area (rectangle of world coordinates)
% Assume `incident_position_cartesian` is on a sphere with radius `r_front`, centered at the
% origin
%
% ## References:
% - Area of a triangle using cross products:
%   - http://mathworld.wolfram.com/TriangleArea.html
%   - http://mathworld.wolfram.com/CrossProduct.html

% See also scatteredInterpolant

nargoutchk(1, 1);
narginchk(4, 5);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
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

% Interpolate image irradiance values and resample on a grid
min_x = min(image_position(:, 1));
max_x = max(image_position(:, 1));
min_y = min(image_position(:, 2));
max_y = max(image_position(:, 2));

x = linspace(min_x, max_x, image_sampling(2));
y = linspace(min_y, max_y, image_sampling(1));
[X,Y] = meshgrid(x,y);

image_interpolant = scatteredInterpolant(image_position, image_irradiance, 'natural');
I = image_interpolant(X,Y);
if verbose
    figure
    surf(X, Y, I, 'EdgeColor', 'none');
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Irradiance');
    c = colorbar;
    c.Label.String = 'Irradiance';
    title('Interpolation of image irradiance values')
end

end

