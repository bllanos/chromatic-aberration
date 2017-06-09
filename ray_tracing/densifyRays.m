function [ I ] = densifyRays( incident_position, incident_position_cartesian, image_position, ray_power, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% Parameters: 
image_sampling = [100, 200]; % Y, X sample counts
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

% Display the input data
if verbose
    figure
    scatter3(...
        image_position(:, 1),...
        image_position(:, 2),...
        ray_power, [], ray_power, 'filled'...
    )
    colorbar
    % axis equal
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    c = colorbar;
    c.Label.String = 'Ray power';
    title('Sparse image points with power values')
end

dt = delaunayTriangulation(image_position);

% Plot the triangulation
% if verbose
%     figure
%     triplot(dt)
%     xlabel('X');
%     ylabel('Y');
% end

% Find triangle areas
triangles = dt.ConnectivityList;
vertices = dt.Points;
u = vertices(triangles(:, 1), :);
v = vertices(triangles(:, 2), :);
w = vertices(triangles(:, 3), :);
uv = v - u;
uw = w - u;
triangles_area = ((uv(:, 1) .* uw(:, 2)) - (uv(:, 2) .* uw(:, 1))) / 2;
centroids = mean(cat(3, u, v, w), 3);

% Interpolate ray power values and resample on a grid
min_x = min(image_position(:, 1));
max_x = max(image_position(:, 1));
min_y = min(image_position(:, 2));
max_y = max(image_position(:, 2));

x = linspace(min_x, max_x, image_sampling(2));
y = linspace(min_y, max_y, image_sampling(1));
[X,Y] = meshgrid(x,y);

power_interpolant = scatteredInterpolant(image_position, ray_power, 'natural');
I_power = power_interpolant(X,Y);
if verbose
    figure
    surf(X, Y, I_power, 'EdgeColor', 'none');
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Ray Power');
    c = colorbar;
    c.Label.String = 'Ray Power';
    title('Interpolation of ray power values')
end

% Interpolate triangle areas and resample on a grid
area_interpolant = scatteredInterpolant(centroids, triangles_area, 'natural');
I_area = area_interpolant(X,Y);
if verbose
    figure
    surf(X, Y, I_area, 'EdgeColor', 'none');
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Area');
    c = colorbar;
    c.Label.String = 'Triangle Area';
    title('Interpolation of Delaunay triangulation triangle areas')
end

% Compute intensity values
I = I_power ./ I_area;
if verbose
    figure
    surf(X, Y, I, 'EdgeColor', 'none');
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Intensity');
    c = colorbar;
    c.Label.String = 'Image Intensity';
    title('Estimated image intensity')
end

% TODO
% Version 1 (done):
%   Just use original points
% Version 2 (done):
%   Delaunay triangulate the image points
%   Find the area of each triangle
% Version 3: Kernel density estimation, with Gaussian kernel width equal to the
% square root of the (convex hull area of the Delaunay triangulation, divided
% by the (number of rays * pi)).
%
% Scattered data interpolation with linear interpolation, using scatteredInterpolant

end

