function [ ellipseValueFun, ellipseBoundaryFun, ellipse_to_world ] = ellipseModel(...
    lambda, theta, t, k, lightness, varargin...
)
% ELLIPSEMODEL  Create functions for evaluating a parametric ellipse
%
% ## Syntax
% ellipseValueFun = ellipseModel(...
%   lambda, theta, t, k, lightness [, image_mode]...
% )
% [ ellipseValueFun, ellipseBoundaryFun ] = ellipseModel(...
%   lambda, theta, t, k, lightness [, image_mode]...
% )
% [ ellipseValueFun, ellipseBoundaryFun, ellipse_to_world ] = ellipseModel(...
%   lambda, theta, t, k, lightness [, image_mode]...
% )
%
% ## Description
% ellipseValueFun = ellipseModel(...
%   lambda, theta, t, k, lightness [, image_mode]...
% )
%   Returns a function giving the value of the ellipse at points in the
%   plane.
%
% [ ellipseValueFun, ellipseBoundaryFun ] = ellipseModel(...
%   lambda, theta, t, k, lightness [, image_mode]...
% )
%   Additionally returns a function giving a discrete indication of a
%   point's location in relation to the ellipse.
%
% [ ellipseValueFun, ellipseBoundaryFun, ellipse_to_world ] = ellipseModel(...
%   lambda, theta, t, k, lightness [, image_mode]...
% )
%   Additionally returns a matrix transformation mapping from the unit disk
%   to the ellipse.
%
% ## Input Arguments
%
% lambda -- Ellipse dimensions
%   A two-element vector containing the major and minor semi-axis lengths
%   of the ellipse, respectively.
%
% theta -- Ellipse orientation
%   The angle in radians from the positive x-axis to the major axis of the
%   ellipse.
%
% t -- Ellipse centre
%   A two-element vector containing the coordinates of the centre of the
%   ellipse.
%
% k -- Ellipse edge width
%   The linear lightness transition between the ellipse and its
%   surroundings extends 1/k units both inside and outside of the geometric
%   boundary of the ellipse.
%
% lightness -- Ellipse values
%   A two-element vector containing the lightness within, and outside the
%   ellipse, respectively.
%
% image_mode -- Image axes conventions
%   If `true`, the coordinate system will be interpreted as image
%   coordinates (reversed y-axis). For instance, the (0, 1) point on the
%   unit circle will map to the point on the ellipse directly above its
%   centre in the image (i.e. at lower pixel y-coordinates).
%
%   Defaults to `false` if not passed.
%
% ## Output Arguments
%
% ellipseValueFun -- Ellipse lightness function
%   A function which accepts an n x 3 array of homogenous 2D coordinates
%   (with the homogenous coordinate scaled to 1), and returns the lightness
%   values of the ellipse model at those coordinates, as a column vector of
%   length 'n'. The input coordinates are in "world space".
%
% ellipseBoundaryFun -- Ellipse boundary function
%   A function which accepts an n x 3 array of homogenous 2D coordinates
%   (with the homogenous coordinate scaled to 1), and returns a column
%   vector of length 'n'. Entries in the vector are -1, 0, or 1 depending
%   on whether the corresponding rows of the input array are positions
%   inside, on the boundary, or outside of the ellipse.
%
% ellipse_to_world -- Ellipse coordinate transformation
%   A matrix transforming 2D points, 'p', on the unit disk, to points on
%   the ellipse. Points on the unit circle will be mapped to points on the
%   boundary of the ellipse. Specifically, a point (x1, y1) and its image
%   (x2, y2) on the ellipse are related by the equation:
%
%     [x2 y2 1] = (ellipse_to_world * [x1 y1 1].').'
%
% ## References
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.
%
% See also plotEllipse

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 20, 2018

nargoutchk(1, 3);
narginchk(5, 6);

if ~isempty(varargin)
    image_mode = varargin{1};
else
    image_mode = false;
end

if image_mode
    theta = -theta;
end
ellipse_to_world = [
        1, 0, t(1);
        0, 1, t(2);
        0, 0,   1
    ] * [
        cos(theta), -sin(theta), 0;
        sin(theta),  cos(theta), 0;
        0, 0, 1
    ] * [
        lambda(1), 0, 0;
        0, lambda(2), 0;
        0, 0, 1
    ];
if image_mode
    ellipse_to_world = ellipse_to_world * [1 0 0; 0 -1 0; 0 0 1];
end
world_to_ellipse = inv(ellipse_to_world);

boundary = [1 - 1/k, 1 + 1/k];
boundary_width = diff(boundary);

    function [ b ] = ellipseBoundary(p)
        p_ellipse = (world_to_ellipse * p.').';
        p_ellipse = p_ellipse(:, 1:2);
        p_mag = sqrt(dot(p_ellipse, p_ellipse, 2));
        b = zeros(size(p, 1), 1);
        b(p_mag < boundary(1)) = -1;
        b(p_mag > boundary(2)) = 1;
    end

    function [ v ] = ellipseValue(p)
        p_ellipse = (world_to_ellipse * p.').';
        p_ellipse = p_ellipse(:, 1:2);
        p_mag = sqrt(dot(p_ellipse, p_ellipse, 2));
        v = zeros(size(p, 1), 1);
        filter_1 = p_mag < boundary(1);
        v(filter_1) = lightness(1);
        filter_2 = p_mag > boundary(2);
        v(filter_2) = lightness(2);
        filter_3 = ~(filter_1 | filter_2);
        interp = (p_mag(filter_3) - boundary(1)) / boundary_width;
        v(filter_3) = (1 - interp) * lightness(1) + interp * lightness(2);
    end

ellipseValueFun = @ellipseValue;
ellipseBoundaryFun = @ellipseBoundary;

end

