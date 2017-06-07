function [x, y, z, theta, phi] = sphereSection( n_phi, n_theta, theta_max, r )
% SPHERESECTION  Sample points on a section of a spherical surface
%
% ## Syntax
% [x, y, z, theta, phi] = sphereSection( n_phi, n_theta, theta_max, r )
%
% ## Description
% [x, y, z, theta, phi] = sphereSection( n_phi, n_theta, theta_max, r )
%   Returns a grid of points on a spherical section, sampled such that they
%   define sections of equal area.
%
%   The number of output arguments can be any number from 1 to 5.
%
% ## Input Arguments
%
% n_phi -- Number of azimuth angle samples
%   The number of segments into which each ring in the sampling pattern is
%   divided.
%
% n_theta -- Number of elevation angle samples
%   The number of rings of sections in the grid of samples.
%
% theta_max -- Maximum elevation angle
%   The lower boundary of the spherical surface, expressed in terms of
%   angular distance from the pole.
%
% r -- Sphere radius
%   The radius of the sphere
%
% ## Output Arguments
%
% Sample point coordinates are output as 2D arrays, of dimension `n_theta +
% 1` x `n_phi + 1`. The first row of the arrays corresponds to the pole.
% The last column of the arrays defines the same positions as the first
% column (as `phi = 0` and `phi = 2 * pi` are equivalent).
%
% The sphere is oriented such that its pole is on the positive z-axis.
%
% The output arguments represent the cartesian (x, y, z) and polar (theta,
% phi, r) coordinates of the sample points.
%
% ## Notes
% - The spherical section can be plotted using the following command:
%   `surf(x, y, z)`.
%
% ## References
% - Uniform sampling of the surface of a sphere:
%   http://mathworld.wolfram.com/SpherePointPicking.html
%
% See also sphere, surf

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

nargoutchk(1, 5);
narginchk(4, 4);

phi = linspace(0, 2 * pi, n_phi + 1);
cos_theta = (cos(theta_max) - 1) / n_theta;
theta = acos(1 + (1:n_theta) * cos_theta);
[phi, theta] = meshgrid(phi, theta);
% Add the pole
phi = [
    zeros(1, size(phi, 2));
    phi
    ];
theta = [
    zeros(1, size(theta, 2));
    theta
    ];
x = r * sin(theta) .* cos(phi);
y = r * sin(theta) .* sin(phi);
z = r * cos(theta);
end