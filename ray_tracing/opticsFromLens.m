function [ imageFn, f_out, f_prime_out, U_out, U_prime_out, P ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
% OPTICSFROMLENS  Obtain imaging quantities from lens geometry
%
% ## Syntax
% imageFn = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
% [ imageFn, f ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
% [ imageFn, f, f_prime ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
% [ imageFn, f, f_prime, U ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
% [ imageFn, f, f_prime, U, U_prime ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
% [ imageFn, f, f_prime, U, U_prime, P ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%
% ## Description
% imageFn = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%   Returns a function for calculating image locations corresponding to
%   object locations.
% [ imageFn, f ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%   Additionally returns the first focal length.
% [ imageFn, f, f_prime ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%   Additionally returns the second focal length.
% [ imageFn, f, f_prime, U ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%   Additionally returns the first principal plane.
% [ imageFn, f, f_prime, U, U_prime ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%   Additionally returns the second principal plane.
% [ imageFn, f, f_prime, U, U_prime, P ] = opticsFromLens( n0, n1, n2, r1, r2, d_lens )
%   Additionally returns the power of the lens
%
% ## Input Arguments
%
% n0 -- First index of refraction
%   A scalar equal to the index of refraction of the medium on the side of
%   the lens with higher z-coordinates.
%
% n1 -- Second index of refraction
%   A scalar equal to the index of refraction of the lens material.
%
% n2 -- Third index of refraction
%   A scalar equal to the index of refraction of the medium on the side of
%   the lens with lower z-coordinates.
%
% r1 -- Front surface radius
%   The radius of curvature of the front surface of the lens (the surface
%   facing the region of higher z-coordinates).
%
% r2 -- Back surface radius
%   The radius of curvature of the back surface of the lens (the surface
%   facing the region of lower z-coordinates).
%
% d_lens -- Lens thickness parameter
%   The distance from the centre of the sphere corresponding to the front
%   of the lens to the centre of the sphere corresponding to the back of
%   the lens. A shift in the positive direction moves the back of the lens
%   in the negative z-direction.
%
% ## Output Arguments
%
% imageFn -- Imaging function
%   A function of a single argument, `X_object`, which is an n x 3 array,
%   where the columns contain the x, y, and z-coordinates, respectively of
%   object points. The object points are assumed to be located at higher
%   z-coordinates than the front surface of the lens.
%
%   `imageFn` returns a n x 3 array, where the columns contain the x, y,
%   and z-coordinates of the image points formed by the lens, corresponding
%   to the object points `X_object`.
%
% f -- First focal length
%   The first focal length of the lens; The distance in the positive
%   z-direction from the first principal plane through which rays parallel
%   to the optical axis incident on the back surface of the lens are
%   refracted.
%
% f_prime -- Second focal length
%   The second focal length of the lens; The distance, in the negative
%   z-direction, from the second principal plane, through which rays
%   parallel to the optical axis incident on the front surface of the lens
%   are refracted.
%
% U -- First principal plane
%   The z-coordinate of the first principal plane of the lens.
%
% U_prime -- Second principal plane
%   The z-coordinate of the second principal plane of the lens.
%
% P -- Lens power
%   The power of the lens (in arbitrary units corresponding to the units of
%   the world coordinate system).
%
% ## Notes
% - Coordinate system:
%   - The front of the lens is spherical, with the centre of the sphere at
%     the origin.
%   - The radii of both faces of the lens are positive when the lens is
%     biconvex.
%   - The positive z-axis points towards the front of the lens, along
%     the optical axis, assuming the front of the lens is convex.
%
% ## References
% - M. Born and E. Wolf, Principles of optics, 4. ed. Oxford: Pergamon Press,
%   1970, isbn: 9780080139876.
%   - I used section 4.4.3, which covers the thick lens, and section 4.3.1,
%     which covers Newton's equation.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 13, 2017

% Focal lengths for light coming from medium 1 to medium 0 through an
% interface of radius r (r > 0 if the surface is convex towards the light)
    function f_s = sphereFocalLength(n0, n1, r)
        f_s = n0 * r / (n1 - n0);
    end

% Focal lengths of the first lens surface (Equation 19)
f0 = sphereFocalLength(n0, n1, r1);
f0_prime = sphereFocalLength(n1, n0, r1);

% Focal lengths of the first lens surface (Equation 20)
r2 = -r2;
f1 = sphereFocalLength(n1, n2, r2);
f1_prime = sphereFocalLength(n2, n1, r2);

% Combined focal lengths
t = r1 + d_lens - r2; % Lens axial thickness
c = t + f0_prime - f1; % (22)
f = - f0 * f1 / c; % (21)
f_prime = f0_prime * f1_prime / c; % (21)

% Lens power
P = n0 / f;

% Assuming the following:
% - The centre of the sphere corresponding to the first lens surface
%   is at the origin.
% - The positive z-axis points towards the first lens surface, if the first
%   lens surface is convex (has a positive radius).
% Then, the locations of the principal planes are:
D = ((n1 - n0) * (n2 - n1) * t) - ...
    n1 * ((n2 - n1) * r1 + (n1 - n0) * r2); % (24)
d = - n0 * (n2 - n1) * r1 * t / D;
d_prime = n2 * (n1 - n0) * r2 * t / D;
U = r1 - d;
U_prime = -d_lens + r2 - d_prime;

% Thick lens law (only holds if f == -f_prime)
%     function z_image = imagePoint(z_object)
%         z1 = z_object - U;
%         z2 = 1 / (1/f - 1/z1);
%         z_image = U_prime + z2;
%     end

    function X_image = imagePoint(X_object)
        n_points = size(X_object, 1);
        z_object = X_object(:, 3);
        z0 = repmat(U, n_points, 1) - z_object;
        % Equation 16 in section 4.3.1
        z1 = -f_prime ./ (1 + (f ./ z0));
        z_image = repmat(U_prime, n_points, 1) - z1;
        % Equation 10 in section 4.3.1, taking into account equation 15
        magnification = f ./ (z0 + f);
        X1 = X_object(:, 1) .* magnification;
        Y1 = X_object(:, 2) .* magnification;
        X_image = [X1, Y1, z_image];
    end

imageFn = @imagePoint;
% Account for differences in the coordinate system relative to that used in
% "Principles of Optics".
f_out = -f;
f_prime_out = -f_prime;
U_out = U;
U_prime_out = U_prime;

end

