function [ v_r, T ] = refract( ior_1, ior_2, n, v_i )
% REFRACT  Refract a 3D ray using Snell's law
%
% ## Syntax
% v_r = refract( ior_1, ior_2, n, v_i )
% [ v_r, T_s, T_p ] = refract( ior_1, ior_2, n, v_i )
%
% ## Description
% v_r = refract( ior_1, ior_2, n, v_i )
%   Returns the direction of the refracted ray.
% [ v_r, T ] = refract( ior_1, ior_2, n, v_i )
%   Additionally returns the transmissivity coefficient for the refracted
%   ray.
%
% ## Input Arguments
%
% ior_1 -- First index of refraction
%   A scalar equal to the index of refraction of the medium containing the
%   incident ray.
%
% ior_2 -- Image colour
%   A scalar equal to the index of refraction of the medium containing the
%   refracted ray.
%
% n -- Surface normal
%   An m x 3 array, where `m` is the number of rays to refract. `n(i, :)`
%   is a unit vector equal to the surface normal at the point of incidence
%   with the medium boundary for the i-th refraction. The surface normals
%   point towards the medium containing the incident rays.
%
% v_i -- Incident ray direction
%   An m x 3 array, where `v_i(i, :)` is a unit vector equal to the
%   direction of the incident ray for the i-th refraction. The incident
%   rays point towards the medium containing the refracted rays.
%
% ## Output Arguments
%
% v_r -- Refracted ray direction
%   An m x 3 array, where `v_r(i, :)` is a unit vector equal to the
%   direction of the refracted ray for the i-th refraction.
%
%   If the i-th angle of refraction exceeds the critical angle (i.e. total
%   internal reflection occurs), `v_r(i, :)` will consist of `NaN` values.
%
% T -- Transmissivity coefficient
%   A scalar between zero and one indicating the fraction of the incident
%   light's power which is transmitted as the refracted ray.
%
%   `T` is the average of the transmissivity coefficients for light
%   polarized in the plane of incidence, and perpendicular to the plane of
%   incidence; We assume the light ray has a uniform distribution of
%   polarizations.
%
% ## References
% - Wikipedia page on Snell's Law:
%   https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
% - M. Born and E. Wolf, Principles of Optics. Cambridge, U.K.;
%   New York: Cambridge University Press, 1999.
%   - Technically, I used section 1.5 in the 1970 edition of the book,
%     which covers the Fresnel equations.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

nargoutchk(1, 2);
narginchk(4, 4);

ior = ior_1 / ior_2;

cosine_i = - dot(n, v_i, 2);

radicand = 1 - (ior ^ 2) * (1 - (cosine_i .^ 2));
radicand(radicand < 0) = NaN; % Total internal reflection
v_r = ior * v_i + ...
    n .* repmat(...
        ior * cosine_i - sqrt(radicand), ...
        1, 3 ...
    );

if nargout > 1
    % Pages 41 - 45 in "Principles of Optics"
    theta_i = acos(cosine_i);
    theta_t = acos(- dot(n, v_r, 2));
    % Equation 35
    T_p = (sin(2 * theta_i) .* sin(2 * theta_t)) ./ ((sin(theta_i + theta_t)) .^ 2);
    T_s = T_p ./ ((cos(theta_i - theta_t)) .^ 2);
    T = (T_p + T_s) / 2;
    T(theta_i < eps) = 4 / (ior * ((1 / ior) + 1) ^ 2);
end

end