function v_r = refract( ior_1, ior_2, n, v_i )
% REFRACT  Refract a 3D ray using Snell's law
%
% ## Syntax
% v_r = refract( ior_1, ior_2, n, v_i )
%
% ## Description
% v_r = refract( ior_1, ior_2, n, v_i )
%   Returns the direction of the refracted ray.
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
% ## References
% - Wikipedia page on Snell's Law:
%   https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

nargoutchk(1, 1);
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

end