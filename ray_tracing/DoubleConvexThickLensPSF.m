%% Ray tracing simulation of chromatic aberration
% Simulate the chromatic point spread function of a thick (biconvex) lens.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
%
% ## References
% - Ray-sphere intersection testing:
%   - http://www.ccs.neu.edu/home/fell/CSU540/programs/RayTracingFormulas.htm
%   - https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
%   - https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
% - Uniform sampling of the surface of a sphere:
%   http://mathworld.wolfram.com/SpherePointPicking.html

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

%% Input data and parameters

% Coordinate system:
% - The front of the lens is spherical, with the centre of the sphere at
%   the origin.
% - The positive z-axis points towards the front of the lens, along
%   the optical axis.
% - 'theta' is the angle between a direction and the positive z-axis.
% - 'phi' is the angle between a direction and the positive x-axis

% ## Raytracing parameters

% Location of light source
s_position = [3, 3, 10];

% Radius of the front of the lens
radius_front = 1.0;

% Exposed angular range of the front of the lens (i.e. aperture)
theta_aperture_front = pi / 4;

% Radius of the back of the lens
radius_back = 2.0;

% Exposed angular range of the back of the lens (i.e. aperture)
theta_aperture_back = pi / 6;

% Separation between the centres of the spheres corresponding to the front
% and back of the lens. (A shift in the positive direction moves the back
% of the lens in the negative z-direction.)
d_lens = -1;

% Number of samples over the front of the lens
n_incident_rays = 100;

% Environment index of refraction
ior_environment = 1.0;

% Lens index of refraction relative to its surroundings
ior_lens = 1.52;

% z-position of photographic film parallel to the xy-plane
film_z = -10;

% ## Visualization parameters

% Sampling resolution for lens surface theta values
n_theta_vis = 10;

% Sampling resolution for lens surface phi values
n_phi_vis = 20;

%% Trace rays through the lens

% Total surface area defined by the aperture
cos_theta_aperture_front = cos(theta_aperture_front);
aperture_area = 2 * pi * (1 - cos_theta_aperture_front);

% Average area for one sample
ray_area = aperture_area / n_incident_rays;

% Point of incidence on the front of the lens, in spherical polar
% coordinates (theta, phi)
% Uniform sampling: See http://mathworld.wolfram.com/SpherePointPicking.html
incident_position = [
    acos(cos_theta_aperture_front + rand(n_incident_rays, 1)*(1 - cos_theta_aperture_front)),...
    2 * pi * rand(n_incident_rays, 1)
    ];

% Surface normal at the point of incidence on the front of the lens
incident_normal = [
    sin(incident_position(:, 1)) .* cos(incident_position(:, 2)),...
    sin(incident_position(:, 1)) .* sin(incident_position(:, 2)),...
    cos(incident_position(:, 1))
    ];

% Point of incidence on the front of the lens
incident_position_cartesian = radius_front * incident_normal;

% Incident ray direction
s_position = repmat(s_position, n_incident_rays, 1);
incident_direction = incident_position_cartesian - s_position;
incident_direction = incident_direction ./ repmat(...
    sqrt(dot(incident_direction, incident_direction, 2)), 1, 3 ...
);

% Cosine of angle with surface normal at the point of incidence
incident_cosine = - dot(incident_normal, incident_direction, 2);

% Cull rays which are occluded by the lens surface before the point of
% incidence
front_occlusion_filter = (incident_cosine <= 0);
front_occlusion_filter = repmat(front_occlusion_filter, 1, 3);
incident_direction(front_occlusion_filter) = NaN;
incident_position_cartesian(front_occlusion_filter) = NaN;

% Power transmitted by the ray
ray_power = ray_area * incident_cosine;

% First refracted ray direction
internal_direction = refract(...
    ior_environment, ior_lens, incident_normal, incident_direction...
    );

% Centre of the sphere representing the back of the lens
back_center = [0, 0, -d_lens];

% Point of emission at the back of the lens
% From
% http://www.ccs.neu.edu/home/fell/CSU540/programs/RayTracingFormulas.htm,
% https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm,
% and https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
back_center_rep = repmat(back_center, n_incident_rays, 1);
b = dot(internal_direction, (incident_position_cartesian - back_center_rep), 2);
c = incident_position_cartesian - back_center_rep;
c = dot(c, c, 2) - repmat(radius_back ^ 2, n_incident_rays, 1);
discriminant = (b .^ 2) - c;
discriminant(discriminant <= 0) = NaN;
t = -b + sqrt(discriminant);
t_minus = -b - sqrt(discriminant);
t(t < 0) = t_minus(t < 0);

emitted_position_cartesian = incident_position_cartesian +...
    repmat(t, 1, 3) .* internal_direction;

% Clip to points inside the back aperture
emitted_normal = emitted_position_cartesian - back_center_rep;
emitted_normal = emitted_normal ./ repmat(...
    sqrt(dot(emitted_normal, emitted_normal, 2)),...
    1, 3 ...
    );
emitted_cos_theta = -emitted_normal(:, 3);
back_aperture_filter = repmat(emitted_cos_theta < cos(theta_aperture_back), 1, 3);

emitted_position_cartesian_culled = emitted_position_cartesian;
emitted_position_cartesian_culled(back_aperture_filter) = NaN;
emitted_normal_culled = emitted_normal;
emitted_normal_culled(back_aperture_filter) = NaN;

% Refraction at the back of the lens
emitted_direction = refract(...
    ior_lens, ior_environment, -emitted_normal_culled, internal_direction...
    );

% Intersection with the film
distance_to_film = film_z - emitted_position_cartesian_culled(:, 3);
steps_to_film = distance_to_film ./ emitted_direction(:, 3);
steps_to_film(steps_to_film < 0) = NaN;
image_position = emitted_position_cartesian_culled +...
    repmat(steps_to_film, 1, 3) .* emitted_direction;

%% 3D visualization

figure
hold on

% Lens
[x, y, z] = sphereSection( n_phi_vis, n_theta_vis, theta_aperture_front, radius_front );
lens_front = surf(x, y, z);
set(lens_front, 'EdgeColor', 'g', 'FaceAlpha', 0.4, 'FaceColor', 'g');
[x, y, z] = sphereSection( n_phi_vis, n_theta_vis, theta_aperture_back, radius_back );
z = -d_lens - z;
lens_back = surf(x, y, z);
set(lens_back, 'EdgeColor', 'r', 'FaceAlpha', 0.4, 'FaceColor', 'r');

% Light Source
scatter3(s_position(:, 1), s_position(:, 2), s_position(:, 3), 'c', 'filled');

% Film
min_x = min(image_position(:, 1));
max_x = max(image_position(:, 1));
min_y = min(image_position(:, 2));
max_y = max(image_position(:, 2));

surf(...
    [min_x, max_x; min_x, max_x],...
    [max_y, max_y; min_y, min_y],...
    [film_z, film_z; film_z, film_z],...
    'EdgeColor', 'k', 'FaceAlpha', 0.4, 'FaceColor', 'y' ...
);

% Colour rays according to their power
ray_power_max = max(ray_power);
ray_power_min = min(ray_power);
ray_power_scaled = (ray_power - ray_power_min) / (ray_power_max - ray_power_min);
colors = [
    ray_power_scaled,...
    zeros(n_incident_rays, 1),...
    1 - ray_power_scaled
];

for i = 1:n_incident_rays

    % Incident Rays
    line(...
        [s_position(i, 1) incident_position_cartesian(i, 1)].',...
        [s_position(i, 2) incident_position_cartesian(i, 2)].',...
        [s_position(i, 3) incident_position_cartesian(i, 3)].',...
        'Color', colors(i, :)...
        );

    % Internal Rays
    line(...
        [incident_position_cartesian(i, 1) emitted_position_cartesian(i, 1)].',...
        [incident_position_cartesian(i, 2) emitted_position_cartesian(i, 2)].',...
        [incident_position_cartesian(i, 3) emitted_position_cartesian(i, 3)].',...
        'Color', colors(i, :)...
        );

    % Emitted Rays
    line(...
        [emitted_position_cartesian_culled(i, 1) image_position(i, 1)].',...
        [emitted_position_cartesian_culled(i, 2) image_position(i, 2)].',...
        [emitted_position_cartesian_culled(i, 3) image_position(i, 3)].',...
        'Color', colors(i, :)...
        );

end

hold off

% General axis properties
axis equal
set(gca, 'Color', 'none');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Raytracing from a point source through a double spherical lens')