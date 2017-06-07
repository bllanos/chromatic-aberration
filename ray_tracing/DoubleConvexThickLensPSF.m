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

% Location of light source
s_position = [0, 0, 10];

% Radius of the front of the lens
radius_front = 1.0;

% Exposed angular range of the front of the lens (i.e. aperture)
theta_aperture_front = pi / 4;

% Radius of the back of the lens
radius_back = 1.0;

% Exposed angular range of the back of the lens (i.e. aperture)
theta_aperture_back = pi / 4;

% Separation between the centres of the spheres corresponding to the front
% and back of the lens. (A shift in the positive direction moves the back
% of the lens in the negative z-direction.)
d_lens = 0;

% Number of samples over the front of the lens
n_incident_rays = 1;

% Environment index of refraction
ior_environment = 1.0;

% Lens index of refraction relative to its surroundings
ior_lens = 1.52;

% z-position of photographic film parallel to the xy-plane
film_z = -10;

% Point of incidence on the front of the lens, in spherical polar
% coordinates (theta, phi)
incident_position = [0, 0];

%% Trace rays through the lens

% Total surface area defined by the aperture
aperture_area = 2 * pi * (1 - cos(theta_aperture_front));

% Average area for one sample
ray_area = aperture_area / n_incident_rays;

% Surface normal at the point of incidence on the front of the lens
incident_normal = [
    sin(incident_position(1)) .* cos(incident_position(2)),...
    sin(incident_position(1)) .* sin(incident_position(2)),...
    cos(incident_position(1))
    ];

% Point of incidence on the front of the lens
incident_position_cartesian = radius_front * incident_normal;

% Incident ray direction
incident_direction = incident_position_cartesian - repmat(s_position, n_incident_rays, 1);
incident_direction = incident_direction ./ repmat(...
    sqrt(dot(incident_direction, incident_direction, 2)), 1, 3 ...
);

% Cosine of angle with surface normal at the point of incidence
incident_cosine = - dot(incident_normal, incident_direction, 2);

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
emitted_theta = acos(-emitted_normal(:, 3));
back_aperture_filter = repmat(emitted_theta <= theta_aperture_back, 1, 3);
emitted_position_cartesian(back_aperture_filter) = NaN;
emitted_normal(back_aperture_filter) = NaN;

% Refraction at the back of the lens
emitted_direction = refract(...
    ior_lens, ior_environment, emitted_normal, internal_direction...
    );

% Intersection with the film
distance_to_film = film_z - emitted_position_cartesian;
steps_to_film = distance_to_film ./ emitted_direction(:, 3);
steps_to_film(steps_to_film < 0) = NaN;
image_position = emitted_position_cartesian +...
    repmat(steps_to_film, 1, 3) .* emitted_direction;