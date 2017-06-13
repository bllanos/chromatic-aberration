function [ ...
    image_position, ray_irradiance, ...
    incident_position, incident_position_cartesian ...
] = doubleSphericalLens( params, varargin )
% DOUBLESPHERICALLENS  Trace rays through a lens with two spherical surfaces
%
% ## Syntax
% [image_position, ray_irradiance] = doubleSphericalLens( params [, verbose] )
%
% ## Description
% [image_position, ray_irradiance] = doubleSphericalLens( params [, verbose] )
%   Sample image positions and intensities by raytracing.
%
% [ ...
%     image_position, ray_irradiance, ...
%     incident_position ...
% ] = doubleSphericalLens( params, varargin )
%   Additionally returns the incidence locations of the rays, in angular
%   coordinates.
%
% [ ...
%     image_position, ray_irradiance, ...
%     incident_position, incident_position_cartesian ...
% ] = doubleSphericalLens( params, varargin )
%   Additionally returns the incidence locations of the rays, in
%   cartesian coordinates.
%
% ## Input Arguments
%
% params -- Parameters structure
%   A structure with the following fields:
%   - source_position: The 3D position of a point light source
%   - radius_front: The radius of curvature of the front surface of the
%     lens
%   - theta_aperture_front: The limit of the front aperture of the lens,
%     expressed as the angle with the optical axis subtended by the
%     aperture. In other words, half the angle subtended by the aperture
%     when viewed from the center of a circle with radius equal to the
%     radius of the front surface.
%   - radius_back: The radius of curvature of the back surface of the
%     lens
%   - theta_aperture_back: Similar to `theta_aperture_front`, but for the
%     back surface of the lens.
%   - d_lens: The distance from the centre of the sphere corresponding to
%     the front of the lens to the centre of the sphere corresponding to
%     the back of the lens. A shift in the positive direction moves the
%     back of the lens in the negative z-direction.
%   - n_incident_rays: Number of samples over the front of the lens. Each
%     sample produces one ray from the point light source through the front
%     surface of the lens. The front aperture is uniformly sampled, but
%     samples are culled if they are occluded by the front lens surface
%     from the perspective of the point light source.
%   - ior_environment: The refractive index of the surrounding medium
%   - ior_lens: The refractive index of the lens
%   - d_film: The distance of the image plane from the centre of the sphere
%     corresponding to the front surface of the lens.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% image_position -- Image coordinates of rays
%   The points of intersection of the light paths with the image plane. The
%   image plane is parallel to the xy-plane, so `image_position` is a
%   two-column array, with the columns containing x, and y coordinates,
%   respectively.
%
% ray_irradiance -- Ray irradiance
%   The incident irradiance at the points of intersection of the light
%   paths with the image plane. Irradiance is determined by:
%   - The foreshortening of the front aperture from the perspective of
%     the light source. (This is necessary to account for the uniform
%     sampling of the front aperture's surface area, which induces a
%     non-uniform sampling of solid angles from the perspective of the
%     light source.)
%   - The Fresnel equations, during refraction (in calls to `refract()`).
%   - The foreshortening of the ray emitted from the lens onto the image
%     plane.
%
%   Ray irradiances are not equal to image intensities, because the change
%   in ray density, between the front aperture and the image plane, is not
%   taken into account. Rather, the image intensity (irradiance) is the sum
%   of the ray irradiances at each image location. (Computation of image
%   irradiance is the purpose of 'densifyRays.m'.)
%
%   A ray which has no forshortening, nor attenuation during refraction,
%   has an irradiance value of unity.
%
%   `ray_irradiance` is a vector, where the i-th element corresponds to the
%   i-th row of `image_position`.
%
% incident_position -- Sampling positions on the front aperture
%   The positions, expressed in terms of the two angular coordinates
%   (theta, phi), of the incident rays on the front aperture of the lens.
%   Each row of `incident_position` corresponds to a row of
%   `image_position`.
%
% incident_position_cartesian -- Cartesian sampling positions on the front aperture
%   The positions, expressed in cartesian coordinates (x, y, z), of the
%   incident rays on the front aperture of the lens. Each row of
%   `incident_position_cartesian` corresponds to a row of `image_position`.
%
% ## Notes
% - For now, both lens surfaces are assumed to be convex (pointing outwards
%   from the centre of the lens).
% - Coordinate system:
%   - The front of the lens is spherical, with the centre of the sphere at
%     the origin.
%   - The positive z-axis points towards the front of the lens, along
%     the optical axis.
%   - In the code, 'theta' is the angle between a direction and the
%     positive z-axis. 'phi' is the angle between a direction and the
%     positive x-axis.
% - The number of output image positions is generally smaller than the
%   number of samples (`params.n_incident_rays`), because rays are culled
%   if they undergo total internal reflection, miss the lens apertures, or
%   are occluded by the front surface of the lens.
%
% ## References
% - Wikipedia page on Snell's Law:
%   https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
%
% See also sphereSection, refract, densifyRays

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 8, 2017

nargoutchk(2, 4);
narginchk(1, 2);

% Parse input arguments

s_position = params.source_position;
radius_front = params.radius_front;
theta_aperture_front = params.theta_aperture_front;
radius_back = params.radius_back;
theta_aperture_back = params.theta_aperture_back;
d_lens = params.d_lens;
n_incident_rays = params.n_incident_rays;
sample_random = params.sample_random;
ior_environment = params.ior_environment;
ior_lens = params.ior_lens;
film_z = -params.d_film;

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

% Visualization parameters

% Sampling resolution for lens surface theta values
n_theta_vis = 10;

% Sampling resolution for lens surface phi values
n_phi_vis = 20;

% Trace rays through the lens

cos_theta_aperture_front = cos(theta_aperture_front);

% Point of incidence on the front of the lens, in spherical polar
% coordinates (theta, phi)
if sample_random
    % Uniform random sampling:
    % See http://mathworld.wolfram.com/SpherePointPicking.html
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

else
    % Deterministic sampling
    sampling_ratio = theta_aperture_front / (2 * pi);
    n_phi = ceil(sqrt(n_incident_rays / sampling_ratio));
    n_theta = ceil(sqrt(n_incident_rays * sampling_ratio));
    [
        incident_normal_x, incident_normal_y, incident_normal_z, ...
        incident_position_theta, incident_position_phi ...
    ] = sphereSection( n_phi, n_theta, theta_aperture_front, 1 );

    incident_normal_x = incident_normal_x(2:end, 1:(end - 1));
    incident_normal_y = incident_normal_y(2:end, 1:(end - 1));
    incident_normal_z = incident_normal_z(2:end, 1:(end - 1));
    incident_position_theta = incident_position_theta(2:end, 1:(end - 1));
    incident_position_phi = incident_position_phi(2:end, 1:(end - 1));
    
    incident_position = [
        0, 0;
        incident_position_theta(:), incident_position_phi(:)
    ];

    % Surface normal at the point of incidence on the front of the lens
    incident_normal = [
        0, 0, 1;
        incident_normal_x(:), incident_normal_y(:), incident_normal_z(:)
    ];

    n_incident_rays = size(incident_normal, 1);
end

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

% Irradiance transmitted by the ray
% Energy is proportional to the solid angle of the patch of area on
% the aperture, from the perspective of the source. The solid angle is
% proportional to the cosine of the foreshortening angle, assuming the
% patch is approximately flat.
ray_irradiance = incident_cosine;

% Display the input data
if verbose
    figure
    scatter3(...
        incident_position_cartesian(:, 1),...
        incident_position_cartesian(:, 2),...
        ray_irradiance, [], ray_irradiance, 'filled'...
    )
    colorbar
    % axis equal
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    c = colorbar;
    c.Label.String = 'Ray irradiance';
    title('Rays incident on the front surface of the lens')
end

% First refracted ray direction
[ internal_direction, T ] = refract(...
    ior_environment, ior_lens, incident_normal, incident_direction...
    );
ray_irradiance = ray_irradiance .* T;

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
[ emitted_direction, T ] = refract(...
    ior_lens, ior_environment, -emitted_normal_culled, internal_direction...
    );
ray_irradiance = ray_irradiance .* T;

% Intersection with the film
distance_to_film = film_z - emitted_position_cartesian_culled(:, 3);
steps_to_film = distance_to_film ./ emitted_direction(:, 3);
steps_to_film(steps_to_film < 0) = NaN;
image_position = emitted_position_cartesian_culled +...
    repmat(steps_to_film, 1, 3) .* emitted_direction;

% Foreshortening at the image
ray_irradiance = ray_irradiance .* -emitted_direction(:, 3);

% 3D visualization
if verbose
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

    % Colour rays according to their irradiance
    ray_irradiance_max = max(ray_irradiance);
    ray_irradiance_min = min(ray_irradiance);
    ray_irradiance_scaled = (ray_irradiance - ray_irradiance_min) / (ray_irradiance_max - ray_irradiance_min);
    colors = [
        ray_irradiance_scaled,...
        zeros(n_incident_rays, 1),...
        1 - ray_irradiance_scaled
    ];
    colors(~isfinite(colors)) = 0.5;

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
end

% Finalize output arguments
output_filter = all(isfinite(image_position), 2);
image_position = image_position(output_filter, 1:2);
ray_irradiance = ray_irradiance(output_filter);
incident_position = incident_position(output_filter, :);
incident_position_cartesian = incident_position_cartesian(output_filter, :);

end