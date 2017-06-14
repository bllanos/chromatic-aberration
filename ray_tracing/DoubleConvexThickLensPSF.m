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

% ## Raytracing parameters
% Refer to the documentation of `doubleSphericalLens` for details

% ### Lens parameters
% Based on
% 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170613_SimpleLenses_EdmundOptics\25mmDiameter40mmFLUncoatedDoubleConvexLens_prnt_45296.pdf'
% and
% 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170613_SimpleLenses_EdmundOptics\25mmDiameter40mmFLUncoatedDoubleConvexLens.pdf'
ray_params.radius_front = 40.42;
lens_radius = 25 / 2;
axial_thickness = 5.30;
edge_thickness = 1.34;
aperture_thickness = (axial_thickness - edge_thickness) / 2;
ray_params.theta_aperture_front = atan2(...
    lens_radius, ray_params.radius_front - aperture_thickness...
);
ray_params.radius_back = ray_params.radius_front;
ray_params.theta_aperture_back = ray_params.theta_aperture_front;
ray_params.d_lens = axial_thickness - ray_params.radius_front - ray_params.radius_back;
ray_params.n_incident_rays = 500;
ray_params.sample_random = false;
ray_params.ior_environment = 1.0;
% The focal length specification wavelength is 587.6
% At this wavelength, N-BK7 glass has a refractive index of 1.51680
% (https://www.pgo-online.com/intl/katalog/BK7.html)
ray_params.ior_lens = 1.51680;

% ### Ray interpolation parameters
% Refer to the documentation of `densifyRays` for details
image_sampling = [400, 400];

% ## Scene setup

% Field of view over which point light sources are distributed
% (Actually, the field of view is larger than this, as this is the angle
% subtended by the midlines of the rectangular grid of lights, not the
% diagonals.)
% Set it to zero to use a single point light source instead.
scene_angle = pi / 4;
single_source = (scene_angle == 0);
if single_source
    % Angle of a single light source from the optical axis
    single_source_angle = pi / 6;
else
    % Number of lights in a grid of lights, along the x-dimension
    n_lights_x = 3;
    % Number of lights in a grid of lights, along the y-dimension
    n_lights_y = 3;
end

% Multiples of the focal length at which the light source plane will be
% positioned
light_distance_factor = 2;

% ## Debugging Flags
verbose_ray_tracing = (single_source | false);
verbose_ray_interpolation = (single_source | false);
display_each_psf = true;

%% Calculate lens imaging properties

[ imageFn, f, f_prime, U, U_prime, P ] = opticsFromLens(...
    ray_params.ior_environment,...
    ray_params.ior_lens,...
    ray_params.ior_environment,...
    ray_params.radius_front, ray_params.radius_back,...
    ray_params.d_lens...
);

%% Create light sources and the image plane

% Light source positions
z_light = (light_distance_factor * abs(f)) + U;
lens_center_z = lens_radius - (axial_thickness / 2);
lens_center_to_light = z_light - lens_center_z;
if single_source
    x_light = tan(single_source_angle) * lens_center_to_light;
    y_light = 0;
else
    scene_half_width = tan(scene_angle) * lens_center_to_light;
    x_light = linspace(-scene_half_width, scene_half_width, n_lights_x);
    y_light = linspace(-scene_half_width, scene_half_width, n_lights_y);
    [x_light,y_light] = meshgrid(x_light,y_light);
end
n_lights = numel(x_light);
X_lights = [
    x_light(:),...
    y_light(:),...
    repmat(z_light, n_lights, 1)
];

% Image plane position
X_image_ideal = imageFn([0, 0, z_light]);
ideal_image_position = X_image_ideal(3);
ray_params.d_film = -ideal_image_position;

% "Theoretical" image positions, on the in-focus plane
X_image_ideal = imageFn(X_lights);
% Adjust magnification to approximately account for the actual image plane
% location
magnification_correction_ideal = ...
    (repmat(U_prime, n_lights, 1) - X_image_ideal(:, 3)) ./ ...
    repmat(U_prime + ray_params.d_film, n_lights, 1);
principal_point = repmat([0, 0, U_prime], n_lights, 1);
X_image_rays = X_image_ideal - principal_point;
X_image_ideal = principal_point + (X_image_rays .* magnification_correction_ideal);

% Find image boundaries
if ~single_source
    X_image_grid_x = reshape(X_image_ideal(:, 1), n_lights_y, n_lights_x);
    X_image_grid_y = reshape(X_image_ideal(:, 2), n_lights_y, n_lights_x);
    left_buffer = max(abs(X_image_grid_x(:, 2) - X_image_grid_x(:, 1)));
    right_buffer = max(abs(X_image_grid_x(:, end) - X_image_grid_x(:, end - 1)));
    bottom_buffer = max(abs(X_image_grid_y(2, :) - X_image_grid_y(1, :)));
    top_buffer = max(abs(X_image_grid_y(end, :) - X_image_grid_y(end - 1, :)));
    image_bounds = [
        min(X_image_ideal(:, 1)) - left_buffer, min(X_image_ideal(:, 2)) - bottom_buffer
    ];
    image_width = max(X_image_ideal(:, 1)) - image_bounds(1) + right_buffer;
    image_height = max(X_image_ideal(:, 2)) - image_bounds(2) + top_buffer;
    image_bounds = [
        image_bounds,...
        image_width image_height...
    ];
end

%% Trace rays through the lens and form rays into an image
I = zeros(image_sampling);
X_image_real = zeros(n_lights, 2);
max_irradiance = zeros(n_lights, 1);
for i = 1:n_lights
    ray_params.source_position = X_lights(i, :);
    [ ...
        image_position, ray_irradiance, ~, incident_position_cartesian ...
    ] = doubleSphericalLens( ray_params, verbose_ray_tracing );

    if single_source
        [ X_image_real(i, :), max_irradiance(i) ] = densifyRays(...
            incident_position_cartesian,...
            ray_params.radius_front,...
            image_position,...
            ray_irradiance,...
            verbose_ray_interpolation ...
        );
    else
        [ X_image_real(i, :), max_irradiance(i), I_i ] = densifyRays(...
            incident_position_cartesian,...
            ray_params.radius_front,...
            image_position,...
            ray_irradiance,...
            image_bounds, image_sampling,...
            verbose_ray_interpolation ...
        );
    
        I = I + (I_i ./ max(max(I_i)));
    end

    if display_each_psf
        figure
        ax = gca;
        imagesc(ax,...
            [image_bounds(1), image_bounds(1) + image_bounds(3)],...
            [image_bounds(2), image_bounds(2) + image_bounds(4)],...
            I_i...
        );
        ax.YDir = 'normal';
        xlabel('X');
        ylabel('Y');
        c = colorbar;
        c.Label.String = 'Irradiance';
        title(...
            sprintf('Estimated PSF for a point source at position [%g, %g, %g]',...
            X_lights(i, 1), X_lights(i, 2), X_lights(i, 3))...
        );
    end
end

%% Visualize the results

figure
ax = gca;
imagesc(...
    [image_bounds(1), image_bounds(1) + image_bounds(3)],...
    [image_bounds(2), image_bounds(2) + image_bounds(4)],...
    I...
);
ax.YDir = 'normal';
xlabel('X');
ylabel('Y');
c = colorbar;
c.Label.String = 'Irradiance';
hold on
scatter(X_image_ideal(:, 1), X_image_ideal(:, 2), [], 'g', 'o');
scatter(X_image_real(:, 1), X_image_real(:, 2), [], 'r', '.');
legend('Thick lens formula', 'Raytracing peaks');
title('Images of point sources seen through a thick lens')
hold off