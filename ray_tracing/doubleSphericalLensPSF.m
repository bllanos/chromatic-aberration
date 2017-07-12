function [...
    stats_real, stats_ideal, X_lights, depth_factors, I, I_color...
] = doubleSphericalLensPSF(...
    lens_params, ray_params, image_params, scene_params,...
    request_spline_smoothing, varargin...
)
% DOUBLESPHERICALLENSPSF  Generate images of point light sources by raytracing
%
% ## Syntax
% [...
%     stats_real, stats_ideal, X_lights, depth_factors, I, I_color...
% ] = doubleSphericalLensPSF(...
%     lens_params, ray_params, image_params, scene_params,...
%     request_spline_smoothing [, verbose]...
% )
%
% ## Description
% [...
%     stats_real, stats_ideal, X_lights, depth_factors, I, I_color...
% ] = doubleSphericalLensPSF(...
%     lens_params, ray_params, image_params, scene_params,...
%     request_spline_smoothing [, verbose]...
% )
%   Simulate the image of grids of point light sources by raytracing. (One
%   to six output arguments can be requested.)
%
% ## Input Arguments
%
% lens_params -- Lens parameters structure
%   A description of a lens formed from two spherical surfaces.
%   Passed as a structure with the following fields:
%   - lens_radius: The radius of the lens (i.e. half the height of the lens
%     when viewed edge-on)
%   - axial_thickness: The thickness of the lens along its optical axis.
%   - radius_front: The radius of curvature of the front surface of the
%     lens
%   - radius_back: The radius of curvature of the back surface of the
%     lens
%   - ior_lens: The refractive indices of the lens, one for each wavelength
%     of the light to be simulated. A row vector of length 'k'.
%   - ior_lens_reference_index: The index into `ior_lens` of the reference
%     index of refraction. The reference index of refraction is the index
%     of refraction used to set the image plane position for the desired
%     focal distance.
%   - wavelengths: The wavelengths of light corresponding to the elements
%     of `ior_lens`. A row vector of length 'k'. This parameter is used for
%     figure legends only, not for calculations.
%   - wavelengths_to_rgb: RGB colours to be used when displaying the
%     images generated by raytracing with the indices of refraction in
%     `ior_lens`. The i-th row of this k x 3 matrix represents the RGB
%     colour corresponding to the i-th wavelength. `wavelengths_to_rgb`
%     allows colour images to be produced by adding together the
%     contributions of each wavelength to the red, green, and blue colour
%     channels.
%
% ray_params -- Raytracing parameters structure
%   A structure with the following fields:
%   - n_incident_rays: The number of rays to sample, over the front
%     aperture of the lens, for each light source in the scene. Each sample
%     produces one ray from the point light source through the front
%     surface of the lens. The front aperture is uniformly sampled, but
%     samples are culled if they are occluded by the front lens surface,
%     from the perspective of the point light source.
%   - sample_random: Whether to sample rays at random over the front
%     aperture of the lens (true) or in a polar grid pattern (false). In
%     either case, the rays will be sampled uniformly per unit area on the
%     front aperture.
%   - ior_environment: The refractive index of the medium surrounding the
%     lens on both sides
%
% image_params -- Image formation parameters
%   A structure with the following fields, describing how to generate
%   images from raytracing results:
%   - image_sampling: A two-element vector containing the image height,
%     and width, respectively, of the final image (in units of pixels).
%   - normalize_psfs_before_combining: Some point spread functions are very
%     sharp, compared to others. To assist with visualization, the combined
%     image can be generated by normalizing each PSF by its maximum
%     intensity before adding the individual images together. However, the
%     final image will be a little misleading. Setting this parameter to
%     `false` will result in an image formed from raw PSF intensities.
%   - normalize_color_images_globally: Colour images are normalized so that
%     colour channel values are no greater than one. Normalization is
%     either performed by dividing by the maximum channel value over all
%     scene depths (true), or the maximum channel value of each scene depth
%     (false). Therefore, if this parameter is true, colour images produced
%     for different depths are directly comparable, but images for
%     out-of-focus depths may be quite dim.
%
% scene_params -- Light source parameters
%   A structure with the following fields, describing the grid or line of
%   light sources illuminating the lens:
%   - theta_max: For the grid of light sources at the reference depth, this
%     is the maximum angle with the optical axis subtended by the grid. The
%     angle is measured between the optical axis, and the rays from the
%     point where the optical axis intersects the first principal plane of
%     the lens to the four sides of the grid. These rays meet the sides of
%     the grid at right angles when viewed along the optical axis. In other
%     words, this is half of the angle subtended by the perpendicular
%     bisectors of the square grid of lights.
%
%     For a line of light sources at the reference depth, this is the angle
%     between the optical axis and the ray from the first principal plane's
%     intersection with the optical axis to the outermost light source on
%     the line.
%
%   - theta_min: Analogous to `theta_max`, but defines a minimum angle with
%     the optical axis.
%
%     Instead of a grid of light sources, a nonzero value of `theta_min`
%     will produce a rectangular frame of light sources, surrounding an
%     empty rectangle with sides defined by `theta_min`. Consequently, the
%     number of light sources will be less than `prod(n_lights)`. However,
%     the output data will have the same dimensions - NaN values will be
%     used for positions within the empty rectangle.
%
%     For a line of light sources, if `theta_min` is nonzero, the line will
%     simply start at a position offset by `theta_min`.
%
%   - n_lights: The number of lights in the scene.
%
%     If `n_lights` is a scalar, it is the number of lights on a line of
%     light sources. The line extends perpendicular to the optical axis,
%     along the positive x-direction.
%
%     If `n_lights` is a two-element vector, it contains the number of
%     lights in the x and y-directions, respectively, in a grid of lights.
%     The grid is centered on the optical axis.
%
%     If `n_lights` is `1`, or `[1 1]`, a single light source is created,
%     and placed on the y = 0 plane, with a positive x-coordinate, such
%     that it corresponds to the angle `theta_max`.
%
%   - light_distance_factor_focused: The reference depth, measured in
%     multiples of the first focal length of the lens. The "reference
%     depth" is the depth of the scene which produces a focused image,
%     according to the thick lens equation. The image plane is positioned
%     so that the lights at the reference depth are in-focus. Depth is
%     measured relative to the first principal plane of the lens.
%   - light_distance_factor_larger: Additional depths at which light
%     sources will be placed, that are greater than the reference depth.
%     Distances will be spaced uniformly in inverse depth space (with depth
%     measured from the first principal plane of the lens). This is a
%     two-element vector, where the first element contains the multiple of
%     lens focal lengths corresponding to the largest depth, and the second
%     element contains the number of depths larger than the reference
%     depth. If the second element is zero, no larger depths will be used
%     in the simulation. In the output arguments of this function, depths
%     are indexed such that higher indices represent larger depths.
%   - light_distance_factor_smaller: Equivalent to
%     `light_distance_factor_larger`, but for depths smaller than the
%     reference depth.
%   - preserve_angle_over_depths: If `true`, lights at depths other than
%     the reference depth will be positioned such that their images are
%     aligned, according to the thick lens formulas. If `false`, the grids
%     of lights at different depths will be aligned so that the lights have
%     the same xy-coordinates regardless of depth. In other words, the
%     lights will lie along rays parallel to the optical axis.
%
% request_spline_smoothing -- Smooth image intensities
%   A flag determining whether image intensities are subject to thin-plate
%   spline smoothing. For low numbers of rays, `request_spline_smoothing`
%   can be set to reduce noise, but for high numbers of rays, spline
%   smoothing will be computationally expensive.
%
%   Presently, if spline smoothing is not requested, but the image output
%   arguments are generated (refer to the notes on Efficiency, below),
%   spline smoothing will be triggered just for the purposes of producing
%   images.
%
% verbose -- Debugging flags
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   All debugging flags default to false if `verbose` is not passed.
%
% ## Output Arguments
%
% stats_real -- Simulated point spread function statistics
%   Point spread function statistics computed for the images produced by
%   each light source, for each wavelength, and for each depth, output as a
%   structure array. `stats_real(i, k, j)` is the `stats` output argument
%   of 'analyzePSF.m' for the i-th light, emitting the k-th wavelength, and
%   positioned at the j-th depth.
%
% stats_ideal -- Theoretical image statistics
%   Point spread function statistics for the images produced by each light
%   source, as predicted by the thick lens equation. `stats_ideal` has
%   the same format as `stats_real`.
%
% X_lights -- Light positions
%   The positions, expressed in cartesian coordinates (x, y, z), of the
%   lights in the scene. `X_lights(i, :, j)` is the 3D position of the i-th
%   light in the grid of lights placed at the j-th depth.
%
% depth_factors -- Depths expressed in focal lengths
%   The depths of the light sources, measured in multiples of the first
%   focal length of the lens from the first principal plane.
%   `depth_factors(j)` corresponds to the z-values in `X_lights(:, :, j)`.
%
% I -- Simulated greyscale images
%   A 4D array containing the images formed by combining all point spread
%   functions for all lights at each depth. `I(:, :, k, j)` is the
%   greyscale image for the k-th wavelength produced by the grid of lights
%   placed at the j-th depth.
%
%   This output argument is only available if there are multiple lights at
%   each depth.
%
% I_color -- Simulated colour images
%   A 4D array containing the images formed by combining all point spread
%   functions for all lights at each depth. `I_color(:, :, :, j)` is the
%   colour image for the produced by the grid of lights placed at the j-th
%   depth. The third dimension of `I_color` indexes Red, Green, and Blue
%   colour channels. `I_color` is produced by combining the images for
%   individual wavelengths in `I` according to the RGB colour values for
%   the different wavelengths in `lens_params.wavelengths_to_rgb`.
%
%   This output argument is only available if there are multiple lights at
%   each depth.
%
% ## Notes
%
% ### Coordinate System
% - The radii of both faces of the lens are positive when the lens is
%   biconvex.
%
% ### Polychromatic light sources
%
% If I wanted to simulate real colour images, I would model point spread
% functions for colour channels, not for individual wavelengths. I could do
% so by adding together the point spread functions for individual
% wavelengths, weighted by the quantum efficiencies of the camera sensor
% for the individual wavelengths, and weighted by the spectral power
% distribution of the light source. Such an approach relies on the
% following assumptions:
%
% - The camera sensor does not perform any on-chip processing to adjust the
%   relative responses of pixels based on assumptions of the light colour.
% - The scene is emitting light with the given spectral power distribution.
%
% I think it is cleaner to assume narrowband sensor responses for now, so
% that a point spread function for an individual wavelength corresponds
% exactly to the image recorded in the nearest colour channel.
%
% Presently, I output colour images by mixing together point spread
% functions for individual wavelengths, but this is purely for
% visualization. All mathematical analysis of chromatic aberration is done
% under narrowband sensor response assumptions. An advantage of this
% segregation is it is possible to estimate the theoretical chromatic
% aberration between individual wavelengths, ignoring the properties of the
% camera sensor. Furthermore, it is possible to run this script on an
% arbitrary number of wavelengths, and obtain as many estimates of
% chromatic aberration (not only three).
%
% ### Efficiency
% - The `I` or `I_color` output arguments are expensive to produce. They
%   are generated if explicitly requested as output, or if any of the
%   following flags are set:
%   - `verbose.display_each_psf`
%   - `verbose.display_all_psf_each_ior`
%   - `verbose.display_all_psf_each_depth`
%
%   Of course, for a single light source, these output arguments are
%   unavailable, as mentioned above. The image boundaries are calculated
%   from the spacing between the image positions of multiple sources,
%   predicted by the thick lens equation. For a single light source, point
%   spread functions can be visualized using `verbose.display_each_psf`,
%   but there is presently no convenient way to define image boundaries.
% - Thin-plate spline smoothing is prohibitively expensive for large
%   numbers of rays (as discussed further in the documentation of
%   `request_spline_smoothing`, above).
% - For a single light source, `verbose.display_each_psf` will presently
%   trigger a scattered data interpolation operation, if
%   `request_spline_smoothing` is `false`.
%
% See also opticsFromLens, doubleSphericalLens, densifyRays, analyzePSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 22, 2017

%% Parse parameters

nargoutchk(1,6);
narginchk(5,6);

ray_params.radius_front = lens_params.radius_front;
ray_params.theta_aperture_front = asin(...
    lens_params.lens_radius / lens_params.radius_front...
);
ray_params.radius_back = lens_params.radius_back;
ray_params.theta_aperture_back = asin(...
    lens_params.lens_radius / lens_params.radius_back...
);
ray_params.d_lens = lens_params.axial_thickness -...
    lens_params.radius_front - lens_params.radius_back;

ior_lens = lens_params.ior_lens;
ior_lens_reference = ior_lens(lens_params.ior_lens_reference_index);
wavelengths = lens_params.wavelengths;
wavelengths_to_rgb = lens_params.wavelengths_to_rgb;

image_sampling = image_params.image_sampling;
normalize_psfs_before_combining = image_params.normalize_psfs_before_combining;
normalize_color_images_globally = image_params.normalize_color_images_globally;

scene_theta_max = scene_params.theta_max;
scene_theta_min = scene_params.theta_min;
if scene_theta_max <= 0
    error('`scene_params.theta_max` must be greater than zero.')
elseif scene_theta_min < 0
    error('`scene_params.theta_min` must be greater than or equal to zero.')
elseif scene_theta_max <= scene_theta_min
    error('`scene_params.theta_max` must be greater than `scene_params.theta_min`.')
end
if length(scene_params.n_lights) == 1
    n_lights_x = scene_params.n_lights;
    n_lights_y = 1;
    radial_lights = true;
elseif length(scene_params.n_lights) == 2
    n_lights_x = scene_params.n_lights(1);
    n_lights_y = scene_params.n_lights(2);
    radial_lights = false;
else
    error('`scene_params.n_lights` must be either a scalar, or a two element vector.');
end
single_source = (n_lights_x == 1 & n_lights_y == 1);

image_output_requested = nargout > 4;
if single_source && image_output_requested
    error('Images (`I` and `I_color`) are not simulated for a single light source.')
end

light_distance_factor_focused = scene_params.light_distance_factor_focused;
light_distance_factor_larger = scene_params.light_distance_factor_larger;
light_distance_factor_smaller = scene_params.light_distance_factor_smaller;
single_depth = (light_distance_factor_larger(2) == 0 &...
    light_distance_factor_smaller(2) == 0);
preserve_angle_over_depths = scene_params.preserve_angle_over_depths;

if ~isempty(varargin)
    verbose = varargin{1};
    plot_light_positions = verbose.plot_light_positions;
    verbose_ray_tracing = verbose.verbose_ray_tracing;
    verbose_ray_interpolation = verbose.verbose_ray_interpolation;
    verbose_psf_analysis = verbose.verbose_psf_analysis;
    display_each_psf = verbose.display_each_psf;
    display_all_psf_each_ior = verbose.display_all_psf_each_ior;
    display_all_psf_each_depth = verbose.display_all_psf_each_depth;
    display_summary = verbose.display_summary;
else
    plot_light_positions = false;
    verbose_ray_tracing = false;
    verbose_ray_interpolation = false;
    verbose_psf_analysis = false;
    display_each_psf = false;
    display_all_psf_each_ior = false;
    display_all_psf_each_depth = false;
    display_summary = false;
end

%% Calculate lens imaging properties

n_ior_lens = length(ior_lens);
[...
    imageFn_reference,...
    f_reference, ~, ...
    U_reference...
] = opticsFromLens(...
    ray_params.ior_environment,...
    ior_lens_reference,...
    ray_params.ior_environment,...
    ray_params.radius_front, ray_params.radius_back,...
    ray_params.d_lens...
);

%% Create light sources and the image plane

% Light source positions on the reference plane
z_light = (light_distance_factor_focused * abs(f_reference)) + U_reference;
principal_point = [0, 0, U_reference];
U_to_light_ref = z_light - U_reference;
if single_source
    x_light = tan(scene_theta_max) * U_to_light_ref;
    y_light = 0;
    scene_theta_min_filter = true;
else
    scene_half_width = tan(scene_theta_max) * U_to_light_ref;
    scene_inner_width = tan(scene_theta_min) * U_to_light_ref;
    if radial_lights
        x_light = linspace(scene_inner_width, scene_half_width, n_lights_x);
        y_light = zeros(1, n_lights_x);
        scene_theta_min_filter = true(n_lights_x, 1);
    else
        x_light = linspace(-scene_half_width, scene_half_width, n_lights_x);
        y_light = linspace(-scene_half_width, scene_half_width, n_lights_y);
        [x_light,y_light] = meshgrid(x_light,y_light);
        scene_theta_min_filter = reshape(...
            (abs(x_light) >= scene_inner_width) |...
            (abs(y_light) >= scene_inner_width),...
            [], 1 ...
        );
    end
end
n_lights = numel(x_light);
X_lights = [
    x_light(:),...
    y_light(:),...
    repmat(z_light, n_lights, 1)
];

% Light source positions at other depths
if single_depth
    depth_ref_index = 1;
    n_depths = 1;
    depth_factors = light_distance_factor_focused;
    X_lights_matrix = X_lights;
else
    depth_ref_inv = 1 / light_distance_factor_focused;
    depth_factors = [
        linspace(1 / light_distance_factor_larger(1), depth_ref_inv, light_distance_factor_larger(2) + 1),...
        linspace(depth_ref_inv, 1 / light_distance_factor_smaller(1), light_distance_factor_smaller(2) + 1)
    ];
    depth_factors = fliplr(1 ./ depth_factors);
    % Remove duplicated reference depth
    depth_ref_index = light_distance_factor_smaller(2) + 1;
    if light_distance_factor_smaller(2) == 0
        depth_factors = depth_factors(2:end);
    else
        depth_factors = [depth_factors(1:depth_ref_index), depth_factors((depth_ref_index + 2):end)];
    end
    depth_factors = depth_factors.';
    n_depths = length(depth_factors);
    
    z_light = reshape((depth_factors * abs(f_reference)) + U_reference, 1, 1, []);
    if preserve_angle_over_depths
        U_to_light_z = z_light - U_reference;
        adjustment_factors = (...
            (U_to_light_z + f_reference) ./ ...
            (U_to_light_ref + f_reference) ...
        ) .* ( ...
            (1 + f_reference ./ U_to_light_ref) ./ ...
            (1 + f_reference ./ U_to_light_z) ...
        );
        adjustment_factors = repmat(adjustment_factors, n_lights, 2, 1);
        X_lights = [
            adjustment_factors .* repmat(X_lights(:, 1:2), 1, 1, n_depths),...
            repmat(z_light, n_lights, 1, 1)
        ];
    else
        X_lights = [repmat(X_lights(:, 1:2), 1, 1, n_depths), repmat(z_light, n_lights, 1, 1)];
    end
    
    X_lights_matrix = reshape(permute(X_lights, [1, 3, 2]), [], 3);
end

n_lights_and_depths = n_lights * n_depths;

% Image plane position
z_film = imageFn_reference([0, 0, z_light(1, 1, depth_ref_index)]);
z_film = z_film(3);
ray_params.d_film = -z_film;

% "Theoretical" image positions, from the thick lens equation
stats_ideal_matrix = analyzePSF([n_lights_and_depths, n_ior_lens]);
for k = 1:n_ior_lens
    [ imageFn, ~, ~, U, U_prime ] = opticsFromLens(...
        ray_params.ior_environment,...
        ior_lens(k),...
        ray_params.ior_environment,...
        ray_params.radius_front, ray_params.radius_back,...
        ray_params.d_lens...
    );
    psfFn = opticsToPSF( imageFn, U, U_prime, lens_params.lens_radius, z_film );

    stats_ideal_matrix(:, k) = psfFn(X_lights_matrix);
    
end

stats_ideal = permute(reshape(stats_ideal_matrix, n_lights, n_depths, n_ior_lens), [1, 3, 2]);

% Find image boundaries
if ~single_source
    X_image_ideal = [stats_ideal.mean_position];
    X_image_ideal = reshape(X_image_ideal, 2, []);
    X_image_ideal = X_image_ideal.';
    X_image_ideal = reshape(X_image_ideal, n_lights, 2, n_ior_lens, n_depths);
    largest_image_abs = max(max(abs(X_image_ideal), [], 4), [], 3);
    % Restore signs
    largest_image = X_image_ideal;
    largest_image(abs(X_image_ideal) ~= repmat(largest_image_abs, 1, 1, n_ior_lens, n_depths)) = 0;
    largest_image = sum(sum(largest_image, 3), 4);
    X_image_grid_x = reshape(largest_image(:, 1), n_lights_y, n_lights_x);
    X_image_grid_y = reshape(largest_image(:, 2), n_lights_y, n_lights_x);
    if n_lights_x > 1
        left_buffer = max(abs(X_image_grid_x(:, 2) - X_image_grid_x(:, 1)));
        right_buffer = max(abs(X_image_grid_x(:, end) - X_image_grid_x(:, end - 1)));
    end
    if n_lights_y > 1
        bottom_buffer = max(abs(X_image_grid_y(2, :) - X_image_grid_y(1, :)));
        top_buffer = max(abs(X_image_grid_y(end, :) - X_image_grid_y(end - 1, :)));
    end
    if n_lights_x == 1
        left_buffer = mean([bottom_buffer, top_buffer]);
        right_buffer = left_buffer;
    end
    if n_lights_y == 1
        bottom_buffer = mean([left_buffer, right_buffer]);
        top_buffer = bottom_buffer;
    end
    image_bounds = [
        min(largest_image(:, 1)) - left_buffer,...
        min(largest_image(:, 2)) - bottom_buffer
    ];
    image_width = max(largest_image(:, 1)) - image_bounds(1) + right_buffer;
    image_height = max(largest_image(:, 2)) - image_bounds(2) + top_buffer;
    image_bounds = [
        image_bounds,...
        image_width image_height...
    ];
end

%% Remove filtered-out light positions
scene_theta_min_filter_rep = repmat(scene_theta_min_filter, n_depths, 1);
stats_ideal = stats_ideal(scene_theta_min_filter, :, :);
stats_ideal_matrix = stats_ideal_matrix(scene_theta_min_filter_rep, :);
X_lights = X_lights(scene_theta_min_filter, :, :);
X_lights_matrix = X_lights_matrix(scene_theta_min_filter_rep, :);
n_lights = size(X_lights, 1);

%% Visualize scene setup
if plot_light_positions
    figure
    scatter3(...
        X_lights_matrix(:, 1),...
        X_lights_matrix(:, 2),...
        X_lights_matrix(:, 3),...
        [], X_lights_matrix(:, 3), 'filled'...
    );
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Light source positions, and lens parameters at the reference wavelength');
    hold on
    
    % First principal plane
    min_x = min(X_lights_matrix(:, 1));
    max_x = max(X_lights_matrix(:, 1));
    min_y = min(X_lights_matrix(:, 2));
    max_y = max(X_lights_matrix(:, 2));

    surf(...
        [min_x, max_x; min_x, max_x],...
        [max_y, max_y; min_y, min_y],...
        [U_reference, U_reference; U_reference, U_reference],...
        'EdgeColor', 'k', 'FaceAlpha', 0.4, 'FaceColor', 'g' ...
    );

    % Focal length
    surf(...
        [min_x, max_x; min_x, max_x],...
        [max_y, max_y; min_y, min_y],...
        [U_reference - f_reference, U_reference - f_reference; U_reference - f_reference, U_reference - f_reference],...
        'EdgeColor', 'k', 'FaceAlpha', 0.4, 'FaceColor', 'y' ...
    );

    % Principal Point
    scatter3(...
        principal_point(1), principal_point(2), principal_point(3),...
        [], 'r', 'filled'...
    );

    % Lens centre
    scatter3(...
        0, 0, lens_params.lens_radius - (lens_params.axial_thickness / 2),...
        [], 'k', 'filled'...
    );
    legend(...
        'Lights', 'First principal plane', 'First focal length',...
        'First principal point', 'Lens centre'...
    );
    set(gca, 'Color', 'none');
    hold off
end

%% Trace rays through the lens and form rays into an image
request_images = ~single_source && (...
    display_each_psf || display_all_psf_each_ior ||...
    display_all_psf_each_depth || image_output_requested ...
);
if request_images
    n_channels = 3;
    I = zeros([image_sampling, n_ior_lens, n_depths]);
    I_color = zeros([image_sampling, n_channels, n_depths]);
else
    I = [];
    I_color = [];
end
stats_real = analyzePSF([n_lights, n_ior_lens, n_depths]);
for j = 1:n_depths
    for k = 1:n_ior_lens
        ray_params.ior_lens = ior_lens(k);
        for i = 1:n_lights
            ray_params.source_position = X_lights(i, :, j);
            [ ...
                image_position, ray_irradiance, ~, incident_position_cartesian ...
                ] = doubleSphericalLens( ray_params, verbose_ray_tracing );
            
            if request_images
                if request_spline_smoothing
                    [ ~, v_adj, image_spline, I_ikj ] = densifyRays(...
                        incident_position_cartesian,...
                        ray_params.radius_front,...
                        image_position,...
                        ray_irradiance,...
                        image_bounds, image_sampling,...
                        verbose_ray_interpolation ...
                    );
                else
                    [ image_values, v_adj, ~, I_ikj ] = densifyRays(...
                        incident_position_cartesian,...
                        ray_params.radius_front,...
                        image_position,...
                        ray_irradiance,...
                        image_bounds, image_sampling,...
                        verbose_ray_interpolation ...
                    );
                end
                
                if normalize_psfs_before_combining
                    I_ikj_scaled = I_ikj ./ max(max(I_ikj));
                else
                    I_ikj_scaled = I_ikj;
                end
                I(:, :, k, j) = I(:, :, k, j) + I_ikj_scaled;
                for c = 1:n_channels
                    I_color(:, :, c, j) = I_color(:, :, c, j) +...
                        (I_ikj_scaled .* wavelengths_to_rgb(k, c));
                end
                
                if display_each_psf
                    figure
                    ax = gca;
                    imagesc(ax,...
                        [image_bounds(1), image_bounds(1) + image_bounds(3)],...
                        [image_bounds(2), image_bounds(2) + image_bounds(4)],...
                        I_ikj...
                        );
                    colormap gray
                    ax.YDir = 'normal';
                    xlabel('X');
                    ylabel('Y');
                    c = colorbar;
                    c.Label.String = 'Irradiance';
                    title(...
                        sprintf('Estimated PSF for a point source at position\n[%g, %g, %g] (%g focal lengths, IOR %g)',...
                        X_lights(i, 1, j), X_lights(i, 2, j), X_lights(i, 3, j),...
                        depth_factors(j), ior_lens(k)...
                        ));
                    axis equal
                end
            else
                if request_spline_smoothing
                    [ ~, v_adj, image_spline ] = densifyRays(...
                        incident_position_cartesian,...
                        ray_params.radius_front,...
                        image_position,...
                        ray_irradiance,...
                        verbose_ray_interpolation ...
                    );
                else
                    [ image_values, v_adj ] = densifyRays(...
                        incident_position_cartesian,...
                        ray_params.radius_front,...
                        image_position,...
                        ray_irradiance,...
                        verbose_ray_interpolation ...
                    );
                end
            end
            
            if single_source && display_each_psf
                figure
                if request_spline_smoothing
                    pts = fnplt(image_spline);
                else
                    plot_resolution = [200 200];
                    x = linspace(...
                        min(image_position(:, 1)),...
                        max(image_position(:, 1)), plot_resolution(2)...
                    );
                    y = linspace(...
                        min(image_position(:, 2)),...
                        max(image_position(:, 2)), plot_resolution(1)...
                    );
                    [pts{1},pts{2}] = meshgrid(x,y);
                    % In the future, this operation could be deduplicated
                    % between this function and `analyzePSF()` (and perhaps
                    % also `densifyRays()`, if applicable).
                    psf_interpolant = scatteredInterpolant(...
                        image_position, image_values...
                    );
                    pts{3} = psf_interpolant(pts{1},pts{2});
                end
                surf(pts{1}, pts{2}, pts{3}, 'EdgeColor', 'none');
                colorbar
                xlabel('X');
                ylabel('Y');
                zlabel('Irradiance');
                colormap summer
                c = colorbar;
                c.Label.String = 'Irradiance';
                title(...
                    sprintf('Estimated PSF for a point source at position\n[%g, %g, %g] (%g focal lengths, IOR %g)',...
                    X_lights(i, 1, j), X_lights(i, 2, j), X_lights(i, 3, j),...
                    depth_factors(j), ior_lens(k)...
                ));
            end
            
            if request_spline_smoothing
                stats_real(i, k, j) = analyzePSF(...
                    image_spline, image_position, v_adj,...
                    verbose_psf_analysis...
                );
            else
                stats_real(i, k, j) = analyzePSF(...
                    image_values, image_position, v_adj,...
                    verbose_psf_analysis...
                );
            end
        end
        
        % Visualize the results, for this wavelength
        if request_images && display_all_psf_each_ior
            figure
            ax = gca;
            imagesc(...
                [image_bounds(1), image_bounds(1) + image_bounds(3)],...
                [image_bounds(2), image_bounds(2) + image_bounds(4)],...
                I(:, :, k, j)...
                );
            colormap gray
            ax.YDir = 'normal';
            xlabel('X');
            ylabel('Y');
            c = colorbar;
            c.Label.String = 'Irradiance';
            hold on
            mean_position_ideal = vertcat(stats_ideal(:, k, j).mean_position);
            mean_position_real = vertcat(stats_real(:, k, j).mean_position);
            scatter(mean_position_ideal(:, 1), mean_position_ideal(:, 2), [], wavelengths_to_rgb(k, :), 'o');
            scatter(mean_position_real(:, 1), mean_position_real(:, 2), [], wavelengths_to_rgb(k, :), '.');
            legend('Thick lens formula', 'Raytracing centroids');
            title(sprintf(...
                'Images of point sources at %g focal lengths, for \\lambda = %g nm',...
                depth_factors(j), wavelengths(k)));
            axis equal
            hold off
        end
    end
end

% Visualize the results, for each depth
if request_images
    if normalize_color_images_globally
        I_color = I_color ./ max(max(max(max(I_color))));
    else
        for j = 1:n_depths
            I_color(:, :, :, j) = I_color(:, :, :, j) ./ max(max(max(I_color(:, :, :, j))));
        end
    end
    
    if display_all_psf_each_depth
        for j = 1:n_depths
            figure
            ax = gca;
            image(...
                [image_bounds(1), image_bounds(1) + image_bounds(3)],...
                [image_bounds(2), image_bounds(2) + image_bounds(4)],...
                I_color(:, :, :, j)...
                );
            ax.YDir = 'normal';
            xlabel('X');
            ylabel('Y');
            hold on
            legend_strings = cell(k * 2, 1);
            for k = 1:n_ior_lens
                mean_position_real = vertcat(stats_real(:, k, j).mean_position);
                mean_position_ideal = vertcat(stats_ideal(:, k, j).mean_position);
                scatter(mean_position_ideal(:, 1), mean_position_ideal(:, 2), [], wavelengths_to_rgb(k, :), 'o');
                scatter(mean_position_real(:, 1), mean_position_real(:, 2), [], wavelengths_to_rgb(k, :), '.');
                legend_strings{k} = sprintf('Thick lens formula, \\lambda = %g nm', wavelengths(k));
                legend_strings{n_ior_lens + k} = sprintf('Raytracing centroids, \\lambda = %g nm', wavelengths(k));
            end
            legend(legend_strings);
            title(sprintf(...
                'Images of point sources at %g focal lengths',...
                depth_factors(j)));
            axis equal
            hold off
        end
    end
end

%% Visualize the results, for all depths
stats_real_matrix = reshape(permute(stats_real, [1, 3, 2]), [], n_ior_lens);

if display_summary
    if single_source
        disp('Point source angle from optical axis:')
        disp(scene_theta_max)
        disp('Distances from first principal plane, in focal lengths:')
        disp(depth_factors)
        disp('Light source position [x, y, z]:')
        disp(X_lights_matrix)
        disp('Image positions:')
        for k = 1:n_ior_lens
            fprintf('Thick lens equation (lambda = %g nm):\n', wavelengths(k))
            disp(vertcat(stats_ideal_matrix(:, k).mean_position))
            fprintf('Raytracing centroids (lambda = %g nm):\n', wavelengths(k))
            disp(vertcat(stats_real_matrix(:, k).mean_position))
        end
    else
        figure
        hold on
        depth_factors_rep = repelem(depth_factors, n_lights);
        legend_strings = cell(n_ior_lens * 2, 1);
        for k = 1:n_ior_lens
            mean_position_ideal = vertcat(stats_ideal_matrix(:, k).mean_position);
            mean_position_real = vertcat(stats_real_matrix(:, k).mean_position);
            scatter3(mean_position_ideal(:, 1), mean_position_ideal(:, 2), depth_factors_rep, [], wavelengths_to_rgb(k, :), 'o');
            scatter3(mean_position_real(:, 1), mean_position_real(:, 2), depth_factors_rep, [], wavelengths_to_rgb(k, :), '.');
            legend_strings{2 * k - 1} = sprintf('Thick lens formula, \\lambda = %g nm', wavelengths(k));
            legend_strings{2 * k} = sprintf('Raytracing centroids, \\lambda = %g nm', wavelengths(k));
        end
        legend(legend_strings);        
        title('Images of point sources seen through a thick lens')
        xlabel('X');
        ylabel('Y');
        zlabel('Light source distance (focal lengths)')
        hold off
    end
end

end

