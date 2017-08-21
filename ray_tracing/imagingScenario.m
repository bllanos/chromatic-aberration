function [...
    X_lights, z_film, lights_filter, depth_factors...
] = imagingScenario(...
    lens_params, ior_environment, scene_params, varargin...
)
% IMAGINGSCENARIO  Generate a set of point light source positions aligned over depths
%
% ## Syntax
% [...
%     X_lights, z_film, lights_filter, depth_factors...
% ] = imagingScenario(...
%     lens_params, ior_environment, scene_params [, verbose]...
% )
%
% ## Description
% [...
%     X_lights, z_film, lights_filter, depth_factors...
% ] = imagingScenario(...
%     lens_params, ior_environment, scene_params [, verbose]...
% )
%   Create a grid of point light source positions, to be used for
%   raytracing. (One to four output arguments can be requested.)
%
% ## Input Arguments
%
% lens_params -- Lens parameters structure
%   A description of a lens formed from two spherical surfaces.
%   Passed as a structure with the following fields:
%   - lens_radius: The radius of the lens (i.e. half the height of the lens
%     when viewed edge-on). This field is used only when `verbose` is
%     `true`.
%   - axial_thickness: The thickness of the lens along its optical axis.
%   - radius_front: The radius of curvature of the front surface of the
%     lens
%   - radius_back: The radius of curvature of the back surface of the
%     lens
%   - ior_lens: The refractive index of the lens. A scalar.
%
% ior_environment -- Environment refractive index
%   The refractive index of the medium surrounding the lens (on both sides)
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
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% X_lights -- Light positions
%   The positions, expressed in cartesian coordinates (x, y, z), of the
%   lights in the scene. `X_lights(i, :, j)` is the 3D position of the i-th
%   light in the grid of lights placed at the j-th depth.
%
% z_film -- Image plane location
%   The z-coordinate of the image plane.
%
% lights_filter -- Light angle filter
%   `lights_filter` is a logical vector of length `size(X_lights, 1)`
%   indicating which light positions are to be ignored, because they are
%   less than `scene_params.theta_min` from the optical axis. Outputting
%   both `X_lights`, and `lights_filter` allows the caller to both use only
%   lights satisfying `scene_params.theta_min`, and use the full grid
%   pattern of the lights (without holes), as desired.
%
% depth_factors -- Depths expressed in focal lengths
%   The depths of the light sources, measured in multiples of the first
%   focal length of the lens from the first principal plane.
%   `depth_factors(j)` corresponds to the z-values in `X_lights(:, :, j)`.
%
% ## Notes
%
% ### Coordinate System
% - The radii of both faces of the lens are positive when the lens is
%   biconvex.
% - The positive z-axis points towards the front of the lens, along
%   the optical axis, assuming the front of the lens is convex.
%
% See also opticsFromLens, doubleSphericalLens, doubleSphericalLensPSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 18, 2017

%% Parse parameters

nargoutchk(1,4);
narginchk(3,4);

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

light_distance_factor_focused = scene_params.light_distance_factor_focused;
light_distance_factor_larger = scene_params.light_distance_factor_larger;
light_distance_factor_smaller = scene_params.light_distance_factor_smaller;
single_depth = (light_distance_factor_larger(2) == 0 &...
    light_distance_factor_smaller(2) == 0);
preserve_angle_over_depths = scene_params.preserve_angle_over_depths;

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

%% Calculate lens imaging properties

ior_lens = lens_params.ior_lens;
d_lens = lens_params.axial_thickness -...
    lens_params.radius_front - lens_params.radius_back;

[...
    imageFn_reference,...
    f_reference, ~, ...
    U_reference...
] = opticsFromLens(...
    ior_environment,...
    ior_lens,...
    ior_environment,...
    lens_params.radius_front, lens_params.radius_back,...
    d_lens...
);

%% Create light sources and the image plane

% Light source positions on the reference plane
z_light = (light_distance_factor_focused * abs(f_reference)) + U_reference;
principal_point = [0, 0, U_reference];
U_to_light_ref = z_light - U_reference;
if single_source
    x_light = tan(scene_theta_max) * U_to_light_ref;
    y_light = 0;
    lights_filter = true;
else
    scene_half_width = tan(scene_theta_max) * U_to_light_ref;
    scene_inner_width = tan(scene_theta_min) * U_to_light_ref;
    if radial_lights
        x_light = linspace(scene_inner_width, scene_half_width, n_lights_x);
        y_light = zeros(1, n_lights_x);
        lights_filter = true(n_lights_x, 1);
    else
        x_light = linspace(-scene_half_width, scene_half_width, n_lights_x);
        y_light = linspace(-scene_half_width, scene_half_width, n_lights_y);
        [x_light,y_light] = meshgrid(x_light,y_light);
        lights_filter = reshape(...
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
    
    z_light = reshape((depth_factors * abs(f_reference)) + U_reference, 1, 1, []);
    n_depths = length(depth_factors);
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

% Image plane position
z_film = imageFn_reference([0, 0, z_light(1, 1, depth_ref_index)]);
z_film = z_film(3);

%% Visualize scene setup
if verbose
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

end

