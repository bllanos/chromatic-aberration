function [...
    X_image_real, X_image_ideal, max_irradiance, X_lights, I, I_color...
] = doubleSphericalLensPSF(...
    lens_params, ray_params, image_params, scene_params, varargin...
)
%% Parse parameters

nargoutchk(1,6);
narginchk(4,5);

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
single_source = (scene_params.n_lights_x == 1 & scene_params.n_lights_x == 1);
if single_source
    single_source_theta = scene_theta_max;
else
    n_lights_x = scene_params.n_lights_x;
    n_lights_y = scene_params.n_lights_y;
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
    display_each_psf = verbose.display_each_psf;
    display_all_psf_each_ior = verbose.display_all_psf_each_ior;
    display_all_psf_each_depth = verbose.display_all_psf_each_depth;
    display_summary = verbose.display_summary;
else
    plot_light_positions = false;
    verbose_ray_tracing = false;
    verbose_ray_interpolation = false;
    display_each_psf = false;
    display_all_psf_each_ior = false;
    display_all_psf_each_depth = false;
    display_summary = false;
end

%% Calculate lens imaging properties

n_ior_lens = length(ior_lens);
[...
    imageFn_reference,...
    f_reference,...
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
lens_center_z = lens_params.lens_radius - (lens_params.axial_thickness / 2);
lens_center = [0, 0, lens_center_z];
lens_center_to_light_ref = z_light - lens_center_z;
if single_source
    x_light = tan(single_source_theta) * lens_center_to_light_ref;
    y_light = 0;
else
    scene_half_width = tan(scene_theta_max) * lens_center_to_light_ref;
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
        lens_center_to_light_relative = (z_light - lens_center_z) ./ lens_center_to_light_ref;
        lens_center_to_light_relative = repmat(lens_center_to_light_relative, n_lights, 3, 1);
        lens_center_rep = repmat(lens_center, n_lights, 1);
        lens_center_to_light =  repmat(X_lights - lens_center_rep, 1, 1, n_depths);
        lens_center_rep = repmat(lens_center_rep, 1, 1, n_depths);
        X_lights = lens_center_rep + (lens_center_to_light .* lens_center_to_light_relative);
    else
        X_lights = [repmat(X_lights(:, 1:2), 1, 1, n_depths), repmat(z_light, n_lights, 1, 1)];
    end
    
    X_lights_matrix = reshape(permute(X_lights, [1, 3, 2]), [], 3);
end

n_lights_and_depths = n_lights * n_depths;

if plot_light_positions
    figure
    scatter3(...
        X_lights_matrix(:, 1), X_lights_matrix(:, 2), X_lights_matrix(:, 3),...
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

    % Lens centre
    scatter3(...
        lens_center(1), lens_center(2), lens_center(3),...
        [], 'k', 'filled'...
    );
    legend('Lights', 'First principal plane', 'First focal length', 'Lens centre');
    set(gca, 'Color', 'none');
    hold off
end

% Image plane position
X_image_ideal_matrix = imageFn_reference([0, 0, z_light(1, 1, depth_ref_index)]);
ideal_image_position = X_image_ideal_matrix(3);
ray_params.d_film = -ideal_image_position;

% "Theoretical" image positions, from the thick lens equation
X_image_ideal_matrix = zeros(n_lights_and_depths, 3, n_ior_lens);
for k = 1:n_ior_lens
    [ imageFn, ~, ~, ~, U_prime ] = opticsFromLens(...
        ray_params.ior_environment,...
        ior_lens(k),...
        ray_params.ior_environment,...
        ray_params.radius_front, ray_params.radius_back,...
        ray_params.d_lens...
    );

    X_image_ideal_matrix(:, :, k) = imageFn(X_lights_matrix);
    % Adjust magnification to approximately account for the actual image plane
    % location
    magnification_correction_ideal = ...
        repmat(U_prime + ray_params.d_film, n_lights_and_depths, 1) ./ ...
        (repmat(U_prime, n_lights_and_depths, 1) - X_image_ideal_matrix(:, 3, k));
    principal_point = repmat([0, 0, U_prime], n_lights_and_depths, 1);
    X_image_rays = X_image_ideal_matrix(:, :, k) - principal_point;
    X_image_ideal_matrix(:, :, k) = principal_point + (X_image_rays .* magnification_correction_ideal);
end

X_image_ideal = permute(reshape(X_image_ideal_matrix, n_lights, n_depths, 3, n_ior_lens), [1, 3, 4, 2]);

% Find image boundaries
if ~single_source
    largest_image_abs = max(max(abs(X_image_ideal), [], 4), [], 3);
    % Restore signs
    largest_image = X_image_ideal;
    largest_image(abs(X_image_ideal) ~= repmat(largest_image_abs, 1, 1, n_ior_lens, n_depths)) = 0;
    largest_image = sum(sum(largest_image, 3), 4);
    X_image_grid_x = reshape(largest_image(:, 1), n_lights_y, n_lights_x);
    X_image_grid_y = reshape(largest_image(:, 2), n_lights_y, n_lights_x);
    left_buffer = max(abs(X_image_grid_x(:, 2) - X_image_grid_x(:, 1)));
    right_buffer = max(abs(X_image_grid_x(:, end) - X_image_grid_x(:, end - 1)));
    bottom_buffer = max(abs(X_image_grid_y(2, :) - X_image_grid_y(1, :)));
    top_buffer = max(abs(X_image_grid_y(end, :) - X_image_grid_y(end - 1, :)));
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

%% Trace rays through the lens and form rays into an image
if single_source
    I = [];
    I_color = [];
else
    n_channels = 3;
    I = zeros([image_sampling, n_ior_lens, n_depths]);
    I_color = zeros([image_sampling, n_channels, n_depths]);
end
X_image_real = zeros(n_lights, 2, n_ior_lens, n_depths);
max_irradiance = zeros(n_lights, 1, n_ior_lens, n_depths);
for j = 1:n_depths
    for k = 1:n_ior_lens
        ray_params.ior_lens = ior_lens(k);
        for i = 1:n_lights
            ray_params.source_position = X_lights(i, :, j);
            [ ...
                image_position, ray_irradiance, ~, incident_position_cartesian ...
                ] = doubleSphericalLens( ray_params, verbose_ray_tracing );
            
            if single_source
                [ X_image_real(i, :, k, j), max_irradiance(i, :, k, j) ] = densifyRays(...
                    incident_position_cartesian,...
                    ray_params.radius_front,...
                    image_position,...
                    ray_irradiance,...
                    verbose_ray_interpolation ...
                    );
            else
                [ X_image_real(i, :, k, j), max_irradiance(i, :, k, j), I_ikj ] = densifyRays(...
                    incident_position_cartesian,...
                    ray_params.radius_front,...
                    image_position,...
                    ray_irradiance,...
                    image_bounds, image_sampling,...
                    verbose_ray_interpolation ...
                    );
                
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
                end
            end
        end
        
        % Visualize the results, for this wavelength
        if ~single_source && display_all_psf_each_ior
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
            scatter(X_image_ideal(:, 1, k, j), X_image_ideal(:, 2, k, j), [], wavelengths_to_rgb(k, :), 'o');
            scatter(X_image_real(:, 1, k, j), X_image_real(:, 2, k, j), [], wavelengths_to_rgb(k, :), '.');
            legend('Thick lens formula', 'Raytracing peaks');
            title(sprintf(...
                'Images of point sources at %g focal lengths, for \\lambda = %g nm',...
                depth_factors(j), wavelengths(k)));
            hold off
        end
    end
end

% Visualize the results, for each depth
if ~single_source && display_all_psf_each_depth
    if normalize_color_images_globally
        I_color = I_color ./ max(max(max(max(I_color))));
    else
        for j = 1:n_depths
            I_color(:, :, :, j) = I_color(:, :, :, j) ./ max(max(max(I_color(:, :, :, j))));
        end
    end
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
            scatter(X_image_ideal(:, 1, k, j), X_image_ideal(:, 2, k, j), [], wavelengths_to_rgb(k, :), 'o');
            scatter(X_image_real(:, 1, k, j), X_image_real(:, 2, k, j), [], wavelengths_to_rgb(k, :), '.');
            legend_strings{k} = sprintf('Thick lens formula, \\lambda = %g nm', wavelengths(k));
            legend_strings{n_ior_lens + k} = sprintf('Raytracing peaks, \\lambda = %g nm', wavelengths(k));
        end
        legend(legend_strings);
        title(sprintf(...
            'Images of point sources at %g focal lengths',...
            depth_factors(j)));
        hold off
    end
end

%% Visualize the results, for all depths
X_image_real_matrix = reshape(permute(X_image_real, [1, 4, 2, 3]), [], 2, n_ior_lens);

if display_summary
    if single_source
        disp('Point source angle from optical axis:')
        disp(single_source_theta)
        disp('Distances from first principal plane, in focal lengths:')
        disp(depth_factors)
        disp('Light source position [x, y, z]:')
        disp(X_lights_matrix)
        disp('Image positions:')
        for k = 1:n_ior_lens
            fprintf('Thick lens equation (lambda = %g nm):\n', wavelengths(k))
            disp(X_image_ideal_matrix(:, 1:2, k))
            fprintf('Raytracing (lambda = %g nm):\n', wavelengths(k))
            disp(X_image_real_matrix(:, :, k))
        end
    else
        figure
        hold on
        depth_factors_rep = repelem(depth_factors, n_lights);
        legend_strings = cell(k * 2, 1);
        for k = 1:n_ior_lens
            scatter3(X_image_ideal_matrix(:, 1, k, :), X_image_ideal_matrix(:, 2, k, :), depth_factors_rep, [], wavelengths_to_rgb(k, :), 'o');
            scatter3(X_image_real_matrix(:, 1, k, :), X_image_real_matrix(:, 2, k, :), depth_factors_rep, [], wavelengths_to_rgb(k, :), '.');
            legend_strings{k} = sprintf('Thick lens formula, \\lambda = %g nm', wavelengths(k));
            legend_strings{n_ior_lens + k} = sprintf('Raytracing peaks, \\lambda = %g nm', wavelengths(k));
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

