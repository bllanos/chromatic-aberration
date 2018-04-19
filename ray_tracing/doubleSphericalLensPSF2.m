function [...
    stats_real_rgb, stats_ideal_rgb, I_color, I...
] = doubleSphericalLensPSF2(...
    lens_params, ray_params, image_params, X_lights, z_film,...
    lights_filter, varargin...
)
% DOUBLESPHERICALLENSPSF2  Simulate pixelated images of point light sources by raytracing
%
% ## Syntax
% [...
%     stats_real, stats_ideal, I_color, I...
% ] = doubleSphericalLensPSF2(...
%     lens_params, ray_params, image_params, X_lights, z_film,...
%     lights_filter [, depth_factors, verbose]...
% )
%
% ## Description
% [...
%     stats_real, stats_ideal, I_color, I...
% ] = doubleSphericalLensPSF2(...
%     lens_params, ray_params, image_params, X_lights, z_film,...
%     lights_filter [, depth_factors, verbose]...
% )
%   Simulate the image of grids of point light sources by raytracing. (One
%   to four output arguments can be requested.)
%
%   'doubleSphericalLensPSF()' specializes in accurate PSF statistics,
%   because its results are not dependent on an image resolution.
%   'doubleSphericalLensPSF2()' specializes in efficient image simulation,
%   and produces PSF statistics which are affected by image resolution.
%   Therefore, 'doubleSphericalLensPSF()' is useful for analyzing the
%   properties of the lens, whereas 'doubleSphericalLensPSF2()' can better
%   simulate the lens and sensor combination.
%
%   Also, the PSF analysis performed by 'doubleSphericalLensPSF()' is
%   per-wavelength, whereas the PSF analysis performed by
%   'doubleSphericalLensPSF2()' is per-colour channel.
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
%   - wavelengths: The wavelengths of light corresponding to the elements
%     of `ior_lens`. A row vector of length 'k'. This parameter is used for
%     figure legends only, not for calculations.
%   - wavelengths_to_rgb: RGB quantum efficiencies for the wavelengths
%     corresponding to the indices of refraction in `ior_lens`. The i-th
%     row of this k x 3 matrix represents the RGB sensitivities
%     corresponding to the i-th wavelength. `wavelengths_to_rgb` allows
%     colour images to be produced by adding together the contributions of
%     each wavelength to the red, green, and blue colour channels.
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
%   - image_bounds: A 4-element vector, in the form of the output argument
%     of 'imageBoundaries()'. (Refer to the documentation of
%     'imageBoundaries.m' for details.) Alternatively, this information can
%     be automatically estimated if there is more than one light source in
%     the scene, in which case the field can be empty (`[]`).
%   - normalize_color_images_globally: Colour images are normalized so that
%     colour channel values are no greater than one. Normalization is
%     either performed by dividing by the maximum channel value over all
%     scene depths (true), or the maximum channel value of each scene depth
%     (false). Therefore, if this parameter is `true`, colour images
%     produced for different depths are directly comparable, but images for
%     out-of-focus depths may be quite dim.
%   - intensity_threshold: When combining point spread functions computed
%     from the thick lens equation to obtain per-colour channel statistics,
%     use this scalar as the `threshold` input argument of
%     'diskStatsToRGB()'. Refer to the documentation of 'diskStatsToRGB.m'
%     for details.
%
% X_lights -- Light positions
%   The positions, expressed in cartesian coordinates (x, y, z), of the
%   lights in the scene. `X_lights(i, :, j)` is the 3D position of the i-th
%   light in the grid of lights placed at the j-th depth.
%
% z_film -- Image plane location
%   The z-coordinate of the image plane.
%
% lights_filter -- Light gaps filter
%   For efficiency, a full grid of lights producing point spread functions
%   over the image plane may not be desirable. `lights_filter` is a logical
%   vector of length `size(X_lights, 1)` indicating which light positions
%   are to be ignored. For instance, lights close to the optical axis may
%   be ignored, as their images show little chromatic aberration.
%
% depth_factors -- Depths expressed in focal lengths
%   The depths of the light sources, measured in multiples of the first
%   focal length of the lens from the first principal plane, for instance.
%   `depth_factors(j)` corresponds to the z-values in `X_lights(:, :, j)`.
%
%   `depth_factors` is a vector of values to substitute for light source
%   depths, for display purposes, and is not needed when
%   debugging/visualization output is disabled.
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
%   each light source, for each colour channel, and for each depth, output
%   as a structure array. `stats_real(i, k, j)` is the `stats` output
%   argument of 'analyzePSFImage.m' for the i-th light and the k-th colour
%   channel, when the light is positioned at the j-th depth.
%
% stats_ideal -- Theoretical image statistics
%   Point spread function statistics for the images produced by each light
%   source, as predicted by the thick lens equation. `stats_ideal` has
%   the same format as `stats_real`. It is produced by summing over
%   statistics for individual wavelengths, weighted by RGB quantum
%   efficiencies.
%
% I_color -- Simulated colour images
%   A 4D array containing the images formed by combining all point spread
%   functions for all lights at each depth. `I_color(:, :, :, j)` is the
%   colour image for the produced by the grid of lights placed at the j-th
%   depth. The third dimension of `I_color` indexes Red, Green, and Blue
%   colour channels. `I_color` is produced by combining the images for
%   individual wavelengths in `I` according to the RGB quantum efficiencies
%   for the different wavelengths in `lens_params.wavelengths_to_rgb`.
%
%   This output argument is only available if there are multiple lights at
%   each depth.
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
% ## Notes
% - Automatic estimation of the image boundaries, for `I_color` and `I`,
%   will yield poor results if all of the lights are on a single line.
%
% ### Coordinate System
% - The radii of both faces of the lens are positive when the lens is
%   biconvex.
% - The positive z-axis points towards the front of the lens, along
%   the optical axis, assuming the front of the lens is convex.
%
% ### Polychromatic light sources
%
% I simulate point spread functions for RGB colour channels by adding
% together the point spread functions for individual wavelengths, weighted
% by the quantum efficiencies of the camera sensor for the individual
% wavelengths. Such an approach relies on the following assumptions:
% - The spectral power distribution of the light source is uniform.
% - The camera sensor does not perform any on-chip processing to adjust the
%   relative responses of pixels based on assumptions of the light colour.
%
% ### Efficiency
% - The `I_color` and `I` output arguments can consume a lot of memory
%   (especially `I`). They are generated if explicitly requested as output.
%
% See also opticsFromLens, doubleSphericalLens, densifyRaysImage, analyzePSFImage, imagingScenario, doubleSphericalLensPSF, diskStatsToRGB, imageBoundaries

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 22, 2017

%% Parse parameters

nargoutchk(1,4);
narginchk(6,8);

ray_params = lensParamsToRayParams(ray_params, lens_params, z_film);

ior_lens = lens_params.ior_lens;
n_ior_lens = length(ior_lens);

wavelengths_to_rgb = lens_params.wavelengths_to_rgb;

image_sampling = image_params.image_sampling;
normalize_color_images_globally = image_params.normalize_color_images_globally;
stats_ideal_intensity_threshold = image_params.intensity_threshold;

n_lights = size(X_lights, 1);
single_source = (n_lights == 1);

I_color_output_requested = nargout > 2;
I_output_requested = nargout > 3;

if ~isempty(varargin)
    if length(varargin) == 2
        depth_factors = varargin{1};
        verbose = varargin{2};
    else
        error('Both `depth_factors`, and `verbose` must be passed together.')
    end
    verbose_ray_tracing = verbose.verbose_ray_tracing;
    verbose_ray_interpolation = verbose.verbose_ray_interpolation;
    verbose_psf_analysis = verbose.verbose_psf_analysis;
    display_each_psf = verbose.display_each_psf;
    display_each_psf_rgb = verbose.display_each_psf_rgb;
    display_all_psf_each_depth = verbose.display_all_psf_each_depth;
    display_summary = verbose.display_summary;
else
    verbose_ray_tracing = false;
    verbose_ray_interpolation = false;
    verbose_psf_analysis = false;
    display_each_psf = false;
    display_each_psf_rgb = false;
    display_all_psf_each_depth = false;
    display_summary = false;
end


%% Imaging setting configuration, and thick lens equation results

n_depths = size(X_lights, 3);
if n_depths == 1
    X_lights_matrix = X_lights;
else
    X_lights_matrix = reshape(permute(X_lights, [1, 3, 2]), [], 3);
end

n_lights_and_depths = n_lights * n_depths;

% "Theoretical" image positions, from the thick lens equation
stats_ideal_matrix_wavelengths = preallocateStats([n_lights_and_depths, n_ior_lens]);
for k = 1:n_ior_lens
    [ imageFn, ~, ~, U, U_prime ] = opticsFromLens(...
        ray_params.ior_environment,...
        ior_lens(k),...
        ray_params.ior_environment,...
        ray_params.radius_front, ray_params.radius_back,...
        ray_params.d_lens...
    );
    psfFn = opticsToPSF( imageFn, U, U_prime, lens_params.lens_radius, z_film );
    stats_ideal_matrix_wavelengths(:, k) = psfFn(X_lights_matrix);
end

% Combine the results into colour channels
n_channels = size(wavelengths_to_rgb, 2);
stats_ideal_wavelengths = permute(reshape(stats_ideal_matrix_wavelengths, n_lights, n_depths, n_ior_lens), [1, 3, 2]);
stats_ideal_rgb = diskStatsToRGB(stats_ideal_wavelengths, wavelengths_to_rgb, stats_ideal_intensity_threshold);
stats_ideal_rgb_matrix = reshape(permute(stats_ideal_rgb, [1, 3, 2]), [], n_channels);

% Find image boundaries
if isempty(image_params.image_bounds) && single_source
    error('Image boundaries cannot be automatically estimated for a single light source.')
end
image_bounds = imageBoundaries( image_params.image_bounds, stats_ideal_rgb );

%% Remove filtered-out light positions

lights_filter_rep = repmat(lights_filter, n_depths, 1);
stats_ideal_rgb = stats_ideal_rgb(lights_filter, :, :);
stats_ideal_rgb_matrix = stats_ideal_rgb_matrix(lights_filter_rep, :);
X_lights = X_lights(lights_filter, :, :);
X_lights_matrix = X_lights_matrix(lights_filter_rep, :);
n_lights = size(X_lights, 1);

%% Trace rays through the lens and form rays into an image
if I_output_requested
    I = zeros([image_sampling, n_ior_lens, n_depths]);
else
    I = [];
end
if I_color_output_requested || display_all_psf_each_depth
    I_color = zeros([image_sampling, n_channels, n_depths]);
else
    I_color = [];
end
stats_real_rgb = preallocateStats([n_lights, n_channels, n_depths]);
color_fmts = {'r', 'g', 'b'};
color_names = {'Red', 'Green', 'Blue'};
for j = 1:n_depths
    for i = 1:n_lights
        I_color_ji = zeros([image_sampling, n_channels]);
        I_color_ji_mask = false(image_sampling);
        ray_params.source_position = X_lights(i, :, j);
        for k = 1:n_ior_lens
            ray_params.ior_lens = ior_lens(k);
            [ ...
                image_position, ray_irradiance ...
            ] = doubleSphericalLens( ray_params, verbose_ray_tracing );
            
            [I_jik, mask ] = densifyRaysImage(...
                image_position, ray_irradiance,...
                image_bounds, image_sampling,...
                verbose_ray_interpolation...
            );
            I_color_ji_mask = I_color_ji_mask | mask;

            if I_output_requested
                I(:, :, k, j) = I(:, :, k, j) + I_jik;
            end
            for c = 1:n_channels
                I_color_ji(:, :, c) = I_color_ji(:, :, c) +...
                    (I_jik .* wavelengths_to_rgb(k, c));
            end

            if display_each_psf
                figure
                ax = gca;
                imagesc(ax,...
                    [image_bounds(1), image_bounds(1) + image_bounds(3)],...
                    [image_bounds(2) + image_bounds(4), image_bounds(2)],...
                    I_jik...
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
        end
        
        for c = 1:n_channels
            stats_real_rgb(i, c, j) = analyzePSFImage(...
                I_color_ji(:, :, c), image_bounds, I_color_ji_mask, verbose_psf_analysis...
            );
        end
        
        if I_color_output_requested || display_all_psf_each_depth
            I_color(:, :, :, j) = I_color(:, :, :, j) + I_color_ji;
        end

        % Visualize the results, for this light source
        if display_each_psf_rgb
            figure
            ax = gca;
            image(...
                [image_bounds(1), image_bounds(1) + image_bounds(3)],...
                [image_bounds(2) + image_bounds(4), image_bounds(2)],...
                I_color_ji / max(max(max(I_color_ji)))...
                );
            ax.YDir = 'normal';
            xlabel('X');
            ylabel('Y');
            hold on
            legend_strings = cell(n_channels * 2, 1);
            for c = 1:n_channels
                mean_position_real = stats_real_rgb(i, c, j).mean_position;
                mean_position_ideal = stats_ideal_rgb(i, c, j).mean_position;
                scatter(mean_position_ideal(1), mean_position_ideal(2), [], color_fmts{c}, 'o');
                scatter(mean_position_real(1), mean_position_real(2), [], color_fmts{c}, '.');
                legend_strings{c} = sprintf('Thick lens formula, %s channel', color_names{c});
                legend_strings{n_channels + c} = sprintf('Raytracing centroids, %s channel', color_names{c});
            end
            legend(legend_strings);
            title(sprintf(...
                'Image of a point source at %g focal lengths',...
                depth_factors(j)));
            axis equal
            hold off
        end
    end
end

% Visualize the results, for each depth
if I_color_output_requested || display_all_psf_each_depth
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
                [image_bounds(2) + image_bounds(4), image_bounds(2)],...
                I_color(:, :, :, j)...
                );
            ax.YDir = 'normal';
            xlabel('X');
            ylabel('Y');
            
            hold on
            legend_strings = cell(n_channels * 2, 1);
            for c = 1:n_channels
                mean_position_real = vertcat(stats_real_rgb(:, c, j).mean_position);
                mean_position_ideal = vertcat(stats_ideal_rgb(:, c, j).mean_position);
                scatter(mean_position_ideal(:, 1), mean_position_ideal(:, 2), [], color_fmts{c}, 'o');
                scatter(mean_position_real(:, 1), mean_position_real(:, 2), [], color_fmts{c}, '.');
                legend_strings{2 * c - 1} = sprintf('Thick lens formula, %s channel', color_names{c});
                legend_strings{2 * c} = sprintf('Raytracing centroids, %s channel', color_names{c});
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
stats_real_rgb_matrix = reshape(permute(stats_real_rgb, [1, 3, 2]), [], n_channels);

if display_summary
    if single_source
        disp('Distances from first principal plane, in focal lengths:')
        disp(depth_factors)
        disp('Light source position [x, y, z]:')
        disp(X_lights_matrix)
        disp('Image positions:')
        for c = 1:n_channels
            fprintf('Thick lens formula, %s channel:\n', color_names{c})
            disp(vertcat(stats_ideal_rgb_matrix(:, c).mean_position))
            fprintf('Raytracing centroids, %s channel:\n', color_names{c})
            disp(vertcat(stats_real_rgb_matrix(:, c).mean_position))
        end
    else
        figure
        hold on
        depth_factors_rep = repelem(depth_factors, n_lights);
        legend_strings = cell(n_channels * 2, 1);
        for c = 1:n_channels
            mean_position_ideal = vertcat(stats_ideal_rgb_matrix(:, c).mean_position);
            mean_position_real = vertcat(stats_real_rgb_matrix(:, c).mean_position);
            scatter3(mean_position_ideal(:, 1), mean_position_ideal(:, 2), depth_factors_rep, [], color_fmts{c}, 'o');
            scatter3(mean_position_real(:, 1), mean_position_real(:, 2), depth_factors_rep, [], color_fmts{c}, '.');
            legend_strings{2 * c - 1} = sprintf('Thick lens formula, %s channel:\n', color_names{c});
            legend_strings{2 * c} = sprintf('Raytracing centroids, %s channel', color_names{c});
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
