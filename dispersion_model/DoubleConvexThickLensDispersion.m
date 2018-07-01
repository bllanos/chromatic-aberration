%% Ray tracing simulation of an image dispersion function
% Obtain a dispersion model corresponding to a thick (biconvex) lens
% projecting an image onto a sensor. Fit the dispersion model to point
% spread function statistics that vary with spatial coordinates and
% wavelength.
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
% ### Graphical output from 'plotXYLambdaPolyfit()'
% - Displayed if `plot_polynomial_model` is `true`.
%
% ### Polynomial fitting results
%
% A '.mat' file containing the following variables:
%
% - 'centers': The independent variables data used for fitting the
%   polynomial model of dispersion. `centers` is a structure array, with
%   one field containing the image positions of the centres of simulated
%   point spread functions. `centers(i, k)` is the centre of the point
%   spread function for the i-th light source position, with the light
%   source emitting light of the k-th wavelength. In this context, "centre"
%   refers to the point spread function statistic named in the
%   `dispersion_fieldname` variable.
% - 'disparity': The dependent variables data used for fitting the
%   polynomial model of dispersion. `disparity` is the first output
%   argument of 'statsToDisparity()', produced when 'statsToDisparity()'
%   was called with `centers` as one of its input arguments. The format of
%   `disparity` is described in the documentation of 'statsToDisparity()'.
%   `disparity` contains the dispersion vectors between the centres of
%   point spread functions for different wavelengths.
% - 'polyfun_data': The polynomial model of dispersion, modeling the
%   mapping from `centers` to `disparity`. `polyfun_data` can be
%   converted to a function form using `polyfun =
%   makePolyfun(polyfun_data)`
% - 'bands': A vector containing the wavelengths at which dispersion can be
%   evaluated to approximate the dispersion between the Red, Green, and
%   Blue colour channels of the sensor. These are the wavelengths at which
%   the colour channels reach their peak quantum efficiencies.
% - 'model_space': A structure describing the range of coordinates over
%   which the polynomial model of dispersion is valid, having the following
%   fields:
%   - 'corners': The first and second rows contain the (x,y) coordinates of
%     the top left and bottom right corners of the region, respectively.
%   - 'pixel_size': The side length of a pixel in the coordinate frame used
%     by the 'corners' field.
%   - 'system': A character vector, 'geometric', indicating that the
%     dispersion model was constructed under geometrical optics coordinate
%     conventions, wherein the y-axis is positive upwards on the image
%     plane, and the origin is the image centre.
%
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## References
% - Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%   single-shot hyperspectral imaging using a prism." ACM Transactions
%   on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%   doi:10.1145/3130800.3130896
% - V. Rudakova and P. Monasse. "Precise Correction of Lateral Chromatic
%   Aberration in Images," Lecture Notes on Computer Science, 8333, pp.
%   12–22, 2014.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 26, 2018

% List of parameters to save with results
parameters_list = {
        'lens_params',...
        'ior_lens_reference_index',...
        'ray_params',...
        'sellmeierConstants',...
        'image_params',...
        'pixel_size',...
        'request_spline_smoothing',...
        'scene_params',...
        'dispersion_fieldname',...
        'max_degree_xy',...
        'max_degree_lambda',...
        'model_from_reference'...
    };

%% Input data and parameters
% Refer to the documentation of `doubleSphericalLensPSF` for details

% ## Raytracing parameters

% ### Lens parameters
% Based on
% '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180226_SmallFLLenses_EdmundOptics/3mmDiameter4dot5mmFLUncoatedDoubleConvexLens_prnt_32022.pdf'
lens_params.lens_radius = 3 / 2;
lens_params.axial_thickness = 2;
lens_params.radius_front = 4.29;
lens_params.radius_back = lens_params.radius_front;

ray_params.n_incident_rays = 125000;
ray_params.sample_random = true;
ray_params.ior_environment = 1.0;

% #### Index of refraction
% The focal length specification wavelength is 587.6 nm
% The lens is made of SCHOTT N-BK7 glass

% Constants for SCHOTT N-BK7 glass retrieved from the SCHOTT glass
% datasheet provided at
% https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
sellmeierConstants.B_1 = 1.03961212;
sellmeierConstants.B_2 = 0.231792344;
sellmeierConstants.B_3 = 1.01046945;
sellmeierConstants.C_1 = 0.00600069867;
sellmeierConstants.C_2 = 0.0200179144;
sellmeierConstants.C_3 = 103.560653;

lens_params.wavelengths = linspace(300, 1100, 100);
lens_params.ior_lens = sellmeierDispersion(lens_params.wavelengths, sellmeierConstants);

% Index of the wavelength/index of refraction to be used to position the
% image plane
[~, ior_lens_reference_index] = min(abs(lens_params.wavelengths - 587.6));

% Used for debugging visualizations only
%
% Obtained using the quantum efficiencies presented in
% '/home/llanos/GoogleDrive/ThesisResearch/Equipment/FLEA3/20170508_FL3_GE_EMVA_Imaging Performance Specification.pdf'
% Image sensor: Sony ICX655, 2/3", Color (page 19)
lens_params.wavelengths_to_rgb = sonyQuantumEfficiency(lens_params.wavelengths);

% Normalize, for improved colour saturation
lens_params.wavelengths_to_rgb = lens_params.wavelengths_to_rgb ./...
    max(max(lens_params.wavelengths_to_rgb));

% ### Ray interpolation parameters
image_params.image_sampling = [2048, 2448];
% Pixel size for the Sony ICX655 sensor is 3.45 micrometres
% (https://www.edmundoptics.com/resources/application-notes/imaging/pixel-sizes-and-optics/)
pixel_size = 3.45e-3; % millimetres
image_min_x = -image_params.image_sampling(2) * pixel_size / 2;
image_min_y = -image_params.image_sampling(1) * pixel_size / 2;
image_width = image_params.image_sampling(2) * pixel_size;
image_height = image_params.image_sampling(1) * pixel_size;
image_params.image_bounds = [image_min_x, image_min_y, image_width, image_height];
image_params.normalize_color_images_globally = false;
image_params.normalize_psfs_before_combining = false;

% ### Use spline interpolation to reduce noise
request_spline_smoothing = false;

% ## Scene setup
scene_params.theta_min = deg2rad(0);
scene_params.theta_max = deg2rad(20);
scene_params.n_lights = [12 12];
scene_params.light_distance_factor_focused = 10;
scene_params.light_distance_factor_larger = [4, 0];
scene_params.light_distance_factor_smaller = [1.5, 0];
scene_params.preserve_angle_over_depths = true;

% ## Dispersion model generation
dispersion_fieldname = 'mean_position';
max_degree_xy = min(12, min(scene_params.n_lights) - 1);
max_degree_lambda = min(12, length(lens_params.wavelengths) - 1);

% If `true`, model dispersion between bands (colour channels or spectral
% bands) as a function of positions in the reference band. If `false`,
% model dispersion as a function of positions in the non-reference bands.
% The first case is useful for warping the other bands to align with the
% reference band, such as when correcting chromatic aberration by image
% warping. The second case is useful for warping an "ideal" image to
% compare it with an observed aberrated image. In both cases, the
% dispersion vectors point from the reference band to the other bands.
model_from_reference = true;

% ## Debugging Flags
plot_light_positions = true;

doubleSphericalLensPSFVerbose.verbose_ray_tracing = false;
doubleSphericalLensPSFVerbose.verbose_ray_interpolation = false;
doubleSphericalLensPSFVerbose.verbose_psf_analysis = false;
doubleSphericalLensPSFVerbose.display_each_psf = false;
doubleSphericalLensPSFVerbose.display_all_psf_each_ior = false;
doubleSphericalLensPSFVerbose.display_all_psf_each_depth = false;
doubleSphericalLensPSFVerbose.display_summary = true;

statsToDisparityVerbose.display_raw_values = true;
statsToDisparityVerbose.display_raw_disparity = true;
statsToDisparityVerbose.filter = struct(...
    'mean_position', true,...
    'mean_value', false,...
    'max_position', true,...
    'max_value', false,...
    'radius', false...
);

xylambdaPolyfitVerbose = true;
plot_polynomial_model = true;
if plot_polynomial_model
    n_lambda_plot = min(20, length(lens_params.wavelengths));
end

%% Create light sources

lens_params_scene = lens_params;
lens_params_scene.ior_lens = lens_params.ior_lens(ior_lens_reference_index);
[...
    X_lights, z_film, lights_filter, depth_factors...
] = imagingScenario(...
    lens_params_scene, ray_params.ior_environment, scene_params, plot_light_positions...
);

%% Run the simulation

[...
    centers...
] = doubleSphericalLensPSF(...
    lens_params, ray_params, image_params, X_lights, z_film, lights_filter,...
    request_spline_smoothing, depth_factors, doubleSphericalLensPSFVerbose...
);

%% Fit a dispersion model to the results

x_fields = struct(...
    'mean_position', 'mean_position',...
    'mean_value', 'mean_position',...
    'max_position', 'max_position',...
    'max_value', 'max_position',...
    'radius', 'mean_position'...
);

disparity = statsToDisparity(...
    centers, ior_lens_reference_index,...
    depth_factors, 0, x_fields, lens_params.wavelengths, lens_params.wavelengths_to_rgb, statsToDisparityVerbose...
);

n_wavelengths = length(lens_params.wavelengths);

if model_from_reference
    centers_for_fitting = repmat(centers(:, ior_lens_reference_index), 1, n_wavelengths);
else
    centers_for_fitting = centers;
end

[ polyfun, polyfun_data ] = xylambdaPolyfit(...
    centers_for_fitting, dispersion_fieldname, max_degree_xy, disparity, dispersion_fieldname,...
    lens_params.wavelengths, max_degree_lambda, xylambdaPolyfitVerbose...
);

%% Visualization

if plot_polynomial_model
    plotXYLambdaPolyfit(...
        centers_for_fitting, dispersion_fieldname, disparity, dispersion_fieldname,...
        lens_params.wavelengths, lens_params.wavelengths(ior_lens_reference_index), n_lambda_plot, polyfun...
    );
end

%% Save results to a file

% Find appropriate wavelengths for approximating RGB chromatic aberration
[~, ind] = max(lens_params.wavelengths_to_rgb, [], 1);
bands = lens_params.wavelengths(ind);

% Indicate where in the image the model is usable
centers_unpacked = permute(reshape([centers.(dispersion_fieldname)], 2, []), [2 1]);
model_space.corners = [
    min(centers_unpacked(:, 1)), max(centers_unpacked(:, 2));
    max(centers_unpacked(:, 1)), min(centers_unpacked(:, 2))
    ];
model_space.pixel_size = pixel_size;
model_space.system = 'geometric';

save_variables_list = [ parameters_list, {...
        'centers',...
        'disparity',...
        'polyfun_data',...
        'bands',...
        'model_space'...
    } ];
uisave(save_variables_list,'DoubleConvexThickLensDispersionResults');