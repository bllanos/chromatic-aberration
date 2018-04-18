%% Ray tracing simulation of an image dispersion function
% Obtain a dispersion model corresponding to a thick (biconvex) lens and an
% image sensor.
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
% TODO
%
% ## References
% - Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%   single-shot hyperspectral imaging using a prism." ACM Transactions
%   on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%   doi:10.1145/3130800.3130896

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 26, 2018

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

ray_params.n_incident_rays = 1000;
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

lens_params.wavelengths = linspace(300, 1100, 20);
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
request_spline_smoothing = true;

% ## Scene setup
scene_params.theta_min = deg2rad(0);
scene_params.theta_max = deg2rad(20);
scene_params.n_lights = [11 11];
scene_params.light_distance_factor_focused = 10;
scene_params.light_distance_factor_larger = [4, 0];
scene_params.light_distance_factor_smaller = [1.5, 0];
scene_params.preserve_angle_over_depths = true;

% ## Dispersion model generation
dispersion_fieldname = 'max_position';
max_degree_xy = min(12, min(scene_params.n_lights) - 1);
max_degree_lambda = min(12, length(lens_params.wavelengths) - 1);

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
    n_lambda_plot = min(5, length(lens_params.wavelengths));
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
    stats_real...
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

disparity_raw_real = statsToDisparity(...
    stats_real, ior_lens_reference_index,...
    depth_factors, 0, x_fields, lens_params.wavelengths, lens_params.wavelengths_to_rgb, statsToDisparityVerbose...
);

polyfun_real = xylambdaPolyfit(...
    stats_real, dispersion_fieldname, max_degree_xy, disparity_raw_real, dispersion_fieldname,...
    lens_params.wavelengths, max_degree_lambda, xylambdaPolyfitVerbose...
);

if plot_polynomial_model
    plotXYLambdaPolyfit(...
        stats_real, dispersion_fieldname, disparity_raw_real, dispersion_fieldname,...
        lens_params.wavelengths, lens_params.wavelengths(ior_lens_reference_index), n_lambda_plot, polyfun_real...
    );
end
