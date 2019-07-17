%% Ray tracing simulation of chromatic aberration
% Simulate the chromatic point spread function of a thick (biconvex) lens.
% Use an image resolution-agnostic simulation method, for highly
% accurate simulations of PSFs for individual wavelengths of light.
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
% Presently, the output arguments of 'radialChromaticAberration.m'.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

%% Input data and parameters
% Refer to the documentation of `doubleSphericalLensPSF` for details

% ## Raytracing parameters

% ### Lens parameters
% Based on
% '20180226_SmallFLLenses_EdmundOptics/3mmDiameter4dot5mmFLUncoatedDoubleConvexLens_prnt_32022.pdf'
lens_params.lens_radius = 3 / 2;
lens_params.axial_thickness = 2;
lens_params.radius_front = 4.29;
lens_params.radius_back = lens_params.radius_front;

ray_params.n_incident_rays = 500000;
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

% Wavelength values corresponding to peak sensitivity in each colour
% channel of the sensor
lens_params.wavelengths = [457, 539, 606];
lens_params.ior_lens = sellmeierDispersion(lens_params.wavelengths, sellmeierConstants);

% Index of the wavelength/index of refraction to be used to position the
% image plane
lens_params.ior_lens_reference_index = 2; % Use the green channel

% Obtained using the quantum efficiencies presented in
% '20170508_FL3_GE_EMVA_Imaging Performance Specification.pdf'
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
scene_params.n_lights = [5 5];
scene_params.light_distance_factor_focused = 2;
scene_params.light_distance_factor_larger = [4, 10];
scene_params.light_distance_factor_smaller = [1.5, 9];
scene_params.preserve_angle_over_depths = true;

% ## Debugging Flags
plot_light_positions = true;

doubleSphericalLensPSFVerbose.verbose_ray_tracing = false;
doubleSphericalLensPSFVerbose.verbose_ray_interpolation = false;
doubleSphericalLensPSFVerbose.verbose_psf_analysis = false;
doubleSphericalLensPSFVerbose.display_each_psf = false;
doubleSphericalLensPSFVerbose.display_all_psf_each_ior = false;
doubleSphericalLensPSFVerbose.display_all_psf_each_depth = false;
doubleSphericalLensPSFVerbose.display_summary = true;

verbose_aberration_ideal = true;
verbose_aberration_real = true;
radialChromaticAberrationVerbose.display_raw_values = true;
radialChromaticAberrationVerbose.display_raw_disparity = true;
radialChromaticAberrationVerbose.display_stats_splines = true;
radialChromaticAberrationVerbose.display_spline_differences = true;
radialChromaticAberrationVerbose.filter = struct(...
    'mean_position', true,...
    'mean_value', true,...
    'max_position', true,...
    'max_value', true,...
    'radius', true...
);

%% Create light sources

lens_params_scene = lens_params;
lens_params_scene.ior_lens = lens_params.ior_lens(lens_params.ior_lens_reference_index);
[...
    X_lights, z_film, lights_filter, depth_factors...
] = imagingScenario(...
    lens_params_scene, ray_params.ior_environment, scene_params, plot_light_positions...
);

%% Run the simulation

[...
%    stats_real, stats_ideal, I, I_color...
    stats_real, stats_ideal...
] = doubleSphericalLensPSF(...
    lens_params, ray_params, image_params, X_lights, z_film, lights_filter,...
    request_spline_smoothing, depth_factors, doubleSphericalLensPSFVerbose...
);

%% Analyze the results

x_fields = struct(...
    'mean_position', 'mean_position',...
    'mean_value', 'mean_position',...
    'max_position', 'max_position',...
    'max_value', 'max_position',...
    'radius', 'mean_position'...
);

if verbose_aberration_ideal
    [...
        disparity_spline_ideal, disparity_ideal, disparity_ideal_radial...
    ] = radialChromaticAberration(...
        stats_ideal, x_fields, lens_params.ior_lens_reference_index,...
        depth_factors, 0,...
        lens_params.wavelengths, lens_params.wavelengths_to_rgb,...
        radialChromaticAberrationVerbose...
    );
else
    [...
        disparity_spline_ideal, disparity_ideal, disparity_ideal_radial...
    ] = radialChromaticAberration(...
        stats_ideal, x_fields, lens_params.ior_lens_reference_index,...
        depth_factors, 0 ...
    );
end

if verbose_aberration_real
    [...
        disparity_spline_real, disparity_real, disparity_real_radial...
    ] = radialChromaticAberration(...
        stats_real, x_fields, lens_params.ior_lens_reference_index,...
        depth_factors, 0,...
        lens_params.wavelengths, lens_params.wavelengths_to_rgb,...
        radialChromaticAberrationVerbose...
    );
else
    [...
        disparity_spline_real, disparity_real, disparity_real_radial...
    ] = radialChromaticAberration(...
        stats_real, x_fields, lens_params.ior_lens_reference_index,...
        depth_factors, 0 ...
    );
end
