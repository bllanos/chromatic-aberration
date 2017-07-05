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
% Presently, the output arguments of 'doubleSphericalLensPSF.m'.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

%% Input data and parameters
% Refer to the documentation of `doubleSphericalLensPSF` for details

% ## Raytracing parameters

% ### Lens parameters
% Based on
% 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170613_SimpleLenses_EdmundOptics\25mmDiameter40mmFLUncoatedDoubleConvexLens_prnt_45296.pdf'
% and
% 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170613_SimpleLenses_EdmundOptics\25mmDiameter40mmFLUncoatedDoubleConvexLens.pdf'
lens_params.radius_front = 40.42;
lens_params.lens_radius = 25 / 2;
lens_params.axial_thickness = 5.30;
lens_params.radius_back = lens_params.radius_front;

ray_params.n_incident_rays = 100;
ray_params.sample_random = false;
ray_params.ior_environment = 1.0;

% #### Index of refraction
% The focal length specification wavelength is 587.6 nm
% At this wavelength, N-BK7 glass has a refractive index of 1.51680
% (https://www.pgo-online.com/intl/katalog/BK7.html)
% The refractive index at 480.0 nm is 1.52283
% The refractive index at 643.8 nm is 1.51472
lens_params.ior_lens = [1.51472, 1.51680, 1.52283]; % Red, Green, Blue

% Index of the wavelength/index of refraction to be used to position the
% image plane
lens_params.ior_lens_reference_index = 2; % Use the green channel

% Wavelength values corresponding to indices of refraction (for display
% purposes only)
lens_params.wavelengths = [643.8, 587.6, 480.0];

% Obtained using the quantum efficiencies presented in
% 'C:\Users\llanos\Google Drive\ThesisResearch\Equipment\FLEA3\20170508_FL3_GE_EMVA_Imaging Performance Specification.pdf'
% Image sensor: Sony ICX655, 2/3", Color (page 19)
lens_params.wavelengths_to_rgb = [
    0.38, 0.04, 0;
    0.12, 0.35, 0;
    0, 0.17, 0.44
    ];
% Normalize, for improved colour saturation
lens_params.wavelengths_to_rgb = lens_params.wavelengths_to_rgb ./...
    max(max(lens_params.wavelengths_to_rgb));

% ### Ray interpolation parameters
image_params.image_sampling = [20, 20];
image_params.normalize_psfs_before_combining = false;
image_params.normalize_color_images_globally = true;

% ## Scene setup
scene_params.theta_max = pi / 12;
scene_params.theta_min = pi / 24;
scene_params.n_lights = [4 4];
scene_params.light_distance_factor_focused = 3;
scene_params.light_distance_factor_larger = [5, 2];
scene_params.light_distance_factor_smaller = [2, 1];
scene_params.preserve_angle_over_depths = true;

% ## Debugging Flags
verbose.plot_light_positions = true;
verbose.verbose_ray_tracing = false;
verbose.verbose_ray_interpolation = false;
verbose.verbose_psf_analysis = false;
verbose.display_each_psf = false;
verbose.display_all_psf_each_ior = false;
verbose.display_all_psf_each_depth = true;
verbose.display_summary = true;

%% Run the simulation

[...
    stats_real, stats_ideal,...
    X_lights, depth_factors, I, I_color...
] = doubleSphericalLensPSF(...
    lens_params, ray_params, image_params, scene_params, verbose...
);

%% Analyze the results

[ disparity_spline_ideal, disparity_ideal ] = radialChromaticAberration(...
    stats_ideal, lens_params.ior_lens_reference_index,...
    depth_factors, 0,...
    lens_params.wavelengths, lens_params.wavelengths_to_rgb...
);

[ disparity_spline_real, disparity_real ] = radialChromaticAberration(...
    stats_real, lens_params.ior_lens_reference_index,...
    depth_factors, 0,...
    lens_params.wavelengths, lens_params.wavelengths_to_rgb...
);