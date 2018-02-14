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
% 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170613_SimpleLenses_EdmundOptics\25mmDiameter40mmFLUncoatedDoubleConvexLens_prnt_45296.pdf'
% and
% 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170613_SimpleLenses_EdmundOptics\25mmDiameter40mmFLUncoatedDoubleConvexLens.pdf'
lens_params.radius_front = 40.42;
lens_params.lens_radius = 25 / 2;
lens_params.axial_thickness = 5.30;
lens_params.radius_back = lens_params.radius_front;

ray_params.n_incident_rays = 250;
ray_params.sample_random = false;
ray_params.ior_environment = 1.0;

% #### Index of refraction
% The focal length specification wavelength is 587.6 nm
% At this wavelength, N-BK7 glass has a refractive index of 1.51680
% (https://www.pgo-online.com/intl/katalog/BK7.html)
% The refractive index at 480.0 nm is 1.52283
% The refractive index at 546.1 nm is 1.51872
% The refractive index at 643.8 nm is 1.51472
lens_params.ior_lens = [1.51472, 1.51680, 1.51872, 1.52283]; % Red, Green, Green, Blue

% Index of the wavelength/index of refraction to be used to position the
% image plane
lens_params.ior_lens_reference_index = 2; % Use the green channel

% Wavelength values corresponding to indices of refraction (for display
% purposes only)
lens_params.wavelengths = [643.8, 587.6, 546.1, 480.0];

% Obtained using the quantum efficiencies presented in
% 'C:\Users\llanos\Google Drive\ThesisResearch\Equipment\FLEA3\20170508_FL3_GE_EMVA_Imaging Performance Specification.pdf'
% Image sensor: Sony ICX655, 2/3", Color (page 19)
lens_params.wavelengths_to_rgb = [
    0.38, 0.04, 0;
    0.12, 0.35, 0;
    0.03, 0.52, 0.04;
    0, 0.17, 0.44
    ];
% Normalize, for improved colour saturation
lens_params.wavelengths_to_rgb = lens_params.wavelengths_to_rgb ./...
    max(max(lens_params.wavelengths_to_rgb));

% ### Ray interpolation parameters
image_params.image_sampling = [600, 400];
image_params.normalize_psfs_before_combining = false;
image_params.normalize_color_images_globally = true;

% ### Use spline interpolation to reduce noise
request_spline_smoothing = false;

% ## Scene setup
scene_params.theta_min = deg2rad(15);
scene_params.theta_max = deg2rad(30);
scene_params.n_lights = [3 2];
scene_params.light_distance_factor_focused = 2;
scene_params.light_distance_factor_larger = [4, 2];
scene_params.light_distance_factor_smaller = [1.5, 2];
scene_params.preserve_angle_over_depths = true;

% ## Debugging Flags
plot_light_positions = true;

doubleSphericalLensPSFVerbose.verbose_ray_tracing = false;
doubleSphericalLensPSFVerbose.verbose_ray_interpolation = false;
doubleSphericalLensPSFVerbose.verbose_psf_analysis = false;
doubleSphericalLensPSFVerbose.display_each_psf = false;
doubleSphericalLensPSFVerbose.display_all_psf_each_ior = false;
doubleSphericalLensPSFVerbose.display_all_psf_each_depth = true;
doubleSphericalLensPSFVerbose.display_summary = true;

verbose_aberration_ideal = true;
verbose_aberration_real = true;
radialChromaticAberrationVerbose.display_raw_values = false;
radialChromaticAberrationVerbose.display_raw_disparity = false;
radialChromaticAberrationVerbose.display_stats_splines = true;
radialChromaticAberrationVerbose.display_spline_differences = true;
radialChromaticAberrationVerbose.filter = struct(...
    'mean_position', true,...
    'mean_value', false,...
    'max_position', false,...
    'max_value', false,...
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
