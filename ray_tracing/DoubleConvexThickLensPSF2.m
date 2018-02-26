%% Ray tracing simulation of chromatic aberration
% Simulate the chromatic point spread function of a thick (biconvex) lens.
% Use an image resolution-dependent simulation method, for realistic
% generation of colour images with chromatic aberration.
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
% File created February 22, 2018

%% Input data and parameters
% Refer to the documentation of `doubleSphericalLensPSF2` for details

% ## Raytracing parameters

% ### Lens parameters
% Based on
% '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20170613_SimpleLenses_EdmundOptics/25mmDiameter40mmFLUncoatedDoubleConvexLens_prnt_45296.pdf'
% and
% '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20170613_SimpleLenses_EdmundOptics/25mmDiameter40mmFLUncoatedDoubleConvexLens.pdf'
lens_params.lens_radius = 25 / 2;
lens_params.axial_thickness = 5.30;
lens_params.radius_front = 40.42;
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

lens_params.wavelengths = linspace(300, 1100, 1000);
lens_params.ior_lens = sellmeierDispersion(lens_params.wavelengths, sellmeierConstants);

% Index of the wavelength/index of refraction to be used to position the
% image plane
[~, ior_lens_reference_index] = min(abs(lens_params.wavelengths - 587.6));

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
image_params.intensity_threshold = 0.1;

% ## Scene setup
scene_params.theta_min = deg2rad(0.1);
scene_params.theta_max = deg2rad(1);
scene_params.n_lights = [9 9];
scene_params.light_distance_factor_focused = 2;
scene_params.light_distance_factor_larger = [4, 2];
scene_params.light_distance_factor_smaller = [1.5, 2];
scene_params.preserve_angle_over_depths = true;

% ## Data analysis parameters
color_names = {'Red', 'Green', 'Blue'};
colors_to_rgb = [
    1, 0, 0;
    0, 1, 0;
    0, 0, 1
    ];
reference_channel_index = 2;

% ## Debugging Flags
plot_light_positions = true;

doubleSphericalLensPSF2Verbose.verbose_ray_tracing = false;
doubleSphericalLensPSF2Verbose.verbose_ray_interpolation = false;
doubleSphericalLensPSF2Verbose.verbose_psf_analysis = false;
doubleSphericalLensPSF2Verbose.display_each_psf = false;
doubleSphericalLensPSF2Verbose.display_each_psf_rgb = true;
doubleSphericalLensPSF2Verbose.display_all_psf_each_depth = true;
doubleSphericalLensPSF2Verbose.display_summary = false;

verbose_aberration_ideal = false;
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
lens_params_scene.ior_lens = lens_params.ior_lens(ior_lens_reference_index);
[...
    X_lights, z_film, lights_filter, depth_factors...
] = imagingScenario(...
    lens_params_scene, ray_params.ior_environment, scene_params, plot_light_positions...
);

%% Run the simulation

[...
    stats_real, stats_ideal,  I_color...
] = doubleSphericalLensPSF2(...
    lens_params, ray_params, image_params, X_lights, z_film, lights_filter,...
    depth_factors, doubleSphericalLensPSF2Verbose...
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
        stats_ideal, x_fields, reference_channel_index,...
        depth_factors, 0,...
        color_names, colors_to_rgb,...
        radialChromaticAberrationVerbose...
    );
else
    [...
        disparity_spline_ideal, disparity_ideal, disparity_ideal_radial...
    ] = radialChromaticAberration(...
        stats_ideal, x_fields, reference_channel_index,...
        depth_factors, 0 ...
    );
end

if verbose_aberration_real
    [...
        disparity_spline_real, disparity_real, disparity_real_radial...
    ] = radialChromaticAberration(...
        stats_real, x_fields, reference_channel_index,...
        depth_factors, 0,...
        color_names, colors_to_rgb,...
        radialChromaticAberrationVerbose...
    );
else
    [...
        disparity_spline_real, disparity_real, disparity_real_radial...
    ] = radialChromaticAberration(...
        stats_real, x_fields, reference_channel_index,...
        depth_factors, 0 ...
    );
end
