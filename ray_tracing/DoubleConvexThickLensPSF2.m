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

% List of parameters to save with results
parameters_list = {
        'lens_params',...
        'ray_params',...
        'sellmeierConstants',...
        'ior_lens_reference_index',...
        'image_params',...
        'pixel_size',...
        'scene_params',...
        'color_names',...
        'colors_to_rgb',...
        'reference_channel_index',...
        'output_directory',...
        'bayer_pattern'...
    };

%% Input data and parameters
% Refer to the documentation of `doubleSphericalLensPSF2` for details

% ## Raytracing parameters

% ### Lens parameters
% Based on
% '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180226_SmallFLLenses_EdmundOptics/3mmDiameter4dot5mmFLUncoatedDoubleConvexLens_prnt_32022.pdf'
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

lens_params.wavelengths = linspace(300, 1100, 100);
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
scene_params.theta_min = deg2rad(0);
scene_params.theta_max = deg2rad(20);
scene_params.n_lights = [6 6];
scene_params.light_distance_factor_focused = 10;
scene_params.light_distance_factor_larger = [4, 0];
scene_params.light_distance_factor_smaller = [1.5, 0];
scene_params.preserve_angle_over_depths = true;

% ## Data analysis parameters
disable_analysis = true;
color_names = {'Red', 'Green', 'Blue'};
colors_to_rgb = [
    1, 0, 0;
    0, 1, 0;
    0, 0, 1
    ];
reference_channel_index = 2;

% ## Output parameters

% Output directory for all images and saved parameters and data
% If empty (`[]`), output will be disabled.
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180523_TestingImageWarping/ray_tracing';

% Colour-filter pattern for RAW image output
bayer_pattern = 'gbrg';

% ## Debugging Flags
plot_light_positions = true;

doubleSphericalLensPSF2Verbose.verbose_ray_tracing = false;
doubleSphericalLensPSF2Verbose.verbose_ray_interpolation = false;
doubleSphericalLensPSF2Verbose.verbose_psf_analysis = false;
doubleSphericalLensPSF2Verbose.display_each_psf = false;
doubleSphericalLensPSF2Verbose.display_each_psf_rgb = false;
doubleSphericalLensPSF2Verbose.display_all_psf_each_depth = false;
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

if disable_analysis
    I_color = doubleSphericalLensPSF2(...
        lens_params, ray_params, image_params, X_lights, z_film, lights_filter,...
        'images_only', depth_factors, doubleSphericalLensPSF2Verbose...
    );
else
    [...
        stats_real, stats_ideal,  I_color...
    ] = doubleSphericalLensPSF2(...
        lens_params, ray_params, image_params, X_lights, z_film, lights_filter,...
        depth_factors, doubleSphericalLensPSF2Verbose...
    );
end

%% Analyze the results
if ~disable_analysis
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
end

%% Save results to a file

if ~isempty(output_directory)
    ext = '.tif';
    cdate = replace(datestr(now, 31), {'-',' ',':'},'_');
    
    I_raw = mosaic(I_color, bayer_pattern);
    for i = 1:length(depth_factors)
        I_color_filename = [ cdate '_color_depth' num2str(depth_factors(i)) ext];
        I_raw_filename = [ cdate '_raw_depth' num2str(depth_factors(i)) ext];
        imwrite(I_color(:, :, :, i), fullfile(output_directory, I_color_filename));
        imwrite(I_raw(:, :, i), fullfile(output_directory, I_raw_filename));
    end

    save_variables_list = [ parameters_list, {...
            'X_lights',...
            'z_film',...
            'lights_filter',...
            'depth_factors'...
        } ];
    if ~disable_analysis
        save_variables_list = [ save_variables_list, {
            'stats_real', 'stats_ideal',...
            'disparity_spline_ideal', 'disparity_ideal', 'disparity_ideal_radial',...
            'disparity_spline_real', 'disparity_real', 'disparity_real_radial'...
        } ];
    end
    save_data_filename = fullfile(output_directory, [cdate '_DoubleConvexThickLensPSF2.mat']);
    save(save_data_filename, save_variables_list{:});
end
