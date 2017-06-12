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
% ## References
% - Ray-sphere intersection testing:
%   - http://www.ccs.neu.edu/home/fell/CSU540/programs/RayTracingFormulas.htm
%   - https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
%   - https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
% - Uniform sampling of the surface of a sphere:
%   http://mathworld.wolfram.com/SpherePointPicking.html

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2017

%% Input data and parameters

% ## Raytracing parameters
% Refer to the documentation of `doubleSpherical
ray_params.source_position = [0, 0, 10];
ray_params.radius_front = 2.0;
ray_params.theta_aperture_front = pi / 2;
ray_params.radius_back = 2.0;
ray_params.theta_aperture_back = pi / 6;
ray_params.d_lens = 3;
ray_params.n_incident_rays = 50;
ray_params.sample_random = false;
ray_params.ior_environment = 1.0;
ray_params.ior_lens = 1.52;
ray_params.d_film = 10;

% Debugging Flags
verbose_ray_tracing = true;
verbose_image_formation = true;

%% Trace rays through the lens

[ ...
    image_position, ray_irradiance, ~, incident_position_cartesian ...
] = doubleSphericalLens( ray_params, verbose_ray_tracing );

%% Form rays into an image

I = densifyRays(...
    incident_position_cartesian,...
    ray_params.radius_front, image_position,...
    ray_irradiance, verbose_image_formation ...
);