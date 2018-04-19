function [ ray_params ] = lensParamsToRayParams(ray_params_in, lens_params, z_film)
% LENSPARAMSTORAYPARAMS  Convert lens parameters to raytracing parameters
%
% ## Syntax
% ray_params = lensParamsToRayParams(ray_params_in, lens_params)
%
% ## Description
% ray_params = lensParamsToRayParams(ray_params_in, lens_params)
%   Convert lens parameters in a format convenient for lens design to lens
%   parameters in a format convenient for raytracing.
%
% ## Input Arguments
%
% ray_params_in -- Raytracing parameters structure
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
% z_film -- Image plane location
%   The z-coordinate of the image plane.
%
% ## Output Arguments
%
% ray_params -- Raytracing parameters structure
%   A structure with the following fields:
%   - radius_front: The radius of curvature of the front surface of the
%     lens
%   - theta_aperture_front: The limit of the front aperture of the lens,
%     expressed as the angle with the optical axis subtended by the
%     aperture. In other words, half the angle subtended by the aperture
%     when viewed from the center of a circle with radius equal to the
%     radius of the front surface.
%   - radius_back: The radius of curvature of the back surface of the
%     lens
%   - theta_aperture_back: Similar to `theta_aperture_front`, but for the
%     back surface of the lens.
%   - d_lens: The distance from the centre of the sphere corresponding to
%     the front of the lens to the centre of the sphere corresponding to
%     the back of the lens. A shift in the positive direction moves the
%     back of the lens in the negative z-direction.
%   - n_incident_rays: Number of samples over the front of the lens. Each
%     sample produces one ray from the point light source through the front
%     surface of the lens. The front aperture is uniformly sampled, but
%     samples are culled if they are occluded by the front lens surface
%     from the perspective of the point light source.
%   - ior_environment: The refractive index of the surrounding medium
%   - d_film: The distance of the image plane from the centre of the sphere
%     corresponding to the front surface of the lens. (A positive value if
%     the image plane is in the negative z-direction relative to the centre
%     of the sphere.)
%
% ## Notes
%
% ### Coordinate System
% - The radii of both faces of the lens are positive when the lens is
%   biconvex.
% - The positive z-axis points towards the front of the lens, along
%   the optical axis, assuming the front of the lens is convex.
%
% See also doubleSphericalLens, doubleSphericalLensPSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 19, 2018

nargoutchk(1,1);
narginchk(3,3);

ray_params = ray_params_in;
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
ray_params.d_film = -z_film;

end

