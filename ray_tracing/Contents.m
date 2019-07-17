% RAY_TRACING
% Version 1.0.0 18-Jul-2019
%
% Raytracing, in particular to simulate images of point light sources at
% different positions relative to a single double-convex lens.
%
% Raytracing simulation
%   DoubleConvexThickLensPSF  - Ray tracing simulation of chromatic aberration
%   DoubleConvexThickLensPSF2 - Ray tracing simulation of chromatic aberration
%   doubleSphericalLens       - Trace rays through a lens with two spherical surfaces
%   doubleSphericalLensPSF    - Generate images of point light sources by raytracing
%   doubleSphericalLensPSF2   - Simulate pixelated images of point light sources by raytracing
%
% Point spread function and dispersion analysis
%   analyzePSF                - Describe a point spread function represented by an interpolant, or by samples
%   analyzePSFImage           - Describe a point spread function represented by an image
%   diskStatsToRGB            - Create the RGB equivalent of monochromatic PSF statistics
%
% Image formation
%   densifyRays               - Model image intensities from discrete samples of ray irradiance
%   densifyRaysImage          - Generate an image from discrete samples of ray irradiance
%
% Configuration
%   imageBoundaries           - Pad a grid of image positions to estimate image boundaries
%   imagingScenario           - Generate a set of point light source positions aligned over depths
%   lensParamsToRayParams     - Convert lens parameters to raytracing parameters
%
% Geometric optics
%   opticsFromLens            - Obtain imaging quantities from lens geometry
%   opticsToPSF               - Theoretical point spread functions from lens geometry
%   refract                   - Refract a 3D ray using Snell's law
%
% Dispersion analysis and calculation
%   radialChromaticAberration - Model chromatic aberration
%   sellmeierDispersion       - Calculate indices of refraction from the Sellmeier dispersion formula
%   statsToDisparity          - Calculate raw chromatic aberration vectors
%   TestSellmeierDispersion   - Test script for the 'sellmeierDispersion()' function
%
% Helper functions
%   neighborVertices          - Compute adjacency relationships in a triangulation
%   neighborVertices2         - Compute adjacency relationships in a triangulation
%   preallocateStats          - Return an empty PSF statistics structure array
%   sphereSection             - Sample points on a section of a spherical surface
