% ABERRATION_DATA
% Version 2.0.0 18-Jul-2019
%
% Image simulation and spectral data processing.
%
% Synthetic image generation
%   BimaterialImages      - Aberrated image generation
%   chirpImage            - Create an image with spatial and spectral linear chirp, and spectral dispersion
%   imageFormation        - Patch-wise conversion of a spectral image to RGB and RAW images
%
% Image noise simulation
%   addNoise              - Add noise to an image
%   noiseFractionToSNR    - Convert relative amounts of noise to a signal-to-noise ratio
%   shotNoise             - Add photon shot noise to a spectral image
%   snrToNoiseFraction    - Convert a signal-to-noise ratio to a relative amount of noise
%
% Spectral data generation
%   ciedIlluminant        - Spectral power distribution of a CIE D-Illuminant
%
% Spectral data manipulation and conversion to colour
%   cieSpectralToColor    - Convert spectral radiance to sRGB using the CIE tristimulus functions
%   reflectanceToColor    - Convert spectral reflectances to RGB colour values
%   reflectanceToRadiance - Convert spectral reflectances to spectral radiances
%
% Image and spectral data analysis
%   segmentColors         - Soft segmentation of image chromaticity values
%   smoothMetamer         - Find a smooth, non-negative metamer
