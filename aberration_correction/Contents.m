% ABERRATION_CORRECTION
% Version 2.0.0 18-Jul-2019
%
% Spectral or colour image reconstruction and chromatic aberration correction.
%
% Image reconstruction or correction scripts
%   CorrectByHyperspectralADMM - Demosaicing and hyperspectral ADMM-based correction of chromatic aberration
%   CorrectByWarping           - Correction of chromatic aberration by image warping
%
% Image transformation
%   channelConversionMatrix    - Create a sparse matrix to convert images between colour spaces
%   mosaicMatrix               - Create a sparse matrix to mosaic an image
%   warpImageSpectral          - Patch-wise warping and conversion to colour of a spectral image
%
% Regularization operators for optimization problems
%   spatialGradient            - Create sparse matrices acting as image spatial gradient operators
%   spatialLaplacian           - Create a sparse matrix acting as an image spatial Laplacian operator
%   spectralGradient           - Create a sparse matrix acting as an image spectral gradient operator
%   TestMatrices               - Test script for matrix generation functions
