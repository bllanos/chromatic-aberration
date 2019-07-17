% DISPERSION_MODEL
% Version 2.0.0 18-Jul-2019
%
% Generation, manipulation, and visualization of models of dispersion
% (approximated as lateral chromatic aberration). Additional visualization
% scripts can be found in '../data_analysis/'.
%
% Dispersion model generation scripts
%   DoubleConvexThickLensDispersion - Ray tracing simulation of an image dispersion function
%   RegistrationDispersion          - Chromatic aberration calibration from captured images
%
% Loading and instantiation
%   dispersionfunToMatrix           - Convert a warp model to a warp matrix, or use it
%   makeDispersionfun               - Create a function to evaluate a model of disparity in three variables
%   modelSpaceTransform             - Create a transformation from pixel to dispersion model coordinates
%   pixelsToWorldTransform          - Create a transformation from pixel to world coordinates
%
% Dispersion model generation
%   registerPatches                 - Find patch-wise displacements between the channels of an image
%   splineKernel2D                  - 2D thin-plate spline kernel function
%   splineKernel3D                  - 3D thin-plate spline kernel function
%   xylambdaPolyfit                 - Fit a polynomial model in three variables
%   xylambdaSplinefit               - Fit a spline model in two or three variables
%
% Dispersion model manipulation
%   mergeDispersionModelROI         - Find the intersection of dispersion model domains
%
% Visualization
%   plotXYLambdaModel               - Plot a model of dispersion in two or three variables
%   plotXYLambdaModel2              - Plot a model of dispersion in two or three variables
%
% Generic helpers
%   matchByVectors                  - Associate structures by a vector-valued field

