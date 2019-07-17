% EVALUATION
% Version 2.0.0 18-Jul-2019
%
% Execution and comparison of image estimation algorithms on datasets of images.
%
% Dataset preprocessing
%   Choi2017Input           - Generate input data for the modified version of Choi et al. 2017's code
%   Choi2017Output          - Convert output data from the modified version of Choi et al. 2017's code
%   SelectWeightsForDataset - Select regularization weights for images in a dataset
%
% Synthetic experiments
%   CorrectChirpImage       - Evaluate hyperspectral ADMM-based reconstruction of a chirp image
%
% General experiments (dataset processing)
%   describeDataset         - Create a structure describing a dataset for evaluation
%   RunOnDataset            - Evaluate demosaicking, spectral reconstruction, and/or chromatic aberration correction
%   SetAlgorithms           - List of algorithms to test
%
% General image similarity evaluation
%   evaluateAndSaveRGB      - Compare colour images and save the comparison results
%   evaluateAndSaveSpectral - Compare spectral images and save the comparison results
%   evaluateRGB             - Compare colour images
%   evaluateSpectral        - Compare spectral images
%   mergeRGBTables          - Combine RGB comparison results across images
%   mergeSpectralTables     - Combine spectral comparison results across images
%   metrics                 - Compute image evaluation metrics
%
% Specialized evaluation
%   EvaluateColorChecker    - Evaluate demosaicking, spectral reconstruction, and chromatic aberration correction
%   PointSpectralEvaluation - Evaluate spectral reconstruction with respect to measured spectra
%
% Dispersion model loading
%   loadDispersionModel     - Load dispersion information from a file
%   makeDispersionForImage  - Create an image-specific dispersion model