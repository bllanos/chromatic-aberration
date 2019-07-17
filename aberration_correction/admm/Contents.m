% ADMM
% Version 2.0.0 18-Jul-2019
%
% Image reconstruction optimization problem setup and resolution.
%
% Alternating Direction Method of Multipliers
%   baek2017Algorithm2LowMemory     - Run ADMM (loosely) as in Algorithm 2 of Baek et al. 2017
%   initBaek2017Algorithm2LowMemory - Allocate memory for 'baek2017Algorithm2LowMemory()'
%   solvePatchesColor               - Run ADMM (loosely) as in Algorithm 2 of Baek et al. 2017, with weight selection and patch-wise decomposition
%   solvePatchesSpectral            - Run ADMM for spectral image estimation
%
% Regularization weight selection for the optimization problem
%   initPenalties                   - Allocate memory for 'penalties()'
%   initWeightsLowMemory            - Allocate memory for 'weightsLowMemory()'
%   penalties                       - Calculate data fitting and regularization errors
%   weightsLowMemory                - Select regularization weights using a grid search