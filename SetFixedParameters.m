%% Set Fixed Parameters
% Set values of parameters for image correction that seldomly need to be
% changed.
%
% ## Usage
%   Modify the parameters in the code below, as desired. This script exists
%   just to deduplicate code, and will be called by other scripts.
%
%   Some parameters can be given multiple values (indexed by row). Some
%   scripts will iterate through all rows, whereas others will just use the
%   first row of a parameter's values.
%
% ## Implementation Notes
% - When modifying this file, remember to update `parameters_list`.
% - Run this script after setting custom parameters in the calling script,
%   in order for the correct value of `parameters_list` to be generated.
%   The calling script must initialize `parameters_list` with its custom
%   parameter variable names.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 27, 2018

%% List of parameters to save with results
if ~exist('parameters_list', 'var')
    error('`parameters_list` should be initialized prior to running SetFixedParameters.m');
end
parameters_list = [parameters_list, {
    'bayer_pattern',...
    'samplingWeightsOptions',...
    'patch_sizes',...
    'paddings',...
    'use_fixed_weights',...
    'solvePatchesADMMOptions'...
    }];

%% Image parameters

% Colour-filter pattern
bayer_pattern = 'gbrg';

%% Spectral resampling parameters

% Options for 'samplingWeights()'. Refer to the documentation of
% 'samplingWeights.m' for more details.

% Integration method to use for colour calculations. If the latent space
% consists of wavelength bands, use this type of numerical integration in
% 'integrationWeights()' within 'samplingWeights()'. (Otherwise,
% 'samplingWeights()' should not even be called.)
samplingWeightsOptions.int_method = 'trap';

samplingWeightsOptions.power_threshold = 0.99;

samplingWeightsOptions.support_threshold = 0.05;

samplingWeightsOptions.bands_padding = 1000;

%% Hyperspectral image estimation parameters

% ## Image estimation options

solvePatchesADMMOptions.admm_options = struct;

% Whether to make the spectral gradient the same size as the image
solvePatchesADMMOptions.admm_options.full_GLambda = false;

% Penalty parameters in ADMM, the `rho` input argument.
% Sample values seem to be in the range 1-10 (see pages 89, 93, and 95 of
% Boyd et al. 2011)
solvePatchesADMMOptions.admm_options.rho = [ 1, 1, 1, 1 ];

% Weights on the prior terms. Baek et al. (2017) used [1e-5, 0.1]. Setting
% elements to zero disables the corresponding regularization term during
% image estimation.
weights = [ 1e-2, 0, 0 ];

% Convergence tolerances in ADMM. Reasonable values for the third element
% are 10^-4 to 10^-3 (page 21 of Boyd et al. 2011).
solvePatchesADMMOptions.admm_options.tol = [ 1e-5, 1e-5, 1e-5 ];

% Maximum number of inner and outer iterations, the `maxit` input argument
solvePatchesADMMOptions.admm_options.maxit = [ 500, 500 ];

% Parameters for adaptively changing the penalty parameters for improved
% convergence speed. (Disable adaptive penalty parameter variation by
% setting this option to an empty array.)
solvePatchesADMMOptions.admm_options.varying_penalty_params = [2, 2, 10];

% Types of norms to use on the prior terms
solvePatchesADMMOptions.admm_options.norms = [false, false, false];

% Whether to apply a non-negativity constraint (in which case, `rho` must
% have four elements)
solvePatchesADMMOptions.admm_options.nonneg = true;

% ## Options for patch-wise image estimation

% Every combination of rows of `patch_sizes` and elements of `paddings`
% will be tested by some image estimation pipelines, and if `patch_sizes`
% is empty only whole image estimation may be performed. Some image
% estimation pipelines only use the first row of `patch_sizes` and the
% first element of `paddings`, via `solvePatchesADMMOptions.patch_options`.
patch_sizes = [ % Each row contains a (number of rows, number of columns) pair
   30 30;
]; 
paddings = 10;

solvePatchesADMMOptions.patch_options = struct;
solvePatchesADMMOptions.patch_options.patch_size = patch_sizes(1, :);
solvePatchesADMMOptions.patch_options.padding = paddings(1);

% ## Options for selecting regularization weights

solvePatchesADMMOptions.reg_options = struct;

solvePatchesADMMOptions.reg_options.enabled = logical(weights(1, :));

use_fixed_weights = false;
if use_fixed_weights
    solvePatchesADMMOptions.reg_options.minimum_weights = weights;
    solvePatchesADMMOptions.reg_options.maximum_weights = weights;
else
    % Minimum values to use for regularization weights (and to use to set
    % the origin of the minimum distance function)
    solvePatchesADMMOptions.reg_options.minimum_weights = eps * ones(1, length(weights));
    % Maximum values to use for regularization weights (and to use to set
    % the origin of the minimum distance function)
    solvePatchesADMMOptions.reg_options.maximum_weights = 1e10 * ones(1, length(weights));
end
solvePatchesADMMOptions.reg_options.low_guess = [1e-3, 1e-3, 1e-3];
solvePatchesADMMOptions.reg_options.high_guess = [100, 100, 100];
solvePatchesADMMOptions.reg_options.tol = 1e-6;

% Maximum and minimum number of grid search iterations
% Song et al. 2016 used a fixed number of 6 iterations, but I don't know
% what range of regularization weights they were searching within.
solvePatchesADMMOptions.reg_options.n_iter = [30, 6];

%% ## Debugging Flags

solvePatchesADMMVerbose = true;
samplingWeightsVerbose = false;