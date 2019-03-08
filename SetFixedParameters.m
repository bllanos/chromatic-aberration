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
%
% ## References
%
% Third-party algorithms:
% - Sun, T., Peng, Y., & Heidrich, W. (2017). "Revisiting cross-channel
%   information transfer for chromatic aberration correction." In 2017 IEEE
%   International Conference on Computer Vision (ICCV) (pp. 3268–3276).
%   doi:10.1109/ICCV.2017.352
% - Krishnan, D., Tay, T. & Fergus, R. (2011). "Blind deconvolution using a
%   normalized sparsity measure." In IEEE Conference on Computer Vision and
%   Pattern Recognition (CVPR) (pp. 233–240).

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 27, 2018

%% List of parameters to save with results
if ~exist('parameters_list', 'var')
    error('`parameters_list` should be initialized prior to running SetFixedParameters.m');
end
parameters_list = [parameters_list, {
    'criteria',...
    'save_all_images',...
    'bayer_pattern',...
    'findSamplingOptions',...
    'patch_sizes',...
    'paddings',...
    'use_fixed_weights',...
    'solvePatchesADMMOptions',...
    'solvePatchesMultiADMMOptions'...
    'krishnan2011Options',...
    'sun2017Options'...
    }];

%% Evaluation parameters

% Enable or disable different methods for selecting regularization weights.
% All enabled methods will be run in 'SelectWeightsForDataset.m'. In
% 'RunOnDataset.m', if regularization weights selected using
% 'SelectWeightsForDataset.m' are used, then this array is ignored.
% Otherwise, if regularization weights are automatically selected, then all
% enabled methods are run.
criteria = [
    false; % Minimum distance criterion
    false; % Similarity with the true image
    true % Similarity with a demosaicing result
    ];
mdc_index = 1;
mse_index = 2;
dm_index = 3;

% Fields used for saving regularization weights in
% 'SelectWeightsForDataset.m'
criteria_fields = {'mdc_weights', 'mse_weights', 'dm_weights'};

% Visualization and file output variables common to multiple scripts
criteria_names = {'Minimum distance criterion', 'Mean square error', 'Demosaic mean square error'};
criteria_abbrev = {'MDC', 'MSE', 'DM'};
criteria_filenames = {'mdc_', 'mse_', 'dm_'};
criteria_colors = eye(3);

% ### Output images

% One of each of the following types of images can be created for each
% input image. The filename of the input image, concatenated with a string
% of parameter information, is represented by '*' below.
% - '*_roi.tif' and '*_roi.mat': A cropped version of the input image
%   (stored in the variable 'I_raw'), containing the portion used as input
%   for ADMM. This region of interest was determined using the
%   `model_space` and `fill` variables saved in an input model of
%   dispersion data file. If these variables were not present, the cropped
%   region is the entire input image. All of the other output images listed
%   below are limited to the region shown in '*_roi.tif'.
% - '*_latent.tif' and '*_latent.mat': The latent image estimated using
%   ADMM (stored in the variable 'I_latent'). The '.tif' image is only
%   output if the latent images are greyscale or 3-channel images.
% - '*_warped.tif' and '*_warped.mat': A version of the latent image
%   (stored in the variable 'I_warped') created by warping the latent image
%   according to the dispersion model. The '.tif' image is only output if
%   the latent images are greyscale or 3-channel images.
% - '*_rgb.tif': A colour image created by converting the latent image to
%   the RGB colour space of the camera.
% - '*_rgb_warped.tif' and '*_rgb_warped.mat': A colour image (stored in
%   the variable 'J_full') created by warping the latent image according to
%   the dispersion model, then converting the image to the RGB colour space
%   of the camera. This output image is, in a sense, a demosaiced version
%   of the input image.
% - '*_reestimated.tif' and '*_reestimated.mat': A simulation (stored in
%   the variable 'J_est') of the input RAW image from the latent image,
%   useful for visually evaluating the convergence of the ADMM algorithm.
%
% Of the above types of images, the following will only be saved if the
% flag below is `true`:
% - '*_roi.tif'
% - '*_warped.tif' and '*_warped.mat'
% - '*_rgb_warped.tif' and '*_rgb_warped.mat'
% - '*_reestimated.tif' and '*_reestimated.mat'
save_all_images = false;

% Not all scripts follow the above guidelines. In particular, image
% generation scripts may output all possible images regardless of this
% flag, especially if it will save time later during image estimation and
% evaluation.

%% Image parameters

% Colour-filter pattern
bayer_pattern = 'gbrg';

%% Spectral resampling parameters

% Options for 'findSampling()'. Refer to the documentation of
% 'findSampling.m' for more details.

% Integration method to use for colour calculations. If the latent space
% consists of wavelength bands, use this type of numerical integration in
% 'integrationWeights()' within 'findSampling()'. (Otherwise,
% 'findSampling()' should not even be called.)
findSamplingOptions.int_method = 'trap';

findSamplingOptions.power_threshold = 0.95;
% As an alternative to automatically determining the number of spectral
% bands, according to `findSamplingOptions.power_threshold`, set it
% explicitly (if the following option is an integer greater than zero).
findSamplingOptions.n_bands = 0;

findSamplingOptions.support_threshold = 0.05;

findSamplingOptions.bands_padding = 1000;

% Interpolation function for estimated spectral data:
%
% `x = 0` is the current interpolation location, and an increment or
% decrement of one unit in `x` represents a shift equal to the spacing
% between samples in the sequence of samples being interpolated. The
% interpolation function `f(x)` returns the weight for a sample at location
% `x` relative to the current  interpolation location.
findSamplingOptions.interpolant = @gaussian;

% Interpolation function for other spectral data, such as sensor spectral
% sensitivities or ground truth spectral radiances
findSamplingOptions.interpolant_ref = @triangle;

% Additional options for 'solvePatchesMultiADMM()'
solvePatchesMultiADMMOptions.sampling_options = findSamplingOptions;

% How to choose spectral resolutions lower than the one given by
% 'findSampling()' based on the above options.
solvePatchesMultiADMMOptions.sampling_options.progression = 'doubling';

% Output the results for the lower spectral resolutions. CAUTION: Not
% recommended when estimating large images, because of memory consumption.
solvePatchesMultiADMMOptions.sampling_options.show_steps = false;

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
% The first element is the tolerance for the conjugate gradients method. MATLAB
% uses a default value of 10^-6
solvePatchesADMMOptions.admm_options.tol = [ 1e-5, 1e-5, 1e-5 ];

% Maximum number of inner and outer iterations, the `maxit` input argument.
% The first element applies to the conjugate gradients method. MATLAB
% uses a default value of 20.
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

% Image initialization method
solvePatchesADMMOptions.admm_options.init = 'zero';

solvePatchesMultiADMMOptions.admm_options = solvePatchesADMMOptions.admm_options;

% ## Options for patch-wise image estimation

% Every combination of rows of `patch_sizes` and elements of `paddings`
% will be tested by some image estimation pipelines, and if `patch_sizes`
% is empty only whole image estimation may be performed.
patch_sizes = [ % Each row contains a (number of rows, number of columns) pair
   30 30;
]; 
paddings = 10;

solvePatchesADMMOptions.patch_options = struct;
solvePatchesADMMOptions.patch_options.patch_size = patch_sizes(1, :);
solvePatchesADMMOptions.patch_options.padding = paddings(1);

solvePatchesMultiADMMOptions.patch_options = solvePatchesADMMOptions.patch_options;

% ## Options for selecting regularization weights

solvePatchesADMMOptions.reg_options = struct;

solvePatchesADMMOptions.reg_options.enabled = logical(weights(1, :));

solvePatchesADMMOptions.reg_options.low_guess = [1e-3, 1e-3, 1e-3];
solvePatchesADMMOptions.reg_options.high_guess = [1e3, 1e3, 1e3];
solvePatchesADMMOptions.reg_options.tol = 1e-6;

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

% Maximum and minimum number of grid search iterations
% Song et al. 2016 used a fixed number of 6 iterations, but I don't know
% what range of regularization weights they were searching within.
solvePatchesADMMOptions.reg_options.n_iter = [24, 6];

% Select regularization weights based on similarity to a demosaicking
% result, instead of using the minimum distance criterion, if no true image
% is provided for regularization weight selection. Scripts which use
% multiple regularization weight selection methods will override this
% option.
solvePatchesADMMOptions.reg_options.demosaic = true;

solvePatchesMultiADMMOptions.reg_options = solvePatchesADMMOptions.reg_options;

%% ## Parameters for third-party algorithms

% ## Parameters for Krishnan et al. 2011

% Kernel sizes must be odd integers. I will set the kernel size based on the
% estimated amount of dispersion.
krishnan2011Options.kernel_sz = 9;

% Value tuned by finding the 'knee' in the error plot output by
% 'TuneSunEtAl2017.m'
krishnan2011Options.min_lambda = 200;

% Window in which to estimate the PSF: (y1, x1, y2, x2) of the top left and
% bottom right corners. (Set it to an empty array to use the entire image.)
krishnan2011Options.kernel_est_win = [];

krishnan2011Options.prescale = 1;
krishnan2011Options.k_reg_wt = 1;
krishnan2011Options.gamma_correct = 1;
krishnan2011Options.k_thresh = 0.0;
krishnan2011Options.kernel_init = 3;
krishnan2011Options.delta = 0.001;
krishnan2011Options.x_in_iter = 2; 
krishnan2011Options.x_out_iter = 2;
krishnan2011Options.xk_iter = 21;
krishnan2011Options.nb_lambda = 3000;
krishnan2011Options.nb_alpha = 1.0;
krishnan2011Options.use_ycbcr = 1;

% ## Parameters for Sun et al. 2017

% PSF estimation window size, as a fraction of the image's largest dimension.
% Optimized using 'TuneSunEtAl2017.m'
sun2017Options.psf_sz = 0.1714;

% CCT implementation window size, in pixels
% Optimized using 'TuneSunEtAl2017.m'
sun2017Options.win_sz = 5;

sun2017Options.alpha = 0.3;
sun2017Options.beta = 0.3;
sun2017Options.iter = 3;

% Index of the reference colour channel (Green)
sun2017Options.reference_channel_index = 2;

%% ## Debugging Flags

findSamplingVerbose = true;
solvePatchesADMMVerbose = true;
solvePatchesMultiADMMVerbose = true;