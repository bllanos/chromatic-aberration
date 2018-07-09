function [ I_3D, image_bounds, varargout ] = baek2017Algorithm2(...
    image_sampling, align, sensitivity,...
    dispersionfun, lambda, J_2D, rho, weights,...
    options, varargin...
    )
% BAEK2017ALGORITHM2  Run ADMM as in Algorithm 2 of Baek et al. 2017
%
% ## Syntax
% I = baek2017Algorithm2(...
%     image_sampling, align, sensitivity,...
%     dispersionfun, lambda, J, rho, weights,...
%     options [, verbose]...
% )
% [ I, image_bounds ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb, J_full ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb, J_full, J_est ] = baek2017Algorithm2(___)
%
% ## Description
% I = baek2017Algorithm2(...
%     image_sampling, align, sensitivity,...
%     dispersionfun, lambda, J, rho, weights,...
%     options [, verbose]...
% )
%   Estimate a latent RGB or hyperspectral image `I` from dispersion in
%   the input RAW image `J`.
%
% [ I, image_bounds ] = baek2017Algorithm2(___)
%   Additionally returns the boundaries of the output image in the
%   coordinate system of the input image.
%
% [ I, image_bounds, I_rgb ] = baek2017Algorithm2(___)
%   Additionally returns the RGB equivalent of the latent image.
%
% [ I, image_bounds, I_rgb, J_full ] = baek2017Algorithm2(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped by the model of chromatic aberration.
%
% [ I, image_bounds, I_rgb, J_full, J_est ] = baek2017Algorithm2(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the height and width, respectively, of
%   the output image `I`.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th input colour channel or spectral
%   band of `I`. `sensitivity` is a matrix mapping colours in `I` to
%   colours in `J`.
%
% dispersionfun -- Model of dispersion
%   A function handle, produced by 'makeDispersionfun()'. `dispersionfun(X)`,
%   where `X` is a three-element row vector (x, y, lambda), returns the
%   dispersion vector for the position (x, y) in `J` corresponding to light
%   with wavelength or colour channel index `lambda`. The dispersion vector
%   corrects for lateral chromatic aberration by pointing from the
%   corresponding position in the reference spectral band or colour channel
%   to position (x, y). This function will negate the dispersion vectors
%   produced by `dispersionfun()` in order to create a warp matrix from `I` to
%   `J`.
%
% lambda -- Wavelength bands
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersionfun`. 'c' is the desired number of spectral bands or colour
%   channels in `I`.
%
% J -- Input RAW image
%   A 2D array containing the raw colour-filter pattern data of an image.
%
% rho -- Penalty parameters
%   A two-element vector containing the penalty parameters, `rho_1` and
%   `rho_2` for the constraints on `Z1` and `Z2`, respectively, in the ADMM
%   optimization problem.
%
% weights -- Regularization weights
%   `weights(1)` is the 'alpha' weight on the regularization of the spatial
%   gradient of the image in the ADMM optimization problem. `weights(2)` is
%   the `beta` weight on the regularization of the spectral gradient of the
%   spatial gradient of the image in the ADMM optimization problem.
%
%   If all elements of `weights` are zero, this function will find a
%   solution to the simpler problem:
%     argmin_i (||M * Omega * Phi * i - j||_2) ^ 2
%   where `M` performs mosaicing, `Omega` converts colours to the RGB
%   colour space of the camera, and `Phi` warps the image according to the
%   dispersion model. `j` is the input RAW image. If the linear system is
%   underdetermined, the function will find the minimum-norm least squares
%   solution. If the problem is sufficiently determined, as it may be in
%   cases where `i` has fewer wavelength bands of colour channels than `j`
%   has colour channels, then the function will find the exact or
%   least-squares solution.
%
%   This function will not terminate if `weights` has a mixture of zero and
%   nonzero elements, because the ADMM algorithm will not converge.
%
% options -- Options and small parameters
%   A structure with the following fields:
%   - 'add_border': A Boolean value indicating whether or not the
%     `image_bounds` input argument of 'dispersionfunToMatrix()' should be empty.
%     If `true`, the output image `I` will be large enough to contain the
%     lateral chromatic aberration-corrected coordinates of all pixels in
%     `J`. If `false`, the output image `I` will be clipped to the region
%     occupied by `J`.
%   - 'full_GLambda': A Boolean value used as the `replicate` input
%     argument of 'spectralGradient()' when creating the spectral gradient
%     matrix needed for regularizing `I`. Refer to the documentation of
%     'spectralGradient.m' for details.
%   - 'int_method': The numerical integration method used for spectral to
%     colour space conversion. `int_method` is passed to
%     'channelConversionMatrix()' as its `int_method` argument. Refer to
%     the documentation of 'channelConversionMatrix.m' for details. If
%     'int_method' is 'none', as should be the case when colour conversion
%     is from a set of colour channels, not a spectral space, numerical
%     integration will not be performed.
%   - 'maxit': Maximum number of iterations: The first element of `maxit`
%     contains the maximum number of iterations to use with MATLAB's
%     'pcg()' function when solving the I-minimization step of the ADMM
%     algorithm. The second element of `maxit` contains the maximum number
%     of ADMM iterations to perform.
%   - 'tol': Convergence tolerances: The first element of `tol` is the
%     tolerance value to use with MATLAB's 'pcg()' function when solving
%     the I-minimization step of the ADMM algorithm. The second and third
%     elements of `tol` are the absolute and relative tolerance values for
%     the ADMM algorithm, as explained in Section 3.3.1 of Boyd et al.
%     2011.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show the progress of the
%   iterative optimization.
%
% ## Output Arguments
%
% I -- Latent image
%   An image_sampling(1) x image_sampling(2) x length(lambda) array,
%   storing the latent image estimated using the ADMM algorithm.
%
% image_bounds -- Latent image coordinate frame
%   The boundaries of `I` expressed in the coordinate frame of `J`.
%   `image_bounds` is an output argument of the same name of
%   'dispersionfunToMatrix()'. If `add_border` is `false`, `image_bounds`
%   will be equal to `[0, 0, size(J, 2), size(J, 1)]`.
%
% I_rgb -- Latent RGB image
%   The RGB equivalent of the latent image, generated using the colour
%   space conversion data in `sensitivity`. An image_sampling(1) x
%   image_sampling(2) x 3 array.
%
% J_full -- Warped latent RGB image
%   An RGB image produced by warping `I` according to the dispersion model,
%   followed by conversion to the RGB colour space of `J`. An size(J, 1) x
%   size(J, 2) x 3 array.
%
% J_est -- Re-estimated input RAW image
%   The mosaiced version of `J_full`, representing the result of passing
%   `I` through the forward model of dispersion and image capture. An array
%   with the same dimensions as `J`.
%
% ## References
%
% This function implements Algorithm 2 in the first set of supplemental
% material of the following article:
%
%   Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%     single-shot hyperspectral imaging using a prism." ACM Transactions
%     on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%     doi:10.1145/3130800.3130896
%
% For more information on ADMM (Alternating Direction Method of
% Multipliers), read:
%
%   Boyd, S, et al.. "Distributed Optimization and Statistical Learning via
%     the Alternating Direction Method of Multipliers." Foundations and
%     Trends in Machine Learning, vol. 3, no. 1, pp. 1-122, 2011.
%     doi:10.1561/2200000016
%
% ## Future Work
%
% There are several modifications and expansions which may improve the
% performance of ADMM:
% - Section 3.4.1 of Boyd et al. 2011, "Varying Penalty Parameter"
% - Section 4.3.2 of Boyd et al. 2011, "Early Termination"
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix,
% spatialGradient, spectralGradient, softThreshold, subproblemI, lsqminnorm

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 27, 2018

nargoutchk(1, 5);
narginchk(9, 10);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

if length(image_sampling) ~= 2
    error('The `image_sampling` input argument must contain an image height and width only.');
end

if isStringScalar(options.int_method) || ischar(options.int_method)
    do_integration = ~strcmp(options.int_method, 'none');
else
    error('`options.int_method` must be a character vector or a string scalar.');
end

J = J_2D(:);

% Create constant matrices
image_sampling_J = size(J_2D);
n_bands = length(lambda);
M = mosaicMatrix(image_sampling_J, align);
if do_integration
    Omega = channelConversionMatrix(image_sampling_J, sensitivity, lambda, options.int_method);
else
    Omega = channelConversionMatrix(image_sampling_J, sensitivity);
end
if options.add_border
    image_bounds = [];
else
    image_bounds = [0, 0, image_sampling_J(2), image_sampling_J(1)];
end
[ Phi, image_bounds ] = dispersionfunToMatrix(...
   dispersionfun, lambda, image_sampling_J, image_sampling, image_bounds, true...
);

if all(weights == 0)
    A = M * Omega * Phi;
    if (size(A, 1) < size(A, 2)) || (rank(A) < size(A, 2))
        % Minimum-norm least squares solution to the non-regularized problem
        I = lsqminnorm(A, J);
    else
        % Unique or least-squares solution
        I = A \ J;
    end  
else
    % Perform ADMM
    G_xy = spatialGradient([image_sampling, n_bands]);
    G_lambda = spectralGradient([image_sampling, n_bands], options.full_GLambda);
    G_lambda_sz1 = size(G_lambda, 1);
    G_lambda_sz2 = size(G_lambda, 2);
    % The product `G_lambda * G_xy` must be defined, so `G_lambda` needs to be
    % replicated to operate on both the x and y-gradients.
    G_lambda = [
        G_lambda, sparse(G_lambda_sz1, G_lambda_sz2);
        sparse(G_lambda_sz1, G_lambda_sz2), G_lambda
        ];

    % Initialization
    J_bilinear = bilinearDemosaic(J_2D, align, [false, true, false]); % Initialize with the Green channel
    if any(image_sampling ~= image_sampling_J)
        I = imresize(J_bilinear, image_sampling, 'bilinear');
    else
        I = J_bilinear;
    end
    I = repmat(I(:), n_bands, 1);
    len_I = length(I);
    Z1 = G_xy * I;
    G_lambda_xy = G_lambda * G_xy;
    Z2 = G_lambda_xy * I;
    len_Z1 = length(Z1);
    U1 = zeros(len_Z1, 1);
    len_Z2 = length(Z2);
    U2 = zeros(len_Z2, 1);

    % Iteration
    [ b, A ] = subproblemI(M, Omega, Phi, G_xy, G_lambda, J, Z1, Z2, U1, U2, rho);
    soft_thresholds = weights ./ rho;
    G_xy_T = G_xy.';
    G_lambda_xy_T = G_lambda_xy.';
    converged = false;
    for iter = 1:options.maxit(2)
        % Optimization
        [ I, flag, relres, iter_pcg ] = pcg(...
            A, b, options.tol(1), options.maxit(1), [], [], I...
        );
        if(verbose)
            fprintf('%d:    PCG (flag = %d, relres = %g, iter = %d)\n',...
                iter, flag, relres, iter_pcg...
                );
        end
        g_xy = G_xy * I;
        g_lambda_xy = G_lambda_xy * I;
        Z1_prev = Z1;
        Z2_prev = Z2;
        Z1 = softThreshold(g_xy + U1, soft_thresholds(1));
        Z2 = softThreshold(g_lambda_xy + U2, soft_thresholds(2));
        R1 = g_xy - Z1;
        R2 = g_lambda_xy - Z2;
        U1 = U1 + R1;
        U2 = U2 + R2;

        % Calculate residuals
        R1_norm = norm(R1);
        R2_norm = norm(R2);
        S1 = rho(1) * G_xy_T * (Z1 - Z1_prev);
        S2 = rho(2) * G_lambda_xy_T * (Z2 - Z2_prev);
        S1_norm = norm(S1);
        S2_norm = norm(S2);

        if(verbose)
            fprintf(' Residuals (R1_norm = %g, R2_norm = %g, S1_norm = %g, S2_norm = %g)\n',...
                R1_norm, R2_norm, S1_norm, S2_norm...
            );
        end

        % Calculate stopping criteria
        % See Section 3.3.1 of Boyd et al. 2011.
        epsilon_pri1 = sqrt(len_Z1) * options.tol(2) +...
            options.tol(3) * max([norm(g_xy), norm(Z1)]);
        epsilon_pri2 = sqrt(len_Z2) * options.tol(2) +...
            options.tol(3) * max([norm(g_lambda_xy), norm(Z2)]);

        Y1 = rho(1) * U1;
        Y2 = rho(2) * U2;
        epsilon_dual1 = sqrt(len_I) * options.tol(2) +...
            options.tol(3) * norm(G_xy_T * Y1);
        epsilon_dual2 = sqrt(len_I) * options.tol(2) +...
            options.tol(3) * norm(G_lambda_xy_T * Y2);

        if(verbose)
            fprintf('Stop Crit. (e_p1 = %g, e_p2 = %g, e_d2 = %g, e_d2 = %g)\n',...
                epsilon_pri1, epsilon_pri2, epsilon_dual1, epsilon_dual2...
            );
        end

        % Check against stopping criteria
        if (R1_norm < epsilon_pri1  && R2_norm < epsilon_pri2 &&...
            S1_norm < epsilon_dual1 && S2_norm < epsilon_dual2)
            converged = true;
            break;
        end

        b = subproblemI(M, Omega, Phi, G_xy, G_lambda, J, Z1, Z2, U1, U2, rho);
    end

    if verbose
        if converged
            fprintf('Convergence after %d iterations.\n', iter);
        else
            fprintf('Maximum number of iterations, %d, reached without convergence.\n', iter);
        end
    end
end

I_3D = reshape(I, image_sampling(1), image_sampling(2), n_bands);

if nargout > 2
    if do_integration
        Omega_I = channelConversionMatrix(image_sampling, sensitivity, lambda, options.int_method);
    else
        Omega_I = channelConversionMatrix(image_sampling, sensitivity);
    end
    I_rgb = Omega_I * I;
    n_channels_rgb = 3;
    varargout{1} = reshape(I_rgb, image_sampling(1), image_sampling(2), n_channels_rgb);
    
    if nargout > 3
        J_full = Omega * Phi * I;
        varargout{2} = reshape(J_full, image_sampling_J(1), image_sampling_J(2), n_channels_rgb);
        
        if nargout > 4
            J_est = M * J_full;
            varargout{3} = reshape(J_est, image_sampling_J);
        end
    end
end
end