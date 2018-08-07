function [ I_3D, image_bounds, varargout ] = baek2017Algorithm2(...
    image_sampling, align, dispersion, sensitivity, lambda, J_2D,...
    rho, weights, options, varargin...
    )
% BAEK2017ALGORITHM2  Run ADMM as in Algorithm 2 of Baek et al. 2017
%
% ## Syntax
% I = baek2017Algorithm2(...
%   image_sampling, align, dispersion, sensitivity, lambda, J,...
%   rho, weights, options [, verbose]...
% )
% [ I, image_bounds ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb, J_full ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb, J_full, J_est ] = baek2017Algorithm2(___)
% [ I, image_bounds, I_rgb, J_full, J_est, I_warped ] = baek2017Algorithm2(___)
%
% ## Description
% I = baek2017Algorithm2(...
%     image_sampling, align, dispersion, sensitivity,...
%     lambda, J, rho, weights, options [, verbose]...
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
%   image, warped according to the model of chromatic aberration.
%
% [ I, image_bounds, I_rgb, J_full, J_est ] = baek2017Algorithm2(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% [ I, image_bounds, I_rgb, J_full, J_est, I_warped ] = baek2017Algorithm2(___)
%   Additionally returns the version of the latent image warped according
%   to the model of chromatic aberration.
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
% dispersion -- Model of dispersion
%   `dispersion` can be empty (`[]`), if there is no model of dispersion.
%   Otherwise, two forms of this argument can be passed:
%
%   `dispersion` can be a function handle, produced by
%   'makeDispersionfun()'. `dispersion(X)`, where `X` is a three-element
%   row vector (x, y, lambda), returns the dispersion vector for the
%   position (x, y) in `J` corresponding to light with wavelength or colour
%   channel index `lambda`. The dispersion vector corrects for lateral
%   chromatic aberration by pointing from the corresponding position in the
%   reference spectral band or colour channel to position (x, y). This
%   function will negate the dispersion vectors of `dispersion` in order to
%   create a warp matrix from `I` to `J`.
%
%   `dispersion` can be a matrix for warping `I` to `J`, where the k-th row
%   contains the weights of pixels in `I` used to re-estimate the k-th
%   pixel in `J`. With a matrix value for `dispersion`, this function will
%   not make use of `lambda` or `options.add_border`, which are only needed
%   for computing a warp matrix.
%
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th input colour channel or spectral
%   band of `I`. `sensitivity` is a matrix mapping colours in `I` to
%   colours in `J`.
%
% lambda -- Wavelength bands
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersion`. 'c' is the desired number of spectral bands or colour
%   channels in `I`.
%
% J -- Input RAW image
%   A 2D array containing the raw colour-filter pattern data of an image.
%
% rho -- Penalty parameters
%   A two or three-element vector containing the penalty parameters,
%   `rho_1` and `rho_2` for the constraints on `Z1` and `Z2`, respectively,
%   in the ADMM optimization problem. The third element is a `rho_3`
%   penalty parameter for a non-negativity constraint on the solution, and
%   is only required if `options.nonneg` is `true`.
%
% weights -- Regularization weights
%   `weights(1)` is the 'alpha' weight on the regularization of the spatial
%   gradient of the image in the ADMM optimization problem. `weights(2)` is
%   the 'beta' weight on the regularization of the spectral gradient of the
%   spatial gradient of the image in the ADMM optimization problem.
%
%   If all elements of `weights` are zero, and `options.nonneg` is `false`,
%   this function will find a solution to the simpler problem:
%     argmin_i (||M * Omega * Phi * i - j||_2) ^ 2
%   where `M` performs mosaicing, `Omega` converts colours to the RGB
%   colour space of the camera, and `Phi` warps the image according to the
%   dispersion model. `j` is the input RAW image. If the linear system is
%   underdetermined, the function will find the minimum-norm least squares
%   solution. If the problem is sufficiently determined, as it may be in
%   cases where `i` has fewer wavelength bands of colour channels than `j`
%   has colour channels, then the function will find the iterative
%   approximation to the exact solution (using MATLAB's 'pcg()' function),
%   or will find a least-squares solution.
%
%   The ADMM iterations will not converge if an element of `weights` is
%   zero, but the corresponding element in `options.norms` is `true`, as
%   ADMM will be trying to optimize a constraint that has no influence.
%
% options -- Options and small parameters
%   A structure with the following fields:
%   - 'add_border': A Boolean value indicating whether or not the
%     `image_bounds` input argument of 'dispersionfunToMatrix()' should be
%     empty. If `true`, the output image `I` will be large enough to
%     contain the lateral chromatic aberration-corrected coordinates of all
%     pixels in `J`. If `false`, the output image `I` will be clipped to
%     the region occupied by `J`. If `dispersion` is empty, or is a matrix,
%     this option is not used, and can be absent.
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
%   - 'norms': Algorithm variants: A two-element logical vector. The first
%     element specifies whether to use an L1 (`true`) or L2 (`false`) norm
%     prior on the latent image's spatial gradient. The second element
%     specifies whether to use an L1 (`true`) or L2 (`false`) norm prior on
%     the spectral gradient of the latent image's spatial gradient. If both
%     elements are `true`, the behaviour of the algorithm matches that
%     described by Baek et al. 2017 (Algorithm 2 in the first set of
%     supplemental material). If only one element is `true`, the ADMM
%     iterations are simplified by eliminating one slack variable. If both
%     elements are `false`, a least-squares solution is obtained
%     (iteratively, using MATLAB's 'pcg()' function).
%   - 'nonneg': A Boolean scalar specifying whether or not to enable a
%     non-negativity constraint on the output latent image. If `true`,
%     `rho` must have three elements.
%   - 'tol': Convergence tolerances: The first element of `tol` is the
%     tolerance value to use with MATLAB's 'pcg()' function, such as when
%     solving the I-minimization step of the ADMM algorithm. The second and
%     third elements of `tol` are the absolute and relative tolerance
%     values for the ADMM algorithm, as explained in Section 3.3.1 of Boyd
%     et al. 2011.
%   - 'varying_penalty_params': If empty, the penalty parameters passed as
%     `rho` will be fixed for all iterations. Otherwise,
%     'varying_penalty_params' is a three-element vector containing the
%     parameters 'tau_incr', 'tau_decr', and 'mu', respectively, in
%     equation 3.13 of Boyd et al. 2011. In this case, the penalty
%     parameters in the ADMM iterations will vary so as to help speed up
%     convergence. Refer to Section 3.4.1 of Boyd et al. 2011 for a full
%     explanation.
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
%   'dispersionfunToMatrix()'. If `options.add_border` is `false`,
%   `image_bounds` will be equal to `[0, 0, size(J, 2), size(J, 1)]`.
%
%   `image_bounds` is empty if `dispersion` is empty, or is a matrix, as
%   'dispersionfunToMatrix()' is not called.
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
% I_warped -- Warped latent image
%   An size(J, 1) x size(J, 2) x length(lambda) array, storing the latent
%   image warped according to the dispersion model.
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
% Depending on the options passed, this function also implements variants
% of the algorithm: L2 priors instead of L1 priors, and a non-negativity
% constraint. I implemented the non-negativity constraint by adding an
% extra term to the ADMM x-minimization step, and an additional
% z-minimization and dual update step. This is different from the
% constrained optimization examples in Boyd et. al. 2011, sections 4.2.5
% and 5.2, but I think it matches the format given at the start of Chapter
% 5.
%
% A non-negativity constraint was used in (among other works):
%
%   Park, J.-I., Lee, M.-H., Grossberg, M. D., & Nayar, S. K. (2007).
%     "Multispectral Imaging Using Multiplexed Illumination." In 2007 IEEE
%     International Conference on Computer Vision (ICCV).
%     doi:10.1109/ICCV.2007.4409090
%
% The initialization method is based on Equation 6 of:
%
%   Sun, T., Peng, Y., & Heidrich, W. (2017). "Revisiting cross-channel
%     information transfer for chromatic aberration correction." In 2017
%     IEEE International Conference on Computer Vision (ICCV) (pp.
%     3268–3276). doi:10.1109/ICCV.2017.352
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
%   - Now activated by a non-empty `options.varying_penalty_params` vector.
%   - Further refinement may be possible by "taking into account the relative
%    magnitudes of [the primal and dual convergence thresholds]."
% - Section 4.3.2 of Boyd et al. 2011, "Early Termination"
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix,
% spatialGradient, spectralGradient, softThreshold, subproblemI,
% subproblemI_L1L2, subproblemI_L2L1, subproblemI_L2L2, lsqminnorm

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 27, 2018

    % Equation 3.13 in Section 3.4.1 of Boyd et al. 2011
    function rho = updatePenalty(R, S, rho)
        if R > mu * S
            rho = rho * tau_incr;
        elseif S > mu * R
            rho = rho / tau_decr;
        end
    end

nargoutchk(1, 6);
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

if any(rho <= 0)
    error('The penalty parameters, `rho`, must be positive numbers.');
end
if options.nonneg && length(rho) < 3
    error('A third penalty parameter must be provided in `rho` when `options.nonneg` is `true`.');
end

vary_penalty_parameters = false;
if ~isempty(options.varying_penalty_params)
    vary_penalty_parameters = true;
    tau_incr = options.varying_penalty_params(1);
    if tau_incr <= 1
        error('The `tau_incr` parameter, `options.varying_penalty_params(1)`, must be greater than one.');
    end
    tau_decr = options.varying_penalty_params(2);
    if tau_decr <= 1
        error('The `tau_decr` parameter, `options.varying_penalty_params(2)`, must be greater than one.');
    end
    mu = options.varying_penalty_params(3);
    if mu <= 1
        error('The `mu` parameter, `options.varying_penalty_params(3)`, must be greater than one.');
    end
end

n_bands = length(lambda);

has_dispersion = ~isempty(dispersion);
if has_dispersion
    dispersion_is_matrix = false;
    if isfloat(dispersion) && ismatrix(dispersion)
        dispersion_is_matrix = true;
        if size(dispersion, 1) ~= (numel(J_2D) * n_bands)
            error('The `dispersion` matrix must have as many rows as there are pixels in `J` times bands.');
        elseif size(dispersion, 2) ~= (prod(image_sampling) * n_bands)
            error('The `dispersion` matrix must have as many columns as there are values in `I`.');
        end
    elseif ~isa(dispersion, 'function_handle')
        error('`dispersion` must be either a floating-point matrix, or a function handle.');
    end
end

J = J_2D(:);

% Create constant matrices
image_sampling_J = size(J_2D);
M = mosaicMatrix(image_sampling_J, align);
if do_integration
    Omega = channelConversionMatrix(image_sampling_J, sensitivity, lambda, options.int_method);
else
    Omega = channelConversionMatrix(image_sampling_J, sensitivity);
end
image_bounds = [];
if has_dispersion
    if dispersion_is_matrix
        Phi = dispersion;
    else
        if ~options.add_border
            image_bounds = [0, 0, image_sampling_J(2), image_sampling_J(1)];
        end
        [ Phi, image_bounds ] = dispersionfunToMatrix(...
           dispersion, lambda, image_sampling_J, image_sampling, image_bounds, true...
        );    
    end
    M_Omega_Phi = M * Omega * Phi;
else
    M_Omega_Phi = M * Omega;
end
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
G_lambda_xy = G_lambda * G_xy;

% Adjust the weights so that they have the same relative importance
% regardless of the differences in the lengths of the vectors whose norms
% are being weighted.
weights(1) = weights(1) * size(M_Omega_Phi, 1) / size(G_xy, 1);
weights(2) = weights(2) * size(M_Omega_Phi, 1) / size(G_lambda_xy, 1);

% Initialization
J_bilinear = bilinearDemosaic(J_2D, align, [false, true, false]); % Initialize with the Green channel
if any(image_sampling ~= image_sampling_J)
    I = imresize(J_bilinear, image_sampling, 'bilinear');
else
    I = J_bilinear;
end
I = repmat(I(:), n_bands, 1);
len_I = length(I);
% Scale to match the mean intensity
J_mean = mean(J);
J_est = M_Omega_Phi * I;
J_est_mean = mean(J_est);
I = I * (J_mean / J_est_mean);

% Select the appropriate algorithm variant
if all(weights == 0) && ~options.nonneg
    A = M_Omega_Phi;
    if (size(A, 1) < size(A, 2)) || (rank(A) < size(A, 2))
        % Minimum-norm least squares solution to the non-regularized problem
        if(verbose)
            fprintf('Computing the minimum-norm least squares solution...\n');
        end
        I = lsqminnorm(A, J);
        if(verbose)
            fprintf('\t...done.\n');
        end
    else
        % Unique or least-squares solution
        if(verbose)
            fprintf(...
                'Computing a non-regularized least squares solution with tolerance %g for up to %d iterations...\n',...
                options.tol(1), options.maxit(1)...
            );
        end
        [ I, flag, relres, iter_pcg ] = pcg(...
            A, J, options.tol(1), options.maxit(1), [], [], I...
        );
        if(verbose)
            fprintf('\tLeast-squares result: PCG (flag = %d, relres = %g, iter = %d)\n',...
                flag, relres, iter_pcg...
            );
        end
    end
elseif all(~options.norms) && ~options.nonneg
    if(verbose)
        fprintf(...
            'Computing a regularized least squares solution with tolerance %g for up to %d iterations...\n',...
            options.tol(1), options.maxit(1)...
        );
    end
    f_Ab = subproblemI_L2L2(M_Omega_Phi, G_xy, G_lambda_xy, J, weights, options.nonneg);
    [ b, A ] = f_Ab({}, {}, []);
    [ I, flag, relres, iter_pcg ] = pcg(...
        A, b, options.tol(1), options.maxit(1), [], [], I...
    );
    if(verbose)
        fprintf('\tLeast-squares result: PCG (flag = %d, relres = %g, iter = %d)\n',...
            flag, relres, iter_pcg...
        );
    end
else
    % Perform ADMM
    if(verbose)
        fprintf('Computing an iterative solution using ADMM...\n');
    end
    
    G{1} = G_xy;
    G{2} = G_lambda_xy;
    
    active_constraints = [options.norms, options.nonneg];
    n_priors = 2;
    n_Z = find(active_constraints, 1, 'last');
    nonneg_ind = 3;
    
    % Initialization
    Z = cell(n_Z, 1);
    len_Z = zeros(n_Z, 1);
    U = cell(n_Z, 1);
    G_T = cell(n_priors, 1);
    g = cell(n_Z, 1);
    Z_prev = cell(n_Z, 1);
    R = cell(n_Z, 1);
    R_norm = zeros(n_Z, 1);
    S = cell(n_Z, 1);
    S_norm = zeros(n_Z, 1);
    epsilon_pri = zeros(n_Z, 1);
    Y = cell(n_Z, 1);
    epsilon_dual = zeros(n_Z, 1);
    
    for z_ind = 1:n_Z
        if active_constraints(z_ind)
            if z_ind == nonneg_ind
                Z{z_ind} =  I;
            else
                Z{z_ind} = G{z_ind} * I;
                G_T{z_ind} = G{z_ind}.';
            end
            len_Z(z_ind) = length(Z{z_ind});
            U{z_ind} = zeros(len_Z(z_ind), 1);
        end
    end

    % Iteration
    if all(options.norms)
        f_Ab = subproblemI(M_Omega_Phi, G_xy, G_lambda_xy, J, options.nonneg);
    elseif options.norms(1)
        f_Ab = subproblemI_L1L2(M_Omega_Phi, G_xy, G_lambda_xy, J, weights(2), options.nonneg);
    elseif options.norms(2)
        f_Ab = subproblemI_L2L1(M_Omega_Phi, G_xy, G_lambda_xy, J, weights(1), options.nonneg);
    else
        f_Ab = subproblemI_L2L2(M_Omega_Phi, G_xy, G_lambda_xy, J, weights, options.nonneg);
    end
    [ b, A ] = f_Ab(Z, U, rho);
    soft_thresholds = weights ./ rho(1:n_priors);

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
        
        converged = true;
        
        for z_ind = 1:n_Z
            if ~active_constraints(z_ind)
                continue;
            end
            if z_ind == nonneg_ind
                g{z_ind} = I;
            else
                g{z_ind} = G{z_ind} * I;
            end
            Z_prev{z_ind} = Z{z_ind};
            if z_ind == nonneg_ind
                % See Section 5.2 of Boyd et al. 2011.
                Z{z_ind} = g{z_ind} + U{z_ind};
                Z{z_ind}(Z{z_ind} < 0) = 0;
            else
                % See Section 6.3 of Boyd et al. 2011.
                Z{z_ind} = softThreshold(...
                    g{z_ind} + U{z_ind}, soft_thresholds(z_ind)...
                );
            end
            % See Section 3.1.1 of Boyd et al. 2011.
            R{z_ind} = g{z_ind} - Z{z_ind};
            U{z_ind} = U{z_ind} + R{z_ind};
            
            % Calculate residuals
            R_norm(z_ind) = norm(R{z_ind});
            if z_ind == nonneg_ind
                S{z_ind} = rho(z_ind) * (Z{z_ind} - Z_prev{z_ind});
            else
                S{z_ind} = rho(z_ind) * G_T{z_ind} * (Z{z_ind} - Z_prev{z_ind});
            end
            S_norm(z_ind) = norm(S{z_ind});
            
            % Calculate stopping criteria
            % See Section 3.3.1 of Boyd et al. 2011.
            epsilon_pri(z_ind) = sqrt(len_Z(z_ind)) * options.tol(2) +...
                options.tol(3) * max([norm(g{z_ind}), norm(Z{z_ind})]);
            Y{z_ind} = rho(z_ind) * U{z_ind};
            if z_ind == nonneg_ind
                epsilon_dual(z_ind) = sqrt(len_I) * options.tol(2) +...
                    options.tol(3) * norm(Y{z_ind});
            else
                epsilon_dual(z_ind) = sqrt(len_I) * options.tol(2) +...
                    options.tol(3) * norm(G_T{z_ind} * Y{z_ind});
            end
            converged = converged &&...
                (R_norm(z_ind) < epsilon_pri(z_ind) && S_norm(z_ind) < epsilon_dual(z_ind));
            
            if(verbose)
                fprintf('%d:    Residuals %d (R1_norm = %g, S1_norm = %g)\n',...
                    iter, z_ind, R_norm(z_ind), S_norm(z_ind)...
                    );
                fprintf('%d:    Stop Crit. %d (e_p1 = %g, e_d1 = %g)\n',...
                    iter, z_ind, epsilon_pri(z_ind), epsilon_dual(z_ind)...
                );
            end
        end

        % Check against stopping criteria
        if converged
            break;
        end

        if vary_penalty_parameters
            for z_ind = 1:n_Z
                if ~active_constraints(z_ind)
                    continue;
                end
                rho(z_ind) = updatePenalty(...
                    R_norm(z_ind), S_norm(z_ind), rho(z_ind)...
                );
                U{z_ind} = Y{z_ind} ./ rho(z_ind);
            end
            soft_thresholds = weights ./ rho(1:n_priors);
            [ b, A ] = f_Ab(Z, U, rho);
        else
            b = f_Ab(Z, U, rho);
        end
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
        if has_dispersion
            I_warped = Phi * I;
        else
            I_warped = I;
        end
        J_full = Omega * I_warped;
        varargout{2} = reshape(J_full, image_sampling_J(1), image_sampling_J(2), n_channels_rgb);
        
        if nargout > 4
            J_est = M * J_full;
            varargout{3} = reshape(J_est, image_sampling_J);
            
            if nargout > 5
                varargout{4} = reshape(I_warped, image_sampling_J(1), image_sampling_J(2), n_bands);
            end
        end
    end
end
end