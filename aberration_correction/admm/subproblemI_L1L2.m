function f = subproblemI_L1L2(M_Omega_Phi, G_xy, G_lambda_xy, J, beta, nonneg)
% SUBPROBLEMIL1L2  I-minimization equation for an L2 prior on the spectral gradient
%
% ## Syntax
% This function generates a function that stores partial computations:
% f = subproblemI_L1L2(M_Omega_Phi, G_xy, G_lambda_xy, J, beta, nonneg);
% b = f(Z, U, rho)
% [ b, A ] = f(Z, U, rho)
%
% ## Description
% b = f(Z, U, rho)
%   Returns an updated left-hand side vector.
%
% [ b, A ] = f(Z, U, rho)
%   Additionally returns an updated right-hand side matrix.
%
% ## Input Arguments
%
% The array dimensions mentioned in this section are defined as follows:
% - `n_px_J` is the number of pixels in the vectorized RAW input image,
%   `J`. `n_px_J = length(J)`.
% - `c` is the number of colour channels or spectral bands in the latent
%   image to be estimated, `I`.
% - `n_px_I` is the number of pixels in the vectorized latent image, `I`.
% - `c_gLambda` is the number of channels in the spectral gradient, as
%   determined by the `replicate` input argument of 'spectralGradient()'.
%
% M_Omega_Phi -- Reprojection matrix
%   The product of three sparse matrices, in the following order:
%   - Mosaicing matrix: M, the output argument of 'mosaicMatrix()'. An
%     (n_px_J)-by-(n_px_J x 3) sparse array.
%   - Colour channel conversion matrix: Omega, the output argument of
%     'channelConversionMatrix()'. An (n_px_J x 3)-by-(n_px_J x c) sparse
%     array which converts the latent image to the RGB colour space of the
%     camera.
%   - Dispersion model warp matrix: Phi, the `W` output argument of
%     'dispersionfunToMatrix()'. An (n_px_J x c)-by-(n_px_I x c) sparse
%     array.
%
% G_xy -- Spatial gradient matrix
%   The output argument of 'spatialGradient()'. An (n_px_I x c x
%   2)-by-(n_px_I x c) sparse array.
%
% G_lambda_xy -- Spectral gradient of the spatial gradient matrix
%   The product of the output argument of 'spectralGradient()', an (n_px_I
%   x c_gLambda)-by-(n_px_I x c) sparse array, with each half of `G_xy`.
%
% J -- Vectorized input RAW image
%   A vectorized form of the input RAW image, where pixels are ordered
%   first by column, then by row. `J` is an n_px_J x 1 vector.
%
% beta -- Spectral gradient regularization weight
%   The weight on the regularization of the spectral gradient of the
%   spatial gradient of the image in the ADMM optimization problem.
%
% nonneg -- Non-negativity constraint
%   A logical scalar indicating whether or not there is a non-negativity
%   constraint on the solution.
%
% Z -- Slack variables
%   A one to three-element cell vector:
%   - The first cell contains the slack variable equal to `G_xy * I` in the
%     ADMM optimization problem. `Z{1}` is an (n_px_I x c x 2) x 1 vector.
%   - The second cell is unused.
%   - The third cell contains the slack variable for the non-negativity
%     constraint on `I`, which is an addition to the method of Baek et al.
%     2017. `Z{3}` is an (n_px_I x c) x 1 vector.
%
%  If `nonneg` is `true`, then `Z` must have three cells.
%
% U -- Scaled dual variables
%   The scaled Lagrange multipliers for the constraints on the
%   corresponding cells of `Z` in the ADMM optimization problem. `U{k}` has
%   the same dimensions as `Z{k}`. `U` must be a cell vector of the same
%   length as `Z`.
%
% rho -- Penalty parameters
%   A one to three-element vector. The first element contains the penalty
%   parameter, `rho_1`, for the constraints on `Z{1}` in the ADMM
%   optimization problem of Baek et al. 2017. The second element is
%   ignored. The third element, needed if `Z` has three cells, contains the
%   `rho_3` penalty parameter for the non-negativity constraint.
%
% ## Output Arguments
%
% b -- Right-hand side vector
%   The right-hand side vector in the equation `A * I = b` for the next
%   estimate of `I` in the I-minimization step of the ADMM optimization
%   problem.
%
% A -- Right-hand side vector
%   The left-hand side matrix in the equation `A * I = b` for the next
%   estimate of `I` in the I-minimization step of the ADMM optimization
%   problem.
%
% ## Notes
% - Requesting `A` as an output argument will trigger more computation than
%   only requesting `b`.
% - `A` is a function of the penalty parameters, `rho`, and the constant
%   matrices `M`, `Omega`, `Phi`, `G_xy`, and `G_lambda`. Therefore, `A`
%   does not need to be recalculated each iteration of ADMM unless the
%   penalty parameter is variable.
% - `b` is a function of both the constant matrices and the intermediate
%   variables in the ADMM optimization problem, and therefore must be
%   recalculated each iteration of ADMM.
%
% ## References
%
% This function calculates the left-hand matrix and constant vector which
% produce the next estimate of `I` in line 2 of Algorithm 2 in the first
% set of supplemental material for the following article:
%
%   Baek, S.-H., Kim, I., Gutierrez, D., & Kim, M. H. (2017). "Compact
%     single-shot hyperspectral imaging using a prism." ACM Transactions
%     on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 217:1–12.
%     doi:10.1145/3130800.3130896
%
% with the variation that there is an L2-norm prior on the spectral
% gradient of the spatial gradient of `I` instead of an L1-norm prior, and
% with the optional variation of a non-negativity constraint on `I`.
%
% For more information on ADMM (Alternating Direction Method of
% Multipliers), read:
%   Boyd, S, et al.. "Distributed Optimization and Statistical Learning via
%     the Alternating Direction Method of Multipliers." Foundations and
%     Trends in Machine Learning, vol. 3, no. 1, pp. 1-122, 2011.
%     doi:10.1561/2200000016
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix,
% spatialGradient, spectralGradient, subproblemI

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 18, 2018

nargoutchk(1, 2);
narginchk(6, 6);

M_Omega_Phi_J = M_Omega_Phi.' * J;
A_const = (M_Omega_Phi.' * M_Omega_Phi) + beta * (G_lambda_xy.' * G_lambda_xy);
G_xy_T = G_xy.';
G_xy2 = G_xy.' * G_xy;
if nonneg
    I_A = speye(size(A_const));
end

    function [ b, A ] = inner(Z, U, rho)
        % Equation 4.1 of Boyd et al. 2011
        b = M_Omega_Phi_J +...
            (rho(1) / 2) * G_xy_T * (Z{1} - U{1});
        if nonneg
            b = b + (rho(3) / 2) * (Z{3} - U{3});
        end
        
        if nargout > 1
            A = A_const + (rho(1) / 2) * G_xy2;
            if nonneg
                A = A + (rho(3) / 2) * I_A;
            end
        else
            A = [];
        end
    end

f = @inner;

end