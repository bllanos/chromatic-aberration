function f = subproblemI(M_Omega_Phi, G_xy, G_lambda_xy, J)
% SUBPROBLEMI  Set up the matrix equation for the I-minimization step of ADMM
%
% ## Syntax
% This function generates a function that stores partial computations:
% f = subproblemI(M_Omega_Phi, G_xy, G_lambda_xy, J);
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
% Z -- Slack variables
%   A two or three-element cell vector:
%   - The first cell contains the slack variable equal to `G_xy * I` in the
%     ADMM optimization problem. `Z{1}` is an (n_px_I x c x 2) x 1 vector.
%   - The second cell contains he slack variable equal to the spectral
%     gradient of `G_xy * I` in the ADMM optimization problem. `Z{2}` is an
%     (n_px_I x c_gLambda x 2) x 1 vector.
%   - The third cell contains the slack variable for the non-negativity
%     constraint on `I`, which is an addition to the method of Baek et al.
%     2017. `Z{3}` is an (n_px_I x c) x 1 vector.
%
%  If `Z` does not have three cells, the non-negativity constraint on `I`
%  will not be applied.
%
% U -- Scaled dual variables
%   The scaled Lagrange multipliers for the constraints on the
%   corresponding cells of `Z` in the ADMM optimization problem. `U{k}` has
%   the same dimensions as `Z{k}`. `U` must be a cell vector of the same
%   length as `Z`.
%
% rho -- Penalty parameters
%   A two or three-element vector. The first two elements contain the
%   penalty parameters, `rho_1` and `rho_2` for the constraints on `Z{1}`
%   and `Z{2}`, respectively, in the ADMM optimization problem of Baek et
%   al. 2017. The third element, needed if `Z` has three cells, contains
%   the `rho_3` penalty parameter for the non-negativity constraint.
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
%   penalty parameters are variable.
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
% spatialGradient, spectralGradient

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 27, 2018

nargoutchk(1, 2);
narginchk(4, 4);

nonneg = (length(Z) > 2);

M_Omega_Phi_J = M_Omega_Phi.' * J;
A_const = (M_Omega_Phi.' * M_Omega_Phi);
G_xy_T = G_xy.';
G_xy2 = G_xy.' * G_xy;
G_lambda_xy_T = G_lambda_xy.';
G_lambda_xy2 = G_lambda_xy.' * G_lambda_xy;
if nonneg
    I_A = speye(size(A));
end

    function [ b, A ] = inner(Z, U, rho)
        b = M_Omega_Phi_J +...
            (rho(1) / 2) * G_xy_T * (Z{1} - U{1}) +...
            (rho(2) / 2) * G_lambda_xy_T * (Z{2} - U{2});
        if nonneg
            b = b + (rho(3) / 2) * (Z{3} - U{3});
        end
        
        if nargout > 1
            A = A_const +...
                (rho(1) / 2) * G_xy2 +...
                (rho(2) / 2) * G_lambda_xy2;
            if nonneg
                A = A + (rho(3) / 2) * I_A;
            end
        else
            A = [];
        end
    end

f = @inner;

end