function f = subproblemI_L2L2(M_Omega_Phi, G_xy, G_lambda_xy, J, weights)
% SUBPROBLEMIL2L2  I-minimization equation for L2 priors on both gradients
%
% ## Syntax
% This function generates a function that stores partial computations:
% f = subproblemI_L2L2(M_Omega_Phi, G_xy, G_lambda_xy, J, weights);
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
% Z -- Slack variable
%   A zero to three-element cell vector:
%   - The first cell is unused.
%   - The second cell is unused.
%   - The third cell contains the slack variable for the non-negativity
%     constraint on `I`, which is an addition to the method of Baek et al.
%     2017. `Z{3}` is an (n_px_I x c) x 1 vector.
%
%  If `Z` does not have three cells, the non-negativity constraint on `I`
%  will not be applied.
%
% U -- Scaled dual variable
%   The scaled Lagrange multipliers for the constraints on the
%   corresponding cells of `Z` in the ADMM optimization problem. `U{k}` has
%   the same dimensions as `Z{k}`. `U` must be a cell vector of the same
%   length as `Z`.
%
% rho -- Penalty parameters
%   A zero to three-element vector. The first and second elements are
%   ignored. The third element, needed if `Z` has three cells, contains the
%   `rho_3` penalty parameter for the non-negativity constraint.
%
% weights -- Regularization weights
%   `weights(1)` is the 'alpha' weight on the regularization of the spatial
%   gradient of the image. `weights(2)` is the 'beta' weight on the
%   regularization of the spectral gradient of the spatial gradient of the
%   image.
%
% ## Output Arguments
%
% b -- Right-hand side vector
%   The right-hand side vector in the equation `A * I = b` for the estimate
%   of `I`.
%
% A -- Right-hand side vector
%   The left-hand side matrix in the equation `A * I = b` for the estimate
%   of `I`.
%
% ## Notes
% - Requesting `A` as an output argument will trigger more computation than
%   only requesting `b`.
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
% with the variation that there are L2-norm priors on the spatial gradient,
% and on the spectral gradient of the spatial gradient of `I` instead of
% L1-norm priors. As such, the solution can be determined by linear least
% squares, as opposed to by iterative optimization. When `Z` has three
% cells, there is a non-negativity constraint on `I`, and so iterative
% optimization is still required.
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix,
% spatialGradient, spectralGradient, subproblemI

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 18, 2018

nargoutchk(1, 2);
narginchk(5, 5);

nonneg = (length(Z) > 2);

b_const = M_Omega_Phi.' * J;
A_const = (M_Omega_Phi.' * M_Omega_Phi) +...
        weights(1) * (G_xy.' * G_xy) +...
        weights(2) * (G_lambda_xy.' * G_lambda_xy);
if nonneg
    I_A = speye(size(A));
end

    function [ b, A ] = inner(Z, U, rho)
        if nonneg
            b = b_const + (rho(3) / 2) * (Z{3} - U{3});
        else
            b = b_const;
        end
        
        if nargout > 1
            if nonneg
                A = A_const + (rho(3) / 2) * I_A;
            else
                A = A_const;
            end
        else
            A = [];
        end
    end

f = @inner;

end