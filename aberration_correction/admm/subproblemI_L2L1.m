function [ b, A ] = subproblemI_L2L1(M, Omega, Phi, G_xy, G_lambda, J, Z2, U2, rho_2, alpha)
% SUBPROBLEMIL2L1  I-minimization equation for an L2 prior on the spatial gradient
%
% ## Syntax
% b = subproblemI_L2L1(M, Omega, Phi, G_xy, G_lambda, J, Z2, U2, rho_2, alpha)
% [ b, A ] = subproblemI_L2L1(M, Omega, Phi, G_xy, G_lambda, J, Z2, U2, rho_2, alpha)
%
% ## Description
% b = subproblemI_L2L1(M, Omega, Phi, G_xy, G_lambda, J, Z2, U2, rho_2, alpha)
%   Returns an updated left-hand side vector.
%
% [ b, A ] = subproblemI_L2L1(M, Omega, Phi, G_xy, G_lambda, J, Z2, U2, rho_2, alpha)
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
% M -- Mosaicing matrix
%   The output argument of 'mosaicMatrix()'. An (n_px_J)-by-(n_px_J x 3)
%   sparse array.
%
% Omega -- Colour channel conversion matrix
%   The output argument of 'channelConversionMatrix()'. An (n_px_J x
%   3)-by-(n_px_J x c) sparse array which converts the latent image to the
%   RGB colour space of the camera.
%
% Phi -- Dispersion model warp matrix
%   The `W` output argument of 'dispersionfunToMatrix()'. An (n_px_J x
%   c)-by-(n_px_I x c) sparse array.
%
% G_xy -- Spatial gradient matrix
%   The output argument of 'spatialGradient()'. An (n_px_I x c x
%   2)-by-(n_px_I x c) sparse array.
%
% G_lambda -- Spectral gradient matrix
%   The output argument of 'spectralGradient()'. An (n_px_I x
%   c_gLambda)-by-(n_px_I x c) sparse array.
%
% J -- Vectorized input RAW image
%   A vectorized form of the input RAW image, where pixels are ordered
%   first by column, then by row. `J` is an n_px_J x 1 vector.
%
% Z2 -- Spectral gradient slack variable
%   The slack variable equal to the spectral gradient of `G_xy * I` in the
%   ADMM optimization problem. `Z2` is an (n_px_I x c_gLambda x 2) x 1
%   vector.
%
% U2 -- Spectral gradient scaled dual variable
%   The scaled Lagrange multiplier for the constraint on `Z2` in the ADMM
%   optimization problem. `U2` is an (n_px_I x c_gLambda x 2) x 1 vector.
%
% rho_2 -- Penalty parameter
%   A scalar containing the penalty parameter, `rho_2` for the constraint
%   on `Z2` in the ADMM optimization problem.
%
% alpha -- Spatial gradient regularization weight
%   The weight on the regularization of the spatial gradient of the image
%   in the ADMM optimization problem.
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
% - `A` is a function of the penalty parameter, `rho_2`, and the constant
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
% with the variation that there is an L2-norm prior on the spatial gradient
% of `I` instead of an L1-norm prior.
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
narginchk(10, 10);

M_Omega_Phi = M * Omega * Phi;
G_lambda_xy = G_lambda * G_xy;

% Equation 4.1 of Boyd et al. 2011
b = (M_Omega_Phi.' * J) +...
    (rho_2 / 2) * G_lambda_xy.' * (Z2 - U2);

if nargout > 1
    A = (M_Omega_Phi.' * M_Omega_Phi) + alpha * (G_xy.' * G_xy) +...
        (rho_2 / 2) * (G_lambda_xy.' * G_lambda_xy);
else
    A = [];
end

end