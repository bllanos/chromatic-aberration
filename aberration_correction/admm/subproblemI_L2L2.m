function [ b, A ] = subproblemI_L2L2(M, Omega, Phi, G_xy, G_lambda, J, alpha, beta)
% SUBPROBLEMIL2L2  I-minimization equation for L2 priors on both gradients
%
% ## Syntax
% b = subproblemI_L2L2(M, Omega, Phi, G_xy, G_lambda, J, alpha, beta)
% [ b, A ] = subproblemI_L2L2(M, Omega, Phi, G_xy, G_lambda, J, alpha, beta)
%
% ## Description
% b = subproblemI_L2L2(M, Omega, Phi, G_xy, G_lambda, J, alpha, beta)
%   Returns an updated left-hand side vector.
%
% [ b, A ] = subproblemI_L2L2(M, Omega, Phi, G_xy, G_lambda, J, alpha, beta)
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
% alpha -- Spatial gradient regularization weight
%   The weight on the regularization of the spatial gradient of the image.
%
% beta -- Spectral gradient regularization weight
%   The weight on the regularization of the spectral gradient of the
%   spatial gradient of the image.
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
% squares, as opposed to by iterative optimization.
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix,
% spatialGradient, spectralGradient, subproblemI

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 18, 2018

nargoutchk(1, 2);
narginchk(8, 8);

M_Omega_Phi = M * Omega * Phi;
G_lambda_xy = G_lambda * G_xy;

b = M_Omega_Phi.' * J;

if nargout > 1
    A = (M_Omega_Phi.' * M_Omega_Phi) +...
        alpha * (G_xy.' * G_xy) +...
        beta * (G_lambda_xy.' * G_lambda_xy);
else
    A = [];
end

end