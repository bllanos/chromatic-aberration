function in = penalties( J_2D, I, J_est_2D, G, norms, in )
% PENALTIES  Calculate data fitting and regularization errors
%
% ## Syntax
% in = penalties( J_2D, I, J_est_2D, G, norms, in )
%
% ## Description
% in = penalties( J_2D, I, J_est_2D, G, norms, in )
%   Returns an updated structure containing penalty values.
%
% ## Input Arguments
%
% J_2D -- Input RAW image
%   A 2D array containing the raw colour-filter pattern data of an image.
%
% I -- Vectorized latent image
%   An vector, storing the vectorized form of the estimated latent image
%   corresponding to `J_2D`.
%
% J_est_2D -- Estimated RAW image
%   A 2D array containing the estimated version of `J_2D` created from `I`.
%
% G -- Regularization operators
%   A cell vector, where `G{i}` is a matrix that can multiply `I` to
%   produce a vector whose norm is a regularization penalty on the solution
%   `I`.
%
% norms -- Penalty norm types
%   A logical vector the same length as `G`. If `norms(i)` is `true`, then
%   the penalty corresponding to `G{i}` is the L1 norm of the result of
%   applying `G{i}` to `I`. Otherwise, the penalty is evaluated by taking
%   the L2 norm.
%
% in -- Preallocated intermediate data and results
%   A structure with the following fields:
%   - 'err': A vector containing the penalty values corresponding to the
%     regularization operators in `G` as its second through to last
%     elements. All penalty values are normalized by the lengths of the
%     vectors that they are the norms of. The first element of 'err' is the
%     data fitting error between the input image and the estimated input
%     image.
%   - 'err_vectors': A cell vector containing the results of applying the
%     regularization operators in `G` to `I`.
%
% ## Output Arguments
%
% in -- Updated intermediate data and results
%   An updated version of the `in` input argument.
%
% ## Notes
% - This function is intended to be used for in-place computation, by
%   having the caller assign the output argument `in` to the same variable
%   input as the input argument `in`.
% - In contrast to 'baek2017Algorithm2()', 'penalties()' does not remove a
%   border region from the image prior to calculating penalty values, for
%   efficiency.
%
% See also baek2017Algorithm2LowMemory, baek2017Algorithm2,
% selectWeightsGrid, solvePatchesADMM

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 7, 2018

nargoutchk(1, 1);
narginchk(6, 6);

in.err(1) = immse( J_2D, J_est_2D );

for w = 1:length(G)
    in.err_vectors{w} = G{w} * I;

    if norms(aw)
        in.err(w + 1) = mean(abs(in.err_vectors{aw}));
    else
        in.err(w + 1) = dot(in.err_vectors{aw}, in.err_vectors{aw}) / length(in.err_vectors{aw});
    end 
end

end