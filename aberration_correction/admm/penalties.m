function in = penalties( J, I, M_Omega_Phi, G, norms, in )
% PENALTIES  Calculate data fitting and regularization errors
%
% ## Syntax
% in = penalties( J, I, M_Omega_Phi, G, norms, in )
%
% ## Description
% in = penalties( J, I, M_Omega_Phi, G, norms, in )
%   Returns an updated structure containing penalty values.
%
% ## Input Arguments
%
% J -- Input RAW image
%   A vector containing the vectorized raw colour-filter pattern data of an
%   image.
%
% I -- Vectorized latent image
%   An vector, storing the vectorized form of the estimated latent image
%   corresponding to `J`.
%
% M_Omega_Phi -- Reprojection matrix
%   A matrix that can multiply `I` to produce the version of `J`
%   corresponding to `I`.
%
% G -- Regularization operators
%   A cell vector, where `G{i}` is a matrix that can multiply `I` to
%   produce a vector whose norm is a regularization penalty on the solution
%   `I`. Cells can also contain empty values (`[]`), indicating that the
%   corresponding regularization penalties are disabled (so no penalty
%   values will be computed for them).
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
%   - 'J_est': A vector containing the estimated version of `J` created
%     from `I`.
%
%   `in` can be initialized by 'initPenalties()'.
%
% ## Output Arguments
%
% in -- Updated intermediate data and results
%   An updated version of the `in` input argument.
%
% ## Notes
% - This function is intended to be used for in-place computation, by
%   having the caller assign the output argument `in` to the same variable
%   input as the input argument `in`. `in` can be created by
%   'initPenalties()'.
% - In contrast to 'baek2017Algorithm2()', 'penalties()' does not remove a
%   border region from the image prior to calculating penalty values, for
%   efficiency.
%
% See also initPenalties, baek2017Algorithm2, weightsLowMemory,
% solvePatchesColor

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 7, 2018

nargoutchk(1, 1);
narginchk(6, 6);

in.J_est = M_Omega_Phi * I;
in.err(1) = immse(J, in.J_est);

for w = 1:length(G)
    if ~isempty(G{w})
        in.err_vectors{w} = G{w} * I;

        if norms(w)
            in.err(w + 1) = mean(abs(in.err_vectors{w}));
        else
            in.err(w + 1) = dot(in.err_vectors{w}, in.err_vectors{w}) / length(in.err_vectors{w});
        end
    end
end

end