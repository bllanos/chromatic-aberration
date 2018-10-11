function out = initPenalties(M_Omega_Phi, G)
% INITPENALTIES  Allocate memory for 'penalties()'
%
% ## Syntax
% out = initPenalties(M_Omega_Phi, G)
%
% ## Description
% out = initPenalties(M_Omega_Phi, G)
%   Returns a structure containing arrays to be used by 'penalties()'
%
% ## Input Arguments
%
% M_Omega_Phi -- Reprojection matrix
%   Refer to the documentation of penalties.m.
%
% G -- Regularization operators
%   Refer to the documentation of penalties.m.
%
% ## Output Arguments
%
% out -- Preallocated arrays and intermediate data
%   The `in` input/output argument of 'penalties()'. Refer to the
%   documentation of penalties.m.
%
% See also penalties

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 10, 2018

nargoutchk(1, 1);
narginchk(2, 2);

n_priors = length(G);
out.err = zeros(1, n_priors + 1);
out.err_vectors = cell(n_priors, 1);
for w = 1:n_priors
    if ~isempty(G{w})
        out.err_vectors{w} = zeros(size(G{w}, 1), 1);
    end
end

out.J_est = zeros(size(M_Omega_Phi, 1), 1);

end