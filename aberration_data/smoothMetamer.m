function I = smoothMetamer(c, color_weights, lambda, tol, verbose)
% SMOOTHMETAMER  Find a smooth, non-negative metamer
% 
% ## Syntax
% I = smoothMetamer(c, color_weights, lambda, [tol, verbose])
%
% ## Description
% I = smoothMetamer(c, color_weights, lambda, [tol, verbose])
%   Returns a non-negative spectral signal that projects to the input colour,
%   and which has the low L2 norm of its spectral first-order derivative.
%
% ## Input Arguments
%
% c -- Colour
%   A vector containing the colour to which the spectral signal must project.
%
% color_weights -- Colour conversion matrix
%   A matrix of dimensions `length(c)` by n which converts spectral information
%   to colour. `color_weights` must account for any numerical intergration that
%   is part of colour conversion.
%
% lambda -- Spectral regularization weight
%   A scalar containing the weight on the norm of the spectral gradient in the
%   following optimization problem:
%
%     find `I` minimizing
%       ||color_weights * I - c||_2^2 + lambda * ||spectralGradient(I)||_2^2
%     such that `I > 0`
%
% tol -- Convergence tolerance
%   An optional 'TolX' convergence tolerance to pass to 'lsqnonneg()'. If `tol`
%   is zero, or is not passed, 'lsqnonneg()' will use its default convergence
%   tolerance.
%
% verbose -- Verbosity flag
%   If `true`, console output will be produced by 'lsqnonneg()' to indicate
%   whether or not the optimization converged.
%
% ## Output Arguments
%
% I -- Metamer
%   A column vector of length `size(color_weights, 2)` containing the spectral
%   signal satisfying the optimization problem described above.
%
% See also lsqnonneg 

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 17, 2019

nargoutchk(1, 1);
narginchk(3, 5);

if nargin < 4
    tol = 0;
end
if nargin > 4 && verbose && tol
    options = optimset('Display', 'final', 'TolX', tol);
elseif nargin > 4 && verbose
    options = optimset('Display', 'final');
elseif tol
    options = optimset('TolX', tol);
else
    options = optimset('Display', 'none');
end

% Set up the optimization problem
G_lambda = full(spectralGradient([1, 1, size(color_weights, 2)], false));
A = [color_weights; sqrt(lambda) * G_lambda];
b = [reshape(c, [], 1); zeros(size(G_lambda, 1), 1)];

I = lsqnonneg(A, b, options);

end