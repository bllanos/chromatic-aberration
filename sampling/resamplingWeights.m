function weights = upsamplingWeights(dst, src, f, padding)
% UPSAMPLINGWEIGHTS Create an upsampling matrix operator
%
% ## Syntax
% weights = upsamplingWeights(dst, src, f, padding)
%
% ## Description
% weights = upsamplingWeights(dst, src, f, padding)
%   Returns a matrix mapping from the sampling space of `src` to `dst`.
%
% ## Input Arguments
%
% dst -- Destination sampling locations
%   A vector containing the locations at which the upsampled signal will be
%   sampled. The values in `dst` are expected to be evenly-spaced.
%
% src -- Original sampling locations
%   A vector containing the locations at which the signal is currently
%   sampled. The values in `src` are expected to be evenly-spaced, and
%   their spacing must be at least as large as that of the values in `dst`.
%
% f -- Interpolation function
%   A function handle used to perform the upsampling. `f(x)` must return
%   the weight of a sample at location `x`, given that the location of the
%   new sample is at the origin (`x = 0`). Changes of one unit in `x`
%   represent shifts by one sampling position in `src`. `f(x)` need not be
%   normalized over the integers (i.e. the integral of `f(x)` from `x =
%   -Inf` to `x = Inf`, where `x` is an integer, need not equal `1`),
%   because this function will renormalize interpolation weights.
%
%   `f` should operate element-wise on an array input.
%
% padding -- Extrapolation conditions
%   A scalar indicating the number of virtual samples outside the domain of
%   `src` which will contribute to the upsampling. `src` will be padded
%   with `padding` samples on each side, taking the values of the signal at
%   the corresponding endpoint of `src`. In other words, if `padding` is
%   a large number, the signal will be extrapolated as though it is equal
%   to its endpoints outside of its domain. If `padding` is zero, the
%   signal will effectively be padded with zeros.
%
% ## Output Arguments
%
% weights -- Upsampling matrix
%   A matrix of dimensions `length(dst)` by `length(src)` defining a
%   mapping from the sampling space of `src` to the sampling space of
%   `dst`. Given a column vector of the same length as `src`, `v`,
%   representing a signal sampled at the locations in `src`, its upsampled
%   version, sampled at the locations in `dst` is `u = weights * v`. `u(i)`
%   is computed by convolution: "u(i) = sum_j (f(src(j) - dst(i) / (src(2)
%   - src(1))) * v(j))", where `src` and `v` have been padded as described
%   in the documentation of `padding` above. In the case where `src` and
%   `dst` are the same, this function ignores `f` and assumes that
%   `weights` should be the identity mapping.
%
% ## Notes
% - In the special case where `src` has only one element, `weights` will be
%   a column vector of ones with the length of `dst`.
% - In the special case where `src` and `dst` are identical, `weights` will
%   be an identity matrix.
%
% ## References
% - The code is loosely based on the "Ideal Bandlimited Interpolation"
%   example in the MATLAB Help page for 'sinc()'.
%
% See also bandlimit, samplingWeights, sinc, interp1

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 5, 2018

% Parse input arguments
narginchk(4, 4);
nargoutchk(1, 1);

if length(src) == 1
    weights = ones(length(dst), 1);
    return;
end

if padding < 0 || padding ~= round(padding)
    error('`padding` must be a nonnegative integer.');
end

if length(dst) > 1
    diff_dst = diff(dst);
    dst_spacing = diff_dst(1);
    if max(abs(diff_dst - dst_spacing)) > 1e-6
        error('`dst` must contain equally-spaced values.')
    end
else
    dst_spacing = Inf;
end

diff_src = diff(src);
src_spacing = diff_src(1);
if max(abs(diff_src - src_spacing)) > 1e-6
    error('`src` must contain equally-spaced values.')
end

if src_spacing < dst_spacing && ~(length(dst) == length(src) && all(dst == src))
    error(['The spacing of the elements of `dst` is greater than that o'
           'f `src`, so `weights` cannot be computed without aliasing']);
end

src_padded = [
    (src(1) - (padding * src_spacing)):src_spacing:(src(1) - src_spacing),...
    reshape(src, 1, []),...
    (src(end) + src_spacing):src_spacing:(src(end) + (padding * src_spacing))
];
[dst_grid, src_grid] = ndgrid(dst, src_padded);
weights = f((dst_grid - src_grid) / src_spacing);
% Adjust the weights of the endpoints so that upsampling assumes the
% value of the signal outside of its domain is equal to its value at
% the nearest endpoint of its domain
weights = [
    sum(weights(:, 1:(padding + 1)), 2),...
    weights(:, (padding + 2):(end - padding - 1)),...
    sum(weights(:, (end - padding):end), 2)
];

% Renormalize
sum_weights = sum(weights, 2);
sum_weights(sum_weights == 0) = 1;
weights = weights ./ repmat(sum_weights, 1, length(src));

end
