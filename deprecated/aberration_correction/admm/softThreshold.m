function [ y ] = softThreshold(x, t)
% SOFTTHRESHOLD  Soft-thresholding operator
%
% ## Syntax
% y = softThreshold(x, t)
%
% ## Description
% y = softThreshold(x, t)
%   Returns the soft-thresholded version of `x` obtained using threshold
%   `t`.
%
% ## Input Arguments
%
% x -- Input array
%   The array to be soft-thresholded.
%
% t -- Soft threshold
%   A scalar used to perform soft-thresholding of `x`.
%
% ## Output Arguments
%
% y -- Output array
%   An array with the same dimensions as `x`, where `y(i)` is:
%   - `x(i) - t`, if `x(i) > t`
%   - `0`, if `abs(x(i)) <= t`
%   - `x(i) + t`, if `x(i) < -t`
%
% ## References
%
% Soft thresholding is the proximity operator of the L1 norm, as discussed
% in section 4.4.3 of:
%   Boyd, S, et al.. "Distributed Optimization and Statistical Learning via
%     the Alternating Direction Method of Multipliers." Foundations and
%     Trends in Machine Learning, vol. 3, no. 1, pp. 1-122, 2011.
%     doi:10.1561/2200000016

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 27, 2018

nargoutchk(1, 1);
narginchk(2, 2);

if ~isscalar(t) || t < 0
    error('The threshold `t` must be a non-negative scalar.');
end

y = zeros(size(x));
y(x > t) = x(x > t) - t;
y(x < -t) = x(x < -t) + t;
end

