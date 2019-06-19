function [s, signals_scaled] = scaleSignals(signals, threshold, weights)
% SCALESIGNALS  Find scaling factors to make a set of signals more uniform
% 
% ## Syntax
% s = scaleSignals(signals, threshold [, weights])
% [s, signals_scaled] = scaleSignals(signals, threshold [, weights])
%
% ## Description
% s = scaleSignals(signals, threshold [, weights])
%   Returns a scaling signal which, when multiplied element-wise by the input
%   signals, produces a set of signals that are closer to unity.
%
% [s, signals_scaled] = scaleSignals(signals, threshold [, weights])
%   Additionally returns scaled versions of the signals, with small values
%   thresholded to zero.
%
% ## Input Arguments
%
% signals -- Set of signals
%   A matrix where the rows index signals, and the columns index points in the
%   domain of the signals.
%
% threshold -- Small value threshold
%   A fraction indicating what proportion of the peak magnitude of a signal the
%   user considers to be effectively zero. Points in the domain where all
%   signals have values below their respective thresholds will be given values
%   of `1` in `s`. Points in the domain where a signal has a value below its
%   threshold will be given values of `0` in `signals_scaled`.
%
% weights -- Signal weights
%   A vector with a length equal to the size of `signals` in its first
%   dimension, containing the relative weights to give to the signals when
%   computing the scaling factors.
%
% ## Output Arguments
%
% s -- Smoothing signal
%   A row vector with the same length as the size of `signals` in its second
%   dimension. `repmat(s, size(signals, 1), 1) .* signals` is close to
%   `ones(size(signals))` in a least-squares sense.
%
% signals_scaled -- Scaled, thresholded signals
%   A matrix equal to `repmat(s, size(signals, 1), 1) .* signals`, followed by
%   replacing values where `signals` was below `threshold` with zeros.
%
% See also smoothMetamer

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 19, 2019

nargoutchk(1, 2);
narginchk(2, 3);

n_signals = size(signals, 1);
if nargin < 3
    weights = ones(n_signals, 1);
else
    if any(weights <= 0)
        error('Elements of `weights` must be nonzero, positive values.');
    elseif length(weights) ~= n_signals
        error('The length of `weights` must equal the number of rows in `signals`.');
    end
    weights = reshape(weights, [], 1);
end

% Filter the signals to the points at which at least some are nonzero
n_points = size(signals, 2);
signals_relative = abs(signals) ./ repmat(max(abs(signals), [], 2), 1, n_points);
filter_neg = (signals_relative < threshold);
filter = ~all(filter_neg, 1);
signals_filtered = signals(:, filter);

signals_filtered_weighted = signals_filtered .* ...
    repmat(weights, 1, size(signals_filtered, 2));
reference = ones(size(signals_filtered));

s_filtered = dot(signals_filtered_weighted, reference, 1) ./ ...
    dot(signals_filtered_weighted, signals_filtered, 1);

s = ones(1, n_points);
s(filter) = s_filtered;

if nargout > 1
    signals_scaled = signals .* repmat(s, n_signals, 1);
    signals_scaled(filter_neg) = 0;
end

end