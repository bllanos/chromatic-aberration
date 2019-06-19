function s = scaleSignals(signals)
% SCALESIGNALS  Find scaling factors to make a set of signals more uniform
% 
% ## Syntax
% s = scaleSignals(signals)
%
% ## Description
% s = scaleSignals(signals)
%   Returns a scaling signal which, when multiplied element-wise by the input
%   signals, produces a set of signals that are closer to unity.
%
% ## Input Arguments
%
% signals -- Set of signals
%   A matrix where the rows index signals, and the columns index points in the
%   domain of the signals.
%
% ## Output Arguments
%
% s -- Smoothing signal
%   A row vector with the same length as the size of `signals` in its second
%   dimension. `repmat(s, size(signals, 1), 1) .* signals` is close to
%   `ones(size(signals))` in a least-squares sense.
%
% See also smoothMetamer

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 19, 2019

nargoutchk(1, 1);
narginchk(1, 1);

% Filter the signals to the points at which at least some are nonzero
filter = any(signals, 1);
signals_filtered = signals(:, filter);
reference = ones(size(signals_filtered));

s_filtered = dot(signals_filtered, reference, 1) ./ ...
    dot(signals_filtered, signals_filtered, 1);

s = ones(1, size(signals, 2));
s(filter) = s_filtered;

end