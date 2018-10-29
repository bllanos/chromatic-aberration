function [limit_freq, limit_rad] = bandlimit(signals, threshold, varargin)
% BANDLIMIT  Find the approximate bandlimits of discrete signals
%
% ## Syntax
% limit_freq = bandlimit(signals, threshold [, verbose])
% [limit_freq, limit_rad] = bandlimit(signals, threshold [, verbose])
%
% ## Description
% limit_freq = bandlimit(signals, threshold [, verbose])
%   Returns the maximum of the bandlimits of the input signals in units of
%   cycles per sampling interval
%
% [limit_freq, limit_rad] = bandlimit(signals, threshold [, verbose])
%   Additionally returns the maximum of the bandlimits as a value in
%   radians per sampling interval
%
% ## Input Arguments
%
% signals -- Time domain signals
%   A 2D array, where `signals(i, :)` is the i-th discrete signal. The
%   signals are assumed to be sampled at evenly-spaced points in a
%   one-dimensional domain.
%
% threshold -- Power spectrum threshold
%   The fraction of the cumulative power spectrum that the user considers
%   to define the bandlimit (because no signal that is time-limited is
%   truly frequency-limited).
%
% verbose -- Debugging flag
%   A Boolean scalar which enables graphical debugging output if `true`.
%   Defaults to `false` if not passed.
%
% ## Output Arguments
%
% limit_freq -- Bandlimit in time units
%   The cumulative power spectrum is computed for each signal. The i-th
%   element of a cumulative power spectrum gives the amount of power in a
%   signal that is represented by frequencies less than or equal to the
%   i-th frequency in the discrete Fourier transform of the signal. The
%   bandlimit for a signal is the frequency at or below which a fraction of
%   at least `threshold` of the total power is contained. `limit_freq` is
%   the maximum of the bandlimits computed for the rows of `signals`.
%
%   `limit_freq` is in units of cycles per interval between adjacent sample
%   points in the discrete representation `signals`.
%
% limit_rad -- Bandlimit in radians
%   `limit_rad` is just the conversion of `limit_freq` to radians by
%   multiplying by `2 * pi` radians per cycle.
%
% ## References
%
% The code is based on the 'FFTOfMatrixRowsExample.mlx' MATLAB example,
% shown on the 'fft' help page.
%
% See also fft

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 26, 2018

narginchk(3, 3);
nargoutchk(1, 2);

verbose = false;
if ~isempty(varargin)
    verbose = varargin{1};
end

if threshold < 0 || threshold > 1
    error('`threshold` must be from 0 to 1.');
end

n = size(signals, 2);
n_pow2 = 2 ^ nextpow2(n);
n_pow2_2 = n_pow2 / 2;
n_pow2_21 = n_pow2_2 + 1;
sampling_freq = 1;
sampling_period = 1 / sampling_freq;
sampling_length = n + 1;
freq_sampling = 0:(sampling_freq/n_pow2):sampling_freq/2;

freq_two_sided = fft(signals, n_pow2, 2);
freq = abs(freq_two_sided / n_pow2);
freq = freq(:, 1:n_pow2_21);
freq(:, 2:(end - 1)) = 2 * freq(:, 2:(end - 1));

% Multiply by 'n_pow2' because of the earlier division by 'n_pow2'
power = (freq .^ 2) * n_pow2;
% Divide by 2 to correct for the earlier multiplication by 2
power(:, 2:(end - 1)) = power(:, 2:(end - 1)) / 2;
power_cum = cumsum(power, 2);
power_cum = power_cum ./ repmat(power_cum(:, end), 1, n_pow2_21);

n_signals = size(signals, 1);
bandlimits_ind = zeros(n_signals, 1);
bandlimits_mask = power_cum >= threshold;
for c = 1:n_signals
    bandlimits_ind(c) = find(bandlimits_mask(c, :), 1);
end

if verbose
    time_sampling = (0:(sampling_length-1)) * sampling_period;
    
    fg = figure;
    for c = 1:n_signals
        subplot(n_signals, 1, c)
        plot(time_sampling, [0, signals(c, :)])
        title(sprintf('Signal %d in the time domain', c));
        xlabel('Time index')
        ylabel('Amplitude')
    end
    
    figure;
    for c=1:n_signals
        subplot(n_signals, 1, c)
        plot(freq_sampling, freq(c, :)...
        )
        title(sprintf('Signal %d in the frequency domain', c));
        xlabel('Frequency [cycles/(time index increment)]')
        ylabel('Amplitude')
    end
    
    figure;
    for c=1:n_signals
        subplot(n_signals, 1, c)
        hold on
        plot(...
            freq_sampling,...
            power(c, :)...
        )
        scatter(freq_sampling(bandlimits_ind(c)), power(c, bandlimits_ind(c)), 'filled')
        hold off
        title(sprintf('Signal %d power spectrum', c));
        xlabel('Frequency [cycles/(time index increment)]')
        ylabel('Power')
        legend('Signal power', 'Bandlimit')
    end
end

max_bandlimit = max(bandlimits_ind);
limit_freq = freq_sampling(max_bandlimit);

limit_rad = 2 * pi * limit_freq;

if verbose
    % Look at the bandlimited versions of the input signals
    freq_bandlimited = freq_two_sided;
    freq_bandlimited((max_bandlimit + 1):(end - max_bandlimit + 1)) = 0;
    signals_bandlimited = ifft(freq_bandlimited, n_pow2, 2, 'symmetric');
    
    figure(fg);
    for c = 1:n_signals
        subplot(n_signals, 1, c)
        hold on
        plot(time_sampling, [0, signals_bandlimited(c, 1:n)], ':')
        hold off
    end
end

end