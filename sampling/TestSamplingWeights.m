%% Test script for the 'samplingWeights()' and 'bandlimit()'

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 2, 2018

verbose = false;

%% Find the bandlimit of cosine waves
% The code is based on the 'FFTOfMatrixRowsExample.mlx' MATLAB example,
% shown on the 'fft' help page.

Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector

x1 = cos(2*pi*50*t);          % First row wave
x2 = cos(2*pi*150*t);         % Second row wave
x3 = cos(2*pi*300*t);         % Third row wave

X = [x1; x2; x3];

power_threshold = 1;
[limit_freq, limit_rad] = bandlimit(X, power_threshold, verbose);
fprintf('With a power threshold of %g, the bandlimit is %g cycles/sampling interval\n', power_threshold, limit_freq);
if limit_freq * 2 * pi ~= limit_rad
    error('`limit_freq` and `limit_rad` are inconsistent');
end

power_threshold = 0.9;
[limit_freq, limit_rad] = bandlimit(X, power_threshold, verbose);
fprintf('With a power threshold of %g, the bandlimit is %g cycles/sampling interval\n', power_threshold, limit_freq);
if limit_freq * 2 * pi ~= limit_rad
    error('`limit_freq` and `limit_rad` are inconsistent');
end

%% Test 'spectralWeights()'

options.int_method = 'trap';
options.power_threshold = 0.99;
options.n_bands = 0;
options.support_threshold = 0.05;
options.bands_padding = 1000;
options.interpolant = @normpdf;

% Create colour channel sensitivities
n_color_bands = 50;
color_bands = linspace(-0.5, 1, n_color_bands);
color_band_spacing = mean(diff(color_bands));
color_band_indices = 0:(n_color_bands - 1);
color_map = [
    ((color_band_indices - n_color_bands / 3) .^ 2 / (n_color_bands ^ 2));
    cos(0.3 * (color_band_indices - 2 * n_color_bands / 3));
    sinc(0.2 * (color_band_indices - n_color_bands / 2))
    ];
n_channels = size(color_map, 1);
color_bands_padded = [
    (color_bands(1) - (options.bands_padding * color_band_spacing)):color_band_spacing:(color_bands(1) - color_band_spacing),...
    color_bands,...
    (color_bands(end) + color_band_spacing):color_band_spacing:(color_bands(end) + (options.bands_padding * color_band_spacing))
];
color_map_padded = [
    repmat(color_map(:, 1), 1, options.bands_padding),...
    color_map,...
    repmat(color_map(:, end), 1, options.bands_padding)...
];

% Create reference spectral signals
n_spectral_bands = [200, 80]; % These sampling rates are higher than those of the colour channels
spectral_bands = {
    linspace(0, 1, n_spectral_bands(1));
    linspace(0.5, 1.5, n_spectral_bands(2))
    };
spectral_signals = {
    sin(spectral_bands{1});
    exp(spectral_bands{2})
};
spectral_bands_spacing = [
    mean(diff(spectral_bands{1}));
    mean(diff(spectral_bands{2}))
    ];

% Determine expected colours, by expressing the colour channels and
% spectral signals in an even higher shared sampling space.
n_signals = length(spectral_bands);
global_band_spacing = min(color_band_spacing, min(spectral_bands_spacing)) / 10;
colors = cell(n_signals, 1);
for i = 1:n_signals
    global_bands = min(spectral_bands{i}(1), color_bands(1)):global_band_spacing:max(spectral_bands{i}(end), color_bands(end));
    colors{i} = zeros(n_channels, 1);
    for c = 1:n_channels
        color_map_resampled = interp1(color_bands, color_map(c, :), global_bands, 'spline', 0);
        color_map_resampled(global_bands < color_bands(1)) = color_map(c, 1);
        color_map_resampled(global_bands > color_bands(end)) = color_map(c, end);
        spectral_signal_resampled = interp1(spectral_bands{i}, spectral_signals{i}, global_bands, 'spline', 0);
        colors{i}(c) = trapz(global_bands, color_map_resampled .* spectral_signal_resampled);
    end
end

% Compare single reference sampling call syntax with multiple reference
% sampling call syntax
[...
  color_weights1, spectral_weights1, bands1, color_weights_reference1...
] = samplingWeights(...
  color_map, color_bands, spectral_bands{1}, options, false...
);
[...
  color_weights, spectral_weights, bands, color_weights_reference...
] = samplingWeights(...
  color_map, color_bands, spectral_bands, options, verbose...
);
if ~all(all(color_weights1 == color_weights))
    error('`color_weights` does not match.');
end
if ~all(all(spectral_weights1 == spectral_weights{1}))
    error('`spectral_weights` does not match.');
end
if ~all(all(bands1 == bands))
    error('`bands` does not match.');
end
if ~all(all(color_weights_reference1 == color_weights_reference{1}))
    error('`color_weights_reference` does not match.');
end

% Construct signals on the output bands
trial_signals = {
    sin(bands);
    exp(bands)
};
band_spacing = mean(diff(bands));
bands_padded = [
    (bands(1) - (options.bands_padding * band_spacing)):band_spacing:(bands(1) - band_spacing),...
    bands,...
    (bands(end) + band_spacing):band_spacing:(bands(end) + (options.bands_padding * band_spacing))
];
trial_signals_padded = {
    [
        repmat(trial_signals{1}(1), 1, options.bands_padding),...
        trial_signals{1},...
        repmat(trial_signals{1}(end), 1, options.bands_padding)...
    ];
    [
        repmat(trial_signals{2}(1), 1, options.bands_padding),...
        trial_signals{2},...
        repmat(trial_signals{2}(end), 1, options.bands_padding)...
    ];
};

% Determine expected colours
trial_colors = cell(n_signals, 1);
for i = 1:n_signals
    global_bands = min(bands_padded(1), color_bands(1)):global_band_spacing:max(bands_padded(end), color_bands(end));
    trial_colors{i} = zeros(n_channels, 1);
    for c = 1:n_channels
        color_map_resampled = interp1(color_bands, color_map(c, :), global_bands, 'spline', 0);
        trial_signal_resampled = interp1(bands_padded, trial_signals_padded{i}, global_bands, 'spline', 0);
        trial_colors{i}(c) = trapz(global_bands, color_map_resampled .* trial_signal_resampled);
    end
end

% Evaluate results
for i = 1:n_signals
   fprintf('Colour of trial signal %d, expected:\n', i);
   disp(trial_colors{i});
   disp('actual:');
   disp(color_weights * reshape(trial_signals{i}, [], 1));
end

for i = 1:n_signals
   fprintf('Colour of reference signal %d, expected:\n', i);
   disp(colors{i});
   disp('actual:');
   disp(color_weights_reference{i} * reshape(spectral_signals{i}, [], 1));
end

for i = 1:n_signals
    figure;
    hold on
    plot(spectral_bands{i}, spectral_signals{i}, 'k');
    plot(bands, trial_signals{i}, 'or');
    plot(spectral_bands{i}, spectral_weights{i} * reshape(trial_signals{i}, [], 1), '-g');
    hold off
    legend('True signal', 'Sampled signal', 'Reconstructed signal');
end

%% Repeat,
% but testing the case where the spectral signals have lower sampling frequencies than the colour channels

% Create colour channel sensitivities
n_color_bands = 1000;
color_bands = linspace(-0.5, 1, n_color_bands);
color_band_spacing = mean(diff(color_bands));
color_band_indices = 0:(n_color_bands - 1);
color_map = [
    ((color_band_indices - n_color_bands / 3) .^ 2 / (n_color_bands ^ 2));
    cos(0.3 * (color_band_indices - 2 * n_color_bands / 3));
    sinc(0.2 * (color_band_indices - n_color_bands / 2))
    ];

options.int_method = 'trap';
options.power_threshold = 0.0;
options.n_bands = 0;
options.support_threshold = 0.05;
options.bands_padding = 1000;

% Create reference spectral signals
n_spectral_bands = [200, 80];
spectral_bands = {
    linspace(0, 1, n_spectral_bands(1));
    linspace(0.5, 1.5, n_spectral_bands(2))
    };
spectral_signals = {
    sin(spectral_bands{1});
    exp(spectral_bands{2})
};
spectral_bands_spacing = [
    mean(diff(spectral_bands{1}));
    mean(diff(spectral_bands{2}))
    ];
spectral_bands_padded = {
    [
        (spectral_bands{1}(1) - (options.bands_padding * spectral_bands_spacing(1))):spectral_bands_spacing(1):(spectral_bands{1}(1) - spectral_bands_spacing(1)),...
        spectral_bands{1},...
        (spectral_bands{1}(end) + spectral_bands_spacing(1)):spectral_bands_spacing(1):(spectral_bands{1}(end) + (options.bands_padding * spectral_bands_spacing(1)))
    ];
    [
        (spectral_bands{2}(1) - (options.bands_padding * spectral_bands_spacing(2))):spectral_bands_spacing(2):(spectral_bands{2}(1) - spectral_bands_spacing(2)),...
        spectral_bands{2},...
        (spectral_bands{2}(end) + spectral_bands_spacing(2)):spectral_bands_spacing(2):(spectral_bands{2}(end) + (options.bands_padding * spectral_bands_spacing(2)))
    ]
};
spectral_signals_padded = {
    [
        repmat(spectral_signals{1}(1), 1, options.bands_padding),...
        spectral_signals{1},...
        repmat(spectral_signals{1}(end), 1, options.bands_padding)...
    ];
    [
        repmat(spectral_signals{2}(1), 1, options.bands_padding),...
        spectral_signals{2},...
        repmat(spectral_signals{2}(end), 1, options.bands_padding)...
    ];
};

% Determine expected colours, by expressing the colour channels and
% spectral signals in an even higher shared sampling space
n_signals = length(spectral_bands);
global_band_spacing = min(color_band_spacing, min(spectral_bands_spacing)) / 10;
colors = cell(n_signals, 1);
for i = 1:n_signals
    global_bands = min(spectral_bands_padded{i}(1), color_bands(1)):global_band_spacing:max(spectral_bands_padded{i}(end), color_bands(end));
    colors{i} = zeros(n_channels, 1);
    for c = 1:n_channels
        color_map_resampled = interp1(color_bands, color_map(c, :), global_bands, 'spline', 0);
        spectral_signal_resampled = interp1(spectral_bands_padded{i}, spectral_signals_padded{i}, global_bands, 'spline', 0);
        spectral_signal_resampled(global_bands < spectral_bands_padded{i}(1)) = spectral_signals_padded{i}(1);
        spectral_signal_resampled(global_bands > spectral_bands_padded{i}(end)) = spectral_signals_padded{i}(end);
        colors{i}(c) = trapz(global_bands, color_map_resampled .* spectral_signal_resampled);
    end
end

[...
  ~, ~, ~, color_weights_reference...
] = samplingWeights(...
  color_map, color_bands, spectral_bands, options, verbose...
);

% Evaluate results
for i = 1:n_signals
   fprintf('Colour of reference signal %d, expected:\n', i);
   disp(colors{i});
   disp('actual:');
   disp(color_weights_reference{i} * reshape(spectral_signals{i}, [], 1));
end

%% Repeat,
% testing the case where the desired sampling has a higher sampling rate
% than the colour channels

% Create colour channel sensitivities
n_color_bands = 100;
color_bands = linspace(-0.5, 1.5, n_color_bands);
color_band_spacing = mean(diff(color_bands));
color_band_indices = 0:(n_color_bands - 1);
color_map = [
    ((color_band_indices - n_color_bands / 3) .^ 2 / (n_color_bands ^ 2));
    cos(0.3 * (color_band_indices - 2 * n_color_bands / 3));
    sinc(0.2 * (color_band_indices - n_color_bands / 2))
    ];
n_channels = size(color_map, 1);

options.int_method = 'trap';
options.power_threshold = 0.99;
options.n_bands = 300;
options.support_threshold = 0.05;
options.bands_padding = 1000;

spectral_bands = {};

[color_weights, ~, bands] = samplingWeights(...
  color_map, color_bands, spectral_bands, options, verbose...
);

% Construct signals on the output bands
trial_signals = {
    sin(bands);
    exp(bands)
};
band_spacing = mean(diff(bands));

% Determine expected colours
trial_colors = cell(n_signals, 1);
for i = 1:n_signals
    global_bands = min(bands_padded(1), color_bands(1)):global_band_spacing:max(bands_padded(end), color_bands(end));
    trial_colors{i} = zeros(n_channels, 1);
    for c = 1:n_channels
        color_map_resampled = interp1(color_bands, color_map(c, :), global_bands, 'spline', 0);
        color_map_resampled(global_bands < color_bands(1)) = color_map(c, 1);
        color_map_resampled(global_bands > color_bands(end)) = color_map(c, end);
        trial_signal_resampled = interp1(bands, trial_signals{i}, global_bands, 'spline', 0);
        trial_colors{i}(c) = trapz(global_bands, color_map_resampled .* trial_signal_resampled);
    end
end

% Evaluate results
for i = 1:n_signals
   fprintf('Colour of trial signal %d, expected:\n', i);
   disp(trial_colors{i});
   disp('actual:');
   disp(color_weights * reshape(trial_signals{i}, [], 1));
end

%% Test that different options give the desired effects

options.power_threshold = 1;
options.support_threshold = 0;
options.n_bands = 0;

[...
  ~, ~, bands_power_threshold_1, ~...
] = samplingWeights(...
  color_map, color_bands, spectral_bands, options, false...
);

if length(bands_power_threshold_1) ~= n_color_bands || any(bands_power_threshold_1 ~= color_bands)
    error('Incorrect bands for `options.power_threshold = 1`');
end

options.power_threshold = 0;
options.n_bands = 5;

[...
  ~, ~, bands_fixed, ~...
] = samplingWeights(...
  color_map, color_bands, spectral_bands, options, false...
);

if length(bands_fixed) ~= options.n_bands
    error('Incorrect number of bands for `options.n_bands = %g`', options.n_bands);
end
fprintf('Maximum different between sampling points and their correct values for %g bands:\n', options.n_bands);
disp(max(abs(bands_fixed - linspace(color_bands(1), color_bands(end), options.n_bands))));