function [...
    color_weights, spectral_weights, bands, color_weights_reference...
] = samplingWeights(color_map, color_bands, spectral_bands, options, varargin)
% SAMPLINGWEIGHTS Find an optimal sampled representation for spectral data
%
% ## About
%
% The purpose of this function is to find the most memory-efficient
% representation for spectral data that is converted to colour, according
% to the colour channel sensitivities in `color_map`. We seek the smallest
% number of samples for the spectral data that still allows us to represent
% changes in the spectral data that would affect its colour. After finding
% such a representation, we then need to specify how we can map data in
% that representation to colours, or to the representation of some
% reference spectral data.
%
% ## Syntax
% color_weights = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
% [color_weights, spectral_weights] = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
% [color_weights, spectral_weights, bands] = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
% [...
%   color_weights, spectral_weights, bands, color_weights_reference...
% ] = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%
% ## Description
% color_weights = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Returns a linear transformation for mapping the sampled representation
%   of spectral data to colour channels.
%
% [color_weights, spectral_weights] = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Additionally returns a linear transformation for mapping the sampled
%   representation of spectral data to the space of comparison spectral
%   data.
%
% [color_weights, spectral_weights, bands] = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Additionally returns the wavelengths at which the sampled
%   representation is defined.
%
% [...
%   color_weights, spectral_weights, bands, color_weights_reference...
% ] = samplingWeights(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Additionally returns a linear transformation for mapping the sampled
%   representation of spectral data corresponding to `spectral_bands` to
%   colour channels.
%
% ## Input Arguments
%
% color_map -- Colour channel spectral sensitivities
%   A 2D array, where `color_map(i, j)` is the sensitivity of the i-th
%   colour channel to the j-th spectral band in `color_bands`.
%
% color_bands -- Wavelength bands for colour channel sensitivities
%   A vector, of length equal to the size of the second dimension of
%   `color_map`, containing the wavelengths at which the sensitivity
%   functions in `color_map` have been sampled. `color_bands(j)` is the
%   wavelength corresponding to `color_map(:, j)`. The values in
%   `color_bands` are expected to be evenly-spaced.
%
% spectral_bands -- Wavelength bands for spectral reference data
%   A vector containing the wavelengths at which data, to which the
%   spectral data is to be compared to, has been sampled. The values in
%   `spectral_bands` are expected to be evenly-spaced.
%
% options -- Options for finding a spectral representation
%   `options` is a structure with the following fields. The fields are
%   described further in the section 'Algorithm' below:
%   - 'int_method': The numerical integration method to use when
%     integrating over the responses of colour channels to compute colour
%     values. `int_method` is passed to `integrationWeights()` as its
%     `method` input argument.
%   - 'power_threshold': In the context of hyperspectral reconstruction, a
%     colour image produced by a camera described by `color_map` can be
%     seen as a convolution of the camera's spectral sensitivity functions
%     with the spectral signals of each pixel. Therefore, reasoning about
%     hyperspectral reconstruction from the perspective of sampling theory,
%     one cannot expect to recover frequencies in the spectral signals of
%     pixels that are outside the bandlimits of the camera's spectral
%     sensitivity functions. `threshold` is the fraction of the cumulative
%     power spectrum that the user considers to define the bandlimit
%     (because no signal that has finite support is frequency-limited).
%   - 'support_threshold': A fraction indicating what proportion of the
%     peak value of a colour channel's sensitivity function the user
%     considers to be effectively zero (i.e. no sensitivity).
%   - 'bands_padding': This function uses bandlimited interpolation to
%     resample spectral signals. (The interpolation is implicitly
%     represented in the output matrices.) Bandlimited interpolation is
%     intended for functions which are defined over an infinite domain. As
%     the spectral signals are defined over finite domains, this function
%     extends them outside of their domains with values equal to their
%     values at the endpoints of their domains. The extension is
%     implemented as an increase in the weight on the values at the
%     endpoints during interpolation. Unfortunately, I do not know of a
%     closed-form expression for the increased weights. Therefore, I have
%     calculated the weights by summing the interpolation weights for
%     'bands_padding' extra positions to each side of the domains of the
%     spectral signals. 'bands_padding' must be a nonnegative scalar.
%
% verbose -- Debugging flag
%   A Boolean scalar which enables graphical debugging output if `true`.
%   Defaults to `false` if not passed.
%
% ## Output Arguments
%
% color_weights -- Colour conversion matrix
%   A matrix of dimensions `size(color_map, 1)` by `length(bands)` defining
%   a mapping from the spectral space of `bands` to the colour channels
%   defined by `color_map`. `color_weights` accounts for the need to
%   integrate over the spectral space when computing colour channel
%   responses, but does not perform division that may be part of conversion
%   to colour. For instance, when converting to XYZ responses, one
%   normalizes by the integral of the 'Y' CIE colour matching function.
%   Such division should be done either by applying it to `color_map`
%   before calling this function, or by applying it to `color_weights`.
%
% spectral_weights -- Spectral resampling matrix
%   A matrix of dimensions `length(spectral_bands)` by `length(bands)`
%   defining a mapping from the spectral space of `bands` to the spectral
%   space of `spectral_bands`. `spectral_weights` can be used to compare
%   data in the spectral space of `bands` with data in the spectral space
%   of `spectral_bands`, such as for quantitative evaluation against
%   reference spectral data.
%
% bands -- Sampling wavelengths
%   A vector containing the wavelengths at which spectral data should be
%   sampled in order to respect the spectral frequency bandlimit, and
%   spectral support, of `color_map`, but not waste memory.
%
% color_weights_reference -- Colour conversion matrix for reference spectral data
%   The equivalent of `color_weights`, but for mapping from spectral data
%   sampled according to `spectral_bands` to colour channels.
%
% ## Algorithm
%
% To determine the wavelengths (`bands`) at which the spectral data will be
% sampled, this function analyzes `color_map`. It finds the endpoints of
% the wavelength range by determining where all colour channels represented
% by `color_map` first exceed, and last exceed, `options.support_threshold`
% times their peak values. It then finds the number of sample points to use
% within this range based on the approximate bandlimits of the colour
% channel sensitivity functions. In detail, it takes the discrete Fourier
% transform of each sensitivity function, and finds the frequency at which
% the cumulative power spectrum of the sensitivity function exceeds a
% proportion of `options.power_threshold` times the total power. This
% frequency is taken to be the bandlimit of the sensitivity function. The
% highest bandlimit, B, among the sensitivity functions, is selected to
% determine the number of samples. The number of samples is set such that
% the spacing between samples is (1 / 2B).
%
% Given the wavelengths at which the spectral data will be sampled, the
% function can construct `spectral_weights`, a matrix for mapping the
% spectral data to the sampling space of some other spectral data, given by
% `spectral_bands`. Each sample in the space of `spectral_bands` is
% computed by centering a sinc() reconstruction function at the given
% wavelength, and weighting all of the samples in the original space
% (`bands`) by the corresponding values of the sinc() function. The weights
% can be computed in advance, as `spectral_weights`, and then applied to
% arbitrary spectral data.
%
% The conversion to colour is similar: The first step is resampling to the
% space defined by `color_bands`. Next, the element-wise product is taken
% between the resampled data and each channel of `color_map`. Lastly, the
% element-wise products are integrated over the domain of `color_bands`,
% according to the numerical integration method given by
% `options.int_method`, to produce colour channel values. If the spectral
% data is the only variable, and all other quantities are fixed, then the
% process can be represented by a matrix transformation of the spectral
% data, `color_weights`.
%
% ## Notes
% - This function requires that `bands` will end up defining a sampling no
%   finer than `spectral_bands`, such that there is no need to eliminate
%   high frequency information before resampling data defined in the space
%   of `bands` to the space of `spectral_bands`.
% - This function does not require that `color_bands` defines a sampling
%   at least as fine as `spectral_bands`, because it can implicitly
%   resample `color_bands` when computing `color_weights_reference`.
%
% ## References
%
% Chapter 2 of the following book describes ideal sampling and
% reconstruction:
%
%   Goodman, J. W. (2005). Introduction to fourier optics (3rd ed.).
%     Englewood, Colorado: Roberts & Company.
%
% The authors of Physically-Based Rendering Tool (https://pbrt.org/)
% describe their approach for resampling spectral curves, conversion from
% spectral data to colour, and how sampling theory is applied in practice
% in computer graphics:
%
%   Pharr, M., & Humphreys, G. (2010). Physically Based Rendering (2nd
%   ed.). Burlington, Massachusetts: Elsevier.
%
% See also bandlimit, sinc, interp1

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 26, 2018

% Parse input arguments
narginchk(4, 5);
nargoutchk(1, 4);

verbose = false;
if ~isempty(varargin)
    verbose = varargin{1};
end

output_reference_weights = nargout > 3;

if options.bands_padding < 0
    error('`options.bands_padding` must be nonnegative.');
end

diff_color_bands = diff(color_bands);
color_bands_spacing = diff_color_bands(1);
if any(diff_color_bands ~= color_bands_spacing)
    error('`color_bands` must contain equally-spaced values.')
end

diff_spectral_bands = diff(spectral_bands);
spectral_bands_spacing = diff_spectral_bands(1);
if any(diff_spectral_bands ~= spectral_bands_spacing)
    error('`spectral_bands` must contain equally-spaced values.')
end

n_channels = size(color_map, 1);
n_color_bands = size(color_map, 2);
if n_color_bands ~= length(color_bands)
    error('`color_bands` and `color_map` have mismatched sizes.');
end

% Find endpoints of domain
color_map_relative = color_map ./ repmat(max(color_map, 2), 1, n_color_bands);
color_map_relative_logical = any(color_map_relative >= options.support_threshold, 1);
start_band_ind = find(color_map_relative_logical, 1, 'first');
start_band = color_bands(start_band_ind);
end_band_ind = find(color_map_relative_logical, 1, 'last');
end_band = color_bands(end_band_ind);

if verbose
    figure;
    plot_colors = jet(n_color_bands);
    hold on
    for c = 1:n_channels
        plot(color_bands, color_map(c, :), 'Color', plot_colors(c, :), 'LineWidth', 2);
    end
    scatter(start_band, 0, [], [0, 0, 0], 'filled');
    scatter(end_band, 0, [], [0, 0, 0], 'filled');
    hold off
    xlabel('Wavelength [nm]')
    ylabel('Relative quantum efficiency')
    title('Colour channel sensitivities with black points marking the ends of the spectral domain')
end

% Find bandlimit, and construct bands
freq = bandlimit(color_map, options.power_threshold, verbose); % Units of cycles per index
freq_wavelengths = freq / diff_color_bands; % Units of cycles per unit change in wavelength
bands_spacing = 1 / (2 * freq_wavelengths); % Nyquist limit sample spacing
bands = start_band:bands_spacing:end_band;
% Add one sample to the end if the full range is not covered
if bands(end) < end_band
    bands(end + 1) = bands(end) + bands_spacing;
end

if spectral_bands_spacing > bands_spacing
    error(['The spacing of `spectral_bands` is greater than that of `bands`',...
        ', so `spectral_weights` cannot be computed without aliasing.']);
end

% Construct a colour conversion matrix
% The following code is based on the "Ideal Bandlimited Interpolation"
% example in the MATLAB Help page for 'sinc()'.
bands_padded = [
    (bands(1) - (options.bands_padding * bands_spacing)):bands_spacing:(bands(1) - bands_spacing),...
    bands,...
    (bands(end) + bands_spacing):bands_spacing:(bands(end) + (options.bands_padding * bands_spacing))
];
[other_bands_grid, bands_grid] = ndgrid(color_bands, bands_padded);
resampling_map = sinc((other_bands_grid - bands_grid) / bands_spacing);
% Adjust the weights of the endpoints to that resampling assumes the value
% of the spectral signal outside of its domain is equal to its value at the
% nearest endpoint of its domain
resampling_map = [
    sum(resampling_map(:, 1:(options.bands_padding + 1)), 2),...
    resampling_map(:, (options.bands_padding + 2):(end - options.bands_padding - 1)),...
    sum(resampling_map(:, (end - options.bands_padding):end), 2)
];
int_weights = integrationWeights(color_bands, options.int_method);
color_weights = color_map * diag(int_weights) * resampling_map;

% Construct a spectral upsampling matrix
[other_bands_grid, bands_grid] = ndgrid(spectral_bands, bands_padded);
spectral_weights = sinc((other_bands_grid - bands_grid) / bands_spacing);
spectral_weights = [
    sum(spectral_weights(:, 1:(options.bands_padding + 1)), 2),...
    spectral_weights(:, (options.bands_padding + 2):(end - options.bands_padding - 1)),...
    sum(spectral_weights(:, (end - options.bands_padding):end), 2)
];

% Construct a reference colour conversion matrix
if output_reference_weights
    if color_bands_spacing > spectral_bands_spacing
        % Upsample the spectral sensitivities
        color_bands_padded = [
            (color_bands(1) - (options.bands_padding * color_bands_spacing)):color_bands_spacing:(color_bands(1) - color_bands_spacing),...
            color_bands,...
            (color_bands(end) + color_bands_spacing):color_bands_spacing:(color_bands(end) + (options.bands_padding * color_bands_spacing))
        ];
        [other_bands_grid, color_bands_grid] = ndgrid(spectral_bands, color_bands_padded);
        resampling_map = sinc((other_bands_grid - color_bands_grid) / color_bands_spacing);
        resampling_map = [
            sum(resampling_map(:, 1:(options.bands_padding + 1)), 2),...
            resampling_map(:, (options.bands_padding + 2):(end - options.bands_padding - 1)),...
            sum(resampling_map(:, (end - options.bands_padding):end), 2)
        ];
        color_map_resampled = (resampling_map * (color_map.')).';
        int_weights = integrationWeights(spectral_bands, options.int_method);
        color_weights_reference = color_map_resampled * diag(int_weights);
    else
        % Upsample the reference spectra
        spectral_bands_padded = [
            (spectral_bands(1) - (options.bands_padding * spectral_bands_spacing)):spectral_bands_spacing:(spectral_bands(1) - spectral_bands_spacing),...
            spectral_bands,...
            (spectral_bands(end) + spectral_bands_spacing):spectral_bands_spacing:(spectral_bands(end) + (options.bands_padding * spectral_bands_spacing))
        ];
        [other_bands_grid, spectral_bands_grid] = ndgrid(color_bands, spectral_bands_padded);
        resampling_map = sinc((other_bands_grid - spectral_bands_grid) / spectral_bands_spacing);
        resampling_map = [
            sum(resampling_map(:, 1:(options.bands_padding + 1)), 2),...
            resampling_map(:, (options.bands_padding + 2):(end - options.bands_padding - 1)),...
            sum(resampling_map(:, (end - options.bands_padding):end), 2)
        ];
        color_weights_reference = color_map * diag(int_weights) * resampling_map;
    end
end

end
