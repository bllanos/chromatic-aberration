function [...
    color_weights, spectral_weights, bands, color_weights_reference...
] = findSampling(color_map, color_bands, spectral_bands, options, varargin)
% FINDSAMPLING Find an optimal sampled representation for spectral data
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
% color_weights = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
% [color_weights, spectral_weights] = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
% [color_weights, spectral_weights, bands] = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
% [...
%   color_weights, spectral_weights, bands, color_weights_reference...
% ] = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%
% ## Description
% color_weights = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Returns a linear transformation for mapping the sampled representation
%   of spectral data to colour channels.
%
% [color_weights, spectral_weights] = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Additionally returns linear transformation(s) for mapping the sampled
%   representation of spectral data to the space of comparison spectral
%   data.
%
% [color_weights, spectral_weights, bands] = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Additionally returns the wavelengths at which the sampled
%   representation is defined.
%
% [...
%   color_weights, spectral_weights, bands, color_weights_reference...
% ] = findSampling(...
%   color_map, color_bands, spectral_bands, options [, verbose]...
% )
%   Additionally returns linear transformation(s) for mapping the sampled
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
%   `spectral_bands` can also be a cell vector of such vectors, if there
%   are multiple reference spectral sampling schemes.
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
%     sensitivity functions. `power_threshold` is the fraction of the
%     cumulative power spectrum that the user considers to define the
%     bandlimit (because no signal that has finite support is
%     frequency-limited). When 'power_threshold' is one, `bands` will have
%     the same sampling interval as `color_bands`. When 'power_threshold'
%     is less than one, `bands` may have a larger sampling interval,
%     reflecting the bandlimit of the spectral sensitivity functions.
%   - 'n_bands': As an alternative to selecting the number of sample points
%     in `bands` using the cumulative power spectrum of the signals in
%     `color_map`, the number of sample points can be given explicitly. If
%     'n_bands' is zero, the number of sample points will be set according
%     to the value of 'power_threshold'. Otherwise, if 'n_bands' is a
%     positive integer, 'n_bands' will be the number of sample points.
%   - 'support_threshold': A fraction indicating what proportion of the
%     peak magnitude of a colour channel's sensitivity function the user
%     considers to be effectively zero (i.e. no sensitivity).
%   - 'bands_padding': This function uses interpolation by convolution to
%     resample spectral signals. (The interpolation is implicitly
%     represented in the output matrices.) Convolution is intended for
%     functions which are defined over an infinite domain. As the spectral
%     signals are defined over finite domains, this function extends them
%     outside of their domains with values equal to their values at the
%     endpoints of their domains. The extension is implemented as an
%     increase in the weight on the values at the endpoints during
%     interpolation. Usually, there is no closed-form expression for the
%     increased weights. Therefore, I have calculated the weights by
%     summing the interpolation weights for 'bands_padding' extra positions
%     to each side of the domains of the spectral signals. 'bands_padding'
%     must be a nonnegative scalar.
%   - 'interpolant': The function convolved with spectral signals, sampled
%     according to `bands`, to interpolate them during upsampling.
%     'interpolant' is passed to 'resamplingWeights()' as its `f` input
%     argument. Refer to the documentation of 'resamplingWeights.m' for
%     more details.
%   - 'interpolant_ref': Similar to 'interpolant', but used for upsampling
%     other spectral data. 'interpolant_ref' is used to generate
%     `color_weights_reference`, for example.
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
%   If `spectral_bands` is a cell vector, `spectral_weights` is a cell
%   vector of the same dimensions, where `spectral_weights{i}` is the
%   mapping matrix corresponding to `spectral_bands{i}`.
%
% bands -- Sampling wavelengths
%   A vector containing the minimal set of wavelengths at which spectral
%   data should be sampled in order to respect the spectral frequency
%   bandlimit (from `options.power_threshold`) of `color_map`, or the
%   desired number of bands (from `options.n_bands`), as appropriate, and
%   the spectral support of `color_map`.
%
% color_weights_reference -- Colour conversion matrix for reference spectral data
%   The equivalent of `color_weights`, but for mapping from spectral data
%   sampled according to `spectral_bands` to colour channels.
%
%   If `spectral_bands` is a cell vector, `color_weights_reference` is a
%   cell vector of the same dimensions, where `color_weights_reference{i}`
%   is the mapping matrix corresponding to `spectral_bands{i}`.
%
% ## Algorithm
%
% To determine the wavelengths (`bands`) at which the spectral data will be
% sampled, this function analyzes `color_map`. It finds the endpoints of
% the wavelength range by determining where any colour channels represented
% by `color_map` first exceed, and last exceed, `options.support_threshold`
% times their peak values.
%
% When `options.n_bands` is zero, it then finds the number of sample points
% to use within this range based on the approximate bandlimits of the
% colour channel sensitivity functions. In detail, it takes the discrete
% Fourier transform of each sensitivity function, and finds the frequency
% at which the cumulative power spectrum (excluding the zero frequency) of
% the sensitivity function exceeds a proportion of
% `options.power_threshold` times the total power. This frequency is taken
% to be the bandlimit of the sensitivity function. The highest bandlimit,
% B, among the sensitivity functions, is selected to determine the number
% of samples. The number of samples is set such that the spacing between
% samples is (1 / 2B).
%
% Given the wavelengths at which the spectral data will be sampled, the
% function can construct `spectral_weights`, a matrix for mapping the
% spectral data to the sampling space of some other spectral data, given by
% `spectral_bands`. Each sample in the space of `spectral_bands` is
% computed by centering a reconstruction function (`options.interpolant`)
% at the given wavelength, and weighting all of the samples in the original
% space (`bands`) by the corresponding values of the reconstruction
% function. The weights can be computed in advance, as `spectral_weights`,
% and then applied to arbitrary spectral data. Using sinc() as the
% reconstruction function corresponds to ideal bandlimited interpolation,
% but will result in ringing, because the samples are distributed over a
% finite domain, and are not subject to periodic extrapolation.
%
% The conversion to colour is similar: The first step is resampling to the
% space defined by `color_bands`. Next, the element-wise product is taken
% between the resampled data, and each channel of `color_map`. Lastly, the
% element-wise products are integrated over the domain of `color_bands`,
% according to the numerical integration method given by
% `options.int_method`, to produce colour channel values. If the spectral
% data is the only variable, and all other quantities are fixed, then the
% process can be represented by a matrix transformation of the spectral
% data, `color_weights`.
%
% ## Notes
% - Preferably, `bands` will end up defining a sampling no finer than (the
%   cells of) `spectral_bands`, such that there is no need
%   to eliminate high frequency information before resampling data defined
%   in the space of `bands` to the space of `spectral_bands`. When this is
%   not the case, this function applies a sinc() filter to eliminate
%   frequencies above the Nyquist limit of `spectral_bands`, to prevent
%   aliasing during downsampling.
% - This function does not require that `color_bands` defines a sampling
%   at least as fine as `spectral_bands`, because it can implicitly
%   resample `color_bands` when computing `color_weights_reference`.
%   However, when doing so, it will assume that the signals in `color_map`
%   take values outside of the domain of `color_bands` that are equal to
%   their values at the endpoints of the domain. (Refer to the
%   documentation of 'options.bands_padding' above, as the same process is
%   followed.)
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
% See also bandlimit, resamplingWeights, sinc, interp1

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

padding = options.bands_padding;
if padding < 0 || padding ~= round(padding)
    error('`options.bands_padding` must be a nonnegative integer.');
end

use_power_threshold = true;
if options.n_bands < 0 || options.n_bands ~= round(options.n_bands)
    error('`options.n_bands` must be a nonnegative integer.');
elseif options.n_bands > 0
    use_power_threshold = false;
end

color_bands = reshape(color_bands, 1, []);
is_cell_spectral_bands = iscell(spectral_bands);
if ~is_cell_spectral_bands
    spectral_bands = {spectral_bands};
end
n_cells = length(spectral_bands);
for i = 1:n_cells
    spectral_bands{i} = reshape(spectral_bands{i}, 1, []);
end

diff_color_bands = diff(color_bands);
color_bands_spacing = diff_color_bands(1);
if max(abs(diff_color_bands - color_bands_spacing)) > 1e-6
    error('`color_bands` must contain equally-spaced values.')
end

for i = 1:n_cells
    diff_spectral_bands = diff(spectral_bands{i});
    if max(abs(diff_spectral_bands - diff_spectral_bands(1))) > 1e-6
        if is_cell_spectral_bands
            error('`spectral_bands{%d}` does not contain equally-spaced values.', i)
        else
            error('`spectral_bands` must contain equally-spaced values.')
        end
    end
end

n_channels = size(color_map, 1);
n_color_bands = size(color_map, 2);
if n_color_bands ~= length(color_bands)
    error('`color_bands` and `color_map` have mismatched sizes.');
end

% Find endpoints of domain
color_map_relative = abs(color_map) ./ repmat(max(abs(color_map), [], 2), 1, n_color_bands);
color_map_relative_logical = any(color_map_relative >= options.support_threshold, 1);
start_domain_ind = find(color_map_relative_logical, 1, 'first');
start_domain = color_bands(start_domain_ind);
end_domain_ind = find(color_map_relative_logical, 1, 'last');
end_domain = color_bands(end_domain_ind);

if verbose
    fg = figure;
    plot_colors = jet(n_channels);
    hold on
    for c = 1:n_channels
        plot(color_bands, color_map(c, :), 'Color', plot_colors(c, :), 'LineWidth', 2);
    end
    scatter(start_domain, 0, [], [0, 0, 0], 'filled');
    scatter(end_domain, 0, [], [0, 0, 0], 'filled');
    hold off
    xlabel('Wavelength [nm]')
    ylabel('Relative quantum efficiency')
    title('Colour channel sensitivities with black points marking the ends of the spectral domain')
end

% Find bandlimit, and construct bands
match_bands = false;
if ~use_power_threshold
    n_bands = options.n_bands;
    bands_spacing = (end_domain - start_domain) / (n_bands - 1);
else
    if start_domain_ind == 1 && end_domain_ind == n_color_bands && options.power_threshold == 1
        bands_spacing = color_bands_spacing;
        n_bands = n_color_bands;
        match_bands = true;
    else
        freq = bandlimit(color_map, options.power_threshold, verbose); % Units of cycles per index
        freq_wavelengths = freq / color_bands_spacing; % Units of cycles per unit change in wavelength
        bands_spacing = 1 / (2 * freq_wavelengths); % Nyquist limit sample spacing
        n_bands = floor((end_domain - start_domain) / bands_spacing) + 1;
    end
end
% Center the bands within the domain
if match_bands
    bands = color_bands;
elseif n_bands > 1
    bands_range = (n_bands - 1) * bands_spacing;
    start_band = (end_domain + start_domain - bands_range) / 2;
    bands = start_band + (bands_spacing * (0:(n_bands - 1)));
else
    bands = (end_domain + start_domain) / 2;
end

if verbose
    figure(fg);
    hold on
    scatter(bands, zeros(n_bands, 1), [], [1, 0, 1], 'filled');
    hold off
end

% Construct a colour conversion matrix
color_weights = colorWeights(color_map, color_bands, bands, options);

% Construct a spectral resampling matrix
spectral_weights = cell(size(spectral_bands));
for i = 1:n_cells
    spectral_weights{i} = resamplingWeights(...
        spectral_bands{i}, bands, options.interpolant, padding...
    );
end
if ~is_cell_spectral_bands
    spectral_weights = spectral_weights{1};
end

% Construct a reference colour conversion matrix
if output_reference_weights        
    color_weights_reference = cell(size(spectral_bands));
    for i = 1:n_cells
        color_weights_reference{i} = colorWeights(color_map, color_bands, spectral_bands{i}, options);
    end
    if ~is_cell_spectral_bands
        color_weights_reference = color_weights_reference{1};
    end
end

end
