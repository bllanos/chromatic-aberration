function color_weights = colorWeights(...
    color_map, color_bands, bands, options...
)
% COLORWEIGHTS Create a colour conversion matrix operator
%
% ## Syntax
% color_weights = colorWeights(color_map, color_bands, bands, options)
%
% ## Description
% color_weights = colorWeights(color_map, color_bands, bands, options)
%   Returns a colour conversion matrix operating on the spectral sampling space
%   of `bands`.
%
% ## Input Arguments
%
% color_map -- Colour channel spectral sensitivities
%   A 2D array, where `color_map(i, j)` is the sensitivity of the i-th colour
%   channel to the j-th spectral band in `color_bands`. `color_map` is not a
%   colour conversion matrix, as it does not perform the desired numerical
%   integration, over the spectrum, that is part of colour conversion.
%
% color_bands -- Wavelength bands for colour channel sensitivities
%   A vector, of length equal to the size of the second dimension of
%   `color_map`, containing the wavelengths at which the sensitivity
%   functions in `color_map` have been sampled. `color_bands(j)` is the
%   wavelength corresponding to `color_map(:, j)`. The values in
%   `color_bands` are expected to be evenly-spaced.
%
% bands -- Input spectral sampling
%   A vector containing the wavelengths at which the spectral images that need
%   to be converted to colour are sampled. The elements of `bands` are expected
%   to be evenly-spaced.
%
% options -- Spectral resampling and integration options
%   `options` is a structure with the following fields:
%   - 'int_method': The numerical integration method to use when
%     integrating over the colour channels sensitivities to compute colour
%     values. `int_method` is passed to `integrationWeights()` as its `method`
%     input argument.
%   - 'bands_padding': A scalar used by 'resamplingWeights()' to control
%     extrapolation outside the limits of `bands` and `color_bands`. Refer to
%     the documentation of the `padding` input parameter of
%     'resamplingWeights()' in 'resamplingWeights.m'.
%   - 'interpolant': A handle to a function used to upsample `bands`, if
%     necessary, by 'resamplingWeights()'. Refer to the documentation of the `f`
%     input parameter of 'resamplingWeights()' in 'resamplingWeights.m'. -
%   - 'interpolant_ref': Similar to 'interpolant', but used for upsampling
%     `color_bands` if necessary.
%
% ## Output Arguments
%
% color_weights -- Colour conversion matrix
%   A matrix of dimensions `size(color_map, 1)` by `length(bands)` which
%   computes the colours, in the colour space of `color_map`, of spectra sampled
%   according to `bands`.
%
% See also resamplingWeights, integrationWeights, findSampling,
% channelConversion

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 1, 2019

narginchk(4, 4);
nargoutchk(1, 1);

if size(color_map, 2) ~= length(color_bands)
    error('`color_bands` and `color_map` have mismatched sizes.');
end

diff_color_bands = diff(color_bands);
color_bands_spacing = diff_color_bands(1);
if max(abs(diff_color_bands - color_bands_spacing)) > 1e-6
    error('`color_bands` must contain equally-spaced values.')
end

diff_bands = diff(bands);
bands_spacing = diff_bands(1);
if max(abs(diff_bands - bands_spacing)) > 1e-6
    error('`bands` must contain equally-spaced values.')
end

if color_bands_spacing > bands_spacing
    % Upsample the spectral sensitivities
    upsampling_map = resamplingWeights(...
        bands, color_bands, options.interpolant_ref, options.bands_padding...
    );
    color_map_upsampled = (upsampling_map * (color_map.')).';
    int_weights = integrationWeights(bands, options.int_method);
    color_weights = color_map_upsampled * diag(int_weights);
    % Also resample the sampled spectra, because the interpolant may not produce
    % an identity mapping
    resampling_map = resamplingWeights(...
        bands, bands, options.interpolant, options.bands_padding...
    );
    color_weights = color_weights * resampling_map;
else
    % Upsample the sampled spectra
    upsampling_map = resamplingWeights(...
        color_bands, bands, options.interpolant, options.bands_padding...
    );
    int_weights = integrationWeights(color_bands, options.int_method);
    color_weights = color_map * diag(int_weights) * upsampling_map;
end

end
