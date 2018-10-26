function [color_weights, spectral_weights, bands] = samplingWeights(color_map, color_bands, spectral_bands, options, varargin)
% SAMPLINGWEIGHTS Find an optimal sampled representation for spectral data
%
% ## About
%
% The purpose of this function is to find the most memory-efficient
% representation for spectral data that is converted to colour according to
% the colour channel sensitivities in `color_map`: We seek the smallest
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
%   wavelength corresponding to `color_map(:, j)`.
%
% spectral_bands -- Wavelength bands for spectral reference data
%   A vector containing the wavelengths at which data to which the spectral
%   data is to be compared to has been sampled.
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
%   sampled to respect the "spectral frequency bandlimit" and spectral
%   support of `color_map`, without wasting memory.
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
% transform of each sensitivity function, and finds the frequency bounding
% `options.power_threshold` of the cumulative powerspectrum of the
% sensitivity function. This frequency is taken to be the bandlimit of the
% sensitivity function. The highest bandlimit, B, among the sensitivity
% functions, is selected to determine the number of samples. The number of
% samples is set such that the spacing between samples is (1 / 2B).
%
% Given the wavelengths at which the spectral data will be sampled, the
% function can construct `spectral_weights`, a matrix for mapping the
% spectral data to the sampling space of some other spectral data, given by
% `spectral_bands`. Each sample in the space of `spectral_bands` is
% computed by centering a sinc() reconstruction function at the given
% wavelength, and weighting all of the samples in the original space
% (`bands`) by the corresponding values of the sinc() function. Lastly, the
% weights are normalized by their sum. The normalized weights can be
% computed in advance, as `spectral_weights`, and then applied to arbitrary
% spectral data.
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
% - Resampling to spectral spaces of larger domains is performed by
%   assuming that the data to be resampled has constant values outside of
%   its domain equal to its values at the endpoints of its domain.
% - This function assumes that `bands` will end up defining a coarser
%   sampling than both `color_bands` and `spectral_bands`, such that there
%   is no need to eliminate high frequency information before resampling
%   data defined in the space of `bands` to the spaces of `color_bands` and
%   `spectral_bands`.
%
% ## References
%
% Chapter 2 of the following book describes ideal sampling and
% reconstruction:
%   Goodman, J. W. (2005). Introduction to fourier optics (3rd ed.).
%     Englewood, Colorado: Roberts & Company.
%
% The authors of Physically-Based Rendering Tool (https://pbrt.org/)
% describe their approach for resampling spectral curves, conversion from
% spectral data to colour, and how sampling theory is applied in practice
% in computer graphics:
%   Pharr, M., & Humphreys, G. (2010). Physically Based Rendering (2nd
%   ed.). Burlington, Massachusetts: Elsevier.
%
% See also interp1

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 26, 2018

% Parse input arguments
narginchk(4, 5);
nargoutchk(1, 3);

end