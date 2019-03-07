function [rgb, XYZ] = reflectanceToColor(lambda_L, L, lambda_Ref, Ref, lambda_C, C, varargin)
% REFLECTANCETOCOLOR Convert spectral reflectances to RGB colour values
%
% ## Syntax
% rgb = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, int_method]...
% )
% lambda_xyzbar, xyzbar...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, int_method]...
% )
%
% ## Description
% rgb = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, int_method]...
% )
%   Returns rgb values for the given reflectances seen under the given
%   illuminant
%
% [rgb, XYZ] = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, int_method]...
% )
%   Additionally returns CIE 1931 tristimulus values for the given
%   reflectances seen under the given illuminant
%
% ## Input Arguments
%
% lambda_L -- Illuminant wavelengths
%   A vector containing wavelengths of light at which the spectral power
%   distribution of the illuminant has been sampled. Units are nanometres.
%
% L -- Illuminant
%   A vector the same length as `lambda_L` containing the spectral power
%   distribution of the illuminant. `L(i)` is the spectral power at the
%   wavelength `lambda_L(i)`. The global scale of the values is not
%   important, as the scale will be compensated for by normalization.
%
% lambda_Ref -- Reflectance wavelengths
%   A vector containing wavelengths of light at which the spectral
%   reflectances have been sampled. Units are nanometres.
%
% Ref -- Spectral reflectances
%   A matrix of size length(lambda_Ref) x n, where 'n' is the number of
%   different spectral reflectances. The rows of `Ref` correspond to the
%   wavelengths in `lambda_Ref`. `Ref(i, j)` is the reflectance of the j-th
%   sample at the i-th wavelength. Reflectances are assumed to be fractions
%   in the range [0, 1].
%
% lambda_C -- Reference wavelength values
%   A vector of wavelengths at which the 'C(lambda)' distributions were
%   sampled. An even sampling of 5 nm, spanning the range from 360 nm
%   (inclusive) to 780 nm (inclusive), is required for this function to
%   conform to the ASTM E308 standard.
%
% C -- CIE 1931 color matching functions
%   A 3-column matrix, where `C(i, j)` contains the value of the j-th CIE
%   tristimulus function evaluated at the wavelength `lambda_C(i)`. The
%   length of `lambda_C` must equal the size of the first dimension of `C`.
%
% int_method -- Numerical integration method
%   The numerical integration method to use when integrating over the
%   tristimulus functions. `int_method` is passed to `integrationWeights()`
%   as its `method` input argument. A value of 'rect' is required for this
%   function to operate according to the ASTM E308 standard.
%
%   Defaults to 'rect' if not passed.
%
% ## Output Arguments
%
% rgb -- sRGB responses
%   An n x 3 matrix, where `rgb(i, :)` is the sRGB response corresponding
%   to the i-th spectral reflectance, given that the imaging system has the
%   spectral response of the CIE 1931 color matching functions. Values are
%   clipped to the range [0, 1].
%
% XYZ -- CIE 1931 tristimulus responses
%   An n x 3 matrix, where `XYZ(i, :)` is the tristimulus (X, Y, and Z)
%   response corresponding to the i-th spectral reflectance, assuming the
%   imaging system has the spectral response of the CIE 1931 color matching
%   functions. Values are clipped to the range [0, 1].
%
% ## References
% - Lindbloom, Bruce J. (2017). Computing XYZ From Spectral Data. Retrieved
%   from http://www.brucelindbloom.com on June 11, 2018.
% - ASTM E308-17 Standard Practice for Computing the Colors of Objects by
%   Using the CIE System, ASTM International, West Conshohocken, PA, 2017,
%   https://doi.org/10.1520/E0308-17
%
% See also reflectanceToRadiance, ciedIlluminant, cieSpectralToColor,
% resampleArrays, integrationWeights

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 12, 2018

nargoutchk(1, 2);
narginchk(6, 7);

if isempty(varargin)
    int_method = 'rect';
else
    int_method = varargin{1};
end

% Calculate the whitepoint by adding an ideal reflectance
Ref_expanded = [Ref, ones(length(lambda_Ref), 1)];

[lambda_resampled, ~, Rad_normalized, lambda_C_resampled, C_resampled] = reflectanceToRadiance(...
    lambda_L, L, lambda_Ref, Ref_expanded, lambda_C, C, 2, int_method...
);

[~, whitepoint] = cieSpectralToColor(...
    lambda_C_resampled, C_resampled,...
    lambda_resampled, Rad_normalized(:, end),...
    [], int_method...
);

[rgb, XYZ] = cieSpectralToColor(...
    lambda_C_resampled, C_resampled,...
    lambda_resampled, Rad_normalized(:, 1:(end - 1)),...
    whitepoint, int_method...
);

end

