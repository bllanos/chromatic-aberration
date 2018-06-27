function [rgb, XYZ] = reflectanceToColor(lambda_L, L, lambda_Ref, Ref, lambda_C, C, varargin)
% REFLECTANCETOCOLOR Convert spectral reflectances to RGB colour values
%
% ## Syntax
% rgb = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, whitepoint]...
% )
% [rgb, XYZ] = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, whitepoint]...
% )
%
% ## Description
% rgb = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, whitepoint]...
% )
%   Returns rgb values for the given reflectances seen under the given
%   illuminant
%
% [rgb, XYZ] = reflectanceToColor(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, whitepoint]...
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
%   sampled. An even sampling of 5 nm is required, spanning the range from
%   360 nm (inclusive) to 780 nm (inclusive).
%
% C -- CIE 1931 color matching functions
%   A 3-column matrix, where `C(i, j)` contains the value of the j-th CIE
%   tristimulus function evaluated at the wavelength `lambda_C(i)`. The
%   length of `lambda_C` must equal the size of the first dimension of `C`.
%
% whitepoint -- Illuminant whitepoint
%   A character vector naming the CIE standard illuminant corresponding to
%   `L`. The 'cieSpectralToColor()' function will use a default value of
%   'd65' if its `whitepoint` input argument is not passed by this
%   function.
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
% resampleArrays

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 12, 2018

nargoutchk(1, 2);
narginchk(6, 7);

[lambda_resampled, ~, Rad_normalized, lambda_C_resampled, C_resampled] = reflectanceToRadiance(...
    lambda_L, L, lambda_Ref, Ref, lambda_C, C, 2 ...
);

if isempty(varargin)
    [rgb, XYZ] = cieSpectralToColor(...
        lambda_resampled, C_resampled,...
        lambda_resampled, Rad_normalized.'...
    );
else
    [rgb, XYZ] = cieSpectralToColor(...
        lambda_C_resampled, C_resampled,...
        lambda_resampled, Rad_normalized.',...
        varargin{1}...
    );
end

end

