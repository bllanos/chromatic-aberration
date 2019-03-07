function [rgb, XYZ] = cieSpectralToColor(lambda_C, C, lambda_R, R, varargin)
% CIESPECTRALTOCOLOR Convert spectral radiance to sRGB using the CIE tristimulus functions
%
% ## Syntax
% rgb = cieSpectralToColor(lambda_C, C, lambda_R, R [, whitepoint, int_method])
% [rgb, XYZ] = cieSpectralToColor(lambda_C, C, lambda_R, R [, whitepoint, int_method])
%
% ## Description
% rgb = cieSpectralToColor(lambda_C, C, lambda_R, R [, whitepoint, int_method])
%   Returns sRGB values corresponding to the input spectral radiances
%
% [rgb, XYZ] = cieSpectralToColor(lambda_C, C, lambda_R, R [, whitepoint, int_method])
%   Additionally returns the CIE 1931 tristimulus values corresponding to
%   the input spectral radiances
%
% ## Input Arguments
%
% lambda_C -- Reference wavelength values
%   A vector of wavelengths at which the 'C(lambda)' distributions were
%   sampled.
%
% C -- CIE 1931 color matching functions
%   A 3-column matrix, where `C(i, j)` contains the value of the j-th CIE
%   tristimulus function evaluated at the wavelength `lambda_C(i)`. The
%   length of `lambda_C` must equal the size of the first dimension of `C`.
%
%   If the range of wavelengths is a sub-interval of 360 to 780 nm, the
%   values of `C(:, 1)` and `C(:, end)` should be equal to the sums of the
%   tristimulus function values from the endpoints of the sub-interval to
%   the endpoints of the 360 to 780 nm interval. This recommendation is in
%   the ASTM E308 standard (Section 7.3.2.2).
%
% lambda_R -- Wavelengths
%   A vector of wavelengths at which the spectral radiances were sampled.
%
% R -- Spectral radiances
%   A matrix of size n x length(lambda_R), where 'n' is the number of
%   different spectral radiances. The columns of `R` correspond to the
%   wavelengths in `lambda_R`. `R(i, j)` is the radiance of the i-th sample
%   at the j-th wavelength.
%
% whitepoint -- Illuminant whitepoint
%   A three-element vector containing the XYZ colour of the illuminant, or a
%   character vector describing the CIE standard illuminant with which the
%   radiances were obtained. `whitepoint` is passed to the MATLAB 'xyz2rgb()'
%   function, and defaults to 'd65' if empty or if not passed.
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
%   to the i-th sample, assuming the imaging system has the spectral
%   response of the CIE 1931 color matching functions. Values are clipped
%   to the range [0, 1].
%
% XYZ -- CIE 1931 tristimulus responses
%   An n x 3 matrix, where `XYZ(i, :)` is the tristimulus (X, Y, and Z) response
%   corresponding to the i-th sample, assuming the imaging system has the
%   spectral response of the CIE 1931 color matching functions. Values are NOT
%   clipped to the range [0, 1], but should be within this range if `R` is
%   normalized (see notes below).
%
% ## Notes
% - Spectral radiances can be obtained from spectral reflectances by
%   multiplying by the spectral power distribution of the illuminant, and
%   normalizing by the "Y" tristimulus value of the illuminant. The
%   'reflectanceToColor()' function can be used for this purpose.
% - This function will resample `R` and `C` if they were sampled at
%   different wavelengths.
%
% ## References
% - Foster, D. H. (2018). Tutorial on Transforming Hyperspectral Images to
%   RGB Colour Images. Retrieved from
%   http://personalpages.manchester.ac.uk/staff/d.h.foster/Tutorial_HSI2RGB/Tutorial_HSI2RGB.html
%   on June 5, 2018.
% - Lindbloom, Bruce J. (2017). Computing XYZ From Spectral Data. Retrieved
%   from http://www.brucelindbloom.com on June 11, 2018.
% - ASTM E308-17 Standard Practice for Computing the Colors of Objects by
%   Using the CIE System, ASTM International, West Conshohocken, PA, 2017,
%   https://doi.org/10.1520/E0308-17
%
% See also xyz2rgb, resampleArrays, reflectanceToColor,
% reflectanceToRadiance, integrationWeights

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2018

nargoutchk(1, 2);
narginchk(4, 6);

whitepoint = [];
int_method = 'rect';
if ~isempty(varargin)
    whitepoint = varargin{1};
    if length(varargin) > 1
        int_method = varargin{2};
    end
end
if isempty(whitepoint)
    whitepoint = 'd65';
end

% Resample the data, if necessary
[C_resampled, R_resampled, lambda] = resampleArrays(...
    lambda_C, C, lambda_R, R, 'spline'...
    );

lambda = reshape(lambda, length(lambda), 1);
weights = repmat(integrationWeights(lambda, int_method).', size(R, 2), 1);

XYZ = (R_resampled.' .* weights) * C_resampled;
% XYZ(XYZ < 0) = 0;
% XYZ(XYZ > 1) = 1;

% Default colour space used by 'xyz2rgb()' is sRGB
rgb = xyz2rgb(XYZ, 'WhitePoint', whitepoint);
rgb(rgb < 0) = 0;
rgb(rgb > 1) = 1;

end

