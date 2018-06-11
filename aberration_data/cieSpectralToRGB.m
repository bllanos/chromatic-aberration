function [rgb] = cieSpectralToRGB(lambda_C, C, lambda_R, R, varargin)
% CIESPECTRALTORGB Convert spectral radiance to sRGB using the CIE tristimulus functions
%
% ## Syntax
% rgb = cieSpectralToRGB(lambda_C, C, lambda_R, R [, whitepoint])
%
% ## Description
% rgb = cieSpectralToRGB(lambda_C, C, lambda_R, R [, whitepoint])
%   Returns sRGB values corresponding to the input spectral radiances
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
%   A character vector describing the CIE standard illuminant with which
%   the radiances were obtained. `whitepoint` is passed to the MATLAB
%   'xyz2rgb()' function, and defaults to 'd65' if not passed.
%
% ## Output Arguments
%
% rgb -- Relative sRGB responses
%   An n x 3 matrix, where `rgb(i, :)` is the sRGB response corresponding
%   to the i-th sample, assuming the imaging system has the spectral
%   response of the CIE 1931 color matching functions. Values are clipped
%   to the range [0, 1].
%
% ## Notes
% - Spectral radiances can be obtained from spectral reflectances by
%   multiplying by the spectral power distribution of the illuminant, and
%   normalizing by the "Y" tristimulus value of the illuminant.
% - The tristimulus values of the spectral radiances will be clipped to the
%   range [0, 1].
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
%
% See also xyz2rgb, resampleArrays

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 7, 2018

nargoutchk(1, 1);
narginchk(4, 5);

if ~isempty(varargin)
    whitepoint = varargin{1};
else
    whitepoint = 'd65';
end

% Resample the data, if necessary
[C_resampled, R_resampled, lambda] = resampleArrays(...
    lambda_C, C, lambda_R, R.', 'spline'...
    );
R_resampled = R_resampled.';

lambda = reshape(lambda, length(lambda), 1);
lambda_diff = diff(lambda);
lambda_diff = [
    lambda_diff;
    lambda_diff(end)
    ];

XYZ = R_resampled * diag(lambda_diff) * C_resampled;
XYZ(XYZ < 0) = 0;
XYZ(XYZ > 1) = 1;

% Default colour space used by 'xyz2rgb()' is sRGB
rgb = xyz2rgb(XYZ, 'WhitePoint', whitepoint);
rgb(rgb < 0) = 0;
rgb(rgb > 1) = 1;

end

