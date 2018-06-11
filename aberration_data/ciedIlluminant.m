function [p, lambda] = ciedIlluminant(T, lambda_S, S, lambda)
% CIEDILLUMINANT Spectral power distribution of a CIE D-Illuminant
%
% ## Syntax
% p = ciedIlluminant(T, lambda_S, S, lambda)
%
% ## Description
% p = ciedIlluminant(T, lambda_S, S, lambda)
%   Returns relative spectral powers at the given wavelengths
%
% ## Input Arguments
%
% T -- Colour temperature
%   The correlated colour temperature of the illuminant, measured in Kelvin.
%
% lambda_S -- Reference wavelength values
%   A vector of wavelengths at which the 'S_n(lambda)' distributions were
%   sampled.
%
% S -- 'S_n(lambda)' distributions
%   The 'S_0(lambda)', 'S_1(lambda)', and 'S_2(lambda)' distributions,
%   evaluated at the wavelengths in `lambda_S`, form the columns of this
%   array.
%
% lambda -- Wavelengths
%   A vector containing wavelengths of light at which to sample the
%   spectral power distribution, measured in nanometres.
%
% ## Output Arguments
%
% p -- Relative spectral powers
%   A column vector the same length as the output argument `lambda`,
%   containing the spectral power distribution values of the CIE
%   D-Illuminant with correlated colour temperature `T` at the wavelengths
%   in `lambda`. `p` is normalized so that its value at approximately 560
%   nm is `1`.
%
% lambda -- Wavelengths
%   The wavelengths corresponding to the values in `p`. `lambda` is the
%   portion of `lambda_S` or the input argument `lambda` in the interval
%   `[max(min(lambda_S), min(lambda), min(max(lambda_S), max(lambda))]`,
%   whichever has more elements in that interval.
%
% ## References
% - Lindbloom, Bruce J. (2017). Spectral Power Distribution of a CIE
%   D-Illuminant. Retrieved from http://www.brucelindbloom.com on June 4,
%   2018.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 6, 2018

nargoutchk(1, 2);
narginchk(4, 4);

if T < 4000
    error('The color temperature, `T`, must be in the range [4000, 25000].');
elseif T <= 7000
    x = (-4.6070e9 / (T^3)) + (2.9678e6 / (T^2)) + (0.09911e3 / T) + 0.244063;
elseif T <= 25000
    x = (-2.0064e9 / (T^3)) + (1.9018e6 / (T^2)) + (0.24748e3 / T) + 0.237040;
else
    error('The color temperature, `T`, must be in the range [4000, 25000].');
end

y = -3 * (x ^ 2) + 2.870 * x - 0.275;

M = 0.0241 + 0.2562 * x - 0.7341 * y;
M1 = (-1.3515 - 1.7703 * x + 5.9114 * y) / M;
M2 = (0.0300 - 31.4424 * x + 30.0717 * y) / M;

[S_resampled, lambda] = resampleArrays(lambda_S, S.', lambda, 'spline');
S_resampled = S_resampled.';

p = S_resampled(:, 1) + M1 * S_resampled(:, 2) + M2 * S_resampled(:, 3);
[~, ind] = min(abs(lambda - 560));
p = p / p(ind);

end

