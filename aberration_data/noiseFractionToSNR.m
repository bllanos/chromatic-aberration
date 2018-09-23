function snr = noiseFractionToSNR(fraction)
% NOISEFRACTIONTOSNR  Convert relative amounts of noise to a signal-to-noise ratio
%
% ## Syntax
% snr = noiseFractionToSNR(fraction)
%
% ## Description
% snr = noiseFractionToSNR(fraction)
%   Returns the signal to noise ratios corresponding to the noise
%   fractions.
%
% ## Input Arguments
%
% fractions -- Noise fractions
%   An array of noise fractions, where `fractions(i)` is the ratio of the
%   standard deviation of the noise to the magnitude of the signal.
%
% ## Output Arguments
%
% snr -- Signal to noise ratios
%   The signal to noise ratios corresponding to `fractions` output as an
%   array with the same dimensions as `fractions`.
%
% ## Notes
% - The signal-to-noise ratio is defined as [1]:
%     SNR = 10 * log_10( ||i||_2^2 / ||e||_2^2 )
%   Where, `i` is the signal and `e` is the noise (zero-mean noise).
%
% ## References
% 1 - Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
%
% See also snrToNoiseFraction, addNoise

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 21, 2018

nargoutchk(1, 1);
narginchk(1, 1);

if any(fraction(:) < 0)
    error('The noise fractions must be nonnegative.')
end

snr = -20 * log10(fraction);

end
