function fraction = snrToNoiseFraction(snr)
% SNRTONOISEFRACTION  Convert a signal-to-noise ratio to a relative amount of noise
%
% ## Syntax
% fraction = snrToNoiseFraction(snr)
%
% ## Description
% fraction = snrToNoiseFraction(snr)
%   Returns the noise fractions corresponding to the signal to noise
%   ratios.
%
% ## Input Arguments
%
% snr -- Signal to noise ratios
%   An array of signal to noise ratios.
%
% ## Output Arguments
%
% fractions -- Noise fractions
%   An array of noise fractions, where `fractions(i)` is the ratio of the
%   standard deviation of the noise to the magnitude of the signal
%   corresponding to `snr(i)`. `fractions` has the same dimensions as
%   `snr`.
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
% See also noiseFractionToSNR, addNoise

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 21, 2018

nargoutchk(1, 1);
narginchk(1, 1);

fraction = (10 ^ -(snr / 20)); % 20 instead of 10 accounts for the square root

end
