function J = addNoise(I, snr)
% ADDNOISE  Add noise to an image
%
% ## Syntax
% J = addNoise(I, snr)
%
% ## Description
% J = addNoise(I, snr)
%   Returns a version of the input image with added noise.
%
% ## Input Arguments
%
% I -- Input image
%   A greyscale or multi-channel image, in the form of a 2D or 3D array.
%
% snr -- Signal to noise ratio
%   The desired signal to noise ratio of the output image `J`.
%
% ## Output Arguments
%
% J -- Noisy image
%   An image produced by adding noise to `I`. One noise distribution will
%   be computed for each pixel, such that the signal to noise ratio is
%   respected for every pixel. Note that this method almost simulates
%   photon shot noise and pixel response non-uniformity noise that scale
%   with image intensity, except that photon shot noise would scale with
%   the square root of the image intensity [3]. A multi-channel image will
%   receive noise which is correlated between colour channels: Each colour
%   channel value at a given pixel will be perturbed by a noise value drawn
%   from the same point in the cumulative probability density function of
%   its noise.
%
% ## Notes
% - The signal-to-noise ratio is defined as [1]:
%     SNR = 10 * log_10( ||i||_2^2 / ||e||_2^2 )
%   Where, in this implementation, `i` is a single element of `I`, as noise
%   is modeled per-channel, per-pixel. `e` is the noise.
% - Usually, image noise is approximated with Gaussian noise, such as in
%   [1-2]. Photon shot noise follows a Poisson distribution, but other
%   sources of noise are approximately Gaussian [3], so this function uses
%   Gaussian noise. Negative corrupted image values are clipped to zero.
% - This function does not add or account for crosstalk. Crosstalk is
%   likely already present in the input image, unless the image is purely
%   synthetic, so adding more crosstalk between colour channels would
%   produce images which are not realistic. Also, crosstalk is generally
%   spatially-varying (and therefore hard to model). See [4] for a method
%   for modelling crosstalk (which requires professional equipment).
%
% ## References
% 1 - Song, Y., Brie, D., Djermoune, E.-H., & Henrot, S.. "Regularization
%     Parameter Estimation for Non-Negative Hyperspectral Image
%     Deconvolution." IEEE Transactions on Image Processing, vol. 25, no.
%     11, pp. 5316-5330, 2016. doi:10.1109/TIP.2016.2601489
% 2 - Tan, H., Zeng, X., Lai, S., Liu, Y., & Zhang, M. (2018). (2018).
%     "Joint demosaicing and denoising of noisy bayer images with ADMM."
%     Proc. International Conference on Image Processing, ICIP, pp.
%     2951-2955. doi:10.1109/ICIP.2017.8296823
% 3 - Martinec, E. (2008). "Noise, dynamic range and bit depth in digital
%     SLR." Retrieved from
%     http://theory.uchicago.edu/∼ejm/pix/20d/tests/noise/
% 4 - Qiu, J., & Xu, H. (2016). "Camera response prediction for various
%     capture settings using the spectral sensitivity and crosstalk model."
%     Applied Optics, 55(25), 6989-6999. doi:10.1364/AO.55.006989
%
% See also PoissonDistribution, bayerMask

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 18, 2018

nargoutchk(1, 1);
narginchk(2, 2);

n_channels = size(I);

if snr <= 0
    error('The signal to noise ratio must be positive.')
elseif isinf(snr)
    J = I;
    return;
end

signal = I .^ 2;
noise_std = sqrt(signal * (10 ^ (-snr / 10)));

J = zeros(size(I));
noise_1 = randn(size(I, 1), size(I, 2));
J(:, :, 1) = I(:, :, 1) + noise_1 .* noise_std(:, :, 1);
for c = 2:n_channels
    J(:, :, c) = I(:, :, c) + noise_1 .* noise_std(:, :, c);
end

J(J < 0) = 0;

end