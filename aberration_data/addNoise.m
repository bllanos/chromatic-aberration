function J = addNoise(I, snr, varargin)
% ADDNOISE  Add noise to an image
%
% ## Syntax
% J = addNoise(I, snr [, align])
%
% ## Description
% J = addNoise(I, snr [, align])
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
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `I`, if `I` is a 2D array representing raw colour
%   filter array data. For example, 'gbrg'. `align` has the same form as
%   the `sensorAlignment` input argument of `demosaic()`.
%
% ## Output Arguments
%
% J -- Noisy image
%   An image produced by adding noise to `I`. If `I` is a multi-channel
%   image, or a colour filter array image, one noise distribution will be
%   computed for each channel, such that the signal-to-noise ratio is the
%   same for each channel. A multi-channel image will receive noise which
%   is correlated between colour channels: Each colour channel value at a
%   given pixel will be perturbed by a noise value drawn from the same
%   point in the cumulative probability density function of its noise.
%
% ## Notes
% - The signal-to-noise ratio is defined as [1]:
%     SNR = 10 * log_10( ||i||_2^2 / ||e||_2^2 )
%   Where, in this implementation, `i` is a single channel of `I`, as noise
%   is modeled per-channel. `e` is the noise.
% - Usually, image noise is approximated with Gaussian noise, such as in
%   [1-2]. Photon shot noise follows a Poisson distribution [3], and a
%   Poisson distribution has the advantage over a Gaussian distribution of
%   being non-negative, so this function uses Poisson noise. However, the
%   noise is shifted so that its mean value is zero, in order to avoid
%   changing the average image intensity. Negative image values are clipped
%   to zero.
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
narginchk(2, 3);

image_height = size(I, 1);
image_width = size(I, 2);
n_px = image_height * image_width;
n_channels = size(I);

is_raw = ~isempty(varargin);
if is_raw && (n_channels > 1)
    error('The `align` argument should not be passed if the image `I` has more than one channel.');
end

if snr <= 0
    error('The signal to noise ratio must be positive.')
elseif isinf(snr)
    J = I;
    return;
end

signal = zeros(n_channels, 1);
if is_raw
    align = varargin{1};
    mask = bayerMask( image_height, image_width, align );
    n_px_c = zeros(n_channels, 1);
    for c = 1:n_channels
        mask_c = mask(:, :, c);
        I_c = I(mask_c);
        n_px_c(c) = length(I_c);
        signal(c) = dot(I_c, I_c) / length(I_c);
    end
else
    for c = 1:n_channels
        signal(c) = sum(sum(I(:, :, c) .^ 2)) / n_px;
    end
end

noise = sqrt(signal * (10 ^ (-snr / 10)));

J = zeros(size(I));
if is_raw
    for c = 1:n_channels
        pd = makedist('Poisson', 'lambda', noise(c));
        noise_c = random(pd, n_px_c(c), 1) - noise(c);
        mask_c = mask(:, :, c);
        J(mask_c) = I(mask_c) + noise_c;
    end
else
    pd = makedist('Poisson', 'lambda', noise(1));
    noise_1 = random(pd, image_height, image_width);
    cdf_1 = cdf(pd, noise_1);
    J(:, :, 1) = I(:, :, 1) + (noise_1 - noise(1));
    for c = 2:n_channels
        pd = makedist('Poisson', 'lambda', noise(c));
        noise_c = icdf(pd, cdf_1) - noise(c);
        J(:, :, c) = I(:, :, c) + noise_c;
    end
end

J(J < 0) = 0;

end