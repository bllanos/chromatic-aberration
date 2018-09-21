function [I_hyper, I_warped] = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
% CHIRPIMAGE  Create an image with spatial and spectral linear chirp, and spectral dispersion
%
% ## Syntax
% I_hyper = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
% [I_hyper, I_warped] = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
%
% ## Description
% I_hyper = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
%   Returns a version of the hyperspectral image without dispersion
%
% [I_hyper, I_warped] = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
%   Additionally returns a version of the hyperspectral image with
%   dispersion
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A three-element vector containing the height, width, and number of
%   wavelength bands, respectively, of the image.
%
% lambda_range -- Range of wavelengths
%   A two element vector containing the minimum and maximum wavelengths
%   captured in the image. If the spacing between the centers of
%   consecutive wavelength bands is '`delta_lambda`, then the center of the
%   first wavelength band is `lambda_range(1) + 0.5 * delta_lambda`, and
%   the center of the last wavelength band (having index
%   `image_sampling(3)`) is `lambda_range(2) - 0.5 * delta_lambda`.
%
% d_fraction -- Amount of dispersion
%   A scalar from 0 to 1, where 0 indicates no dispersion, and 1 indicates
%   maximum dispersion. The amount of dispersion scales linearly with
%   `d_fraction`. Under the maximum dispersion, the hypothetical image for
%   wavelength `lambda_range(2)` is a shifted to the right relative to the
%   hypothetical image for wavelength `lambda_range(1)` by an amount
%   corresponding to the distance from the left border of the latter image
%   to the x-coordinate of its first intensity minimum.
%
% n_samples -- Number of samples per pixel
%   As far as I know, there is no analytical form for the integral of the
%   image intensity function. This integral would need to be taken over the
%   x, y, and spectral domain of a pixel in order to compute the pixel's
%   value. This function approximates the integral with the average of a
%   grid of samples in the pixel's domain. `n_samples` is a three-element
%   vector containing the number of grid voxels along the x, y, and
%   spectral dimensions, respectively. If `n_samples` is empty, or is a
%   vector of ones, then each pixel is given the value of the image
%   intensity function evaluated at its centroid.
%
% ## Output Arguments
%
% I_hyper -- Chirp image
%   An image produced by evaluating a linear chirp function. As the
%   x-coordinate increases, the amplitude of the spectral signal oscillates
%   increasingly rapidly, until reaching a maximum rate of one cycle per
%   two pixels at the right edge of the image. As the y-coordinate
%   increases, the frequency of oscillation in the spectral dimension
%   increases, until reaching a maximum rate of one cycle across two
%   spectral bands at the bottom edge of the image. (The maximum
%   frequencies are intended to respect the Whittaker-Shannon sampling
%   theorem described by Goodman 2005.) Consequently, this image contains a
%   range of both spatial and spectral frequencies.
%
% I_warped -- Dispersed image
%   A version of `I_hyper` altered by spectral dispersion, with the
%   magnitude of the dispersion given by `d_fraction`. The dispersion is a
%   shift in the positive x-direction with wavelength. Note that the
%   dispersed image is synthesized directly, rather than as a warped
%   version of `I_hyper`.
%
% ## References
% The chirp function is described in Chapter 2 of:
%
%   Goodman, J. W. (2005). Introduction to fourier optics (3rd ed.).
%     Englewood, Colorado: Roberts & Company.
%
% See also imageFormation

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 20, 2018

nargoutchk(1, 2);
narginchk(24, 4);

output_warped = nargout > 1;

image_height = image_sampling(1);
image_width = image_sampling(2);
n_bands = image_sampling(3);

alpha = pi / image_width;
lambda_0 = lambda_range(1);
lambda_1 = lambda_range(2);
lambda_span = diff(lambda_range);
delta_lambda = lambda_span / n_bands;
beta = pi / (image_height * delta_lambda);
d_max = sqrt(image_width) / lambda_span;
d_absolute = d_fraction * d_max;

    function intensity = imageIntensity(x, y, lambda, d)
        lambda_rel = lambda - lambda_0;
        intensity = (cos(alpha .* (x - d .* lambda_rel) .^ 2) + 1) .* ...
            (cos(beta .* y .* lambda_rel) + 1) ./ 4;
    end

if isempty(n_samples) || all(n_samples == 1)
    lambda = linspace(lambda_0 + delta_lambda / 2, lambda_1 - delta_lambda / 2, n_bands);
    [Y, X, Lambda] = ndgrid(0.5:1:image_height, 0.5:1:image_width, lambda);
    I_hyper = imageIntensity(X, Y, Lambda, 0);
    if output_warped
        I_warped = imageIntensity(X, Y, Lambda, d_absolute);
    end
else
    n_dims = 3;
    if n_dims ~= length(n_samples)
        error('`n_samples` should be a %d-element vector.', n_dims);
    end
    offsets = cell(n_dims, 1);
    delta = [1, 1, delta_lambda];
    for d = 1:n_dims
        offsets{d} = linspace(...
            delta / (2 * n_samples(d)), delta - (delta / (2 * n_samples(d))),...
            n_samples(d)...
        );
    end
    n_samples_all = prod(n_samples);
    offsets_all = zeros(n_samples_all, n_dims);
    for d = 1:n_dims
        offsets_all(:, d) = repmat(...
            repelem(offsets{d}, prod(n_samples((d + 1):end))),...
            prod(n_samples(1:(d - 1))), 1 ...
        );
    end
    lambda = linspace(lambda_0, lambda_1 - delta_lambda, n_bands);
    [Y, X, Lambda] = ndgrid(0:(image_height - 1), 0:(image_width - 1), lambda);
    I_hyper = zeros(image_sampling);
    if output_warped
        I_warped = zeros(image_sampling);
    end
    for s = 1:n_samples_all
        Y_s = Y + offsets_all(s, 2);
        X_s = X + offsets_all(s, 1);
        Lambda_s = Lambda + offsets_all(s, 3);
        I_hyper = I_hyper + imageIntensity(X_s, Y_s, Lambda_s, 0);
        if output_warped
            I_warped = I_warped + imageIntensity(X_s, Y_s, Lambda_s, d_absolute);
        end
    end
    n_px_all = prod(image_sampling);
    I_hyper = I_hyper / n_px_all;
    if output_warped
        I_warped = I_warped / n_px_all;
    end
end

end
