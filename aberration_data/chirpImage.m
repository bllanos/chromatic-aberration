function varargout = chirpImage(image_sampling, lambda_range, d_fraction, n_samples, varargin)
% CHIRPIMAGE  Create an image with spatial and spectral linear chirp, and spectral dispersion
%
% ## Syntax
% I_hyper = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
% [I_hyper, I_warped] = chirpImage(...
%   image_sampling, lambda_range, d_fraction, n_samples...
% )
% [I_hyper, I_warped, dispersionfun] = chirpImage(...
%   image_sampling, lambda_range, d_fraction, n_samples...
% )
% [delta_lambda, dispersion_max, bands] = chirpImage(...
%   image_sampling, lambda_range, d_fraction, n_samples, 'params'...
% )
%
% ## Description
% I_hyper = chirpImage(image_sampling, lambda_range, d_fraction, n_samples)
%   Returns a version of the hyperspectral image without dispersion
%
% [I_hyper, I_warped] = chirpImage(...
%   image_sampling, lambda_range, d_fraction, n_samples...
% )
%   Additionally returns a version of the hyperspectral image with
%   dispersion
%
% [I_hyper, I_warped, dispersionfun] = chirpImage(...
%   image_sampling, lambda_range, d_fraction, n_samples...
% )
%   Additionally returns a model of the dispersion in the second output
%   image.
%
% [delta_lambda, dispersion_max, bands] = chirpImage(...
%   image_sampling, lambda_range, d_fraction, n_samples, 'params'...
% )
%   Returns intermediate parameters used for image generation.
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
%   wavelength `lambda_range(2)` is shifted to the right relative to the
%   hypothetical image for wavelength `lambda_range(1)` by an amount
%   corresponding to the distance from the left border of the latter image
%   to the x-coordinate of its first intensity minimum.
%
% n_samples -- Number of samples per pixel
%   As far as I know, there is no analytical form for the integral of the
%   image intensity function. This integral would need to be taken over the
%   x, y, and spectral domain of a pixel in order to compute the pixel's
%   value. This function approximates the integral with the average of
%   `n_samples` random samples in the pixel's domain. If `n_samples` is
%   empty, or is one, then each pixel is given the value of the image
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
% dispersionfun -- Dispersion model
%   A function which takes an input 2D array 'xylambda', where the three
%   columns represent two spatial coordinates and a wavelength, (x, y,
%   lambda). The output of the function is a 2D array, 'disparity', where
%   the columns represent the x and y-components of dispersion vectors.
%
%   If `d_fraction` is zero, this output argument is empty (`[]`).
%
% delta_lambda -- Wavelength spacing
%   The distance between the centre wavelengths of consecutive spectral
%   bands in the output images, measured in nanometres (assuming
%   `lambda_range` is measured in nanometres).
%
% dispersion_max -- Magnitude of dispersion
%   The size of the offset between the hypothetical image for wavelength
%   `lambda_range(2)` and the hypothetical image for wavelength
%   `lambda_range(1)`, measured in pixels.
%
% bands -- Centre wavelengths
%   A column vector containing the centre wavelengths of the spectral bands
%   in the output images.
%
% ## References
% The chirp function is described in Chapter 2 of:
%
%   Goodman, J. W. (2005). Introduction to fourier optics (3rd ed.).
%     Englewood, Colorado: Roberts & Company.
%
% See also imageFormation, makeDispersionfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 20, 2018


narginchk(4, 5);
nargoutchk(1, 3);

output_params = ~isempty(varargin) && strcmp(varargin{1}, 'params');
if ~output_params
    output_warped = (nargout > 1);
end

image_height = image_sampling(1);
image_width = image_sampling(2);
n_bands = image_sampling(3);

alpha = pi / image_width;
lambda_0 = lambda_range(1);
lambda_1 = lambda_range(2);
lambda_span = diff(lambda_range);
delta_lambda = lambda_span / n_bands;
beta = pi / (image_height * delta_lambda);
dispersion_max = d_fraction * sqrt(image_width);
d_scaled = dispersion_max / lambda_span;
bands = linspace(lambda_0 + delta_lambda / 2, lambda_1 - delta_lambda / 2, n_bands).';

    function intensity = imageIntensity(x, y, lambda, d)
        lambda_rel = lambda - lambda_0;
        intensity = (cos(alpha .* (x - d .* lambda_rel) .^ 2) + 1) .* ...
            (cos(beta .* y .* lambda_rel) + 1) ./ 4;
    end

    function disparity = dispersionfun(xylambda)
        disparity = [...
            d_scaled * (xylambda(3) - lambda_0),...
            zeros(size(xylambda, 1), 1)...
        ];
    end

if output_params
    varargout{1} = delta_lambda;
    if nargout > 1
        varargout{2} = dispersion_max;
        if nargout > 2
            varargout{3} = bands;
        end
    end
else
    if isempty(n_samples) || n_samples == 1
        [Y, X, Lambda] = ndgrid(0.5:1:image_height, 0.5:1:image_width, bands);
        varargout{1} = imageIntensity(X, Y, Lambda, 0);
        if output_warped
            varargout{2} = imageIntensity(X, Y, Lambda, d_scaled);
        end
    else
        % Reference points at the corners of each pixel
        bands = linspace(lambda_0, lambda_1 - delta_lambda, n_bands);
        [Y, X, Lambda] = ndgrid(0:(image_height - 1), 0:(image_width - 1), bands);
        I_hyper = zeros(image_sampling);
        if output_warped
            I_warped = zeros(image_sampling);
        end
        for s = 1:n_samples
            Y_s = Y + rand(image_sampling);
            X_s = X + rand(image_sampling);
            Lambda_s = Lambda + delta_lambda * rand(image_sampling);
            I_hyper = I_hyper + imageIntensity(X_s, Y_s, Lambda_s, 0);
            if output_warped
                I_warped = I_warped + imageIntensity(X_s, Y_s, Lambda_s, d_scaled);
            end
        end
        varargout{1} = I_hyper / n_samples;
        if output_warped
            varargout{2} = I_warped / n_samples;
        end
    end
    
    if nargout > 2
        if d_fraction > 0
            varargout{3} = dispersionfun;
        else
            varargout{3} = [];
        end
    end
end

end
