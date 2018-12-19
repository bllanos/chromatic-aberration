function varargout = chirpImage(image_sampling, max_freq, lambda, varargin)
% CHIRPIMAGE  Create an image with spatial and spectral linear chirp, and spectral dispersion
%
% ## Syntax
% I_hyper = chirpImage(...
%   image_sampling, max_freq, lambda, dispersion_mag, n_samples, type...
% )
% [I_hyper, I_warped] = chirpImage(...
%   image_sampling, max_freq, lambda, dispersion_mag, n_samples, type...
% )
% [I_hyper, I_warped, dispersionfun] = chirpImage(...
%   image_sampling, max_freq, lambda, dispersion_mag, n_samples, type...
% )
% dispersion_max = chirpImage(...
%   image_sampling, max_freq, lambda, 'params'...
% )
%
% ## Description
% I_hyper = chirpImage(...
%   image_sampling, max_freq, lambda, dispersion_mag, n_samples, type...
% )
%   Returns a version of the hyperspectral image without dispersion
%
% [I_hyper, I_warped] = chirpImage(...
%   image_sampling, max_freq, lambda, dispersion_mag, n_samples, type...
% )
%   Additionally returns a version of the hyperspectral image with
%   dispersion
%
% [I_hyper, I_warped, dispersionfun] = chirpImage(...
%   image_sampling, max_freq, lambda, dispersion_mag, n_samples, type...
% )
%   Additionally returns a function model of the dispersion in the second
%   output image.
%
% dispersion_max = chirpImage(...
%   image_sampling, max_freq, lambda, 'params'...
% )
%   Returns intermediate parameters used for image generation.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the height and width, respectively, of
%   the image.
%
% max_freq -- Maximum spatial and spectral frequencies
%   A two-element vector containing the maximum frequencies of variation in
%   the spatial and spectral dimensions of the image content. The first
%   element is the maximum number of cycles in the spatial variation of
%   image content per displacement of one pixel (in any direction). If this
%   value is greater than `0.5`, aliasing will occur. The second element is
%   the maximum number of cycles in the spectral variation of image content
%   per change in wavelength of one unit (the same unit as used for the
%   values in `lambda`). If this value is greater than `0.5`, divided by
%   the spacing between values in `lambda`, aliasing will occur.
%
%   The suggested maximum values based on `0.5` are intended to respect the
%   Whittaker-Shannon sampling theorem (as described by Goodman 2005, among
%   other references).
%
% lambda -- Sampling wavelengths
%   A vector containing the wavelengths captured in the image. The
%   wavelengths must be evenly-spaced.
%
% dispersion_mag -- Amount of dispersion
%   The amount of dispersion to add to `I_warped` relative to `I_hyper`,
%   measured in pixels. If the spacing between consecutive elements of
%   `lambda` is `delta_lambda`, then `dispersion_mag` is the offset between
%   the hypothetical images for wavelength `lambda(1) - delta_lambda / 2`
%   and wavelength `lambda(end) + delta_lambda / 2`. The image for the
%   latter wavelength is shifted to the right with dispersion. Note that
%   the sign of `dispersion_mag` is ignored.
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
% type -- Subtype of image
%   A character vector describing how the image should change in the
%   x-dimension:
%   - 'phase': The phase of the spectral signal at each pixel advances
%     quadratically with the image x-coordinate. As such, for a
%     sufficiently high spectral frequency, the intensity of image pixels
%     (taken across all spectral bands) is constant with x, even though the
%     amplitude varies within a given spectral band.
%   - 'blend-black': The spectral signal at each pixel has a constant phase
%     with respect to x, but is blended with black. The blending weight is
%     a quadratic-phase cosine function of x. As such, the intensity of
%     image pixels is a chirp function of x.
%   - 'blend-negative': Similar to 'blend-black', but the spectral signal
%     is blended with a version of itself which has the opposite phase. As
%     such, the intensity of image pixels is a chirp function of x, but is
%     more uniform than for 'blend-black'. The image appears to alternate
%     between two colours as x changes.
%   - 'blend-black-sharp': Similar to 'blend-black', but the blending
%     factor is a binarized version of the chirp function of x. As such,
%     the image intensity is a series of vertical sharp-edged bands.
%   - 'blend-negative-sharp': The equivalent of 'blend-black-sharp' for
%     blending with an opposite-phase signal, as in 'blend-negative'.
%
% ## Output Arguments
%
% I_hyper -- Chirp image
%   An image produced by evaluating a linear chirp function. The spectral
%   curve at each pixel is a cosine function. As the x-coordinate
%   increases, either the phase of the spectral signal increases, or the
%   blending factor with black or with an opposite phase cosine function
%   changes, as described in the documentation of `type` above. The rate at
%   which the phase or blending factor changes increases linearly with the
%   x-coordinate, until the rate reaches a maximum at the right edge of the
%   image (given by `max_freq(1)`).
%
%   As the y-coordinate increases, the frequency of oscillation in the
%   spectral dimension increases. The frequency reaches a maximum at the
%   bottom edge of the image (given by `max_freq(2)`). The caller is
%   responsible for ensuring that this maximum frequency is within the
%   effective bandlimit of the camera spectral response functions that will
%   be used to convert the image to colour.
%
% I_warped -- Dispersed image
%   A version of `I_hyper` altered by spectral dispersion, with the
%   magnitude of the dispersion given by `dispersion_mag`. The dispersion
%   is a shift in the positive x-direction with wavelength. Note that the
%   dispersed image is synthesized directly, rather than created by warping
%   `I_hyper`.
%
% dispersionfun -- Function model of dispersion
%   A function which takes an input 2D array 'xylambda', where the three
%   columns represent two spatial coordinates and a wavelength, (x, y,
%   lambda). The output of the function is a 2D array, 'disparity', where
%   the columns represent the x and y-components of dispersion vectors.
%
%   If `dispersion_mag` is zero, this output argument is empty (`[]`).
%
% dispersion_max -- Magnitude of maximum dispersion
%   The suggested maximum value of `dispersion_mag`. If the spacing between
%   consecutive elements of `lambda` is `delta_lambda`, then
%   `dispersion_max` is the offset between the hypothetical images for
%   wavelength `lambda(1) - delta_lambda / 2` and wavelength `lambda(end) +
%   delta_lambda / 2`. `dispersion_max` is equal to the distance from the
%   left border of the former image to the x-coordinate where it reaches a
%   half-cycle change in phase. As such, the images for the endpoints of
%   the spectral range are misaligned by half the width of the central
%   maximum (which extends to both sides of the image's left border).
%
% ## Notes
% - This function uses its own random number generator state, then restores
%   the state of the random number generator before returning. This allows
%   multiple calls to this function, using the same parameters, to produce
%   exactly the same images.
%
% ## References
% The chirp function is described in Chapter 2 of:
%
%   Goodman, J. W. (2005). Introduction to fourier optics (3rd ed.).
%     Englewood, Colorado: Roberts & Company.
%
% See also findSampling, bandlimit, imageFormation, makeDispersionfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 20, 2018

narginchk(4, 6);

output_params = strcmp(varargin{1}, 'params');
if output_params
    nargoutchk(1, 1);
    narginchk(4, 4);
else
    nargoutchk(1, 3);
    narginchk(6, 6);
    dispersion_mag = varargin{1};
    dispersion_mag = abs(dispersion_mag);
    n_samples = varargin{2};
    type = varargin{3};
    output_warped = (nargout > 1);
    do_blend = true;
    do_blend_black = false;
    do_sharp_edges = false;
    if strcmp(type, 'phase')
        do_blend = false;
    elseif strcmp(type, 'blend-black')
        do_blend_black = true;
    elseif strcmp(type, 'blend-negative')
    elseif strcmp(type, 'blend-black-sharp')
        do_blend_black = true;
        do_sharp_edges = true;
    elseif strcmp(type, 'blend-negative-sharp')
        do_sharp_edges = true;
    else
        error('Unrecognized value "%s" of the `type` input argument.', type);
    end
end

diff_lambda = diff(lambda);
delta_lambda = diff_lambda(1);
if max(abs(diff_lambda - delta_lambda)) > 1e-6
    error('`lambda` must contain equally-spaced values.')
end

image_height = image_sampling(1);
image_width = image_sampling(2);
n_bands = length(lambda);

alpha = (pi * max_freq(1)) / image_width;
lambda_0 = lambda(1) - (delta_lambda / 2);
lambda_1 = lambda(end) + (delta_lambda / 2);
lambda_span = lambda_1 - lambda_0;

if ~output_params
    d_scaled = dispersion_mag / lambda_span;
    beta = (2 * pi * max_freq(2)) / image_height;
end

    function intensity = imageIntensity(x, y, lambda, d)
        lambda_rel = lambda - lambda_0;
        if do_blend
            blend_factor = cos(alpha .* (x - d .* lambda_rel) .^ 2);
            if do_sharp_edges
                blend_factor = double(blend_factor > 0);
            else
                blend_factor = (blend_factor + 1) / 2;
            end
            intensity = cos(beta .* y .* lambda_rel);
            if do_blend_black
                intensity = blend_factor .* (intensity + 1) / 2;
            else
                intensity = (1 - intensity + 2 * blend_factor .* intensity) / 2;
            end
        else
            intensity = (cos(...
                    (beta .* y .* lambda_rel) +...
                    (alpha .* (x - d .* lambda_rel) .^ 2)...
                ) + 1) ./ 2;
        end
    end

    function disparity = dispersionfun(xylambda)
        disparity = [...
            d_scaled * (xylambda(:, 3) - lambda_0),...
            zeros(size(xylambda, 1), 1)...
        ];
    end

if output_params
    varargout{1} = sqrt(2 * image_width);
else
    if isempty(n_samples) || n_samples == 1
        [Y, X, Lambda] = ndgrid(0.5:1:image_height, 0.5:1:image_width, lambda);
        varargout{1} = imageIntensity(X, Y, Lambda, 0);
        if output_warped
            varargout{2} = imageIntensity(X, Y, Lambda, d_scaled);
        end
    else
        % Reference points at the corners of each pixel
        lambda_corners = linspace(lambda_0, lambda_1 - delta_lambda, n_bands);
        [Y, X, Lambda] = ndgrid(0:(image_height - 1), 0:(image_width - 1), lambda_corners);
        I_hyper = zeros(image_sampling);
        if output_warped
            I_warped = zeros(image_sampling);
        end
        
        % Use reproducible random numbers, allowing this image to be
        % duplicated
        random_saved = rng;
        rng(1, 'twister');
        for s = 1:n_samples
            Y_s = Y + rand(image_sampling);
            X_s = X + rand(image_sampling);
            Lambda_s = Lambda + delta_lambda * rand(image_sampling);
            I_hyper = I_hyper + imageIntensity(X_s, Y_s, Lambda_s, 0);
            if output_warped
                I_warped = I_warped + imageIntensity(X_s, Y_s, Lambda_s, d_scaled);
            end
        end
        rng(random_saved);
        
        varargout{1} = I_hyper / n_samples;
        if output_warped
            varargout{2} = I_warped / n_samples;
        end
    end
    
    if nargout > 2
        if dispersion_mag == 0
            varargout{3} = [];
        else
            varargout{3} = @dispersionfun;
        end
    end
end

end
