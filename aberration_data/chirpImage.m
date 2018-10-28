function varargout = chirpImage(image_sampling, lambda_range, varargin)
% CHIRPIMAGE  Create an image with spatial and spectral linear chirp, and spectral dispersion
%
% ## Syntax
% I_hyper = chirpImage(...
%   image_sampling, lambda_range, dispersion_mag, n_samples,...
%   type [, sensitivity, threshold, verbose]...
% )
% [I_hyper, I_warped] = chirpImage(...
%   image_sampling, lambda_range, dispersion_mag, n_samples,...
%   type [, sensitivity, threshold, verbose]...
% )
% [I_hyper, I_warped, dispersionfun] = chirpImage(...
%   image_sampling, lambda_range, dispersion_mag, n_samples,...
%   type [, sensitivity, threshold, verbose]...
% )
% [delta_lambda, dispersion_max, bands] = chirpImage(...
%   image_sampling, lambda_range, 'params'...
% )
%
% ## Description
% I_hyper = chirpImage(...
%   image_sampling, lambda_range, dispersion_mag, n_samples,...
%   type [, sensitivity, threshold, verbose]...
% )
%   Returns a version of the hyperspectral image without dispersion
%
% [I_hyper, I_warped] = chirpImage(...
%   image_sampling, lambda_range, dispersion_mag, n_samples,...
%   type [, sensitivity, threshold, verbose]...
% )
%   Additionally returns a version of the hyperspectral image with
%   dispersion
%
% [I_hyper, I_warped, dispersionfun] = chirpImage(...
%   image_sampling, lambda_range, dispersion_mag, n_samples,...
%   type [, sensitivity, threshold, verbose]...
% )
%   Additionally returns a function model of the dispersion in the second
%   output image.
%
% [delta_lambda, dispersion_max, bands] = chirpImage(...
%   image_sampling, lambda_range, 'params'...
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
% dispersion_mag -- Amount of dispersion
%   The amount of dispersion to add to `I_warped` relative to `I_hyper`,
%   measured in pixels. `dispersion_mag` is the offset between the
%   hypothetical images for wavelength `lambda_range(2)` and wavelength
%   `lambda_range(1)`. The image for `lambda_range(2)` is shifted to the
%   right with dispersion. Note that the sign of `dispersion_mag` is
%   ignored.
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
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel to the j-th spectral band. `sensitivity` is a matrix
%   with `image_sampling(3)` columns mapping the hyperspectral
%   representation of the image to a colour image. `sensitivity` and
%   `threshold` must both be passed together, and are used to determine the
%   maximum spectral frequency in the image.
%
% threshold -- Power spectrum threshold
%   In the context of hyperspectral reconstruction, a colour image produced
%   by a camera described by `sensitivity` can be seen as a convolution of
%   the camera's spectral sensitivity (`sensitivity`) with the spectral
%   signals of each pixel. Such an interpretation is valid especially for
%   the 'phase' value of the `type` input argument, as the change in phase
%   of the spectral signals over the image in the horizontal direction can
%   be seen as sampling the convolution at all possible wavelengths.
%   Therefore, reasoning about hyperspectral reconstruction from the
%   perspective of sampling theory, one cannot expect to recover
%   frequencies in the spectral signals of pixels which are outside the
%   bandlimits of the camera's spectral sensitivity functions.
%
%   `threshold` is the fraction of the cumulative power spectrum that the
%   user considers to define the bandlimit (because no signal that is
%   time-limited is frequency-limited). The output image will be
%   constrained such that its maximum frequency in the spectral domain (the
%   third dimension of the output image) is limited to the bandlimit of the
%   spectral sensitivity function with the highest bandlimit.
%
% verbose -- Debugging flag
%   A Boolean scalar which enables graphical debugging output if `true`.
%   Defaults to `false` if not passed.
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
%   x-coordinate, until the rate reaches a maximum of one cycle per two
%   pixels at the right edge of the image.
%
%   As the y-coordinate increases, the frequency of oscillation in the
%   spectral dimension increases. If `sensitivity` and `threshold` are not
%   passed, the frequency reaches a maximum of one cycle per two spectral
%   bands at the bottom edge of the image. Otherwise, if `sensitivity` and
%   `threshold` are passed, and if the bandlimit of the spectral
%   sensitivity function in `sensitivity` with the highest bandlimit is
%   smaller than one cycle per two spectral bands, then this bandlimit is
%   the maximum frequency in the spectral dimension.
%
%   The maximum frequencies of spatial and spectral variation are intended
%   to respect the Whittaker-Shannon sampling theorem (as described by
%   Goodman 2005, among other references).
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
% delta_lambda -- Wavelength spacing
%   The distance between the centre wavelengths of consecutive spectral
%   bands in the output images, measured in nanometres (assuming
%   `lambda_range` is measured in nanometres).
%
% dispersion_max -- Magnitude of maximum dispersion
%   The size of the offset between the hypothetical images for wavelength
%   `lambda_range(2)` and wavelength `lambda_range(1)`, measured in pixels,
%   at the suggested maximum dispersion. Under the maximum dispersion, the
%   hypothetical image for wavelength `lambda_range(2)` is shifted to the
%   right relative to the hypothetical image for wavelength
%   `lambda_range(1)`, by an amount equal to the distance from the left
%   border of the latter image to the x-coordinate where it reaches a
%   half-cycle change in phase. As such, the images for the endpoints of
%   the spectral range are misaligned by half the width of the central
%   maximum (which extends to both sides of the image's left border).
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

narginchk(3, 8);
nargoutchk(1, 3);

verbose = false;
output_params = strcmp(varargin{1}, 'params');
if output_params
    narginchk(3, 3);
else
    narginchk(5, 8);
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
    clip_spectral_frequency = false;
    if length(varargin) > 3
        clip_spectral_frequency = true;
        sensitivity = varargin{4};
        if size(sensitivity, 2) ~= image_sampling(3)
            error('The number of spectral bands must be consistent between `image_sampling` and `sensitivity`.');
        end
        if length(varargin) < 5
            error('Both `sensitivity` and `threshold` must be passed, or neither must be passed.');
        end
        threshold = varargin{5};
        if length(varargin) > 5
            verbose = varargin{6};
        end
    end
end

image_height = image_sampling(1);
image_width = image_sampling(2);
n_bands = image_sampling(3);

alpha = pi / (2 * image_width);
lambda_0 = lambda_range(1);
lambda_1 = lambda_range(2);
lambda_span = diff(lambda_range);
delta_lambda = lambda_span / n_bands;

bands = linspace(lambda_0 + delta_lambda / 2, lambda_1 - delta_lambda / 2, n_bands).';

if ~output_params
    d_scaled = dispersion_mag / lambda_span;
    if clip_spectral_frequency
        [limits_freq, delta_theta_max] = bandlimit(sensitivity, threshold, verbose);
        if verbose
            fprintf('Maximum spectral frequency will be %g [cycles/(band index increment)].\n', limits_freq);
        end
    else
        delta_theta_max = Inf;
    end
    beta = min(pi, delta_theta_max) / (delta_lambda * image_height);
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
    varargout{1} = delta_lambda;
    if nargout > 1
        varargout{2} = sqrt(2 * image_width);
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
        if dispersion_mag == 0
            varargout{3} = [];
        else
            varargout{3} = @dispersionfun;
        end
    end
end

end
