function [e_rgb, varargout] = evaluateRGB(...
    I_rgb, R_rgb, options...
)
% EVALUATERGB  Compare colour images
%
% ## Syntax
% e_rgb = evaluateRGB(I_rgb, R_rgb, options)
% [e_rgb, fg_rgb] = evaluateRGB(I_rgb, R_rgb, options)
%
% ## Description
% e_rgb = evaluateRGB(I_rgb, R_rgb, options)
%   Returns a structure containing the values of different evaluation
%   metrics comparing the two colour images, and generates graphical
%   output, if requested.
%
% [e_rgb, fg_rgb] = evaluateRGB(I_rgb, R_rgb, options)
%   Additionally returns a structure of figure handles from the colour
%   image evaluation graphical output.
%
% ## Input Arguments
%
% I_rgb -- Test colour image
%   An h x w x 3 array containing an estimated colour image.
%
% R_rgb -- Reference colour image
%   An h x w x 3 array containing the ideal/true colour image.
%
% options -- Graphical output options for colour images
%   A structure which enables or disables graphical output to figures
%   relating to the input colour images. `options` has the following
%   fields:
%   - 'error_map': If `true`, then this function will produce figures
%     showing the absolute error between the two colour images. One figure
%     will be produced per channel.
%
% ## Output Arguments
%
% e_rgb -- Colour error statistics
%   A structure with the following fields:
%   - 'mse': The mean square error between the two images. 'mse' is a
%     three-element vector, with one element per colour channel.
%   - 'psnr': The peak signal-to-noise ratio between the two images, in the
%     same format as 'mse'.
%   - 'cpsnr': A scalar storing the CPSNR value. CPSNR is an extension of
%     PSNR to colour images, computed by taking the mean square error over
%     all colour channels instead of one.
%   - 'ssim': The Structural Similarity Index Measure computed between the
%     two images. 'ssim' is a four-element vector, where the first three
%     elements are the SSIM values for individual colour channels, and the
%     last is their average.
%   - 'mi_within': A 3 x 2 array, where the first column contains the
%     mutual information between each pair of channels in the reference
%     image, and the second column contains the mutual information between
%     each pair of channels in the test image. The ordering of channel
%     pairs is Red-Green, Green-Blue, and Red-Blue.
%   - 'mi_betweeen': A 3-element vector, containing the mutual information
%     between the channnels of the test image and the reference image.
%
% fg_rgb -- Colour error evaluation figures
%   A structure with the following fields, all of which store figure
%   handles:
%   - 'error_map': A vector of three figure handles corresponding to the
%     output triggered by `options.error_map`.
%
% ## Notes
% - A border of `border` pixels (a local variable defined in the code) is
%   excluded from the images when calculating image-wide statistics, to
%   ignore artifacts in image estimation arising from the image borders.
% - Figures will not be generated if the corresponding fields of
%   `options` are missing.
% - Figures are produced with titles and axis labels, but without legends.
% - Images can be input in integer or floating-point formats. For peak
%   signal-to-noise ratio calculations, the peak value will be 1.0 for
%   floating-point formats, and the maximum possible positive integer value
%   for integer formats.
%
% ## References
% - Image borders are excluded from image similarity measurements in the
%   image demosaicking literature, such as in:
%
%   Monno, Y., Kiku, D., Tanaka, M., & Okutomi, M. (2017). "Adaptive
%     residual interpolation for color and multispectral image
%     demosaicking." Sensors (Switzerland), 17(12). doi:10.3390/s17122787
%
% - The code for calculating mutual information was retrieved from MATLAB
%   Central,
%   https://www.mathworks.com/matlabcentral/fileexchange/36538-very-fast-mutual-information-betweentwo-images
%   The function 'third_party/MI_GG/MI_GG.m' was written by Generoso
%   Giangregorio, and corresponds to the following article:
%
%   M. Ceccarelli, M. di Bisceglie, C. Galdi, G. Giangregorio, S.L. Ullo,
%     "Image Registration Using Non–Linear Diffusion", IGARSS 2008.
%
% - The idea of using mutual information to evaluate image alignment is
%   mentioned in, among other articles,
%
%   Brauers, J., Schulte, B., & Aach, T. (2008). "Multispectral
%     Filter-Wheel Cameras: Geometric Distortion Model and Compensation
%     Algorithms." IEEE Transactions on Image Processing, 17(12),
%     2368-2380. doi:10.1109/TIP.2008.2006605
%
% See also immse, psnr, ssim, evaluateSpectral

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 14, 2018

    function psnr = peakSignalToNoiseRatio(mse, peak)
        psnr = 10 * log10((peak ^ 2) / mse);
    end

narginchk(3, 3);
nargoutchk(1, 2);

border = 10;

class_rgb = class(I_rgb);
if ~isa(R_rgb, class_rgb)
    error('The two colour images do not have the same datatype.')
end
if isinteger(I_rgb)
    peak_rgb = intmax(class_rgb);
elseif isfloat(I_rgb)
    peak_rgb = 1;
else
    error('The two colour images are of an unexpected datatype.')
end

I_rgb_clipped = I_rgb((border + 1):(end - border), (border + 1):(end - border), :);
R_rgb_clipped = R_rgb((border + 1):(end - border), (border + 1):(end - border), :);

n_channels_rgb = 3;
e_rgb.mse = zeros(n_channels_rgb, 1);
e_rgb.psnr = zeros(n_channels_rgb, 1);
e_rgb.ssim = zeros(n_channels_rgb + 1, 1);
for c = 1:n_channels_rgb
    e_rgb.mse(c) = immse(I_rgb_clipped(:, :, c), R_rgb_clipped(:, :, c));
    e_rgb.psnr(c) = peakSignalToNoiseRatio(e_rgb.mse(c), peak_rgb);
    e_rgb.ssim(c) = ssim(I_rgb_clipped(:, :, c), R_rgb_clipped(:, :, c));
end
e_rgb.cpsnr = peakSignalToNoiseRatio(mean(e_rgb.mse), peak_rgb);
e_rgb.ssim(end) = mean(e_rgb.ssim(1:(end - 1)));

mi_class = 'uint8';
if ~isa(I_rgb, mi_class)
    peak_int = double(intmax(mi_class));
    I_rgb_clipped_int = uint8(I_rgb_clipped * peak_int / peak_rgb);
    R_rgb_clipped_int = uint8(R_rgb_clipped * peak_int / peak_rgb);
else
    I_rgb_clipped_int = I_rgb_clipped;
    R_rgb_clipped_int = R_rgb_clipped;
end

e_rgb.mi_within = zeros(n_channels_rgb, 2);
e_rgb.mi_between = zeros(n_channels_rgb, 1);
for c = 1:n_channels_rgb
    e_rgb.mi_within(c, 1) = MI_GG(...
        R_rgb_clipped_int(:, :, c),...
        R_rgb_clipped_int(:, :, mod(c, n_channels_rgb) + 1)...
    );
    e_rgb.mi_within(c, 2) = MI_GG(...
        I_rgb_clipped_int(:, :, c),...
        I_rgb_clipped_int(:, :, mod(c, n_channels_rgb) + 1)...
    );
    e_rgb.mi_between(c) = MI_GG(...
        I_rgb_clipped_int(:, :, c),...
        R_rgb_clipped_int(:, :, c)...
    );
end

fg_rgb = struct;
if isfield(options, 'error_map') && options.error_map
    diff_rgb = I_rgb - R_rgb;
    for c = 1:n_channels_rgb
        fg_rgb.error_map(c) = figure;
        imagesc(diff_rgb(:, :, c));
        colorbar;
        title(sprintf('Difference image for channel %d', c));
    end
end

if nargout > 1
    varargout{1} = fg_rgb;
end
    
end
