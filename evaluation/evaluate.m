function [e_rgb, varargout] = evaluate(...
    I_rgb, R_rgb, options_rgb, varargin...
)
% EVALUATE  Compare colour and spectral images
%
% ## Syntax
% e_rgb = evaluate(I_rgb, R_rgb, options_rgb)
% [e_rgb, fg_rgb] = evaluate(I_rgb, R_rgb, options_rgb)
% [e_rgb, fg_rgb, e_spectral] = evaluate(...
%     I_rgb, R_rgb, options_rgb,...
%     I_spectral, R_spectral, lambda, options_spectral...
% )
% [e_rgb, fg_rgb, e_spectral, fg_spectral] = evaluate(...
%     I_rgb, R_rgb, options_rgb,...
%     I_spectral, R_spectral, lambda, options_spectral...
% )
%
% ## Description
% e_rgb = evaluate(I_rgb, R_rgb, options_rgb)
%   Returns a structure containing the values of different evaluation
%   metrics comparing the two colour images, and generates graphical
%   output, if requested.
%
% [e_rgb, fg_rgb] = evaluate(I_rgb, R_rgb, options_rgb)
%   Additionally returns a structure of figure handles from the colour
%   image evaluation graphical output.
%
% [e_rgb, e_spectral] = evaluate(...
%     I_rgb, R_rgb, options_rgb,...
%     I_spectral, R_spectral, lambda, options_spectral...
% )
%   Additionally returns a structure containing the values of different
%   evaluation metrics comparing the two spectral images, and generates
%   corresponding graphical output, if requested.
%
% [e_rgb, fg_rgb, e_spectral, fg_spectral] = evaluate(...
%     I_rgb, R_rgb, options_rgb,...
%     I_spectral, R_spectral, lambda, options_spectral...
% )
%   Additionally returns a structure of figure handles from the spectral
%   image evaluation graphical output.
%
% ## Input Arguments
%
% I_rgb -- Test colour image
%   An h1 x w1 x 3 array containing an estimated colour image.
%
% R_rgb -- Reference colour image
%   An h1 x w1 x 3 array containing the ideal/true colour image.
%
% options_rgb -- Graphical output options for colour images
%   A structure which enables or disables graphical output to figures
%   relating to the input colour images. `options_rgb` has the following
%   fields:
%   - 'error_map': If `true`, then this function will produce figures
%     showing the absolute error between the two colour images. One figure
%     will be produced per channel.
%
% I_spectral -- Test spectral image
%   An h2 x w2 x c array containing an estimated spectral image.
%
% R_spectral -- Reference spectral image
%   An h2 x w2 x c array containing the ideal/true spectral image.
%
% lambda -- Wavelengths
%   A vector of length 'c' containing the wavelengths corresponding to the
%   third dimension of `I_spectral` and `R_spectral`.
%
% options_spectral -- Output options for spectral images
%   A structure which controls graphical output to figures relating to
%   the spectral images, and also determines some of the output in
%   `e_spectral`. `options_spectral` has the following fields:
%   - 'error_map': If `true`, then this function will produce a figure
%     showing the absolute error between the two images. The figure will be
%     produced for the band having the highest mean squared error.
%   - 'radiance': A cell vector, where each vector is a four-element vector
%     describing an image patch (center pixel x-coordinate, center pixel
%     y-coordinate, width, and height). For each image patch, a figure will
%     be generated containing plots of spectral information from the two
%     images, averaged within the image patch.
%   - 'radiance_fg': A vector of figure handles. The output produced for
%     'radiance' will be added to these figures, instead of being shown in
%     new figures. It is assumed that the existing figures already have
%     plotlines for the reference image, so only the test image's
%     information will be added.
%   - 'scanlines': A cell vector, where each vector is a four-element vector
%     describing a line segment (start pixel x-coordinate, start pixel
%     y-coordinate, end pixel x-coordinate, and end pixel y-coordinate).
%     For each line segment, two figures will be generated, containing a
%     plot of spectral root mean square error along the line, and a plot of
%     spectral goodness-of-fit coefficient along the line, respectively.
%   - 'scanlines_fg': A structure vector where each element has fields
%     'rmse' and 'gof' storing the root mean square error and
%     goodness-of-fit coefficient figure handles, respectively.
%     'scanlines_fg' has the same role with respecct to 'scanlines' as
%     'radiance_fg' has for 'radiance'.
%   - 'mi_bands': A two-element vector containing the indices of the
%     spectral bands to compute mutual information between.
%   - 'bands_diff': A two-element vector containing the indices of the
%     spectral bands to create difference images between (to show
%     dispersion between wavelengths). The band with the first index will
%     be subtracted from the band with the second index. Two figures will
%     be created, one for the reference image, and one for the test image.
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
%     all colour channels.
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
%     output triggered by `options_rgb.error_map`.
%
% e_spectral -- Spectral error statistics
%   A structure with the following fields:
%   - 'mse': The mean square error between the two images. 'mse' is a
%     structure with the following fields:
%     - 'max': The mean square error for the band with the highest mean
%       square error.
%     - 'mean': The average mean square error over all bands
%     - 'median': The mean square error for the band with the median mean
%       square error.
%     - 'raw': A c-element vector, containing the mean square errors for
%       each band.
%   - 'psnr': The peak signal-to-noise ratio between the two images, in the
%     same format as 'mse', except with 'max' replaced with 'min' (as PSNR
%     decreases when mean square error increases).
%   - 'ssim': The Structural Similarity Index Measure computed between the
%     two images. 'ssim' is a structure with the following fields:
%     - 'min': The SSIM value for the band with the lowest SSIM value.
%     - 'mean': The average SSIM value over all bands
%     - 'median': The SSIM value for the band with the median mean SSIM value.
%     - 'raw': A c-element vector, containing the SSIM values for each band.
%   - 'mi_within': A two element vector, containing the mutual information
%     between two bands in the reference image, and between two bands in
%     the test image, respectively. The indices of the bands are given in
%     `options_spectral.mi_bands`.
%   - 'mi_between': A c-element vector, containing the mutual information
%     between corresponding bands in the reference image and the test
%     image.
%   - 'rmse': The root mean square error between the spectral information
%     at pixels in the test image and pixels in the reference image. 'rmse'
%     is a structure with the following fields:
%     - 'max': The highest RMSE value in the image
%     - 'mean': The average RMSE value over all pixels
%     - 'median': The median RMSE value over all pixels
%   - 'gof': The goodness-of-fit coefficient between the spectral
%     information of pixels in the test image and in the reference image.
%     'gof' is a structure with the same format as 'rmse', but with 'max'
%     replaced with 'min'.
%   - 'radiance': A structure vector, with the same length as
%     `options_spectral.radiance`. Each element contains metrics describing
%     the image patch given by the corresponding element of
%     `options_spectral.radiance`, and has the following fields:
%     - 'rmse': The root mean square error of the average spectral radiance
%       in the test patch relative to the average spectral radiance in the
%       reference patch
%     - 'gof': The goodness-of-fit coefficient of the average spectral
%       radiance in the test patch relative to the average spectral
%       radiance in the reference patch.
%
% fg_spectral -- Spectral error evaluation figures
%   A structure with the following fields, all of which store figure
%   handles:
%   - 'error_map': A figure handle corresponding to the output triggered
%     by `options_spectral.error_map`.
%   - 'radiance': A vector of figure handles corresponding to the output
%     triggered by `options_spectral.radiance`.
%   - 'scanlines': A structure vector corresponding to the output triggered
%     by `options_spectral.scanlines`. Each structure has fields 'rmse' and
%     'gof' storing the root mean square error and goodness-of-fit
%     coefficient figure handles, respectively.
%   - 'bands_diff': A two-element vector, containing figure handles for the
%     difference images between spectral bands. The first element is for
%     the reference image, whereas the second element is for the test
%     image. The spectral bands being subtracted are given by
%     `options_spectral.bands_diff`.
%
% ## Notes
% - A border of `border` pixels (a local variable defined in the code) is
%   excluded from the images when calculating image-wide statistics, to
%   ignore artifacts in image estimation arising from the image borders.
% - Figures are produced with titles and axis labels, but without legends.
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
% - The goodness-of-fit coefficient is used as a spectral error measure
%   (Equation 15) in:
%
%   Nguyen R.M.H., Prasad D.K., Brown M.S. (2014) Training-Based Spectral
%     Reconstruction from a Single RGB Image. In: Fleet D., Pajdla T.,
%     Schiele B., Tuytelaars T. (eds) Computer Vision – ECCV 2014. ECCV
%     2014. Lecture Notes in Computer Science, vol 8695. Springer, Cham
%
% See also immse, psnr, ssim

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 14, 2018

narginchk(3, 7);
if nargin == 3 && nargout > 2
    error('Spectral evaluation cannot be requested as an output argument if spectral input arguments are not provided.');
elseif nargin < 7 && nargout > 2
    error('Spectral evaluation cannot be requested as an output argument if some spectral input arguments are not provided.');
end
    
end
