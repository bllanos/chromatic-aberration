function [M_Omega_Phi, image_bounds, Phi, Omega, M] = projectionMatrix(...
    image_sampling_I, align, dispersion, sensitivity, lambda,...
    image_sampling_J, int_method, add_border...
)
% PROJECTIONMATRIX  Create matrices for latent image to RAW image conversion
%
% ## Syntax
% [M_Omega_Phi, image_bounds, Phi, Omega, M] = projectionMatrix(...
%     image_sampling_I, align, dispersion, sensitivity, lambda,...
%     image_sampling_J, int_method, add_border...
% )
%
% ## Description
% [M_Omega_Phi, image_bounds, Phi, Omega, M] = projectionMatrix(...
%     image_sampling_I, align, dispersion, sensitivity, lambda,...
%     image_sampling_J, int_method, add_border...
% )
%   Returns matrices for simulating the image formation process. One to
%   five output arguments can be requested.
%
% ## Input Arguments
%
% image_sampling_I -- Latent image dimensions
%   A two-element vector containing the height and width, respectively, of
%   the latent image `I`.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the RAW image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% dispersion -- Model of dispersion
%   `dispersion` can be empty (`[]`), if there is no model of dispersion.
%   Otherwise, two forms of this argument can be passed:
%
%   `dispersion` can be a function handle, produced by
%   'makeDispersionfun()'. `dispersion(X)`, where `X` is a three-element
%   row vector (x, y, lambda), returns the dispersion vector for the
%   position (x, y) in `J` corresponding to light with wavelength or colour
%   channel index `lambda`. The dispersion vector corrects for lateral
%   chromatic aberration by pointing from the corresponding position in the
%   reference spectral band or colour channel to position (x, y). This
%   function will negate the dispersion vectors of `dispersion` in order to
%   create a warp matrix from `I` to `J`.
%
%   `dispersion` can be a matrix for warping `I` to `J`, where the k-th row
%   contains the weights of pixels in `I` used to re-estimate the k-th
%   pixel in `J`. With a matrix value for `dispersion`, this function will
%   not make use of `lambda` or `add_border`, which are only needed for
%   computing a warp matrix.
%
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th colour channel or spectral band of
%   `I`. `sensitivity` is a matrix mapping colours in `I` to colours in
%   `J`.
%
% lambda -- Wavelength bands
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersion`. 'c' is the desired number of spectral bands or colour
%   channels in `I`.
%
% image_sampling_J -- RAW image dimensions
%   A two-element vector containing the height and width, respectively, of
%   the RAW image `J`.
%
% int_method -- Colour conversion integration method
%   The numerical integration method used for latent to colour space
%   conversion. `int_method` is passed to 'channelConversionMatrix()' as
%   its `int_method` argument. Refer to the documentation of
%   'channelConversionMatrix.m' for details. If 'int_method' is 'none', as
%   should be the case when colour conversion is from a set of colour
%   channels, not a spectral space, numerical integration will not be
%   performed.
%
% add_border -- Padding flag for latent image
%   A Boolean value indicating whether or not the `image_bounds` input
%   argument of 'dispersionfunToMatrix()' should be empty. If `true`, the
%   latent image `I` will be large enough to contain the
%   dispersion-corrected coordinates of all pixels in `J`. If `false`, the
%   latent image `I` will be clipped to the region occupied by `J`. If
%   `dispersion` is empty, or is a matrix, this argument is not used, and
%   can be empty.
%
% ## Output Arguments
%
% The array dimensions mentioned in this section are defined as follows:
% - `n_px_J` is the number of pixels in the RAW input image, `J`.
%   `n_px_J = prod(image_sampling_J)`.
% - `c` is the number of colour channels or spectral bands in the latent
%   image, `I`. `c = length(lambda)`.
% - `n_px_I` is the number of pixels in the latent image, `I`.
%   `n_px_I = prod(image_sampling_I)`
%
% M_Omega_Phi -- Reprojection matrix
%   The product of three sparse matrices, `M`, `Omega`, and `Phi`.
%   `M_Omega_Phi` converts the latent image `I` to the RAW image `J`.
%
% image_bounds -- Latent image coordinate frame
%   The boundaries of `I` expressed in the coordinate frame of `J`.
%   `image_bounds` is an output argument of the same name of
%   'dispersionfunToMatrix()'. If `add_border` is `false`, `image_bounds`
%   will be equal to `[0, 0, image_sampling_J(2), image_sampling_J(1)]`.
%
%   `image_bounds` is empty if `dispersion` is empty, or is a matrix, as
%   'dispersionfunToMatrix()' is not called.
%
% M -- Mosaicing matrix
%   The output argument of 'mosaicMatrix()'. An (n_px_J)-by-(n_px_J x 3)
%   sparse array.
%
% Omega -- Colour channel conversion matrix
%   The output argument of 'channelConversionMatrix()'. An (n_px_J x
%   3)-by-(n_px_J x c) sparse array which converts the latent image to the
%   RGB colour space of the camera.
%
% Phi -- Dispersion model warp matrix
%   The `W` output argument of 'dispersionfunToMatrix()'. An (n_px_J x
%   c)-by-(n_px_I x c) sparse array. `Phi` is empty if there is no model of
%   dispersion.
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 21, 2018

nargoutchk(1, 5);
narginchk(8, 8);

if length(image_sampling_I) ~= 2
    error('The `image_sampling_I` input argument must contain an image height and width only.');
end

if isStringScalar(int_method) || ischar(int_method)
    do_integration = ~strcmp(int_method, 'none');
else
    error('`int_method` must be a character vector or a string scalar.');
end

n_bands = length(lambda);

has_dispersion = ~isempty(dispersion);
if has_dispersion
    dispersion_is_matrix = false;
    if isfloat(dispersion) && ismatrix(dispersion)
        dispersion_is_matrix = true;
        if size(dispersion, 1) ~= (prod(image_sampling_J) * n_bands)
            error('The `dispersion` matrix must have as many rows as there are pixels in `J` times bands.');
        elseif size(dispersion, 2) ~= (prod(image_sampling_I) * n_bands)
            error('The `dispersion` matrix must have as many columns as there are values in `I`.');
        end
    elseif ~isa(dispersion, 'function_handle')
        error('`dispersion` must be either a floating-point matrix, or a function handle.');
    end
end

M = mosaicMatrix(image_sampling_J, align);
if do_integration
    Omega = channelConversionMatrix(image_sampling_J, sensitivity, lambda, int_method);
else
    Omega = channelConversionMatrix(image_sampling_J, sensitivity);
end
image_bounds = [];
if has_dispersion
    if dispersion_is_matrix
        Phi = dispersion;
    else
        if ~add_border
            image_bounds = [0, 0, image_sampling_J(2), image_sampling_J(1)];
        end
        [ Phi, image_bounds ] = dispersionfunToMatrix(...
           dispersion, lambda, image_sampling_J, image_sampling_I, image_bounds, true...
        );    
    end
    Omega_Phi = Omega * Phi;
else
    Phi = [];
    Omega_Phi = Omega;
end
M_Omega_Phi = M * Omega_Phi;

end