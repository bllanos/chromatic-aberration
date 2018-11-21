function varargout = imageFormation(...
    I_hyper, sensitivity, options, varargin...
    )
% IMAGEFORMATION  Patch-wise conversion of a spectral image to RGB and RAW images
%
% ## Syntax
% I_rgb = imageFormation(...
%   I_hyper, sensitivity, options [, target_patch]...
% )
% [I_rgb, J_full] = imageFormation(...
%   I_hyper, sensitivity, options, dispersionfun, lambda [, target_patch]...
% )
% [I_rgb, J_full, J_est] = imageFormation(...
%   I_hyper, sensitivity, options, dispersionfun, lambda, align [, target_patch]...
% )
% [I_rgb, J_full, J_est, I_warped] = imageFormation(...
%   I_hyper, sensitivity, options, dispersionfun, lambda, align [, target_patch]...
% )
%
% ## Description
% I_rgb = imageFormation(...
%   I_hyper, sensitivity, options [, target_patch]...
% )
%   Returns the colour version of the spectral image.
%
% [I_rgb, J_full] = imageFormation(...
%   I_hyper, sensitivity, options, dispersionfun, lambda [, target_patch]...
% )
%   Additionally returns a version of the colour image equivalent of the
%   spectral image, warped according to a dispersion model.
%
% [I_rgb, J_full, J_est] = imageFormation(...
%   I_hyper, sensitivity, options, dispersionfun, lambda, align [, target_patch]...
% )
%   Additionally returns RAW image corresponding to the warped colour
%   image.
%
% [I_rgb, J_full, J_est, I_warped] = imageFormation(...
%   I_hyper, sensitivity, options, dispersionfun, lambda, align [, target_patch]...
% )
%   Additionally returns a version of the spectral image, warped according
%   to a dispersion model.
%
% ## Input Arguments
%
% I_hyper -- Input image
%   A 2D or 3D array containing the spectral image.
%
% sensitivity -- Colour space conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `I_rgb` to the j-th input colour channel or spectral
%   band of `I_hyper`. `sensitivity` is a matrix mapping colours in
%   `I_hyper` to colours in `I_rgb`. `sensitivity` must account for any
%   numerical intergration that is part of colour conversion.
%
% options -- Options and small parameters
%   A structure with the following fields:
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of the image patches to be estimated. 'patch_size' does
%     not include padding used to eliminate artifacts from the patch-wise
%     estimation.
%   - 'padding': A scalar containing the pixel width of the border
%     surrounding each image patch. 'padding' is needed to prevent edges
%     from appearing between patches when the image patches are warped
%     using `dispersionfun`.
%
%   'patch_size' and 'padding' are used internally by
%   'solvePatchesAlignedSerial()', the documentation of which describes
%   these parameters in more detail.
%
% dispersionfun -- Model of dispersion
%   A function handle, produced by 'makeDispersionfun()'.
%   `dispersionfun(X)`, where `X` is a three-element row vector (x, y,
%   lambda), returns the dispersion vector for the position (x, y) in
%   `J_full` corresponding to light with wavelength or colour channel index
%   `lambda`. The dispersion vector corrects for lateral chromatic
%   aberration by pointing from the corresponding position in the reference
%   spectral band or colour channel to position (x, y). This function will
%   negate the dispersion vectors produced by `dispersionfun()` in order to
%   create a warp matrix from `I_hyper` to `J_full`.
%
%   `dispersionfun` can be empty (`[]`), if there is no model of
%   dispersion.
%
% lambda -- Wavelength bands
%   A vector of length `size(sensitivity, 2)` containing the wavelengths or
%   colour channel indices of the spectral bands or colour channels in
%   `I_hyper`. `lambda` is needed to evaluate the dispersion model in
%   `dispersionfun`. `lambda` can be empty if `dispersionfun` is empty.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% target_patch -- Single patch coordinates
%   A two-element vector containing the row and column, respectively, of
%   the top-left corner of the image patch to be estimated. When
%   `target_patch` is passed, all output arguments are calculated for a
%   single image patch, rather than for the entire image. The border around
%   the patch, with a width given by `options.padding`, will not be
%   included in the output, but will have been used to limit artifacts from
%   image warping when applying `dispersionfun`.
%
% ## Output Arguments
%
% I_rgb -- Colour image
%   The colour equivalent of the input image, generated using the colour
%   space conversion data in `sensitivity`. An size(I_hyper, 1) x
%   size(I_hyper, 2) x size(sensitivity, 1) array.
%
%   If `target_patch` is passed, then `I_rgb` has the same first two
%   dimensions as `options.patch_size` (unless it has been clipped to the
%   image borders).
%
% J_full -- Warped colour image
%   A colour image produced by warping `I_hyper` according to the
%   dispersion model, followed by conversion to the colour space of
%   `I_rgb`. A size(I_hyper, 1) x size(I_hyper, 2) x size(sensitivity, 1)
%   array.
%
%   If `target_patch` is passed, then `J_full` has the same spatial
%   dimensions as `I_rgb`.
%
% J_est -- Estimated RAW image
%   The mosaiced version of `J_full`, representing the result of passing
%   `I_hyper` through the forward model of dispersion and image capture. A
%   size(I_hyper, 1) x size(I_hyper, 2) array.
%
%   If `target_patch` is passed, then `J_est` is a 2D array with the same
%   sizes in its first two dimensions as `I_rgb`.
%
% I_warped -- Warped latent image
%   An size(I_hyper, 1) x size(I_hyper, 2) x length(lambda) array, storing
%   the latent image warped according to the dispersion model.
%
%   If `target_patch` is passed, then `I_warped` has the same spatial
%   dimensions as `I_rgb`.
%
% See also solvePatchesAlignedSerial, mosaicMatrix, channelConversionMatrix,
% dispersionfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 7, 2018

nargoutchk(1, 4);
single_patch = false;
output_rgb_warped = (nargout > 1);
output_raw = (nargout > 2);
if nargout == 1
    narginchk(3, 4);
    if ~isempty(varargin)
        target_patch = varargin{1};
        single_patch = true;
    end 
elseif nargout == 2
    narginchk(5, 6);
    dispersionfun = varargin{1};
    lambda = varargin{2};
    if length(varargin) > 2
        target_patch = varargin{3};
        single_patch = true;
    end
elseif nargout > 2
    narginchk(6, 7);
    dispersionfun = varargin{1};
    lambda = varargin{2};
    align = varargin{3};
    if length(varargin) > 3
        target_patch = varargin{4};
        single_patch = true;
    end
end

if output_rgb_warped && ~isempty(dispersionfun) && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

varargout = cell(1, nargout);

    function I_in = identity(~, ~, ~, ~, I_in)
    end

if single_patch
    if output_raw
        [~, ~, varargout{:}] = solvePatchesAlignedSerial(...
            I_hyper, align, dispersionfun, sensitivity, lambda, options,...
            @identity, {}, target_patch...
        );
    elseif output_rgb_warped
        [~, ~, varargout{:}] = solvePatchesAlignedSerial(...
            I_hyper, [], dispersionfun, sensitivity, lambda, options,...
            @identity, {}, target_patch...
        );
    else
        [~, ~, varargout{:}] = solvePatchesAlignedSerial(...
            I_hyper, [], [], sensitivity, [], options,...
            @identity, {}, target_patch...
        );
    end
else
    if output_raw
        [~, ~, varargout{:}] = solvePatchesAlignedSerial(...
            I_hyper, align, dispersionfun, sensitivity, lambda, options,...
            @identity, {}...
        );
    elseif output_rgb_warped
        [~, ~, varargout{:}] = solvePatchesAlignedSerial(...
            I_hyper, [], dispersionfun, sensitivity, lambda, options,...
            @identity, {}...
        );
    else
        [~, ~, varargout{:}] = solvePatchesAlignedSerial(...
            I_hyper, [], [], sensitivity, [], options,...
            @identity, {}...
        );
    end
end
    
end