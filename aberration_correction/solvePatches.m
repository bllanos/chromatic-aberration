function [ I_3D, image_bounds, varargout ] = solvePatches(...
    image_sampling, J_2D, align, dispersionfun, sensitivity,...
    lambda, options, f, f_args, varargin...
    )
% SOLVEPATCHES  Run an image estimation algorithm on image patches
%
% ## Syntax
% I = solvePatches(...
%   image_sampling, J, align, dispersionfun, sensitivity,...
%   lambda, options, f, f_args...
% )
% [ I, image_bounds ] = solvePatches(___)
% [ I, image_bounds, I_rgb ] = solvePatches(___)
% [ I, image_bounds, I_rgb, J_full ] = solvePatches(___)
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatches(___)
% [...
%   I, image_bounds, I_rgb, J_full, J_est, varargout...
% ] = solvePatches(___, target_patch)
%
% ## Description
% I = solvePatches(...
%   image_sampling, J, align, dispersionfun, sensitivity,...
%   lambda, options, f, f_args...
% )
%   Run an image estimation algorithm on image patches and stitch together
%   the results from all patches.
%
% [ I, image_bounds ] = solvePatches(___)
%   Additionally returns the boundaries of the output image in the
%   coordinate system of the input image.
%
% [ I, image_bounds, I_rgb ] = solvePatches(___)
%   Additionally returns the RGB equivalent of the output image.
%
% [ I, image_bounds, I_rgb, J_full ] = solvePatches(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped by the dispersion model.
%
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatches(___)
%   Additionally returns the forward model estimate of the input RAW image
%   `J`.
%
% [...
%   I, image_bounds, I_rgb, J_full, J_est, varargout...
% ] = solvePatches(___, target_patch)
%   Run an image estimation algorithm on a single image patch, and return
%   its additional output arguments.
%
% ## Input Arguments
%
% image_sampling -- Image dimensions
%   A two-element vector containing the height and width, respectively, of
%   the output image `I`, in pixels.   
%
% J -- Input RAW image
%   A 2D array containing the raw colour-filter pattern data of an image.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% dispersionfun -- Model of dispersion
%   A function handle, produced by 'makeDispersionfun()'.
%   `dispersionfun(X)`, where `X` is a three-element row vector (x, y,
%   lambda), returns the dispersion vector for the position (x, y) in `J`
%   corresponding to light with wavelength or colour channel index
%   `lambda`. The dispersion vector corrects for lateral chromatic
%   aberration by pointing from the corresponding position in the reference
%   spectral band or colour channel to position (x, y). This function will
%   negate the dispersion vectors produced by `dispersionfun()` in order to
%   create a warp matrix from `I` to `J`.
%
% sensitivity -- Spectral band conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel of `J` to the j-th input colour channel or spectral
%   band of `I`. `sensitivity` is a matrix mapping colours in `I` to
%   colours in `J`.
%
% lambda -- Wavelength bands
%   A vector of length 'c' containing the wavelengths or colour channel
%   indices at which to evaluate the dispersion model encapsulated by
%   `dispersionfun`. 'c' is the desired number of spectral bands or colour
%   channels in `I`.
%
% f -- Image estimation algorithm
%   The handle to a function implementing an image estimation algorithm.
%   The first output argument of the function must be the output latent
%   image to be estimated. The first few input arguments of the function
%   must be the following:
%   - image_sampling: A two-element vector containing the height and width,
%     respectively, of the output image patch.
%   - align: A four-character character vector, specifying the Bayer tile
%     pattern of the input image patch from `J`.
%   - dispersion_matrix: A warp matrix from coordinates in the output image
%     patch to coordinates in the input image patch.
%   - sensitivity: A copy of this function's argument of the same name.
%   - lambda: A copy of this function's argument of the same name.
%   - J_patch: A patch of the input image `J`.
%
% f_args -- Additional image estimation algorithm parameters
%   A cell vector of input arguments to `f`, beyond those listed above.
%   `f_args` can be an empty cell array, `{}`.
%
% options -- Options and small parameters
%   A structure with the following fields:
%   - 'add_border': A Boolean value indicating whether or not the
%     `image_bounds` input argument of 'dispersionfunToMatrix()' should be
%     empty. If `true`, the output image `I` will be large enough to
%     contain the un-warped coordinates of all pixels in `J`. If `false`,
%     the output image `I` will be clipped to the region occupied by `J`.
%   - 'int_method': The numerical integration method used for spectral to
%     colour space conversion. 'int_method' is passed to
%     'channelConversionMatrix()' as its 'int_method' argument. Refer to
%     the documentation of 'channelConversionMatrix.m' for details. If
%     'int_method' is 'none' or empty (`[]`), as should be the case when
%     colour conversion is from a set of colour channels, not from a
%     spectral space, numerical integration will not be performed.
%     'int_method' is unused if none of `I_rgb`, `J_full`, or `J_est` are
%     requested as output arguments.
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of the image patches to be estimated. `patch_size` does
%     not include padding used to eliminate artifacts from the patch-wise
%     estimation.
%   - 'padding': A scalar containing the pixel width of the border
%     surrounding each image patch. When running the image estimation
%     algorithm, this function passes the algorithm a patch of the input
%     image. The input patch is large enough for estimating a patch of
%     dimensions `patch_size` in the output image, as well as a border of
%     width `padding` around the output patch. (Note that the patch of the
%     input image accommodates for the warping given by `dispersionfun`.)
%     The final output image, `I`, is assembled from only the central
%     (non-padding) regions of patches, so this function discards the
%     padding regions estimated by the algorithm.
%
% target_patch -- Single patch coordinates
%   A two-element vector containing the row and column, respectively, of
%   the top-left corner of the image patch to be estimated. When
%   `target_patch` is passed, all output arguments are calculated for a
%   single image patch, rather than for the entire image. While a border
%   around `I` will have been estimated, with a width given by
%   `options.padding`, it will not be included in the output. The border
%   region will not be used when calculating `I_rgb`, `J_full`, and
%   `J_est`.
%
% ## Output Arguments
%
% I -- Latent image
%   An image_sampling(1) x image_sampling(2) x length(lambda) array,
%   storing the latent image estimated on a patch-wise basis using the
%   function `f`.
%
%   If `target_patch` is passed, then `I` is at most an
%   options.patch_size(1) x options.patch_size(2) x length(lambda) array
%   containing the patch of the latent image having its top left corner at
%   position `target_patch`.
%
% image_bounds -- Latent image coordinate frame
%   The boundaries of the latent image expressed in the coordinate frame of
%   `J`. `image_bounds` has the form of the output argument of the same
%   name of 'dispersionfunToMatrix()'. If `options.add_border` is `false`,
%   `image_bounds` will be equal to `[0, 0, size(J, 2), size(J, 1)]`. Note
%   that `image_bounds` pertains to the entire latent image, even if only a
%   patch of the latent image is requested (using `target_patch`).
%
% I_rgb -- Latent RGB image
%   The RGB equivalent of the latent image, generated using the colour
%   space conversion data in `sensitivity`. An image_sampling(1) x
%   image_sampling(2) x 3 array.
%
%   If `target_patch` is passed, then `I_rgb` has the same first two
%   dimensions as `I`, and a length of 3 in its third dimension.
%
% J_full -- Warped latent RGB image
%   An RGB image produced by warping `I` according to the dispersion model,
%   followed by conversion to the RGB colour space of `J`. An size(J, 1) x
%   size(J, 2) x 3 array.
%
%   If `target_patch` is passed, then `J_full` has spatial dimensions which
%   depend on the nature of the warp between `I` and `J`, and on the
%   padding region surrounding `I`.
%
% J_est -- Re-estimated input RAW image
%   The mosaiced version of `J_full`, representing the result of passing
%   `I` through the forward model of dispersion and image capture. An array
%   with the same dimensions as `J`.
%
%   If `target_patch` is passed, then `J_est` is a 2D array with the same
%   sizes in its first two dimensions as `J_full`.
%
% varargout -- Additional output arguments of the image estimation algorithm
%   When `target_patch` is passed, this function can return additional
%   output arguments from `f`, as there is no concern over how to combine
%   them from multiple image patches.
%
% ## References
%
% This function combines image patches by placing them side-by-side and
% discarding the overlapping regions. There are other methods for combining
% overlapping image patches, such as by averaging them. The following
% article contains a brief discussion of patch-wise solutions, in the
% context of chromatic aberration correction:
%
% Sun, T., Peng, Y., and Heidrich, W. (2017). "Revisiting cross-channel
%   information transfer for chromatic aberration correction." In 2017 IEEE
%   International Conference on Computer Vision (ICCV) (pp. 3268–3276).
%   doi:10.1109/ICCV.2017.352
%
% See also mosaicMatrix, channelConversionMatrix, dispersionfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 12, 2018

narginchk(9, 10);

single_patch = false;
if ~isempty(varargin)
    target_patch = varargin{1};
    single_patch = true;
end

if length(image_sampling) ~= 2
    error('The `image_sampling` input argument must contain an image height and width only.');
end

if ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

varargout = cell(nargout - 2, 1);

% Create the dispersion matrix
image_sampling_J = size(J_2D);
if options.add_border
    image_bounds = [];
else
    image_bounds = [0, 0, image_sampling_J(2), image_sampling_J(1)];
end
fprintf('Calculating the dispersion matrix...\n');
[ dispersion_matrix, image_bounds ] = dispersionfunToMatrix(...
   dispersionfun, lambda, image_sampling_J, image_sampling, image_bounds, true...
);
fprintf('\t...done\n');    

if single_patch
    [...
        I_3D, image_sampling_J_patch, dispersion_matrix_patch, varargout{4:end}...
    ] = solveOnePatch(...
        image_sampling, J_2D, align, dispersion_matrix, sensitivity,...
        lambda, options.patch_size, options.padding, f, f_args, target_patch...
    );
else
    I_3D = zeros([image_sampling, length(lambda)]);
    i_vector = 1:options.patch_size(1):image_sampling(1);
    n_i_vector = length(i_vector);
    j_vector = 1:options.patch_size(2):image_sampling(2);
    n_j_vector = length(j_vector);
    patches = cell(length(i_vector), length(j_vector));
    parfor i = 1:n_i_vector
        % This line avoids a variable classification error with `patches`;
        % "If you use a nested for-loop to index into a sliced array, you
        % cannot use that array elsewhere in the parfor-loop."
        patches_i = cell(1, n_j_vector);
        for j = 1:n_j_vector
            patches_i{j} = solveOnePatch(...
                image_sampling, J_2D, align, dispersion_matrix, sensitivity,...
                lambda, options.patch_size, options.padding, f, f_args,...
                [i_vector(i), j_vector(j)]...
            );
        end
        patches{i, :} = patches_i;
    end
    for i = 1:n_i_vector
        end_i = min(i_vector(i) + options.patch_size(1) - 1, image_sampling(1));
        for j = 1:n_j_vector
            end_j = min(j_vector(j) + options.patch_size(2) - 1, image_sampling(2));
            I_3D(...
                i_vector(i):end_i,...
                j_vector(j):end_j, :...
            ) = patches{i, j};
        end
    end
end

if nargout > 2
    I = I_3D(:);
    if single_patch
        image_sampling_I_patch = [size(I_3D, 1), size(I_3D, 2)];
    else
        image_sampling_I_patch = image_sampling;
        image_sampling_J_patch = image_sampling_J;
        dispersion_matrix_patch = dispersion_matrix;
        target_patch = [1, 1];
    end
    do_integration = ~(...
        isempty(options.int_method) || strcmp(options.int_method, 'none')...
    );
    if do_integration
        Omega_I = channelConversionMatrix(...
            image_sampling_I_patch, sensitivity, lambda, options.int_method...
        );
    else
        Omega_I = channelConversionMatrix(image_sampling_I_patch, sensitivity);
    end
    I_rgb = Omega_I * I;
    n_channels_rgb = 3;
    varargout{1} = reshape(I_rgb, image_sampling_I_patch(1), image_sampling_I_patch(2), n_channels_rgb);
    
    if nargout > 3
        if do_integration
            Omega_J = channelConversionMatrix(...
                image_sampling_J_patch, sensitivity, lambda, options.int_method...
            );
        else
            Omega_J = channelConversionMatrix(image_sampling_J_patch, sensitivity);
        end
        J_full = Omega_J * dispersion_matrix_patch * I;
        varargout{2} = reshape(...
            J_full,...
            image_sampling_J_patch(1), image_sampling_J_patch(2), n_channels_rgb...
        );
        
        if nargout > 4
            M = mosaicMatrix(...
                image_sampling_J_patch,...
                offsetBayerPattern(target_patch, align)...
            );
            J_est = M * J_full;
            varargout{3} = reshape(J_est, image_sampling_J_patch);
        end
    end
end
    
end