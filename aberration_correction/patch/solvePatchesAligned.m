function [ I_3D, image_bounds, varargout ] = solvePatchesAligned(...
    J_2D, align, dispersionfun, sensitivity,...
    lambda, options, f, f_args, varargin...
    )
% SOLVEPATCHESALIGNED  Run an image estimation algorithm on image patches
%
% ## Usage
%
% This function is a lower-memory version of 'solvePatches()', for the case
% where the output image has the same spatial resolution as the input
% image, and where the boundaries of the two images are made to coincide in
% space. Individual image patches are estimated under these same
% constraints, such that a larger padding parameter (`options.padding`) may
% be necessary to avoid artifacts from patch-wise image estimation in the
% presence of large dispersion. Specifically, if the dispersion vectors
% modelled by `dispersionfun` are large, then an input patch may not
% contain all the pixels which, because of dispersion, influence the
% corresponding output patch. By increasing the amount of padding around
% each output patch, these missing pixels will only have a significant
% effect within the padding area, which is discarded when assembling the
% full output image.
%
% In contrast with 'solvePatches()', because of the restrictions on the
% relationship between the input and output images, this function does not
% use an `add_border` field in its `options` structure. Note that both
% 'solvePatches()' and 'solvePatchesAligned()' calculate per-patch
% dispersion matrices as input for the image estimation algorithms that
% they are managing. Consequently, the algorithms cannot make use of an
% `add_border` option when called on individual patches, as they are not in
% control of generating dispersion matrices.
%
% 'solvePatchesAligned()' is lower-memory than 'solvePatches()' because all
% matrices involved in image estimation are estimated per-patch, rather
% than globally for the whole image.
%
% ## Syntax
% I = solvePatchesAligned(...
%   J, align, dispersionfun, sensitivity, lambda, options, f, f_args...
% )
% [ I, image_bounds ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb, J_full ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatchesAligned(___)
% [...
%   I, image_bounds, I_rgb, J_full, J_est, varargout...
% ] = solvePatchesAligned(___, target_patch)
%
% ## Description
% I = solvePatchesAligned(...
%   J, align, dispersionfun, sensitivity, lambda, options, f, f_args...
% )
%   Run an image estimation algorithm on image patches and stitch together
%   the results from all patches.
%
% [ I, image_bounds ] = solvePatchesAligned(___)
%   Additionally returns the boundaries of the output image in the
%   coordinate system of the input image.
%
% [ I, image_bounds, I_rgb ] = solvePatchesAligned(___)
%   Additionally returns the RGB equivalent of the output image.
%
% [ I, image_bounds, I_rgb, J_full ] = solvePatchesAligned(___)
%   Additionally returns a version of the RGB equivalent of the latent
%   image, warped by the dispersion model.
%
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatchesAligned(___)
%   Additionally returns the forward model estimate of the RAW image.
%
% [...
%   I, image_bounds, I_rgb, J_full, J_est, varargout...
% ] = solvePatchesAligned(___, target_patch)
%   Run an image estimation algorithm on a single image patch, and return
%   its additional output arguments.
%
% ## Input Arguments
%
% J -- Input image
%   A 2D or 3D array containing the input image.
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
%   `dispersionfun` can be empty (`[]`), if there is no model of
%   dispersion.
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
%   - dispersion_matrix: Either a warp matrix from coordinates in the
%     output image patch to coordinates in the input image patch, or an
%     empty array (`[]`), if there is no model of dispersion.
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
%     input image does not accommodate for the warping given by
%     `dispersionfun`.) The final output image, `I`, is assembled from only
%     the central (non-padding) regions of patches, so this function
%     discards the padding regions estimated by the algorithm.
%
% target_patch -- Single patch coordinates
%   A two-element vector containing the row and column, respectively, of
%   the top-left corner of the image patch to be estimated. When
%   `target_patch` is passed, all output arguments are calculated for a
%   single image patch, rather than for the entire image. While a border
%   around `I` will have been estimated, with a width given by
%   `options.padding`, it will not be included in the output. The border
%   region will be used when calculating `I_rgb`, `J_full`, and `J_est`, to
%   limit artifacts from image warping, before being stripped off.
%
% ## Output Arguments
%
% I -- Latent image
%   A size(J, 1) x size(J, 2) x length(lambda) array, storing the latent
%   image estimated on a patch-wise basis using the function `f`.
%
%   If `target_patch` is passed, then `I` is at most an
%   options.patch_size(1) x options.patch_size(2) x length(lambda) array
%   containing the patch of the latent image having its top left corner at
%   position `target_patch`.
%
% image_bounds -- Latent image coordinate frame
%   The boundaries of the latent image expressed in the coordinate frame of
%   `J`. `image_bounds` has the form of the output argument of the same
%   name of 'dispersionfunToMatrix()', and will be equal to `[0, 0, size(J,
%   2), size(J, 1)]`. Note that `image_bounds` pertains to the entire
%   latent image, even if only a patch of the latent image is requested
%   (using `target_patch`).
%
% I_rgb -- Latent RGB image
%   The RGB equivalent of the latent image, generated using the colour
%   space conversion data in `sensitivity`. An size(J, 1) x size(J, 2) x 3
%   array.
%
%   If `target_patch` is passed, then `I_rgb` has the same first two
%   dimensions as `I`, and a length of 3 in its third dimension.
%
% J_full -- Warped latent RGB image
%   An RGB image produced by warping `I` according to the dispersion model,
%   followed by conversion to the RGB colour space of `J`. A size(J, 1) x
%   size(J, 2) x 3 array.
%
%   If `target_patch` is passed, then `J_full` has the same spatial
%   dimensions as `I`.
%
% J_est -- Estimated RAW image
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
% See also solvePatches, mosaicMatrix, channelConversionMatrix,
% dispersionfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 6, 2018

narginchk(8, 9);

single_patch = false;
if ~isempty(varargin)
    target_patch = varargin{1};
    single_patch = true;
end

if ~isempty(dispersionfun) && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

varargout = cell(1, nargout - 2);
n_auxiliary_images = min(nargout - 2, 3);
image_sampling = [size(J_2D, 1), size(J_2D, 2)];
image_bounds = [0, 0, image_sampling(2), image_sampling(1)];

if single_patch
    [...
        I_3D, dispersion_matrix_patch, padding_filter, I_size,...
        varargout{4:end}...
    ] = solveOnePatchAligned(...
        J_2D, align, dispersionfun, sensitivity,...
        lambda, options.patch_size, options.padding, f, f_args, target_patch...
    );

    [I_3D, varargout(1:n_auxiliary_images)] = estimateAuxiliaryImages(...
            I_3D, dispersion_matrix_patch, padding_filter, I_size,...
            sensitivity, lambda, options.int_method,...
            align, target_patch, n_auxiliary_images...
        );

else
    n_channels_rgb = 3;
    i_vector = 1:options.patch_size(1):image_sampling(1);
    n_i_vector = length(i_vector);
    j_vector = 1:options.patch_size(2):image_sampling(2);
    n_j_vector = length(j_vector);
    patches = cell(length(i_vector), length(j_vector));
    patches_auxiliary = cell(length(i_vector), length(j_vector), n_auxiliary_images);
    
    parfor i = 1:n_i_vector
        % These lines avoid a variable classification error with `patches*`;
        % "If you use a nested for-loop to index into a sliced array, you
        % cannot use that array elsewhere in the parfor-loop."
        patches_i = cell(1, n_j_vector);
        patches_auxiliary_i = cell(1, n_j_vector, n_auxiliary_images);
        for j = 1:n_j_vector
            corner_ij = [i_vector(i), j_vector(j)];
            [...
                patches_i{j}, dispersion_matrix_patch, padding_filter, I_size...
            ] = solveOnePatchAligned(...
                J_2D, align, dispersionfun, sensitivity,...
                lambda, options.patch_size, options.padding, f, f_args,...
                corner_ij...
            );
        
            [...
                patches_i{j}, patches_auxiliary_i(1, j, :)...
            ] = estimateAuxiliaryImages(...
                patches_i{j}, dispersion_matrix_patch,...
                padding_filter, I_size,...
                sensitivity, lambda, options.int_method,...
                align, corner_ij, n_auxiliary_images...
            );
        end
        patches(i, :) = patches_i;
        patches_auxiliary(i, :, :) = patches_auxiliary_i;
    end
    
    I_3D = zeros([image_sampling, length(lambda)]);
    output_rgb = n_auxiliary_images > 0;
    if output_rgb
        I_rgb = zeros([image_sampling, n_channels_rgb]);
    end
    output_rgb_warped = n_auxiliary_images > 1;
    if output_rgb_warped
        J_full = zeros([image_sampling, n_channels_rgb]);
    end
    output_raw = n_auxiliary_images > 2;
    if output_raw
        J_est = zeros(image_sampling);
    end
    for i = 1:n_i_vector
        end_i = min(i_vector(i) + options.patch_size(1) - 1, image_sampling(1));
        for j = 1:n_j_vector
            end_j = min(j_vector(j) + options.patch_size(2) - 1, image_sampling(2));
            I_3D(...
                i_vector(i):end_i,...
                j_vector(j):end_j, :...
            ) = patches{i, j};
            if output_rgb
                I_rgb(...
                    i_vector(i):end_i,...
                    j_vector(j):end_j, :...
                ) = patches_auxiliary{i, j, 1};
            end
            if output_rgb_warped
                J_full(...
                    i_vector(i):end_i,...
                    j_vector(j):end_j, :...
                ) = patches_auxiliary{i, j, 2};
            end
            if output_raw
                J_est(...
                    i_vector(i):end_i,...
                    j_vector(j):end_j, :...
                ) = patches_auxiliary{i, j, 3};
            end
        end
    end
    
    if output_rgb
        varargout{1} = I_rgb;
    end
    if output_rgb_warped
        varargout{2} = J_full;
    end
    if output_raw
        varargout{3} = J_est;
    end
    
end
    
end