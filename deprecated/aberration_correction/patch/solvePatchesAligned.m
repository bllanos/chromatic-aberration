function [ I_3D, image_bounds, varargout ] = solvePatchesAligned(...
    J, align, dispersionfun, sensitivity,...
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
% than globally for the whole image. Moreover, the input image is not
% copied to each parallel worker.
%
% ## Syntax
% I = solvePatchesAligned(...
%   J, align, dispersionfun, sensitivity, lambda, options, f, f_args...
% )
% [ I, image_bounds ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb, J_full ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatchesAligned(___)
% [ I, image_bounds, I_rgb, J_full, J_est, I_warped ] = solvePatchesAligned(___)
% [...
%   I, image_bounds, I_rgb, J_full, J_est, I_warped, varargout...
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
%   image, warped according to the dispersion model.
%
% [ I, image_bounds, I_rgb, J_full, J_est ] = solvePatchesAligned(___)
%   Additionally returns the forward model estimate of the RAW image.
%
% [ I, image_bounds, I_rgb, J_full, J_est, I_warped ] = solvePatchesAligned(___)
%   Additionally returns a version of the latent image warped according to
%   the dispersion model.
%
% [...
%   I, image_bounds, I_rgb, J_full, J_est, I_warped, varargout...
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
%   as the `sensorAlignment` input argument of `demosaic()`. `align` can
%   also be empty, indicating that the input image is not mosaiced.
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
%   channels in `I`. If `dispersionfun` is empty, `lambda` can also be
%   empty.
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
%   - J_patch: A patch of the input image `J`.
%
% f_args -- Additional image estimation algorithm parameters
%   A cell vector of input arguments to `f`, beyond those listed above.
%   `f_args` can be an empty cell array, `{}`.
%
% options -- Options and small parameters
%   A structure with the following fields:
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
% I_warped -- Warped latent image
%   An size(J, 1) x size(J, 2) x length(lambda) array, storing the latent
%   image warped according to the dispersion model.
%
%   If `target_patch` is passed, then `I_warped` has the same spatial
%   dimensions as `I`.
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

verbose = true;
if verbose
    tic
end

single_patch = false;
if ~isempty(varargin)
    target_patch = varargin{1};
    single_patch = true;
end

has_dispersion = ~isempty(dispersionfun);
if has_dispersion && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

varargout = cell(1, nargout - 2);
n_auxiliary_images = min(nargout - 2, 4);
image_sampling = [size(J, 1), size(J, 2)];
image_bounds = [0, 0, image_sampling(2), image_sampling(1)];
patch_size = options.patch_size;
n_bands = size(sensitivity, 2);

if single_patch
    [ patch_lim, trim ] = patchBoundaries(...
        image_sampling, patch_size, options.padding, target_patch...
    );

    % Construct arguments for the image estimation algorithm
    if isempty(align)
        align_f = [];
    else
        align_f = offsetBayerPattern(patch_lim(1, :), align);
    end
    image_sampling_f = diff(patch_lim, 1, 1) + 1;
    if has_dispersion
        dispersion_matrix_patch = dispersionfunToMatrix(...
            dispersionfun, lambda, image_sampling_f, image_sampling_f,...
            [0, 0, image_sampling_f(2), image_sampling_f(1)], true, flip(target_patch) - 1 ...
        );
    else
        dispersion_matrix_patch = [];
    end
    J_f = J(patch_lim(1, 1):patch_lim(2, 1), patch_lim(1, 2):patch_lim(2, 2), :);

    % Solve for the output patch
    [I_3D, varargout{(n_auxiliary_images + 1):end}] = f(...
        image_sampling_f, align_f, dispersion_matrix_patch, sensitivity,...
        J_f, f_args{:}...
    );

    padding_filter = false(image_sampling_f);
    padding_filter((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2))) = true;
    [I_3D, varargout(1:n_auxiliary_images)] = estimateAuxiliaryImages(...
            I_3D, dispersion_matrix_patch, padding_filter, diff(trim, 1, 1) + 1,...
            sensitivity, align_f, n_auxiliary_images...
        );

else
    if verbose
        disp('Splitting the input image into patches...');
    end
    i_vector = 1:patch_size(1):image_sampling(1);
    n_i_vector = length(i_vector);
    j_vector = 1:patch_size(2):image_sampling(2);
    n_j_vector = length(j_vector);
    n_patches = n_i_vector * n_j_vector;
    patches_J = cell(1, n_j_vector);
    patches_I = cell(1, n_j_vector);
    patches_auxiliary = cell(n_i_vector, n_j_vector, n_auxiliary_images);
    patch_limits = zeros(n_i_vector, n_j_vector, 4);
    patch_trim = zeros(n_i_vector, n_j_vector, 4);
    corners = zeros(n_i_vector, n_j_vector, 2);
    
    % Divide the input image into patches to be sent to individual parallel
    % workers
    for j = 1:n_j_vector
        for i = 1:n_i_vector
            corners(i, j, :) = [i_vector(i), j_vector(j)];
            [ patch_lim, trim ] = patchBoundaries(...
                image_sampling, patch_size, options.padding, corners(i, j, :)...
            );
            patch_trim(i, j, :) = reshape(trim, 1, 1, 4);
            patch_limits(i, j, :) = reshape(patch_lim, 1, 1, 4);
            
            if i == 1
                patches_J{j} = J(:, patch_lim(1, 2):patch_lim(2, 2), :);
            end
        end
    end
    
    if verbose
        fprintf('\tDone.\n');
        disp('Parallel processing of patches...');
    end
    
    % Process each patch
    parfor j = 1:n_j_vector
        patch_limits_j = patch_limits(:, j, :);
        patch_trim_j = patch_trim(:, j, :);
        corners_j = corners(:, j, :);
        patches_J_j = patches_J{j};
        patches_I_j = zeros(...
            size(patches_J_j, 1),...
            patch_trim_j(1, 1, 4) - patch_trim_j(1, 1, 3) + 1,...
            n_bands...
        );
        patches_auxiliary_j = cell(n_i_vector, 1, n_auxiliary_images);
        for i = 1:n_i_vector
            patch_lim = reshape(patch_limits_j(i, 1, :), 2, 2);
            if isempty(align)
                align_f = [];
            else
                align_f = offsetBayerPattern(patch_lim(1, :), align);
            end
            image_sampling_f = diff(patch_lim, 1, 1) + 1;
            trim = reshape(patch_trim_j(i, 1, :), 2, 2);
            if has_dispersion
                dispersion_matrix_patch = dispersionfunToMatrix(...
                    dispersionfun, lambda, image_sampling_f, image_sampling_f,...
                    [0, 0, image_sampling_f(2), image_sampling_f(1)], true,...
                    [corners_j(i, 1, 2), corners_j(i, 1, 1)] - 1 ...
                );
            else
                dispersion_matrix_patch = [];
            end

            % Solve for the output patch
            patches_I_ij = f(...
                image_sampling_f, align_f, dispersion_matrix_patch, sensitivity,...
                patches_J_j(patch_lim(1, 1):patch_lim(2, 1), :, :), f_args{:}...
            );

            padding_filter = false(image_sampling_f);
            padding_filter((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2))) = true;
            patch_trimmed_size = diff(trim, 1, 1) + 1;
            [patches_I_ij, patches_auxiliary_j(i, 1, :)] = estimateAuxiliaryImages(...
                    patches_I_ij, dispersion_matrix_patch, padding_filter,...
                    patch_trimmed_size,...
                    sensitivity, align_f, n_auxiliary_images...
            );
            patches_I_j(...
                    ((i - 1) * patch_size(1) + 1):((i - 1) * patch_size(1) + patch_trimmed_size(1)),...
                    :, :...
                ) = patches_I_ij;
            if verbose
                fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i_vector, n_patches);
            end
        end
        patches_I{j} = patches_I_j;
        patches_auxiliary(:, j, :) = patches_auxiliary_j;
    end
    
    if verbose
        fprintf('\tDone.\n');
        disp('Recombining results from patches...');
    end
    
    % Recombine patches
    I_3D = cell2mat(patches_I);
    for im = 1:n_auxiliary_images
        varargout{im} = cell2mat(patches_auxiliary(:, :, im));
    end
    
    if verbose
        fprintf('\tDone.\n');
        toc
    end
end
    
end