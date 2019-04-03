function varargout = imageFormation(...
    I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
    dispersionfun, align...
)
% IMAGEFORMATION  Patch-wise conversion of a spectral image to RGB and RAW images
%
% ## Syntax
% I_rgb = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options...
% )
% [I_rgb, J_full] = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
%   dispersionfun...
% )
% [I_rgb, J_full, J_est] = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
%   dispersionfun, align...
% )
% [I_rgb, J_full, J_est, I_warped] = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
%   dispersionfun, align...
% )
%
% ## Description
% I_rgb = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options...
% )
%   Returns the colour version of the spectral image.
%
% [I_rgb, J_full] = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
%   dispersionfun...
% )
%   Additionally returns a version of the colour image equivalent of the
%   spectral image, warped according to a dispersion model.
%
% [I_rgb, J_full, J_est] = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
%   dispersionfun, align...
% )
%   Additionally returns RAW image corresponding to the warped colour
%   image.
%
% [I_rgb, J_full, J_est, I_warped] = imageFormation(...
%   I_hyper, bands, color_map, color_bands, sampling_options, patch_options,...
%   dispersionfun, align...
% )
%   Additionally returns a version of the spectral image, warped according
%   to a dispersion model.
%
% ## Input Arguments
%
% I_hyper -- Input image
%   A 2D or 3D array containing the spectral image.
%
% bands -- Input spectral sampling
%   A vector containing the wavelengths corresponding to the third dimension of
%   `I_hyper`. The elements of `bands` are expected to be evenly-spaced.
%
% color_map -- Colour channel spectral sensitivities
%   A 2D array, where `color_map(i, j)` is the sensitivity of the i-th colour
%   channel to the j-th spectral band in `color_bands`. `color_map` is not a
%   colour conversion matrix, as it does not perform the desired numerical
%   integration, over the spectrum, that is part of colour conversion.
%
% color_bands -- Wavelength bands for colour channel sensitivities
%   A vector, of length equal to the size of the second dimension of
%   `color_map`, containing the wavelengths at which the sensitivity
%   functions in `color_map` have been sampled. `color_bands(j)` is the
%   wavelength corresponding to `color_map(:, j)`. The values in
%   `color_bands` are expected to be evenly-spaced.
%
% sampling_options -- Spectral sampling options
%   `sampling_options` is a structure with the following fields used by
%   'colorWeights()' and 'dispersionfunToMatrix()':
%   - 'resolution': A non-negative scalar providing the desired approximate
%     spacing, in pixels, between the images for consecutive wavelengths at
%     which the dispersion is to be sampled. If 'resolution' is zero or is
%     missing, the dispersion function will be evaluated at the same spectral
%     sampling as `I_hyper`.
%   - 'int_method': The numerical integration method to use when
%     integrating over the responses of colour channels to compute colour
%     values. `int_method` is used by `colorWeights()`.
%   - 'support_threshold': A fraction indicating what proportion of the
%     peak magnitude of a colour channel's sensitivity function the user
%     considers to be effectively zero (i.e. no sensitivity).
%   - 'bands_padding': An option controlling how spectral signals are
%     extrapolated.
%   - 'interpolant': The function convolved with spectral signals to
%     interpolate them from the sampling space of `I_hyper` to another sampling
%     space. 'interpolant' is passed to 'resamplingWeights()' as its `f` input
%     argument. Refer to the documentation of 'resamplingWeights.m' for more
%     details.
%   - 'interpolant_ref': Similar to 'interpolant', but used for resampling
%     spectral signals in the sampling space of `color_bands`.
%
%   The 'resolution' and 'support_threshold' fields are mandatory only if
%   `dispersionfun` is passed.
%
% patch_options -- Options for patch-wise image estimation
%   A structure containing the following fields:
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of image patches to be converted individually. `patch_size`
%     does not include padding used to eliminate patch border artifacts. Patches
%     along the bottom and right edges of the image may be made smaller to fit
%     within the image's borders.
%   - 'padding': A scalar containing the pixel width of the border
%     surrounding each image patch. The image patches actually processed are of
%     size `patch_size + padding`, but a border of width 'padding' is stripped
%     when combining patches to form the output image. Note that patches along
%     the edges of the image are not padded to extend outside the image's
%     borders, and so will only have padding towards the interior of the image.
%     'padding' helps mitigate border artifacts, and should be at least as large
%     as the amount of dispersion in the image formation model.
%   - 'target_patch': An optional field. If it exists, 'target_patch' is a
%     two-element vector containing the row and column, respectively, of
%     the top-left corner of the single image patch to be processed. For
%     simplicity, so that it is not necessary to handle edge cases, wherein
%     the patch may not be a valid colour-filter array image, the elements
%     of 'target_patch' must be odd integers. When `target_patch` is
%     passed, all output arguments are calculated for a single image patch,
%     rather than for the entire image. While a border around the patch
%     will have been estimated, with a width given by 'padding', it will
%     not be included in the output. Prior to its removal, the border
%     region will be used to calculate output images aside from `I`, to
%     limit artifacts from image warping.
%
% dispersionfun -- Model of dispersion
%   `dispersionfun` can be empty (`[]`), if there is no model of dispersion.
%   Otherwise, `dispersionfun` must be a function handle, such as produced by
%   'makeDispersionfun()'. `dispersionfun(X)`, where `X` is a three-element row
%   vector (x, y, l), returns the dispersion vector for the position (x, y) in
%   `I_warped` corresponding to light with wavelength `l`. The dispersion vector
%   points from the corresponding position in the reference spectral band to
%   position (x, y). This function will negate dispersion vectors in order to
%   create a warp matrix from `I` to `I_warped`.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the input image `J`. For example, 'gbrg'. `align` has the same form
%   as the `sensorAlignment` input argument of `demosaic()`.
%
% ## Output Arguments
%
% I_rgb -- Colour image
%   The colour equivalent of the input image, generated using the colour channel
%   sensitivity data in `color_map`. An size(I_hyper, 1) x size(I_hyper, 2) x
%   size(color_map, 1) array.
%
%   If `target_patch` is passed, then `I_rgb` has the same first two dimensions
%   as `options.patch_size` (unless it has been clipped to the image borders).
%
% J_full -- Warped colour image
%   A colour image produced by warping `I_hyper` according to the dispersion
%   model, followed by conversion to the colour space of `I_rgb`. A
%   size(I_hyper, 1) x size(I_hyper, 2) x size(color_map, 1) array.
%
%   If `target_patch` is passed, then `J_full` has the same spatial dimensions
%   as `I_rgb`.
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
%   An size(I_hyper, 1) x size(I_hyper, 2) x length(bands) array, storing the
%   latent image warped according to the dispersion model.
%
%   If `target_patch` is passed, then `I_warped` has the same spatial dimensions
%   as `I_rgb`.
%
% See also solvePatchesSpectral, mosaicMatrix, channelConversionMatrix,
% dispersionfunToMatrix

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 7, 2018

narginchk(6, 8);
nargoutchk(1, 4);

n_output_images = nargout;

verbose = true;

if verbose
    tic
end

% Input argument parsing
do_single_patch = isfield(patch_options, 'target_patch');

has_dispersion = ~isempty(dispersionfun);
if has_dispersion && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

image_sampling = [size(I_hyper, 1), size(I_hyper, 2)];
patch_size = patch_options.patch_size;
padding = patch_options.padding;
n_channels_rgb = size(color_map, 1);
n_channels_rgb1 = n_channels_rgb - 1;

if any(mod(image_sampling, 2) ~= 0)
    error('The input image must have dimensions which are even integers, or it cannot represent a colour-filter array image.');
end
if any(mod(patch_size, 2) ~= 0)
    error('`patch_options.patch_size` must be even integers to produce patches which are valid colour-filter array images.');
end
if mod(padding, 2) ~= 0
    error('`patch_options.padding` must be an even integer to produce patches which are valid colour-filter array images.');
end

if length(color_bands) ~= size(color_map, 2)
    error('The number of columns of `color_map` must equal the length of `color_bands`.');
end
n_bands = length(bands);
n_bands1 = n_bands - 1;
if n_bands ~= size(I_hyper, 3)
    error('The length of `bands` must equal the size of `I_hyper` in its third dimension.');
end

if do_single_patch
    target_patch = patch_options.target_patch;
    if target_patch(1) < 1 || target_patch(1) > image_sampling(1) ||...
       target_patch(2) < 1 || target_patch(2) > image_sampling(2)
        error('The target patch is outside of the image bounds.');
    end
    if any(mod(target_patch, 2) ~= 1)
        error(['To prevent having an odd-sized patch because of clippin',...
            'g at the image boundaries, the target patch corner must hav',...
            'e odd integer coordinates.']);
    end
end

color_weights = colorWeights(color_map, color_bands, bands, sampling_options);

if has_dispersion
    dispersion_options_I_warped = struct(...
        'bands_in', bands,...
        'resolution', sampling_options.resolution,...
        'support_threshold', sampling_options.support_threshold,...
        'bands_padding', sampling_options.bands_padding,...
        'interpolant', sampling_options.interpolant_ref,...
        'interpolant_ref', sampling_options.interpolant_ref...
    );
    dispersion_options = dispersion_options_I_warped;
    dispersion_options.bands_out = color_bands;
    dispersion_options.color_map = color_map;
    dispersion_options.int_method = sampling_options.int_method;
else
    dispersion_options_I_warped = [];
    dispersion_options = [];
end

if verbose
    disp('Splitting the input image into columns...');
end

% Channel indices in the input concatenation of images
n_channels_in = 0;
channels_in.I_in = [ 1, size(I_hyper, 3) ];
n_channels_in = n_channels_in + channels_in.I_in(2);

% Channel indices in the output concatenation of images
n_channels_out = 0;
if n_output_images > 0
    channels_out.I_rgb = n_channels_out + [1, n_channels_rgb];
    n_channels_out = n_channels_out + channels_out.I_rgb(2);
    if n_output_images > 1
        channels_out.J_full = n_channels_out + [1, n_channels_rgb];
        n_channels_out = n_channels_out + channels_out.J_full(2);
        if n_output_images > 2
            channels_out.J_est = n_channels_out + [1, 1];
            n_channels_out = n_channels_out + channels_out.J_est(2);
            if n_output_images > 3
                channels_out.I_warped = n_channels_out + [1, n_bands];
                n_channels_out = n_channels_out + channels_out.I_warped(2);
            end
        end
    end
end

% Divide the input images into columns which will be sent to individual
% parallel workers
if do_single_patch
    n_i = 1;
    n_j = 1;
    patch_offset = target_patch - 1;
else
    n_i = ceil(image_sampling(1) / patch_size(1));
    n_j = ceil(image_sampling(2) / patch_size(2));
    patch_offset = [0, 0];
end
n_patches = n_i * n_j;
columns_in = cell(1, n_j);
for j = 1:n_j
    cols_ind_in = [
        max((j - 1) * patch_size(2) + 1 - padding + patch_offset(2), 1);
        min(j * patch_size(2) + padding + patch_offset(2), image_sampling(2))
    ];
    columns_in{j} = zeros(image_sampling(1), diff(cols_ind_in) + 1, n_channels_in);
    columns_in{j}(:, :, channels_in.I_in(1):channels_in.I_in(2)) = I_hyper(:, cols_ind_in(1):cols_ind_in(2), :);
end

if verbose
    fprintf('\tDone.\n');
    disp('Parallel processing of columns...');
end

% Process each column
columns_out = cell(1, n_j);
parfor j = 1:n_j
    column_in_j = columns_in{j};
    image_sampling_p = [0, size(column_in_j, 2)];
    corner = [0, (j - 1) * patch_size(2) + 1 + patch_offset(2)];
    cols_trim_out = [ min(padding + 1, corner(2)), 0 ];
    cols_trim_out(2) = cols_trim_out(1) + min(patch_size(2) - 1, image_sampling(2) - corner(2));
    col_out_width = diff(cols_trim_out) + 1;
    column_out_j = zeros(size(column_in_j, 1), col_out_width, n_channels_out);
    
    % Process each patch within the column
    for i = 1:n_i
        corner(1) = (i - 1) * patch_size(1) + 1 + patch_offset(1);
        rows_trim_out = [ min(padding + 1, corner(1)), 0 ];
        rows_trim_out(2) = rows_trim_out(1) + min(patch_size(1) - 1, image_sampling(1) - corner(1));
        patch_lim_rows = [
            max(corner(1) - padding, 1);
            min(corner(1) + patch_size(1) + padding - 1, image_sampling(1));
        ];
        patch_end_row = min(corner(1) + patch_size(1) - 1, image_sampling(1));
        image_sampling_p(1) = diff(patch_lim_rows) + 1;

        % Generate the output patch
        patches_I_ij = column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.I_in(1):channels_in.I_in(2));
        patches_I_ij_3D = reshape(patches_I_ij, [image_sampling_p n_bands]);

        if n_output_images > 0
            patches_I_rgb_ij_3D = channelConversion(patches_I_ij_3D, color_weights);
            column_out_j(...
                corner(1):patch_end_row, :, channels_out.I_rgb(1):(channels_out.I_rgb(1) + n_channels_rgb1)...
            ) = patches_I_rgb_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);

            if n_output_images > 1
                if has_dispersion
                    patches_J_full_ij_3D = dispersionfunToMatrix(...
                        dispersionfun, dispersion_options, patches_I_ij_3D, true,...
                        flip(corner) - 1 ...
                        );
                    column_out_j(...
                        corner(1):patch_end_row, :, channels_out.J_full(1):(channels_out.J_full(1) + n_channels_rgb1)...
                    ) = patches_J_full_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                else
                    patches_J_full_ij_3D = patches_I_rgb_ij_3D;
                    column_out_j(...
                        corner(1):patch_end_row, :, channels_out.J_full(1):(channels_out.J_full(1) + n_channels_rgb1)...
                    ) = column_out_j(...
                        corner(1):patch_end_row, :, channels_out.I_rgb(1):(channels_out.I_rgb(1) + n_channels_rgb1)...
                    );
                end

                if n_output_images > 2
                    patches_J_est_ij_3D = mosaic(patches_J_full_ij_3D, align);
                    column_out_j(...
                        corner(1):patch_end_row, :, channels_out.J_est(1)...
                    ) = patches_J_est_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);

                    if n_output_images > 3
                        if has_dispersion
                            patches_I_warped_ij_3D = dispersionfunToMatrix(...
                                dispersionfun, dispersion_options_I_warped, patches_I_ij_3D, true,...
                                flip(corner) - 1 ...
                            );
                            column_out_j(...
                                corner(1):patch_end_row, :, channels_out.I_warped(1):(channels_out.I_warped(1) + n_bands1)...
                            ) = patches_I_warped_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                        else
                            column_out_j(...
                                corner(1):patch_end_row, :, channels_out.I_warped(1):(channels_out.I_warped(1) + n_bands1)...
                            ) = patches_I_ij_3D(rows_trim_out(1):rows_trim_out(2), cols_trim_out(1):cols_trim_out(2), :);
                        end
                    end
                end
            end
        end
        
        if verbose
            fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i, n_patches);
        end
    end
    columns_out{j} = column_out_j;
end

if verbose
    fprintf('\tDone.\n');
    disp('Recombining results from patches...');
end

% Recombine patches
if do_single_patch
    images_out = columns_out{1}(target_patch(1):min(target_patch(1) + patch_size(1) - 1, image_sampling(1)), :, :);
else
    images_out = zeros(image_sampling(1), image_sampling(2), n_channels_out);
    offset = 0;
    for j = 1:n_j
        width_j = size(columns_out{j}, 2);
        images_out(:, (offset + 1):(offset + width_j), :) = columns_out{j};
        offset = offset + width_j;
    end
end

if n_output_images > 0
    varargout = cell(1, n_output_images);
    varargout{1} = images_out(:, :, channels_out.I_rgb(1):channels_out.I_rgb(2));
    if n_output_images > 1
        varargout{2} = images_out(:, :, channels_out.J_full(1):channels_out.J_full(2));
        if n_output_images > 2
            varargout{3} = images_out(:, :, channels_out.J_est(1):channels_out.J_est(2));
            if n_output_images > 3
                varargout{4} = images_out(:, :, channels_out.I_warped(1):channels_out.I_warped(2));
            end
        end
    end
end

if verbose
    fprintf('\tDone.\n');
    toc
end
    
end