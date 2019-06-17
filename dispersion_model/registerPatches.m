function [ centers ] = registerPatches(...
    I, reference_index, patch_options, optimizer, metric,...
    pyramid_levels, varargin...
)
% REGISTERPATCHES  Find patch-wise displacements between the channels of an image
%
% ## Syntax
% centers = registerPatches(...
%   I, reference_index, patch_options, optimizer, metric,...
%   pyramid_levels [, verbose]...
% )
%
% centers = registerPatches(...
%   I, reference_index, patch_options, optimizer, metric,...
%   pyramid_levels [, bands, dispersionfun, verbose]...
% )
%
% ## Description
% centers = registerPatches(...
%   I, reference_index, patch_options, optimizer, metric,...
%   pyramid_levels [, verbose]...
% )
%   Returns the locations of patches in the colour channels or spectral
%   bands of the image following patch-wise registration of the
%   channels/bands.
%
% centers = registerPatches(...
%   I, reference_index, patch_options, optimizer, metric,...
%   pyramid_levels [, bands, dispersionfun, verbose]...
% )
%   Same as the previous syntax, but provides an initial guess for the
%   image registration transformations.
%
% ## Input Arguments
%
% I -- Input image
%   A 3D array containing the image.
%
% reference_index -- Reference channel/band index
%   The index of the colour channel or spectral band of `I` which is to be
%   used as the reference image for image registration.
%
% patch_options -- Options for patch-wise image decomposition
%   A structure containing the following fields:
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of image regions for which individual channel-wise
%     registration transformations are to be defined. Regions along the
%     bottom and right edges of the image may be made smaller to fit within
%     the image's borders. 'patch_size' defines the resolution of the
%     patch-wise image registration, whereas the combination of
%     'patch_size' and 'padding' defines the size of the patches used to
%     compute registration transformations.
%   - 'padding': A scalar containing the pixel width of the border
%     surrounding each image region. The image patches actually registered
%     are of size `patch_size + 2 * padding`. Note that patches along the edges
%     of the image are not padded to extend outside the image's borders,
%     and so will only have padding towards the interior of the image.
%   - 'target_patch': An optional field. If it exists, 'target_patch' is a
%     two-element vector containing the row and column, respectively, of
%     the top-left corner of the single image region to be processed. When
%     `target_patch` is passed, output arguments are calculated for a
%     single image region, rather than for the entire image.
%
% optimizer -- Image registration optimizer
%   The `optimizer` input argument of 'imregtform()'.
%
% metric -- Image registration similarity metric
%   The `metric` input argument of 'imregtform()'.
%
% pyramid_levels -- Number of image pyramid levels
%   The `'PyramidLevels'` input argument of 'imregtform()'.
%
% bands -- Input spectral sampling
%   A vector containing the wavelengths or colour channel indices
%   corresponding to the third dimension of `I`, used to evaluate
%   `dispersionfun`. Either both `bands` and `dispersionfun` must be
%   passed, or neither must be passed.
%
% dispersionfun -- Initial guess for image registration
%   `dispersionfun` must be a function handle, such as produced by
%   'makeDispersionfun()'. `dispersionfun(X)`, where `X` is a three-element
%   row vector (x, y, l), returns the dispersion vector for the position
%   (x, y) in the image for light with wavelength or colour channel index
%   `l`. The dispersion vector originates from position (x, y) in a
%   reference spectral band or colour channel. This function does not
%   assume that the reference channel/band is the same as that referred to
%   by `reference_index`.
%
%   This function uses `dispersionfun` to initialize image registration
%   transformations.
%
% verbose -- Verbosity flag
%   If `true`, console output will be displayed to show progress
%   information.
%
% ## Output Arguments
%
% centers -- Registered region centres
%   A cell vector of structure vectors, where each element of each
%   structure vector has a field, 'center', storing a two-element vector
%   of the x and y-coordinates of the centre of an image region.
%
%   `centers` has length 'm', and each cell contains a structure vector of
%   length 'n'. 'n' is the number of regions, determined by
%   `patch_options.patch_size`. 'm' is equal to the size of `I` in its
%   third dimension. `centers{reference_index}(k)` is the position of the
%   centre of the k-th image region in the input image, whereas
%   `centers{j}(k)` is the position of the k-th image region in the j-th
%   colour channel or spectral band of the input image, when this
%   channel/band is registered with the reference channel/band.
%
% ## Future Work
% - It may be possible to extend this function to process colour-filter
%   array images without demosaicing, to find dispersion between colour
%   channels. A similar extension would be to process 3D images created by
%   stacking colour-filter array images captured under bandpass-filtered
%   illuminations, to find dispersion resulting from different filtered
%   illuminations. In either case, the idea is to imitate
%   RAWDiskDispersion.m. Unfortunately, these extensions would require
%   implementing a custom image registration subroutine.
%
% See also imregtform, imregconfig, findAndFitDisks

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 17, 2019

narginchk(6, 9);
nargoutchk(1, 1);

verbose = false;
has_dispersion = false;
dispersionfun = [];
bands = [];
if ~isempty(varargin)
    if islogical(varargin{1})
        verbose = varargin{1};
        if length(varargin) > 1
            error('`verbose` must be the last optional input argument.');
        end
    elseif length(varargin) == 1 || islogical(varargin{2})
        error('Both `bands` and `disperisonfun` must be passed together.');
    else
        has_dispersion = true;
        bands = varargin{1};
        dispersionfun = varargin{2};
        if length(varargin) > 2
            if islogical(varargin{3})
                verbose = varargin{3};
            else
                error('`verbose` must be the last optional input argument.');
            end
        end
    end
end

if verbose
    tic
end

% Input argument parsing
do_single_patch = isfield(patch_options, 'target_patch');

n_bands = size(I, 3);
if n_bands == 1
    error('The input image must have more than one channel/band.');
end
if has_dispersion
    if ~isa(dispersionfun, 'function_handle')
        error('`dispersionfun` must be a function handle.');
    end
    if ~isvector(bands)
        error('`bands` must be a vector.');
    end
    if n_bands ~= length(bands) 
        error('The length of `bands` must equal the size of `I` in its third dimension.');
    end
end

image_sampling = [size(I, 1), size(I, 2)];
patch_size = patch_options.patch_size;
padding = patch_options.padding;

if do_single_patch
    target_patch = patch_options.target_patch;
    if target_patch(1) < 1 || target_patch(1) > image_sampling(1) ||...
       target_patch(2) < 1 || target_patch(2) > image_sampling(2)
        error('The target patch is outside of the image bounds.');
    end
end

if verbose
    disp('Splitting the input image into columns...');
end

% Channel indices in the input concatenation of images
channels_in.I_in = [ 1, size(I, 3) ];
n_channels_in = channels_in.I_in(2);

% Divide the input images into columns which will be sent to individual
% parallel workers
if do_single_patch
    n_i = 1;
    n_j = 1;
    patch_offset = target_patch - 1;
    display_optimization = verbose;
else
    n_i = ceil(image_sampling(1) / patch_size(1));
    n_j = ceil(image_sampling(2) / patch_size(2));
    patch_offset = [0, 0];
    display_optimization = false;
end
n_patches = n_i * n_j;
columns_in = cell(1, n_j);
for j = 1:n_j
    cols_ind_in = [
        max((j - 1) * patch_size(2) + 1 - padding + patch_offset(2), 1);
        min(j * patch_size(2) + padding + patch_offset(2), image_sampling(2))
    ];
    columns_in{j} = zeros(image_sampling(1), diff(cols_ind_in) + 1, n_channels_in);
    columns_in{j}(:, :, channels_in.I_in(1):channels_in.I_in(2)) = I(:, cols_ind_in(1):cols_ind_in(2), :);
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
    column_out_j = zeros(n_i, n_bands * 2);
    
    % Process each patch within the column
    for i = 1:n_i
        corner(1) = (i - 1) * patch_size(1) + 1 + patch_offset(1);
        patch_lim_rows = [
            max(corner(1) - padding, 1);
            min(corner(1) + patch_size(1) + padding - 1, image_sampling(1));
        ];
        image_sampling_p(1) = diff(patch_lim_rows) + 1;

        % Generate the patch
        patches_I_ij = column_in_j(patch_lim_rows(1):patch_lim_rows(2), :, channels_in.I_in(1):channels_in.I_in(2));
        patches_I_ij_3D = reshape(patches_I_ij, [image_sampling_p n_bands]);
        
        center_reference = [
            corner(2) + (patch_size(2) - 1) / 2, corner(1) + (patch_size(1) - 1) / 2
        ];
    
        if has_dispersion
            d = dispersionfun([repmat(center_reference, n_bands, 1), bands]);
            d = d - repmat(d(reference_index, :), n_bands, 1);
        else
            d = zeros(n_bands, 2);
        end
        
        % Patch registration
        for b = 1:n_bands
            if b ~= reference_index
                tform = imregtform(...
                    patches_I_ij_3D(:, :, reference_index),...
                    patches_I_ij_3D(:, :, b),...
                    'translation', optimizer, metric,...
                    'DisplayOptimization', display_optimization,...
                    'InitialTransformation', affine2d([
                            1, 0, 0;
                            0, 1, 0;
                            d(b, 1), d(b, 2), 1
                        ]),...
                    'PyramidLevels', pyramid_levels...
                );
                column_out_j(i, 2 * b - 1) = tform.T(end, 1) + center_reference(1);
                column_out_j(i, 2 * b) = tform.T(end, 2) + center_reference(2);
            else
                column_out_j(i, 2 * b - 1) = center_reference(1);
                column_out_j(i, 2 * b) = center_reference(2);
            end
        end
        
        if verbose
            fprintf('\tProcessed region %d of %d\n', i + (j-1) * n_i, n_patches);
        end
    end
    columns_out{j} = column_out_j;
end

if verbose
    fprintf('\tDone.\n');
    disp('Recombining results from regions...');
end

% Recombine results
centers = cell(n_bands, 1);
for b = 1:n_bands
    centers{b} = struct('center', num2cell(zeros(n_patches, 2), 2));
    for j = 1:n_j
        centers_cell = num2cell(columns_out{j}(:, (2 * b - 1):(2 * b)), 2);
        [centers{b}((1 + n_i * (j - 1)):(n_i * j)).center] = centers_cell{:};
    end
end

if verbose
    fprintf('\tDone.\n');
    toc
end
    
end