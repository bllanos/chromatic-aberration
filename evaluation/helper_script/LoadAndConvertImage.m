%% Image loading and colour conversion helper script
% Common code initially extracted from 'RunOnDataset.m' and
% 'SelectWeightsForDataset.m'

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 10, 2018

% Make sure dispersion models give the same region of interest
if has_dispersion_rgb
    td_image = mergeDispersionModelROI(td_rgb_reverse, td_rgb_forward);
end
if has_dispersion_spectral && has_dispersion_rgb
    td_image = mergeDispersionModelROI(td_spectral_reverse, td_spectral_forward, td_image);
elseif has_dispersion_spectral
    td_image = mergeDispersionModelROI(td_spectral_reverse, td_spectral_forward);
end

if has_spectral
    I_spectral_gt = loadImage(spectral_filenames{i}, dp.spectral_images_variable);
    image_sampling_full = [size(I_spectral_gt, 1), size(I_spectral_gt, 2)];

    % Convert to radiance images, if required
    if dp.spectral_reflectances
        I_spectral_gt = channelConversion(I_spectral_gt, radiance_normalized_weights, 3);
    end

    if has_dispersion_spectral
        [df_spectral_reverse, I_spectral_gt] = makeDispersionForImage(...
            dd_spectral_reverse, I_spectral_gt, td_image, true...
        );
        df_spectral_forward = makeDispersionForImage(...
            dd_spectral_forward, I_spectral_gt, td_image, true...
        );
    else
        df_spectral_reverse = [];
        df_spectral_forward = [];
    end
    image_sampling = [size(I_spectral_gt, 1), size(I_spectral_gt, 2)];
elseif has_color_map && has_dispersion_spectral
    df_spectral_reverse = makeDispersionForImage(dd_spectral_reverse);
    df_spectral_forward = makeDispersionForImage(dd_spectral_forward);
elseif has_color_map
    df_spectral_reverse = [];
    df_spectral_forward = [];
end

if has_spectral && has_color_map
    if dp.is_aberrated
        df_spectral_reverse_imageFormation = [];
    else
        df_spectral_reverse_imageFormation = df_spectral_reverse;
    end
    [...
        I_rgb_gt_simulated, I_rgb_gt_warped, I_raw_gt_simulated,...
        I_spectral_gt_warped...
    ] = imageFormation(...
        I_spectral_gt, bands_spectral, sensor_map, bands_color,...
        imageFormationSamplingOptions, imageFormationPatchOptions,...
        df_spectral_reverse_imageFormation, bayer_pattern...
    );
else
    I_rgb_gt_warped = [];
    I_spectral_gt_warped = [];
end

if has_rgb
    I_rgb_gt = loadImage(rgb_filenames{i}, dp.rgb_images_variable);
    if has_spectral
        if any([size(I_rgb_gt, 1), size(I_rgb_gt, 2)] ~= image_sampling_full)
            error('The colour version of %s has different spatial dimensions from its spectral version.',...
                names{i}...
            );
        end
    else
        image_sampling_full = [size(I_rgb_gt, 1), size(I_rgb_gt, 2)];
    end
    if has_dispersion_rgb
        [df_rgb_reverse, I_rgb_gt] = makeDispersionForImage(...
                dd_rgb_reverse, I_rgb_gt, td_image, true...
            );
        df_rgb_forward = makeDispersionForImage(...
                dd_rgb_forward, I_rgb_gt, td_image, true...
            );
    else
        df_rgb_reverse = [];
        df_rgb_forward = [];
    end
    if ~has_spectral
        image_sampling = [size(I_rgb_gt, 1), size(I_rgb_gt, 2)];
    end
elseif has_spectral && has_color_map
    I_rgb_gt = I_rgb_gt_simulated;
    I_raw_gt = I_raw_gt_simulated;
end
if ~has_rgb && has_dispersion_rgb
    df_rgb_reverse = makeDispersionForImage(dd_rgb_reverse);
    df_rgb_forward = makeDispersionForImage(dd_rgb_forward);
elseif ~has_dispersion_rgb
    df_rgb_reverse = [];
    df_rgb_forward = [];
end

if has_raw
    I_raw_gt = loadImage(raw_filenames{i}, dp.raw_images_variable);
    if ~ismatrix(I_raw_gt)
        error('Expected a RAW image, represented as a 2D array, not a higher-dimensional array.');
    end
    if any([size(I_raw_gt, 1), size(I_raw_gt, 2)] ~= image_sampling_full)
        error('The RAW version of %s has different spatial dimensions from its other versions.',...
            names{i}...
            );
    end
end

% Convert image patches and lines to the new coordinate system
roi = [];
if has_dispersion_rgb || has_dispersion_spectral
    roi = modelSpaceTransform(...
        image_sampling_full, td_image.model_space, td_image.fill, true...
        );
    if isfield(dp, 'params_patches') && isfield(dp.params_patches, names{i})
        params_patches_i = dp.params_patches.(names{i});
        params_patches_i(:, 1) = params_patches_i(:, 1) - roi(3) + 1;
        params_patches_i(:, 2) = params_patches_i(:, 2) - roi(1) + 1;
        dp.params_patches.(names{i}) = params_patches_i;
        % Checks for invalid patch locations will be done in SelectWeightsForDataset.m
    end
    if isfield(dp.evaluation.custom_spectral, names{i})
        eval_options = dp.evaluation.custom_spectral.(names{i});
        if isfield(eval_options, 'radiance')
            eval_options_radiance = eval_options.radiance;
            eval_options_radiance(:, 1) = eval_options_radiance(:, 1) - roi(3) + 1;
            eval_options_radiance(:, 2) = eval_options_radiance(:, 2) - roi(1) + 1;
            eval_options_radiance_half_sizes = (eval_options_radiance(:, 3:4) - 1) / 2;
            if any(...
                    ((eval_options_radiance(:, 1) - eval_options_radiance_half_sizes(:, 1)) < 1) | ...
                    ((eval_options_radiance(:, 1) + eval_options_radiance_half_sizes(:, 1)) > image_sampling(2)) | ...
                    ((eval_options_radiance(:, 2) - eval_options_radiance_half_sizes(:, 2)) < 1) | ...
                    ((eval_options_radiance(:, 2) + eval_options_radiance_half_sizes(:, 2)) > image_sampling(1))...
                )
                error('The evaluation patches for image %s are outside the region of valid dispersion.',...
                    names{i}...
                    );
            end
            dp.evaluation.custom_spectral.(names{i}).radiance = eval_options_radiance;
        end
        if isfield(eval_options, 'scanlines')
            eval_options_scanlines = eval_options.scanlines;
            eval_options_scanlines(:, [1 3]) = eval_options_scanlines(:, [1 3]) - roi(3) + 1;
            eval_options_scanlines(:, [2 4]) = eval_options_scanlines(:, [2 4]) - roi(1) + 1;
            if any(any(eval_options_scanlines < 1)) || any(any(...
                    (eval_options_scanlines(:, [1 3]) > image_sampling(2)) |...
                    (eval_options_scanlines(:, [2 4]) > image_sampling(1))...
                    ))
                error('The evaluation lines for image %s are outside the region of valid dispersion.',...
                    names{i}...
                    );
            end
            dp.evaluation.custom_spectral.(names{i}).scanlines = eval_options_scanlines;
        end
        if isfield(eval_options, 'reference_patch')
            reference_patch = eval_options.reference_patch;
            reference_patch(1) = reference_patch(1) - roi(3) + 1;
            reference_patch(2) = reference_patch(2) - roi(1) + 1;
            reference_patch_half_size = (reference_patch(3:4) - 1) / 2;
            if any((reference_patch(:, 1:2) - reference_patch_half_size) < 1) || ...
               ((reference_patch(1) + reference_patch_half_size(1)) > image_sampling(2)) || ...
               ((reference_patch(2) + reference_patch_half_size(2)) > image_sampling(1))
                error('The reference patch for image %s is outside the region of valid dispersion.',...
                    names{i}...
                    );
            end
            dp.evaluation.custom_spectral.(names{i}).reference_patch = reference_patch;
        end
    end
end

if has_raw
    % Crop to the region of valid dispersion
    if ~isempty(roi)
        I_raw_gt = I_raw_gt(roi(1):roi(2), roi(3):roi(4), :);
    end
else
    if has_rgb
        if has_dispersion_rgb && ~dp.is_aberrated
            if verbose
                fprintf('[image %d] Warping the true color image...\n', i);
            end
            dispersion_options = struct('bands_in', bands_rgb);
            I_rgb_warped = dispersionfunToMatrix(...
                df_rgb_reverse, dispersion_options, I_rgb_gt, true...
            );
            if verbose
                fprintf('\t...done\n');
            end
        else
            I_rgb_warped = I_rgb_gt;
        end
        I_raw_gt = mosaic(I_rgb_warped, bayer_pattern);
    end
end