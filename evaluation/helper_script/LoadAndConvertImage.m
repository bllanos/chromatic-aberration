%% Image loading and colour conversion helper script
% Common code initially extracted from 'RunOnDataset.m' and
% 'SelectWeightsForDataset.m'

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 10, 2018

if has_spectral
    I_spectral_gt = loadImage(spectral_filenames{i}, dp.spectral_images_variable);

    % Convert to radiance images, if required
    if dp.spectral_reflectances
        I_spectral_gt = channelConversion(I_spectral_gt, radiance_normalized_weights, 3);
    end

    if has_dispersion_spectral
        [df_spectral_reverse, I_spectral_gt] = makeDispersionForImage(...
            dd_spectral_reverse, I_spectral_gt, td_spectral_reverse...
        );
    else
        df_spectral_reverse = [];
    end
    image_sampling = [size(I_spectral_gt, 1), size(I_spectral_gt, 2)];
elseif has_color_map && has_dispersion_spectral
    df_spectral_reverse = makeDispersionForImage(dd_spectral_reverse);
elseif has_color_map
    df_spectral_reverse = [];
end

if has_spectral && has_color_map
    [...
        I_rgb_gt_simulated, I_rgb_gt_warped, I_raw_gt_simulated,...
        I_spectral_gt_warped...
    ] = imageFormation(...
        I_spectral_gt, color_weights_reference, imageFormationOptions,...
        df_spectral_reverse, bands_spectral, bayer_pattern...
    );
else
    I_rgb_gt_warped = [];
    I_spectral_gt_warped = [];
end

if has_rgb
    I_rgb_gt = loadImage(rgb_filenames{i}, dp.rgb_images_variable);
    if has_dispersion_rgb
        [df_rgb_reverse, I_rgb_gt] = makeDispersionForImage(...
                dd_rgb_reverse, I_rgb_gt, td_rgb_reverse...
            );
        df_rgb_forward = makeDispersionForImage(...
                dd_rgb_forward, I_rgb_gt, td_rgb_forward...
            );
    else
        df_rgb_reverse = [];
        df_rgb_forward = [];
    end
    if has_spectral
        if any([size(I_rgb_gt, 1), size(I_rgb_gt, 2)] ~= image_sampling)
            error('The colour version of %s has different spatial dimensions from its spectral version.',...
                names{i}...
            );
        end
    else
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

    % Crop to the region of valid dispersion
    roi = [];
    if has_spectral && has_dispersion_spectral
        roi = modelSpaceTransform(...
            [size(I_raw_gt, 1), size(I_raw_gt, 2)], td_spectral_reverse.model_space, td_spectral_reverse.fill...
            );
    elseif has_rgb && has_dispersion_rgb
        roi = modelSpaceTransform(...
            [size(I_raw_gt, 1), size(I_raw_gt, 2)], td_rgb_reverse.model_space, td_rgb_reverse.fill...
            );
    end
    if ~isempty(roi)
        I_raw_gt = I_raw_gt(roi(1):roi(2), roi(3):roi(4), :);
    end

    if any([size(I_raw_gt, 1), size(I_raw_gt, 2)] ~= image_sampling)
        error('The RAW version of %s has different spatial dimensions from its other versions.',...
            names{i}...
            );
    end
else
    if has_rgb
        if has_dispersion_rgb
            if verbose
                fprintf('[image %d] Calculating the reverse colour dispersion matrix...\n', i);
            end
            W_reverse = dispersionfunToMatrix(...
                df_rgb_reverse, bands_rgb, image_sampling, image_sampling,...
                [0, 0, image_sampling(2),  image_sampling(1)], true...
                );
            if verbose
                fprintf('\t...done\n');
            end
            I_rgb_warped = warpImage(I_rgb_gt, W_reverse, image_sampling);
        else
            I_rgb_warped = I_rgb_gt;
        end
        I_raw_gt = mosaic(I_rgb_warped, bayer_pattern);
    end
end