%% Dataset preprocessing helper script
% Common code initially extracted from 'RunOnDataset.m' and
% 'SelectWeightsForDataset.m'

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created September 10, 2018

%% List of parameters to save with results
if ~exist('parameters_list', 'var')
    error('`parameters_list` should be initialized prior to running PreprocessDataset.m');
end
parameters_list = [parameters_list, {
    'illuminant_filename',...
    'illuminant_temperature',...
    'normalization_channel',...
}];

% ## Parameters for creating radiance images

% CIE D-illuminant
illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
illuminant_temperature = 6504; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D

% Colour channel to use for radiance normalization
normalization_channel = 2;

%% RGB colour space

n_channels_rgb = 3;
bands_rgb = 1:n_channels_rgb;

%% Find and/or prepare to generate the dataset images

has_raw = ~isempty(dp.raw_images_wildcard);
has_rgb = ~isempty(dp.rgb_images_wildcard);
has_spectral = ~isempty(dp.spectral_images_wildcard);
has_color_map = ~isempty(dp.color_map);

if ~has_spectral && ~has_rgb
    error('Ground truth (colour or spectral images) must be provided for evaluation.')
end
if ~has_raw && ~has_rgb && (~has_spectral || (has_spectral && ~has_color_map))
    error('RAW images are not available, and cannot be generated.')
end
if has_spectral && ~has_color_map
    error('Spectral images are provided, but not the camera''s spectral sensitivity data.');
end

if has_spectral
    spectral_filenames = listFiles(dp.spectral_images_wildcard);
    n_images_spectral = length(spectral_filenames);
    n_images = n_images_spectral;
    names = trimCommon(spectral_filenames, [false, true]);
    
    bands_spectral = loadVariables(dp.wavelengths, 'bands');
end
if has_rgb
    rgb_filenames = listFiles(dp.rgb_images_wildcard);
    n_images_rgb = length(rgb_filenames);
    rgb_names = trimCommon(rgb_filenames, [false, true]);
    if has_spectral
        if (n_images_spectral ~= n_images_rgb)
            error('Mismatched number of spectral and colour images.');
        end
        if n_images_spectral > 1
            for i = 1:n_images_rgb
                if ~strcmp(rgb_names{i}, names{i})
                    error('Not all spectral image filenames and colour image filenames match.');
                end
            end
        end
    else
        names = rgb_names;
    end
    n_images = n_images_rgb;
end
if has_raw
    raw_filenames = listFiles(dp.raw_images_wildcard);
    n_images_raw = length(raw_filenames);
    raw_names = trimCommon(raw_filenames, [false, true]);
    if (n_images ~= n_images_raw)
        error('Mismatched number of spectral/colour and RAW images.');
    end
    if n_images > 1
        for i = 1:n_images
            if ~strcmp(raw_names{i}, names{i})
                error('Not all spectral/colour image filenames and RAW image filenames match.');
            end
        end
    end
end

if has_color_map
    % Load colour conversion data
    [sensor_map, channel_mode, bands_color] = loadColorMap(dp.color_map);
    
    if channel_mode
        if (length(bands_spectral) ~= length(bands_color)) ||...
                any(bands_spectral ~= bands_color)
            error(['For a colour channel-based colour conversion, the true l'...
                'atent images must be defined at the same colour channels a'...
                's the colour conversion data.']);
        end
        spectral_weights = eye(length(bands_color));
    else
        if isfield(dp, 'fix_bands') && dp.fix_bands
            findSamplingOptions.power_threshold = 1;
            findSamplingOptions.n_bands = 0;
            findSamplingOptions.support_threshold = 0;
            solvePatchesSpectralOptions.sampling_options.power_threshold = findSamplingOptions.power_threshold;
            solvePatchesSpectralOptions.sampling_options.n_bands = findSamplingOptions.n_bands;
            solvePatchesSpectralOptions.sampling_options.support_threshold = findSamplingOptions.support_threshold;
        end
        if has_spectral
            [~, spectral_weights, bands] = findSampling(...
              sensor_map, bands_color, bands_spectral, findSamplingOptions, findSamplingVerbose...
            );
        else
            [~, ~, bands] = findSampling(...
              sensor_map, bands_color, bands_color, findSamplingOptions, findSamplingVerbose...
            );
        end
    end
    n_bands = length(bands);
end

%% Load dispersion models

has_dispersion_rgb = ~isempty(dp.dispersion_rgb_forward) && ~isempty(dp.dispersion_rgb_reverse);
has_dispersion_spectral = ~isempty(dp.dispersion_spectral_forward) && ~isempty(dp.dispersion_spectral_reverse);
evaluate_aberrated_rgb = (has_dispersion_rgb || has_dispersion_spectral) && dp.is_aberrated;
evaluate_aberrated_spectral = has_dispersion_spectral && dp.is_aberrated;

% For regularization weight selection
use_warped_rgb = criteria(mse_index) && has_rgb && has_dispersion_rgb && dp.is_aberrated;
use_warped_spectral = criteria(mse_index) && has_spectral && has_dispersion_spectral && dp.is_aberrated;

if has_dispersion_rgb
    [...
        dd_rgb_forward, bands_dispersionfun, td_rgb_forward...
    ] = loadDispersionModel(dp.dispersion_rgb_forward, true);
    if (n_channels_rgb ~= length(bands_dispersionfun)) ||...
       any(bands_rgb ~= bands_dispersionfun)
        error('The forward model of colour dispersion does not have the right colour channels.');
    end
    [...
        dd_rgb_reverse, bands_dispersionfun, td_rgb_reverse...
    ] = loadDispersionModel(dp.dispersion_rgb_reverse, false);
    if (n_channels_rgb ~= length(bands_dispersionfun)) ||...
       any(bands_rgb ~= bands_dispersionfun)
        error('The reverse model of colour dispersion does not have the right colour channels.');
    end
end
if has_dispersion_spectral
    [...
        dd_spectral_forward, bands_dispersionfun, td_spectral_forward...
    ] = loadDispersionModel(dp.dispersion_spectral_forward, true);
    if channel_mode && ...
       ((length(bands_color) ~= length(bands_dispersionfun)) ||...
       any(bands_color(:) ~= bands_dispersionfun(:)))
        error('When estimating colour images, the same colour channels must be used by the forward model of dispersion.');
    end
    [...
        dd_spectral_reverse, bands_dispersionfun, td_spectral_reverse...
    ] = loadDispersionModel(dp.dispersion_spectral_reverse, false);
    if channel_mode && ...
       ((length(bands_color) ~= length(bands_dispersionfun)) ||...
       any(bands_color(:) ~= bands_dispersionfun(:)))
        error('When estimating colour images, the same colour channels must be used by the reverse model of dispersion.');
    end
end

%% Load the illuminant and create a conversion matrix from reflectances to radiances

if has_spectral && dp.spectral_reflectances
    illuminant_data = csvread(illuminant_filename);
    bands_illuminant = illuminant_data(:, 1);
    S_illuminant = illuminant_data(:, 2:end);
    spd_illuminant = ciedIlluminant(...
        illuminant_temperature, bands_illuminant, S_illuminant, bands_illuminant...
    );

    reflectances = eye(length(bands_spectral));

    [bands_spectral, ~, radiance_normalized_weights] = reflectanceToRadiance(...
        bands_illuminant, spd_illuminant,...
        bands_spectral, reflectances,...
        bands_color, sensor_map.',...
        normalization_channel, findSamplingOptions.int_method...
    );
end

%% Check for comparison methods

has_choi_rgb = isfield(dp, 'choi_rgb_wildcard');
has_choi_spectral = isfield(dp, 'choi_spectral_wildcard');

if has_choi_rgb
    image_type = 'colour';
    choi_rgb_filenames = listFiles(dp.choi_rgb_wildcard);
    n_choi = length(choi_rgb_filenames);
    choi_names = trimCommon(choi_rgb_filenames, [false, true]);
    if (n_images ~= n_choi)
        error('Expected %d %s images from Choi et al. 2017, not %d.', n_images, image_type, n_choi);
    end
    if n_images > 1
        for i = 1:n_images
            str_found = strfind(choi_names{i}, names{i});
            if length(str_found) ~= 1 || str_found ~= 1
                error('The %s image from Choi et al. 2017,\n"%s", does not have the name of the %d-th dataset image,\n"%s".',...
                    image_type, choi_rgb_filenames{i}, i, names{i});
            end
        end
    end
end

if has_choi_spectral
    image_type = 'spectral';
    choi_spectral_filenames = listFiles(dp.choi_spectral_wildcard);
    n_choi = length(choi_spectral_filenames);
    choi_names = trimCommon(choi_spectral_filenames, [false, true]);
    if (n_images ~= n_choi)
        error('Expected %d %s images from Choi et al. 2017, not %d.', n_images, image_type, n_choi);
    end
    if n_images > 1
        for i = 1:n_images
            str_found = strfind(choi_names{i}, names{i});
            if length(str_found) ~= 1 || str_found ~= 1
                error('The %s image from Choi et al. 2017,\n"%s", does not have the name of the %d-th dataset image,\n"".',...
                    image_type, choi_spectral_filenames{i}, names{i});
            end
        end
    end
end

%% Miscellaneous

solvePatchesColorOptions.patch_options.patch_size = dp.patch_size;
solvePatchesColorOptions.patch_options.padding = dp.padding;
solvePatchesSpectralOptions.patch_options.patch_size = dp.patch_size;
solvePatchesSpectralOptions.patch_options.padding = dp.padding;
imageFormationPatchOptions.patch_size = dp.patch_size;
imageFormationPatchOptions.padding = dp.padding;