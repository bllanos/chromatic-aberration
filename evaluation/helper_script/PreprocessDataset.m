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
    'illuminant_name',...
    'normalization_channel',...
}];

% ## Parameters for creating radiance images

% CIE D-illuminant
illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
illuminant_temperature = 6504; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
illuminant_name = 'd65';

% Colour channel to use for radiance normalization
normalization_channel = 2;

% ## Operational parameters

% Patch size and padding to use for converting spectral images to colour
% and RAW images
imageFormationOptions.patch_size = [100, 100];
imageFormationOptions.padding = 10;

%% RGB colour space

n_channels_rgb = 3;
bands_rgb = 1:n_channels_rgb;
sensor_map_rgb = eye(n_channels_rgb);

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

bands = [];
color_weights_reference = [];
bands_variable = 'bands';
if has_spectral
    spectral_filenames = listFiles(dp.spectral_images_wildcard);
    n_images_spectral = length(spectral_filenames);
    n_images = n_images_spectral;
    names = trimCommon(spectral_filenames);
    
    load(dp.wavelengths, bands_variable);
    if isempty(bands)
        error('No wavelength band information is associated with the spectral images.')
    end
    bands_spectral = bands;
end
if has_rgb
    rgb_filenames = listFiles(dp.rgb_images_wildcard);
    n_images_rgb = length(rgb_filenames);
    rgb_names = trimCommon(rgb_filenames);
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
    raw_names = trimCommon(raw_filenames);
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

bands = [];
if has_color_map
    % Load colour conversion data
    
    model_variables_required = { 'sensor_map', 'channel_mode' };
    load(dp.color_map, model_variables_required{:}, bands_variable);
    if ~all(ismember(model_variables_required, who))
        error('One or more of the required colour space conversion variables is not loaded.')
    end
    if isempty(bands)
        error('No wavelength band information is associated with the colour conversion data.')
    end
    bands_color = bands;
    
    if channel_mode
        if (length(bands_spectral) ~= length(bands_color)) ||...
                any(bands_spectral ~= bands_color)
            error(['For a colour channel-based colour conversion, the true l'...
                'atent images must be defined at the same colour channels a'...
                's the colour conversion data.']);
        end
        color_weights = sensor_map;
        spectral_weights = eye(length(bands_color));
        color_weights_reference = color_weights;
    elseif has_spectral
        [color_weights, spectral_weights, bands, color_weights_reference] = samplingWeights(...
          sensor_map, bands_color, bands_spectral, samplingWeightsOptions, samplingWeightsVerbose...
        );
    else
        [color_weights, ~, bands] = samplingWeights(...
          sensor_map, bands_color, bands_color, samplingWeightsOptions, samplingWeightsVerbose...
        );
    end
    n_bands = length(bands);
end

%% Load dispersion models

has_dispersion_rgb = ~isempty(dp.dispersion_rgb_forward) && ~isempty(dp.dispersion_rgb_reverse);
has_dispersion_spectral = ~isempty(dp.dispersion_spectral_reverse);

if has_dispersion_rgb
    [...
        dd_rgb_forward, bands_dispersionfun, td_rgb_forward...
    ] = loadDispersionModel(dp.dispersion_rgb_forward, true, false);
    if (n_channels_rgb ~= length(bands_dispersionfun)) ||...
       any(bands_rgb ~= bands_dispersionfun)
        error('The forward model of colour dispersion does not have the right colour channels.');
    end
    [...
        dd_rgb_reverse, bands_dispersionfun, td_rgb_reverse...
    ] = loadDispersionModel(dp.dispersion_rgb_reverse, false, false);
    if (n_channels_rgb ~= length(bands_dispersionfun)) ||...
       any(bands_rgb ~= bands_dispersionfun)
        error('The reverse model of colour dispersion does not have the right colour channels.');
    end
end
if has_dispersion_spectral
    [...
        dd_spectral_reverse, ~, td_spectral_reverse...
    ] = loadDispersionModel(dp.dispersion_spectral_reverse, false, false);
    if channel_mode && ...
       ((length(bands_color) ~= length(bands_dispersionfun)) ||...
       any(bands_color ~= bands_dispersionfun))
        error('When estimating colour images, the same colour channels must be used by the model of dispersion.');
    end
end

%% Load the illuminant and create a conversion matrix from reflectances to radiances

if dp.spectral_reflectances
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
        normalization_channel, samplingWeightsOptions.int_method...
    );
end

%% Check for comparison methods

has_choi_rgb = isfield(dp, 'choi_rgb_wildcard');
has_choi_spectral = isfield(dp, 'choi_spectral_wildcard');

if has_choi_rgb
    image_type = 'colour';
    choi_rgb_filenames = listFiles(dp.choi_rgb_wildcard);
    n_choi = length(choi_rgb_filenames);
    choi_names = trimCommon(choi_rgb_filenames);
    if (n_images ~= n_choi)
        error('Expected %d %s images from Choi et al. 2017, not %d.', n_images, image_type, n_choi);
    end
    if n_images > 1
        for i = 1:n_images
            str_found = strfind(choi_names{i}, names{i});
            if length(str_found) ~= 1 || str_found ~= 1
                error('The %s image from Choi et al. 2017,\n"%s", does not have the name of the %d-th dataset image,\n"".',...
                    image_type, choi_rgb_filenames{i}, names{i});
            end
        end
    end
end

if has_choi_spectral
    image_type = 'spectral';
    choi_spectral_filenames = listFiles(dp.choi_spectral_wildcard);
    n_choi = length(choi_spectral_filenames);
    choi_names = trimCommon(choi_spectral_filenames);
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