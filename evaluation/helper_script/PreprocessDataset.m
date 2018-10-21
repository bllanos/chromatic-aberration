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
imageFormationOptions.int_method = int_method;

%% RGB colour space

n_channels_rgb = 3;
bands_rgb = 1:n_channels_rgb;
sensor_map_rgb = eye(n_channels_rgb);

%% Load the illuminant

illuminant_data = csvread(illuminant_filename);
bands_illuminant = illuminant_data(:, 1);
S_illuminant = illuminant_data(:, 2:end);
spd_illuminant = ciedIlluminant(...
    illuminant_temperature, bands_illuminant, S_illuminant, bands_illuminant...
);

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

bands_color = [];
if has_color_map
    % Load colour conversion data
    bands_script = bands;
    bands = [];
    
    optional_variable = 'bands';
    model_variables_required = { 'sensor_map', 'channel_mode' };
    load(dp.color_map, model_variables_required{:}, optional_variable);
    if ~all(ismember(model_variables_required, who))
        error('One or more of the required colour space conversion variables is not loaded.')
    end
    bands_color = bands;
    bands = bands_script;
    
    % Compare with colour space conversion data
    n_bands = length(bands);
    n_bands_sensor_map = size(sensor_map, 2);
    resample_bands = false;
    if ~isempty(bands_color)
        bands_for_interp = bands_color;
        if n_bands ~= length(bands_color) || any(bands ~= bands_color)
            % Resampling is needed
            resample_bands = true;
        end
    elseif n_bands_sensor_map ~= n_bands
        % Resampling is needed, but will be "blind"
        resample_bands = true;
        bands_for_interp = linspace(bands(1), bands(end), n_bands_sensor_map);
    end
    % Resample colour space conversion data if necessary
    if resample_bands
        [sensor_map_resampled, bands] = resampleArrays(...
            bands_for_interp, sensor_map.', bands,...
            bands_interp_method...
            );
        n_bands = length(bands);
        sensor_map_resampled = sensor_map_resampled.';
    else
        sensor_map_resampled = sensor_map;
    end
end

bands_spectral = [];
sensor_map_spectral = [];
if has_spectral
    spectral_filenames = listFiles(dp.spectral_images_wildcard);
    n_images_spectral = length(spectral_filenames);
    n_images = n_images_spectral;
    names = trimCommon(spectral_filenames);
    
    bands = [];
    load(dp.wavelengths, optional_variable);
    if isempty(bands)
        error('No wavelength band information is associated with the spectral images.')
    end
    bands_spectral = bands;
    bands = bands_script;
    
    if has_color_map
        % Check if quantitative evaluation of spectral images is possible
        can_evaluate_spectral = (length(bands_spectral) == n_bands) && all(bands_spectral == bands); 
        if can_evaluate_spectral
            sensor_map_spectral = sensor_map_resampled;
        else
            warning(['Quantitative evaluation of latent images is not possibl'...
                'e because they will be produced at different wavelength'...
                's from the true latent images.']...
            );
        
            % Allow for conversion to colour images
            [sensor_map_spectral, bands_spectral] = resampleArrays(...
                bands_for_interp, sensor_map.', bands_spectral,...
                bands_interp_method...
                );
            sensor_map_spectral = sensor_map_spectral.';
        end
    else
        can_evaluate_spectral = false;
    end
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

%% Load dispersion models

has_dispersion_rgb = ~isempty(dp.dispersion_rgb_forward) && ~isempty(dp.dispersion_rgb_reverse);
has_dispersion_spectral = ~isempty(dp.dispersion_spectral_reverse);

if has_dispersion_rgb
    [...
        dd_rgb_forward, ~, td_rgb_forward...
    ] = loadDispersionModel(dp.dispersion_rgb_forward, true, false);
    [...
        dd_rgb_reverse, ~, td_rgb_reverse...
    ] = loadDispersionModel(dp.dispersion_rgb_reverse, false, false);
end
if has_dispersion_spectral
    [...
        dd_spectral_reverse, ~, td_spectral_reverse...
    ] = loadDispersionModel(dp.dispersion_spectral_reverse, false, false);
end