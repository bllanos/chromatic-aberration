%% Aberrated image generation
% Create ideal, and warped hyperspectral images and their RGB equivalents,
% from distributions of pairs of object spectral reflectances.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Input images
%
% #### Chromaticity maps
%
% A filepath wildcard pointing to either monochromatic or RGB images. Each
% image will be used to generate a chromaticity map by soft clustering of
% its colours. A monochromatic image will be used directly as the input for
% clustering, whereas an RGB image will be converted to the CIE 1976 L*a*b*
% colour space, after which its chromaticity values (a* and b*) will be
% reduced to a single monochromatic image using principal components
% analysis.
%
% A chromaticity map describes the blending between two different
% reflectances in the output images.
%
% #### Illumination maps
% A directory containing either monochromatic or RGB images. Monochromatic
% images will be used directly as illumination maps, whereas the CIE 1976
% L*a*b* luminance channels of RGB images will be used as illumination
% maps.
%
% Chromaticity maps and illumination maps will be paired based on filename.
%
% ### Model of dispersion
%
% There are several types of data that can be provided:
%
% #### Polynomial model
% A '.mat' file containing several variables, which is the output of
% 'RAWDiskDispersion.m', for example. The following variables are required:
% - 'polyfun_data': A polynomial model of chromatic aberration, modeling
%   the warping from the reference wavelength band to the other wavelength
%   bands. `polyfun_data` can be converted to a function form using
%   `polyfun = makePolyfun(polyfun_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `false`.
% The following variables are sometimes required:
% - 'bands': A vector containing the wavelengths to use as the `lambda`
%   input argument of 'polyfunToMatrix()' (to evaluate a dispersion model),
%   and to form hyperspectral images. This variable is required only if not
%   provided in the colour space conversion data file, or directly in this
%   script (see below).
%
% The following variables are optional. If they are present, they are
% assumed to define a conversion between a geometrical optics coordinate
% system in which the polynomial model of chromatic aberration was
% constructed, and the image coordinate system:
% - 'image_params': A structure with an 'image_sampling' field, which is a
%   two-element vector containing the pixel height and width of the image.
% - 'pixel_size': A scalar containing the side length of a pixel.
%
% #### Function model
% A function for evaluating a model of dispersion in terms of three
% variables, X, Y, and lambda (wavelength). The function must have the same
% usage as the `polyfun` output argument of 'makePolyfun()'.
%
% ### Colour space conversion data
% A '.mat' file containing several variables, which is the output of
% 'SonyColorMap.m', for example. The following variables are required:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel of the output RGB images to the j-th
%   spectral band of the synthetic hyperspectral images.
% The following variables are sometimes required:
% - 'bands': A vector containing the wavelengths to use as the `lambda`
%   input argument of 'polyfunToMatrix()' (to evaluate a dispersion model),
%   and to form hyperspectral images. This variable takes precedence over
%   the same variable provided with the dispersion model (see above), but
%   can be overridden by a version of the variable provided directly in
%   this script (see below).
%
% ### Discrete spectral space
%
% The 'bands' variable defined in the parameters section of the code below
% can be empty (`[]`), in which case it is loaded from the dispersion model
% data or from the colour space conversion data (see above).
%
% The final value of 'bands' is determined according to the following list,
% in order by decreasing priority:
% - 'bands' defined in this script
% - 'bands' loaded from with colour space conversion data
% - 'bands' loaded from a polynomial model of dispersion
%
% The final 'bands' vector is clipped to the intervals defined by the other
% vectors of wavelengths, to avoid extrapolation when resampling data to
% conform to the final value of 'bands'. Resampling is needed to express
% all spectral quantities in the same discrete space of wavelengths.
%
% If no variable 'bands' was found in the colour space conversion data, but
% the final value of 'bands' has a length equal to the size of the second
% dimension of 'sensor_map', 'bands' is assumed to be compatible with
% 'sensor_map'. Otherwise, 'sensor_map' is resampled along its second
% dimension, by assuming that its first and last columns correspond to the
% first and last elements of 'bands', and that its remaining columns
% correspond to equally-spaced wavelengths in-between the first and last
% elements of 'bands'.
%
% ### Illuminant
% The following two sources of data can be provided:
%
% #### CIE D-illuminant data file
% A '.csv' file containing the following columns (refer to the webpage by
% Bruce Lindbloom cited below):
% - Wavelength, in nanometres
% - The 'S_0' function value at the corresponding wavelength
% - The 'S_1' function value at the corresponding wavelength
% - The 'S_2' function value at the corresponding wavelength
%
% The Kelvin temperature of the illuminant must be provided in the script
% parameters section below.
%
% The data will be resampled so that it is compatible with the 'bands'
% variable described above.
%
% #### Arbitrary illuminant
% The name of a function of a vector of wavelengths that returns a vector
% of relative spectral intensities.
%
% ### Spectral reflectances
%
% Each output image will be a blend of two reflectance spectra. If the
% input chromaticity map was a monochromatic image, the reflectance spectra
% will be selected randomly from the 24 colour patches of the 30-chart
% averaged ColorChecker data provided by BabelColor (cited below). If the
% input chromaticity map was an RGB image, the reflectance spectra from
% BabelColor with the closest RGB values to the RGB clusters of the image
% will be selected.
%
% ## Output
%
% ### Simulated images
%
% One of the following types of images is created for each chromaticity and
% illumination map pair. The filename of the input maps is represented by
% '*' below.
% - '*_hyperspectral.mat': A hyperspectral image (stored in the variable
%   'I_hyper') produced by blending two reflectance spectra according to
%   the per-pixel weights in the soft segmentation of the chromaticity map,
%   and then by illuminating the spectra according to the illumination
%   weights in the illumination map.
% - '*_rgb.tif': A colour image created by converting the hyperspectral
%   image to the RGB colour space of the camera.
% - '*_hyperspectral_warped.mat': A warped version of the hyperspectral
%   image (stored in the variable 'I_hyper_warped') created by applying the
%   dispersion model to the image.
% - '*_rgb_warped.tif': A colour image created by converting the
%   warped hyperspectral image to the RGB colour space of the camera.
% - '*_raw_warped.tif': A colour-filter array image produced by mosaicing
%   the warped RGB image according to the colour-filter pattern of the
%   camera.
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'bands': The final value of the 'bands' variable, determined as
%   discussed under the section "Discrete spectral space" above.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference. 'bands_color' is empty if no
%   such variable was found in the data file, or if the value of the loaded
%   variable was empty.
% - 'bands_polyfun': The 'bands' variable loaded from the dispersion model
%   data file, for reference. 'bands_polyfun' is empty if no such variable
%   was found in the data file, or if the value of the loaded variable was
%   empty.
% - 'chromaticity_filenames': A cell vector containing the input
%   chromaticity map filenames retrieved based on the wildcard provided in
%   the parameters section of the script.
% - 'illumination_filenames': A cell vector containing the input
%   illumination map filenames retrieved based on the wildcard provided in
%   the parameters section of the script.
% - 'sensor_map_resampled': The resampled version of the 'sensor_map'
%   variable, determined as discussed under the section "Discrete spectral
%   space" above.
% 
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)
%
% ## References
% - Foster, D. H. (2018). Tutorial on Transforming Hyperspectral Images to
%   RGB Colour Images. Retrieved from
%   http://personalpages.manchester.ac.uk/staff/d.h.foster/Tutorial_HSI2RGB/Tutorial_HSI2RGB.html
%   on June 5, 2018.
% - Lindbloom, Bruce J. (2017). Spectral Power Distribution of a CIE
%   D-Illuminant. Retrieved from http://www.brucelindbloom.com on June 4,
%   2018.
% - Pascale, Danny (2016). The ColorChecker Pages. Retrieved from
%   http://www.babelcolor.com/colorchecker.htm on June 4, 2018.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 6, 2018

% List of parameters to save with results
parameters_list = {
        'data_dispersion_model',...
        'polynomial_model_filename',...
        'dispersonFun_name',...
        'color_map_filename',...
        'bands_interp_method',...
        'use_cie_illuminant',...
        'illuminant_filename',...
        'illuminant_temperature',...
        'illuminant_name',...
        'illuminant_function_name',...
        'xyzbar_filename',...
        'lambda_xyzbar',...
        'reflectances_filename',...
        'n_colors',...
        'bayer_pattern',...
        'output_directory'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the chromaticity maps.
input_chromaticity_maps_wildcard = '';

% Wildcard for 'ls()' to find the illumination maps.
input_illumination_maps_wildcard = '';

data_dispersion_model = true;
if data_dispersion_model
    % Polynomial model of dispersion
    polynomial_model_filename = '';
    dispersonFun_name = [];
else
    % Function model of dispersion
    dispersonFun_name = 'sin';
    polynomial_model_filename = [];
end

% Colour space conversion data
color_map_filename = '';

% Override the wavelengths at which to evaluate the model of dispersion, if
% desired.
bands = [];
% Interpolation method used when resampling spectral data
bands_interp_method = 'linear';

use_cie_illuminant = true;
if use_cie_illuminant
    % CIE D-illuminant
    illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
    illuminant_temperature = 5003; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
    illuminant_name = 'd50';
    illuminant_function_name = [];
else
    % Arbitrary illuminant function
    illuminant_function_name = 'sqrt';
    illuminant_filename = [];
    illuminant_temperature = [];
    illuminant_name = [];
end

% CIE tristimulus functions
xyzbar_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180605_HyperspectralToSRGB_DHFoster/Tutorial_HSI2RGB/xyzbar.mat';
% Wavelengths at which the tristimulus functions were sampled
lambda_xyzbar = (400:10:720).';

% ColorChecker spectral reflectances
reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180604_ColorCheckerSpectralData_BabelColor/ColorChecker_spectra_reformatted_llanos.csv';

% Number of colours to blend in each output image
n_colors = 2;

% Colour-filter pattern
bayer_pattern = 'gbrg';

% Output directory for all images and saved data
output_directory = '';

% ## Debugging Flags
segmentColorsVerbose = true;

%% Load ColorChecker spectral reflectances and calculate their RGB values

% ColorChecker spectral reflectances
colorChecker_table = readtable(reflectances_filename);
variable_names = colorChecker_table.Properties.VariableNames;
lambda_colorChecker = colorChecker_table.(variable_names{1});
reflectances = colorChecker_table{:, 2:end};
n_patches = length(variable_names) - 1;

% Illuminant spectral power distribution
if use_cie_illuminant
    illuminant_data = csvread(illuminant_filename);
    lambda_illuminant = illuminant_data(:, 1);
    S_illuminant = illuminant_data(:, 2:end);
    spd_illuminant = ciedIlluminant(...
        illuminant_temperature, lambda_illuminant, S_illuminant, lambda_illuminant...
    );
else
    lambda_illuminant = lambda_colorChecker;
    spd_illuminant = feval(illuminant_function_name, lambda_illuminant);
end

% Find patch colours
variables_required = { 'xyzbar' };
load(xyzbar_filename, variables_required{:});
if ~all(ismember(variables_required, who))
    error('One or more of the CIE tristimulus functions variables is not loaded.')
end

colorChecker_rgb = reflectanceToColor(...
    lambda_illuminant, spd_illuminant,...
    lambda_colorChecker, reflectances,...
    lambda_xyzbar, xyzbar,...
    illuminant_name...
    );

%% Load dispersion model

bands_script = bands;
bands = [];
optional_variable = 'bands';

if data_dispersion_model
    model_variables_required = { 'polyfun_data', 'model_from_reference' };
    model_variables_transform = { 'image_params', 'pixel_size' };
    load(...
        polynomial_model_filename,...
        model_variables_required{:}, model_variables_transform{:},...
        optional_variable...
        );
    if ~all(ismember(model_variables_required, who))
        error('One or more of the dispersion model variables is not loaded.')
    end
    if model_from_reference
        error('Dispersion model is in the wrong frame of reference.')
    end

    if all(ismember(model_variables_transform, who))
        T = pixelsToWorldTransform(image_params.image_sampling, pixel_size);
        dispersonFun = makePolyfun(polyfun_data, T);
    else
        dispersonFun = makePolyfun(polyfun_data);
    end
else
    dispersonFun = str2func(dispersonFun_name);
end

bands_polyfun = bands;
bands = [];

%% Load colour space conversion data

model_variables_required = { 'sensor_map' };
load(color_map_filename, model_variables_required{:}, optional_variable);
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end

bands_color = bands;

%% Select a value of 'bands' and resample spectral arrays

% Select the highest-priority value of `bands`
if ~isempty(bands_script)
    bands = bands_script;
elseif ~isempty(bands_color)
    bands = bands_color;
elseif ~isempty(bands_polyfun)
    bands = bands_polyfun;
else
    error('The variable `bands` is not defined, or is empty');
end

% Use the intersection of all values of `bands` corresponding to arrays
if ~isempty(bands_color)
    bands = bands(bands >= min(bands_color) & bands <= max(bands_color));
end
bands = bands(bands >= min(lambda_colorChecker) & bands <= max(lambda_colorChecker));

% Compare with colour space conversion data
n_bands = length(bands);
n_bands_sensor_map = size(sensor_map, 2);
resample_bands = false;
if ~isempty(bands_color)
    if n_bands ~= length(bands_color) || any(bands ~= bands_color)
        % Resampling is needed
        resample_bands = true;
        bands_for_interp = bands_color;
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
    sensor_map_resampled = sensor_map_resampled.';
else
    sensor_map_resampled = sensor_map;
end

% Resample ColorChecker reflectances
reflectances_resampled = resampleArrays(...
    lambda_colorChecker, reflectances, bands,...
    bands_interp_method...
    );

% Resample the illuminant
if use_cie_illuminant
    spd_illuminant_resampled = ciedIlluminant(...
        illuminant_temperature, lambda_illuminant, S_illuminant, bands...
    );
else
    spd_illuminant_resampled = feval(illuminant_function_name, bands);
end

% Resample the CIE 'y-bar' function
xyzbar_resampled = resampleArrays(...
    lambda_xyzbar, xyzbar, bands,...
    bands_interp_method...
    );

n_bands = length(bands);

%% Compute the normalization constant for radiance

bands_diff = diff(bands);
bands_diff = [
    bands_diff;
    bands_diff(end)
    ];
radiance_normalization_constant = sum(spd_illuminant_resampled .* xyzbar_resampled(:, 2) .* bands_diff);

%% Find the chromaticity and illumination maps

chromaticity_filenames = listFiles(input_chromaticity_maps_wildcard);
n_images = length(chromaticity_filenames);

illumination_filenames = listFiles(input_illumination_maps_wildcard);
if n_images ~= length(chromaticity_filenames)
    error('There are different numbers of chromaticity maps and illumination maps.');
end

% Associate by filename
[~, chromaticity_names] = cellfun(fileparts, chromaticity_filenames, 'UniformOutput', false);
[~, illumination_names] = cellfun(fileparts, illumination_filenames, 'UniformOutput', false);
names_match = cellfun(strcmp, chromaticity_names, illumination_names, 'UniformOutput', true);
if ~all(names_match)
    error('Not all chromaticity map filenames and illumination map filenames match.');
end

%% Process the images

n_channels_rgb = 3;
ext = '.tif';

for i = 1:n_images
    % Cluster the colours in the chromaticity map
    C_map = imread(chromaticity_filenames{i});
    image_height = size(C_map, 1);
    image_width = size(C_map, 2);
    image_sampling = [image_height, image_width];
    n_channels = size(C_map, 3);
    n_px = image_height * image_width;
    [C_centers, ~, C_softSegmentation] = segmentColors(C_map, n_colors, segmentColorsVerbose);    
    
    if n_channels == n_channels_rgb
        % Select matching colours to blend
        color_indices = zeros(n_colors, 1);
        for k = 1:color_indices
            distances = repmat(C_centers(k, :), n_patches, 1) - colorChecker_rgb;
            distances = dot(distances, distances, 2);
            [ ~, color_indices(k) ] = min(distances);
            % Avoid repeating colours
            while sum(color_indices(k) == color_indices(1:(k - 1))) ~= 0
                distances(color_indices(k)) = Inf;
                [ ~, color_indices(k) ] = min(distances);
            end
        end
    elseif n_channels == 1
        % Select random colours to blend
        color_indices = randperm(n_patches, n_colors).';
    end
    
    % Blend spectral reflectances together: Sum weighted reflectances
    reflectance_per_pixel = sum(repmat(...
        permute(reflectances_resampled(:, color_indices), [3, 1, 2]),...
        n_px, 1, 1 ...
    ) * repmat(...
        permute(reshape(C_softSegmentation, n_px, n_colors, 1), [1, 3, 2]),...
        1, n_bands, 1 ...
    ), 3);

    % Compute radiances
    Rad = reflectance_per_pixel .* repmat(...
        spd_illuminant_resampled.', n_px, 1 ...
        ) ./ radiance_normalization_constant;
    
    % Load the illumination map
    I_map = imread(illumination_filenames{i});
    % Check that image size is compatible
    if (size(I_map, 1) ~= image_height) || (size(I_map, 2) ~= image_width)
        error('Illumination map dimensions do not match chromaticity map dimensions, filename "%s".', illumination_filenames{i});
    end
    I_map = im2double(I_map);
    if size(I_map, 3) == n_channels_rgb
        I_map = rgb2lab(I_map);
        I_map = I_map(:, :, 1);
        I_map_min = min(I_map(:));
        I_map = (I_map - I_map_min) ./ (max(I_map(:)) - I_map_min);
    elseif size(I_map, 3) ~= 1
        error('Illumination map does not have one or three colour channels, filename "%s".', illumination_filenames{i});
    end
    
    % Multiply by the illumination map and reshape into a hyperspectral image
    H = reshape(Rad .* repmat(reshape(I_map, n_px, 1), 1, n_bands), [], 1);
    I_hyper = reshape(H, image_height, image_width, n_bands);
    
    % Compute the equivalent RGB image
    Omega = channelConversionMatrix(image_sampling, sensor_map_resampled);
    rgb = Omega * H;
    rgb_3D = reshape(rgb, image_height, image_width, n_channels_rgb);

    % Simulate dispersion
    W = polyfunToMatrix(...
            dispersonFun, bands,...
            image_sampling, image_sampling,...
            [0, 0, image_width,  image_height], true...
        );
    H_warped = W * H;
    I_hyper_warped = reshape(H_warped, image_height, image_width, n_channels_rgb);
    
    % Compute the equivalent RGB image
    rgb_warped = Omega * H_warped;
    rgb_warped_3D = reshape(rgb_warped, image_height, image_width, n_bands);
    
    % Compute the RAW image
    M = mosaicMatrix(image_sampling, bayer_pattern);
    raw = M * rgb_warped;
    raw_2D = reshape(raw, image_sampling);
            
    % Save the results
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_hyperspectral.mat']);
    save(output_filename, 'I_hyper');
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_rgb' ext]);
    imwrite(rgb_3D, output_filename);
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_hyperspectral_warped.mat']);
    save(output_filename, 'I_hyper_warped');
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_rgb_warped' ext]);
    imwrite(rgb_warped_3D, output_filename);
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_raw_warped' ext]);
    imwrite(raw_2D, output_filename);
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'chromaticity_filenames',...
        'illumination_filenames',...
        'bands_polyfun',...
        'bands_color',...
        'bands',...
        'sensor_map_resampled'...
    } ];
save_data_filename = fullfile(output_directory, 'BilateralImages.mat');
save(save_data_filename, save_variables_list{:});