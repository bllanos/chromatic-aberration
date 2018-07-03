%% Aberrated image generation
% Create ideal, and warped hyperspectral images and their 3-channel
% equivalents, from distributions of pairs of object spectral reflectances.
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
% - 'polyfun_data': A polynomial model of dispersion, modeling the warping
%   from the reference wavelength band to the other wavelength bands.
%   `polyfun_data` can be converted to a function form using `polyfun =
%   makePolyfun(polyfun_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `false`.
% - 'model_space': A structure with same form as the `model_space` input
%   argument of 'modelSpaceTransform()', required to stretch the valid
%   domain of the dispersion model to match the output images.
%
% ### Colour space conversion data
% A '.mat' file containing several variables, which is the output of
% 'SonyColorMap.m', for example. The following variables are required:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel of the output sensor response images to the
%   j-th spectral band of the synthetic hyperspectral images.
% - 'bands': A vector containing the wavelengths to use as the `lambda`
%   input argument of 'polyfunToMatrix()' (to evaluate a dispersion model),
%   and to form hyperspectral images. This variable can be overridden by a
%   version of the variable provided directly in this script (see below).
%   However, it is still needed from this file to normalize spectral
%   radiances, even when overridden in its other functions.
% - 'channel_mode': A Boolean value indicating whether the input colour
%   space is a set of colour channels (true) or a set of spectral bands
%   (false). A value of `false` is required.
%
% ### Discrete spectral space
%
% The 'bands' variable defined in the parameters section of the code below
% can be empty (`[]`), in which case it is loaded from the colour space
% conversion data (see above).
%
% The final 'bands' vector is clipped to the intervals defined by the other
% vectors of wavelengths, to avoid extrapolation when resampling data to
% conform to the final value of 'bands'. Resampling is needed to express
% all spectral quantities in the same discrete space of wavelengths.
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
%   weights in the illumination map. Image values are normalized radiances.
% - '*_3.tif' and '*_3.mat': A colour image created by converting the
%   hyperspectral image to the raw colour space of the camera.
% - '*_hyperspectral_warped.mat': A warped version of the hyperspectral
%   image (stored in the variable 'I_hyper_warped') created by applying the
%   dispersion model to the image.
% - '*_3_warped.tif' and '*_3_warped.mat': A colour image created by
%   converting the warped hyperspectral image to the raw colour space of
%   the camera.
% - '*_raw_warped.tif' and '*_raw_warped.mat': A colour-filter array image
%   produced by mosaicing the warped sensor response image according to the
%   colour-filter pattern of the camera.
%
% The raw colour space of the camera is determined by the colour space
% conversion data provided as input to this script. A camera may apply
% additional operations to convert sensor responses from the raw colour
% space to, for example, sRGB colours.
%
% Note that both '.mat' and '.tif' files are output for monochromatic or
% three-channel images to provide both easy display ('.tif' files) and
% lossless storage ('.mat' files).
%
% ### Data file output
%
% A '.mat' file containing the following variables:
%
% - 'bands': The final value of the 'bands' variable, determined as
%   discussed under the section "Discrete spectral space" above.
% - 'bands_color': The 'bands' variable loaded from the colour space
%   conversion data file, for reference.
% - 'chromaticity_filenames': A cell vector containing the input
%   chromaticity map filenames retrieved based on the wildcard provided in
%   the parameters section of the script.
% - 'illumination_filenames': A cell vector containing the input
%   illumination map filenames retrieved based on the wildcard provided in
%   the parameters section of the script.
% - 'sensor_map_resampled': The resampled version of the 'sensor_map'
%   variable, determined as discussed under the section "Discrete spectral
%   space" above.
% - 'polyfun_data': A copy of the same variable from the input polynomial
%   model of dispersion
% - 'model_from_reference': A copy of the same variable from the input
%   polynomial model of dispersion
% - 'model_space': A copy of the same variable from the input polynomial
%   model of dispersion
% - 'fill': A Boolean value of `true` indicating that the valid domain of
%   the model of dispersion is to be stretched to match the domain of each
%   output image.
%
% The dispersion model variables are needed as input for scripts evaluating
% different methods for correcting chromatic aberration.
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
        'polynomial_model_filename',...
        'color_map_filename',...
        'normalization_channel',...
        'int_method',...
        'bands_interp_method',...
        'use_cie_illuminant',...
        'illuminant_filename',...
        'illuminant_temperature',...
        'illuminant_name',...
        'illuminant_function_name',...
        'xyzbar_filename',...
        'reflectances_filename',...
        'n_colors',...
        'bayer_pattern',...
        'output_directory'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the chromaticity maps.
input_chromaticity_maps_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180629_TestingBimaterialImages/*.jpg';

% Wildcard for 'ls()' to find the illumination maps.
input_illumination_maps_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180629_TestingBimaterialImages/*.jpg';

% Polynomial model of dispersion
polynomial_model_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180524_PolynomialDispersion_AnalyzePSF_accurate/DoubleConvexThickLensDispersionResults_modelFromReference_false_withModelSpace.mat';

% Colour space conversion data
color_map_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180629_TestingBimaterialImages/SonyColorMapData.mat';
% Colour channel to use for radiance normalization
normalization_channel = 2;
% Integration method to use for colour calculations
int_method = 'trap';

% Override the wavelengths at which to evaluate the model of dispersion, if
% desired.
bands = 430:10:650;
% Interpolation method used when resampling spectral data
bands_interp_method = 'linear';

use_cie_illuminant = true;
if use_cie_illuminant
    % CIE D-illuminant
    illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
    illuminant_temperature = 5003; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
    illuminant_name = 'd50';
    illuminant_function_name = []; % Unused
else
    % Arbitrary illuminant function
    illuminant_function_name = 'sqrt';
    illuminant_filename = []; % Unused
    illuminant_temperature = []; % Unused
    illuminant_name = []; % Unused
end

% CIE tristimulus functions
xyzbar_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180614_ASTM_E308/Table1_CIE1931_2DegStandardObserver.csv';

% ColorChecker spectral reflectances
reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180604_ColorCheckerSpectralData_BabelColor/ColorChecker_spectra_reformatted_llanos.csv';

% Number of colours to blend in each output image
n_colors = 2;

% Colour-filter pattern
bayer_pattern = 'gbrg';

% Output directory for all images and saved data
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180629_TestingBimaterialImages';

% ## Debugging Flags
segmentColorsVerbose = false;

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
xyzbar_table = readtable(xyzbar_filename);
lambda_xyzbar = xyzbar_table{:, 1};
xyzbar = xyzbar_table{:, 2:end};

colorChecker_rgb = reflectanceToColor(...
    lambda_illuminant, spd_illuminant,...
    lambda_colorChecker, reflectances,...
    lambda_xyzbar, xyzbar,...
    illuminant_name, int_method...
    );

%% Load dispersion model

model_variables_required = {...
    'polyfun_data', 'model_from_reference', 'model_space'...
};
load(polynomial_model_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the dispersion model variables is not loaded.')
end
if model_from_reference
    error('Dispersion model is in the wrong frame of reference.')
end

%% Load colour space conversion data

bands_script = bands;
bands = [];

optional_variable = 'bands';
model_variables_required = { 'sensor_map', 'channel_mode', optional_variable };
load(color_map_filename, model_variables_required{:});
if ~all(ismember(model_variables_required, who))
    error('One or more of the required colour space conversion variables is not loaded.')
end
if channel_mode
    error('The input space of the colour conversion data must be a spectral space, not a space of colour channels.')
end

bands_color = bands;

%% Select a value of 'bands' and resample spectral arrays

% Select the highest-priority value of `bands`
if ~isempty(bands_script)
    bands = bands_script;
else
    bands = bands_color;
end

% Use the intersection of all values of `bands` corresponding to arrays
bands = bands(bands >= min(bands_color) & bands <= max(bands_color));
bands = bands(bands >= min(lambda_colorChecker) & bands <= max(lambda_colorChecker));

% Resample colour space conversion data if necessary
if length(bands) ~= length(bands_color) || any(bands ~= bands_color)
    [sensor_map_resampled, bands] = resampleArrays(...
        bands_color, sensor_map.', bands,...
        bands_interp_method...
        );
    sensor_map_resampled = sensor_map_resampled.';
else
    sensor_map_resampled = sensor_map;
end
n_bands = length(bands);

%% Calculate spectral radiances

[lambda_Rad, ~, Rad_normalized] = reflectanceToRadiance(...
    lambda_illuminant, spd_illuminant,...
    lambda_colorChecker, reflectances,...
    bands_color, sensor_map.',...
    normalization_channel, int_method...
);

% Resample radiances
Rad_normalized_resampled = resampleArrays(...
    lambda_Rad, Rad_normalized, bands,...
    bands_interp_method...
    );

%% Find the chromaticity and illumination maps

chromaticity_filenames = listFiles(input_chromaticity_maps_wildcard);
n_images = length(chromaticity_filenames);

illumination_filenames = listFiles(input_illumination_maps_wildcard);
if n_images ~= length(chromaticity_filenames)
    error('There are different numbers of chromaticity maps and illumination maps.');
end

% Associate by filename
chromaticity_names = cell(n_images, 1);
illumination_names = cell(n_images, 1);
for i = 1:n_images
    [~, chromaticity_names{i}] = fileparts(chromaticity_filenames{i});
    [~, illumination_names{i}] = fileparts(illumination_filenames{i});
    names_match = strcmp(chromaticity_names{i}, illumination_names{i});
    if ~names_match
        error('Not all chromaticity map filenames and illumination map filenames match.');
    end
end

%% Process the images

n_channels_rgb = 3;
n_channels_raw = 3;
ext = '.tif';
fill = 'true';

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
        for k = 1:n_colors
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
    
    % Blend radiances together: Sum weighted radiances
    Rad_per_pixel = sum(repmat(...
        permute(Rad_normalized_resampled(:, color_indices), [3, 1, 2]),...
        n_px, 1, 1 ...
    ) .* repmat(...
        permute(reshape(C_softSegmentation, n_px, n_colors, 1), [1, 3, 2]),...
        1, n_bands, 1 ...
    ), 3);
    
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
    H = reshape(Rad_per_pixel .* repmat(reshape(I_map, n_px, 1), 1, n_bands), [], 1);
    I_hyper = reshape(H, image_height, image_width, n_bands);
    
    % Compute the equivalent sensor response image
    Omega = channelConversionMatrix(image_sampling, sensor_map_resampled, bands, int_method);
    raw_full = Omega * H;    
    raw_full_3D = reshape(raw_full, image_height, image_width, n_channels_raw);
    
    % Simulate dispersion
    [~, T_roi] = modelSpaceTransform(...
        [image_height, image_width], model_space, fill...
    );
    dispersonFun = makePolyfun(polyfun_data, T_roi);
    W = polyfunToMatrix(...
            dispersonFun, bands,...
            image_sampling, image_sampling,...
            [0, 0, image_width,  image_height], true...
        );
    H_warped = W * H;
    I_hyper_warped = reshape(H_warped, image_height, image_width, n_bands);
    
    % Compute the equivalent sensor response image
    raw_warped = Omega * H_warped;
    raw_warped_3D = reshape(raw_warped, image_height, image_width, n_channels_raw);
    
    % Compute the RAW image
    M = mosaicMatrix(image_sampling, bayer_pattern);
    raw = M * raw_warped;
    raw_2D = reshape(raw, image_sampling);
            
    % Save the results
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_hyperspectral.mat']);
    save(output_filename, 'I_hyper');
    
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_3.mat']);
    save(output_filename, 'raw_full_3D');
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_3' ext]);
    imwrite(raw_full_3D, output_filename);
    
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_hyperspectral_warped.mat']);
    save(output_filename, 'I_hyper_warped');
    
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_3_warped.mat']);
    save(output_filename, 'raw_warped_3D');
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_3_warped' ext]);
    imwrite(raw_warped_3D, output_filename);
    
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_raw_warped.mat']);
    save(output_filename, 'raw_2D');
    output_filename = fullfile(output_directory, [chromaticity_names{i} '_raw_warped' ext]);
    imwrite(raw_2D, output_filename);
end

%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
        'chromaticity_filenames',...
        'illumination_filenames',...
        'bands_color',...
        'bands',...
        'sensor_map_resampled',...
        'polyfun_data',...
        'model_from_reference',...
        'model_space',...
        'fill'...
    } ];
save_data_filename = fullfile(output_directory, 'BilateralImages.mat');
save(save_data_filename, save_variables_list{:});