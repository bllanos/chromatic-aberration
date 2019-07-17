%% Display June 26, 2018 averaged spectral data
% Plot spectral reflectances and display sRGB values for the spectral data
% collected on June 26, 2018
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### CIE D-illuminant data file
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
% ### Spectral reflectances
% A '.csv' file containing a header row, a first column for wavelength
% values, and remaining columns for relative spectral reflectances of each
% sample.
%
% ### CIE tristimulus functions
% A '.mat' file containing a variable 'xyzbar', which can be used as the
% `C` input argument of 'cieSpectralToColor()'.
%
% ## Output
%
% ### Graphical output
% - Plots of the spectral reflectances of the samples, with plotlines
%   coloured by their sRGB colours.
%
% ### Console output
% - sRGB values of the samples
%
% ## References
% - Foster, D. H. (2018). Tutorial on Transforming Hyperspectral Images to
%   RGB Colour Images. Retrieved from
%   http://personalpages.manchester.ac.uk/staff/d.h.foster/Tutorial_HSI2RGB/Tutorial_HSI2RGB.html
%   on June 5, 2018.
% - Lindbloom, Bruce J. (2017). Computing XYZ From Spectral Data. Retrieved
%   from http://www.brucelindbloom.com on June 11, 2018.
% - Lindbloom, Bruce J. (2017). Spectral Power Distribution of a CIE
%   D-Illuminant. Retrieved from http://www.brucelindbloom.com on June 4,
%   2018.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 26, 2018

%% Input data and parameters

% CIE D-illuminant
illuminant_filename = '${FILEPATH}';
illuminant_temperature = 5003; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D

% CIE tristimulus functions
xyzbar_filename = '${FILEPATH}';

% Sample spectral reflectances
reflectances_filename = fullfile('.', 'demo_data', 'spectral_data', 'spectra_averaged.csv');

% Categorization of the samples
indices.bandpass_filters = 2:8;
indices.large_chequerboard = 9:10;
indices.paper = 11:12;
indices.patches = 13:36;
indices.small_chequerboard = 37:38;

% Names for plots
set_names.bandpass_filters = 'Bandpass filters';
set_names.large_chequerboard = 'Large chequerboard';
set_names.paper = 'Laser printer paper';
set_names.patches = 'ColorChecker Classic';
set_names.small_chequerboard = 'Small chequerboard';

% Wavelength limits for plots
x_limits.bandpass_filters = [350 750];
x_limits.large_chequerboard = [];
x_limits.paper = [];
x_limits.patches = [];
x_limits.small_chequerboard = [];

%% Load and display the illuminant spectral power distribution

illuminant_data = csvread(illuminant_filename);
lambda_illuminant = illuminant_data(:, 1);
S_illuminant = illuminant_data(:, 2:end);
spd_illuminant = ciedIlluminant(...
    illuminant_temperature, lambda_illuminant, S_illuminant, lambda_illuminant...
);

figure;
plot(lambda_illuminant, spd_illuminant);
title(sprintf(...
    'CIE D-illuminant spectral power distribution for CCT = %g Kelvin',...
    illuminant_temperature...
    ));
xlabel('\lambda [nm]');
ylabel('Relative power');

%% Load and display sample spectra

sample_table = readtable(reflectances_filename);
variable_names = sample_table.Properties.VariableNames;
lambda_samples = sample_table.(variable_names{1});
reflectances = sample_table{:, :};

% Find sample colours
xyzbar_table = readtable(xyzbar_filename);
lambda_xyzbar = xyzbar_table{:, 1};
xyzbar = xyzbar_table{:, 2:end};

rgb = reflectanceToColor(...
    lambda_illuminant, spd_illuminant,...
    lambda_samples, reflectances,...
    lambda_xyzbar, xyzbar...
    );
rgb_integer = floor(256 * rgb);

% Visualization
sets = fieldnames(indices);
for s = 1:length(sets)
    current_indices = indices.(sets{s});
    set_name = set_names.(sets{s});

    figure;
    hold on
    names_legend = cell(length(current_indices), 1);
    fprintf('%s sRGB colours under a CCT = %g Kelvin illuminant:\n', set_name, illuminant_temperature);
    for i = 1:length(current_indices)
        j = current_indices(i);
        plot(...
            lambda_samples, reflectances(:, j),...
            'Color', rgb(j, :), 'LineWidth', 2, 'Marker', 'none'...
        );
        % Recover original variable names, which contained spaces
        names_legend{i} = strsplit(sample_table.Properties.VariableDescriptions{j}, ':');
        names_legend{i} = names_legend{i}{end};
        if isempty(names_legend{i})
            names_legend{i} = sample_table.Properties.VariableNames{j};
        end

        fprintf('\t%s: %d, %d, %d\n', names_legend{i}, rgb_integer(j, 1), rgb_integer(j, 2), rgb_integer(j, 3));
    end
    hold off
    title(sprintf('%s spectral signals', set_name))
    if ~isempty(x_limits.(sets{s})) 
        xlim(x_limits.(sets{s}));
    end
    ylim([0, 1.3]);
    xlabel('\lambda [nm]')
    ylabel('Relative spectral signal')
    legend(names_legend);
    ax = gca;
    ax.Color = [0.5 0.5 0.5];
end