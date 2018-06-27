%% Display June 26, 2018 spectral data
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
% ### Reference spectral reflectances
% A '.csv' file containing a header row, a first column for wavelength
% values, and remaining columns for relative spectral reflectances of each
% patch in the ColorChecker chart, as collected by Danny Pascale (cited
% below).
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
% - Pascale, Danny (2016). The ColorChecker Pages. Retrieved from
%   http://www.babelcolor.com/colorchecker.htm on June 4, 2018.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 26, 2018

%% Input data and parameters

% CIE D-illuminant
illuminant_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180604_Spectral power distributions_BruceLindbloom/DIlluminants.csv';
illuminant_temperature = 5003; % From https://en.wikipedia.org/wiki/Standard_illuminant#Illuminant_series_D
illuminant_name = 'd50';

% CIE tristimulus functions
xyzbar_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180614_ASTM_E308/Table1_CIE1931_2DegStandardObserver.csv';

% Sample spectral reflectances
reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180626_SpectralCharacterizationOfSetup/spectra_named.csv';

% Comparison data, collected prior to November 2014, by Danny Pascale
% (cited above)
reference_reflectances_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180604_ColorCheckerSpectralData_BabelColor/ColorChecker_spectra_reformatted_llanos.csv';

% Categorization of the samples
indices.bandpass_filters = 2:9;
indices.large_chequerboard = 10:13;
indices.paper = 14:15;
indices.patches = 16:39;
indices.projector_instability = 40:41;
indices.bandpass_filters_trial = 42:43;
indices.small_chequerboard = 44:47;

% Names for plots
set_names.bandpass_filters = 'Bandpass filters';
set_names.large_chequerboard = 'Large chequerboard';
set_names.paper = 'Laser printer paper';
set_names.patches = 'ColorChecker Classic';
set_names.projector_instability = 'Drift of projector light';
set_names.bandpass_filters_trial = 'Bandpass filters (initial trial)';
set_names.small_chequerboard = 'Small chequerboard';

% Wavelength limits for plots
x_limits.bandpass_filters = [350 750];
x_limits.large_chequerboard = [];
x_limits.paper = [];
x_limits.patches = [];
x_limits.projector_instability = [];
x_limits.bandpass_filters_trial = [350 750];
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
    lambda_xyzbar, xyzbar,...
    illuminant_name...
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
    fprintf('%s sRGB colours under a %s illuminant:\n', set_name, illuminant_name);
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

%% Compare with Danny Pascale's data

reference_table = readtable(reference_reflectances_filename);
variable_names = reference_table.Properties.VariableNames;
lambda_references = reference_table.(variable_names{1});
reference_reflectances = reference_table{:, :};

% Find reference colours
reference_rgb = reflectanceToColor(...
    lambda_illuminant, spd_illuminant,...
    lambda_references, reference_reflectances,...
    lambda_xyzbar, xyzbar,...
    illuminant_name...
    );

current_indices = indices.patches;
set_name = set_names.patches;

figure;
hold on
names_legend = cell(length(current_indices), 1);
for i = 1:length(current_indices)
    j = current_indices(i);
    plot(...
        lambda_samples, reflectances(:, j),...
        'Color', rgb(j, :), 'LineWidth', 2, 'Marker', 'none'...
    );
    plot(...
        lambda_references, reference_reflectances(:, i+1),...
        'Color', reference_rgb(i+1, :),...
        'LineWidth', 2, 'LineStyle', ':', 'Marker', 'none'...
    );
    % Recover original variable names, which contained spaces
    names_legend{i} = strsplit(sample_table.Properties.VariableDescriptions{j}, ':');
    names_legend{i} = names_legend{i}{end};
    if isempty(names_legend{i})
        names_legend{i} = sample_table.Properties.VariableNames{j};
    end
end
hold off
title(sprintf('%s spectral signals compared with Danny Pascale''s data (dotted lines)', set_name))
xlim([min(lambda_references), max(lambda_references)]);
ylim([0, 1.3]);
xlabel('\lambda [nm]')
ylabel('Relative spectral reflectance')
legend(repelem(names_legend, 2));
ax = gca;
ax.Color = [0.5 0.5 0.5];