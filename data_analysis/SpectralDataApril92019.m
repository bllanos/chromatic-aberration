%% Display April 2019 spectral data
% Plot the spectral data collected on April 7 and 9, 2019.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 9, 2019

%% Input data and parameters

% Sample spectra
radiance_filenames = {
    '${FILEPATH}';
    '${FILEPATH}';
    '${FILEPATH}';
    '${FILEPATH}'
};

% Names for plots
set_names = {
    'ColorChecker Classic (Afternooon)';
    'ColorChecker Classic (Morning)';
    'Points in "book" scene';
    'Projector light with and without filters'
};

% Wavelength limits for plots
x_limits = {
    [];
    [];
    [];
    []
};

% Reference set index and sample index in that set to use as the divisor for
% normalization
references = {
    [1, 19];
    [2, 19];
    [3, 12];
    [4, 1];
};

%% Load and display sample spectra

n_sets = length(radiance_filenames);
radiances = cell(n_sets, 1);
for s = 1:n_sets
    sample_table = readtable(radiance_filenames{s});
    radiances{s} = sample_table{:, 2:end};
end

for s = 1:n_sets
    sample_table = readtable(radiance_filenames{s});
    variable_names = sample_table.Properties.VariableNames{2:end};
    lambda_samples = sample_table.(sample_table.Properties.VariableNames{1});
    radiances_s = sample_table{:, 2:end};
    n_radiances = size(radiances_s, 2);    

    % Visualization
    set_name = set_names{s};

    % Absolute radiances
    figure;
    hold on
    names_legend = cell(n_radiances, 1);
    names_legend_escaped = cell(n_radiances, 1);
    for i = 1:n_radiances
        plot(...
            lambda_samples, radiances_s(:, i),...
            'LineWidth', 2, 'Marker', 'none'...
        );
        % Recover original variable names, which contained spaces
        if ~isempty(sample_table.Properties.VariableDescriptions) && ~isempty(sample_table.Properties.VariableDescriptions{i  + 1})
            names_legend{i} = strsplit(sample_table.Properties.VariableDescriptions{i  + 1}, ':');
            names_legend{i} = names_legend{i}{end};
        else
            names_legend{i} = sample_table.Properties.VariableNames{i + 1};
        end
        names_legend_escaped{i} = strrep(names_legend{i}, '_', '\_');
    end
    hold off
    title(sprintf('%s absolute spectral signals', set_name))
    if ~isempty(x_limits{s}) 
        xlim(x_limits{s});
    end
    xlabel('\lambda [nm]')
    ylabel('Spectral signal')
    legend(names_legend_escaped);
    ax = gca;
    ax.Color = [0.5 0.5 0.5];
    
    % Relative radiances
    reference_radiance = radiances{references{s}(1)}(:, references{s}(2));
    figure;
    hold on
    for i = 1:n_radiances
        plot(...
            lambda_samples, radiances_s(:, i) ./ reference_radiance,...
            'LineWidth', 2, 'Marker', 'none'...
        );
    end
    hold off
    title(sprintf('%s relative spectral signals', set_name))
    if ~isempty(x_limits{s}) 
        xlim(x_limits{s});
    end
    ylim([0, 2]);
    xlabel('\lambda [nm]')
    ylabel('Relative spectral signal')
    legend(names_legend_escaped);
    ax = gca;
    ax.Color = [0.5 0.5 0.5];
end