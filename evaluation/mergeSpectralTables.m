function t = mergeSpectralTables(ts)
% MERGESPECTRALTABLES  Combine spectral comparison results across images
%
% ## Syntax
% t = mergeSpectralTables(ts)
%
% ## Description
% t = mergeSpectralTables(ts)
%   Returns a table aggregating the data in the input tables.
%
% ## Input Arguments
%
% ts -- Spectral comparison results
%   A cell vector containing tables of the form of the `e_spectral_table`
%   output argument of 'evaluateAndSaveSpectral()'. Each element contains
%   evaluation results for one image in a dataset.
%
% ## Output Arguments
%
% t -- Aggregated spectral error statistics
%   A table summarizing the data in `ts` across images. Variables
%   representing maximum, minimum, median, and mean values will be
%   summarized by taking the maximum, minimum, median, and mean,
%   respectively.
%
% See also evaluateAndSaveSpectral, mergeRGBTables, writetable

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 20, 2018

narginchk(1, 1);
nargoutchk(1, 1);

% Strip image-specific variables
n_tables = length(ts);
variables = ts{1}.Properties.VariableNames;
for i = 2:n_tables
    variables = intersect(variables, ts{i}.Properties.VariableNames, 'stable');
end

% Merge tables
t_all = ts{1}(:, variables);
t_all.Image = ones(size(t_all, 1), 1);
for i = 2:n_tables
    t_i = ts{i}(:, variables);
    t_i.Image = i * ones(size(t_i, 1), 1);
    t_all = union(t_all, t_i);
end

% Compute summary statistics
max_variables = {'MSE_max', 'RMSE_max'};
t_max = varfun(...
    @max, t_all, 'GroupingVariables', 'Algorithm', 'InputVariables',...
    max_variables...
);
min_variables = {'PSNR_min', 'SSIM_min', 'GOF_min'};
t_min = varfun(...
    @min, t_all, 'GroupingVariables', 'Algorithm', 'InputVariables',...
    min_variables...
);
median_variables = {'MSE_median', 'PSNR_median', 'SSIM_median', 'RMSE_median', 'GOF_median'};
t_median = varfun(...
    @median, t_all, 'GroupingVariables', 'Algorithm', 'InputVariables',...
    median_variables...
);
mean_variables = setdiff(variables, max_variables);
mean_variables = setdiff(mean_variables, min_variables);
mean_variables = setdiff(mean_variables, median_variables);
mean_variables = setdiff(mean_variables, {'Image', 'Algorithm'});
t_mean = varfun(...
    @mean, t_all, 'GroupingVariables', 'Algorithm',...
    'InputVariables', mean_variables...
);

% Combine summary statistics
t = t_mean;
t(:, t_max.Properties.VariableNames) = t_max(:, :);
t(:, t_min.Properties.VariableNames) = t_min(:, :);
t(:, t_median.Properties.VariableNames) = t_median(:, :);

% Restore original variable names
variables_new = t.Properties.VariableNames;
search_strings = {'max_', 'min_', 'mean_', 'median_'};
for i = 1:length(variables_new)
    for j = 1:length(search_strings)
        if strcmp(variables_new{i}(1:length(search_strings{j})), search_strings{j})
            t.Properties.VariableNames{variables_new{i}} =...
                variables_new{i}((length(search_strings{j}) + 1):end);
        end
    end
end

% Use the original ordering of variables
t = t(:, variables);

end
