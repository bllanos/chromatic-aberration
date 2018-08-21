function t = mergeRGBTables(ts)
% MERGERGBTABLES  Combine RGB comparison results across images
%
% ## Syntax
% t = mergeRGBTables(ts)
%
% ## Description
% t = mergeRGBTables(ts)
%   Returns a table aggregating the data in the input tables.
%
% ## Input Arguments
%
% ts -- RGB comparison results
%   A cell vector containing tables of the form of the `e_rgb_table`
%   output argument of 'evaluateAndSaveRGB()'. Each element contains
%   evaluation results for one image in a dataset.
%
% ## Output Arguments
%
% t -- Aggregated RGB error statistics
%   A table summarizing the data in `ts` across images. Variables
%   representing maximum, minimum, median, and mean values will be
%   summarized by taking the maximum, minimum, median, and mean,
%   respectively.
%
% See also evaluateAndSaveRGB, mergeSpectralTables, writetable

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 21, 2018

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
mean_variables = setdiff(variables, {'Image', 'Algorithm'});
t_mean = varfun(...
    @mean, t_all, 'GroupingVariables', 'Algorithm',...
    'InputVariables', mean_variables...
);

% Combine summary statistics
t = t_mean;

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
