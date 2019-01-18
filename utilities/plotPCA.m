function fg = plotPCA(x, coeff, mu, label_str, legend_str, title_str)
% PLOTPCA  Plot a point cloud and its PCA approximation
%
% ## Syntax
% plotPCA(x, coeff [, mu, label_str, legend_str, title_str])
% fg = plotPCA(x, coeff [, mu, label_str, legend_str, title_str])
%
% ## Description
% plotPCA(x, coeff [, mu, label_str, legend_str, title_str])
%   Creates a figure showing the original data, PCA vectors, and optionally PCA
%   lines (when `mu` is passed, and is not empty).
%
% fg = plotPCA(x, coeff [, mu, label_str, legend_str, title_str])
%   Returns the figure in which the data has been visualized.
%
% ## Input Arguments
%
% x -- Point cloud data
%   A matrix, where the rows index points, and the columns index variables. `x`
%   can have two or three columns.
%
% coeff -- PCA vectors
%   A matrix, with as many rows as there are columns in `x`, containing the PCA
%   coefficients. Each column of `coeff` is one PCA coefficient to be plotted.
%   The PCA coefficients will be plotted as dashed lines through the bounding
%   box of the point cloud that start at the origin (which will not be visible
%   depending upon the bounding box).
%
% mu -- Point cloud centroid
%   The mean of the dataset in `x`. When `mu` is passed, and is not empty
%   (`[]`), lines will be drawn through the point `mu`, parallel to the
%   directions defined by `coeff`.
%
% label_str -- Axes labels
%   A cell array of character vectors containing labels for the plot axes. There
%   should be as many elements as there are columns of `x`.
%
% legend_str -- Legend labels
%   A cell array of character vectors containing labels for the point cloud
%   data, and for the PCA coefficient vectors. There should be as many elements
%   as there are columns of `coeff`, plus one.
%
% title_str -- Plot title
%   A character vector containing the title to add to the plot.
%
% See also pca

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 20, 2018

nargoutchk(0, 1);
narginchk(2, 6);

n_dims = size(x, 2);
if all(n_dims ~= [2, 3])
    error('`x` must have two or three columns, not %d.', n_dims);
end

% Subsample the points
max_points = 2000;
if size(x, 1) > max_points
    x_ind = randperm(size(x, 1), max_points);
else
    x_ind = 1:size(x, 1);
end

% Plot the point cloud
fg = figure;
if n_dims == 2
    scatter(x(x_ind, 1), x(x_ind, 2), 'Marker', '.', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
elseif n_dims == 3
    scatter3(x(x_ind, 1), x(x_ind, 2), x(x_ind, 3), 'Marker', '.', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
else
    error('Unexpected value of `n_dims`, %d.', n_dims);
end

hold on

% Preprocessing
points_xlim = [min(x(:, 1)), max(x(:, 1))];
points_ylim = [min(x(:, 2)), max(x(:, 2))];
if n_dims == 3
    points_zlim = [min(x(:, 3)), max(x(:, 3))];
else
    points_zlim = [-Inf, Inf];
end

n_coeff = size(coeff, 2);
coeff_rows = coeff.';
if n_dims == 2
    coeff_rows = [coeff_rows, zeros(n_coeff, 1)];
end

% Plot PCA vectors, showing them as lines passing through the bounding box of
% the point cloud
line_colors = {[0, 0, 0], [0.5, 0.5, 0.5], [0.75, 0.75, 0.75]};
if size(coeff, 1) ~= n_dims
    error('`coeff` should have the same number of rows as there are columns in `x`.');
end
origin = zeros(1, 3);
for i = 1:n_coeff
    coeff_i = coeff_rows(i, :);
    [...
        intersections, ~, box_filter...
    ] = lineBoxIntersections(...
        origin, coeff_i + origin,...
        points_xlim, points_ylim, points_zlim...
    );
    intersections = intersections(box_filter, :);
    if n_dims == 2
        line(...
            intersections(:, 1), intersections(:, 2),...
            'LineWidth', 2, 'Color', line_colors{i}, 'LineStyle', '--'...
        )
    elseif n_dims == 3
        line(...
            intersections(:, 1), intersections(:, 2), intersections(:, 3),...
            'LineWidth', 2, 'Color', line_colors{i}, 'LineStyle', '--'...
        )
    else
        error('Unexpected value of `n_dims`, %d.', n_dims);
    end
end

% Plot PCA lines, scaled to cover the bounding box of the data
mu_passed = nargin > 2 && ~isempty(mu);
if mu_passed
    mu = reshape(mu, 1, []);
    if n_dims == 2
        mu = [mu, 0];
    end
    for i = 1:n_coeff
        coeff_i = coeff_rows(i, :);
        [...
            intersections, ~, box_filter...
        ] = lineBoxIntersections(...
            mu, coeff_i + mu,...
            points_xlim, points_ylim, points_zlim...
        );
        intersections = intersections(box_filter, :);
        if n_dims == 2
            line(...
                intersections(:, 1), intersections(:, 2),...
                'LineWidth', 2, 'Color', line_colors{i}...
            )
        elseif n_dims == 3
            line(...
                intersections(:, 1), intersections(:, 2), intersections(:, 3),...
                'LineWidth', 2, 'Color', line_colors{i}...
            )
        else
            error('Unexpected value of `n_dims`, %d.', n_dims);
        end
    end
end

hold off

if nargin > 3
    xlabel(label_str{1});
    ylabel(label_str{2});
    if n_dims == 3
        zlabel(label_str{3});
    end
end

if nargin > 4
    if mu_passed
        legend_str = [legend_str(1), repmat(reshape(legend_str(2:end), 1, []), 1, 2)];
    end
    legend(legend_str);
end

if nargin > 5
    title(title_str);
end
end