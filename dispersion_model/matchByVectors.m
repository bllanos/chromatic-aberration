function [vectors_matrix] = matchByVectors(vectors_cell, vector_field, reference_index, varargin)
% MATCHBYVECTORS  Associate structures by a vector-valued field
%
% ## Syntax
% vectors_matrix = matchByVectors(vectors_cell, vector_field, reference_index [, threshold])
%
% ## Description
% vectors_matrix = matchByVectors(vectors_cell, vector_field, reference_index [, threshold])
%   Organize structures into an array of rows of matched structures
%
% ## Input Arguments
%
% vectors_cell -- Sets of series of structures
%   A cell vector, where each cell contains a cell vector of length 'n'.
%   Each cell of these inner cell vectors contains a column vector of
%   structures. The column vectors may have different lengths.
%
%   Each top-level cell represents a domain within which matches can be
%   made. In other words, structures in different top-level cells cannot be
%   matched. The inner cells represent different data series. The purpose
%   of this function is to associate structures between data series.
%
% vector_field -- Structure field name
%   A character vector containing the name of the field by which structures
%   are to be matched. The value of this field must be a row vector of the
%   same length for all structures.
%
% reference_index -- Reference data series index
%   The index into the inner cell vectors in `vectors_cell` which is to be
%   treated as the reference index. Matches will be made between the
%   structures in the inner cells of this index and the closest structures
%   in the other inner cells, where proximity is measured by the Euclidean
%   distance between the values of the `vector_field` field.
%
%   A reference index is required, because, for efficiency, this function
%   will not assess matches for n-way reciprocity (i.e. whether the
%   structures associated together are all mutual nearest neighbours
%   between inner cells).
%
% threshold -- Outlier detection threshold
%   The number of sample standard deviations larger than the mean Euclidean
%   distance between matched structures beyond which matches are rejected. The
%   mean and sample standard deviation are calculated for each of the inner
%   cells (other than the reference inner cell), rather than globally over all
%   matches made within each cell of `vectors_cell`.
%
%   Outlier detection is disabled if `threshold` is zero, or is not passed.
%
% ## Output Arguments
%
% vectors_matrix -- Matched structures
%   A 2D structure array, created by calling 'cell2mat()' on a new version
%   of `vectors_cell`, `x`, where the structures within each cell, `x(i)`,
%   are rearranged and filtered relative to those in `vectors_cell(i)`.
%
%   Each 2D structure array `x{i}` is created by iterating through the
%   structures in the reference inner cell's structure vector,
%   `vectors_cell{i}{reference_index}`. `x{i}(m, reference_index)` is a
%   copy of `vectors_cell{i}{reference_index}(k)`, whereas `x{i}(m, j)` is
%   the closest element of `vectors_cell{i}{j}` to `x{i}(m,
%   reference_index)`. Proximity is measured as the Euclidean distance
%   between the vector `x{i}(m, reference_index).(vector_field)` and the
%   vectors `vectors_cell{i}{j}.(vector_field)`.
%
%   The relationship between the indices 'm' and 'k' in the above
%   description is determined as follows:
%   - The elements of `x{i}(:, reference_index)` are generally in the order
%     that they were in `vectors_cell{i}{reference_index}`.
%   - Any element `vectors_cell{i}{reference_index}(k)` having a closest
%     match in another series which is closer still to a different element
%     `vectors_cell{i}{reference_index}(p)` is not output in `x`.
%   - After the preceding test, any element
%     `vectors_cell{i}{reference_index}(k)` having a match in another series
%     which is further away than `threshold` sample standard deviations from the
%     mean distance to matches in that series is not output in `x`.
%
%   Consequently, only the elements of `vectors_cell{i}{q}`, `q ~=
%   reference_index` which are the closest to some element of
%   `vectors_cell{i}{reference_index}` are output in `x`.
%
% See also cell2mat

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 27, 2018

narginchk(3, 4);
nargoutchk(1, 1);

threshold = 0;
if ~isempty(varargin)
    threshold = varargin{1};
    if threshold < 0
        error('`threshold` must be a nonnegative number.');
    end
end

n_cells = length(vectors_cell);
n_channels = length(vectors_cell{1});
vectors_cell_matched = cell(n_cells, 1);

names = fieldnames(vectors_cell{1}{1});
n_names = length(names);
stats_args = cell(n_names * 2, 1);

for i = 1:n_cells
    % Find the nearest neighbours of the reference series elements
    vectors_reference_i = vertcat(vectors_cell{i}{reference_index}.(vector_field));
    n_reference_i = size(vectors_reference_i, 1);
    neighbours_in_other_series = zeros(n_reference_i, n_channels - 1);
    distances_to_other_series = zeros(n_reference_i, n_channels - 1);
    nearest_neighbour_from_other_series = false(n_reference_i, n_channels - 1);
    channel_ind = 1;
    for j = 1:n_channels
        if j == reference_index
            continue;
        end
        vectors_channel_j = vertcat(vectors_cell{i}{j}.(vector_field));
        n_channel_j = size(vectors_channel_j, 1);
        % Outgoing nearest neighbours
        for k = 1:n_reference_i
            vector_k_rep = repmat(vectors_reference_i(k, :), n_channel_j, 1);
            distances = vectors_channel_j - vector_k_rep;
            distances = dot(distances, distances, 2);
            [...
                distances_to_other_series(k, channel_ind),...
                neighbours_in_other_series(k, channel_ind)...
            ] = min(distances);
        end
        % Incoming nearest neighbours
        for k = 1:n_channel_j
            vector_k_rep = repmat(vectors_channel_j(k, :), n_reference_i, 1);
            distances = vectors_reference_i - vector_k_rep;
            distances = dot(distances, distances, 2);
            [~, min_ind] = min(distances);
            nearest_neighbour_from_other_series(min_ind, channel_ind) = true;
        end
        channel_ind = channel_ind + 1;
    end
    
    % Filter using a mutual consistency constraint centered on the
    % reference series
    filter_i = all(nearest_neighbour_from_other_series, 2);
    n_reference_filtered = sum(filter_i);
    neighbours_in_other_series = neighbours_in_other_series(filter_i, :);
    distances_to_other_series = distances_to_other_series(filter_i, :);
    
    % Filter outlier distances
    % Reference: MATLAB documentation page on "Inconsistent Data"
    if threshold > 0 && n_reference_filtered > 2
        channel_ind = 1;
        inliers_filter = true(n_reference_filtered, 1);
        for j = 1:n_channels
            if j == reference_index
                continue;
            end
            sigma_distances = std(distances_to_other_series(:, channel_ind));
            if sigma_distances > 0
                mu_distances = mean(distances_to_other_series(:, channel_ind));
                mu_distances_rep = repmat(mu_distances, n_reference_filtered, 1);
                sigma_distances_rep = repmat(sigma_distances, n_reference_filtered, 1);
                inliers_filter = inliers_filter &...
                    ((distances_to_other_series(:, channel_ind) - mu_distances_rep) < (threshold * sigma_distances_rep));
            end
            channel_ind = channel_ind + 1;
        end
        n_reference_filtered = sum(inliers_filter);
        filter_i(filter_i) = inliers_filter;
        neighbours_in_other_series = neighbours_in_other_series(inliers_filter, :);
    end
    
    % Preallocate the output structure array
    for m = 1:length(stats_args)
        if mod(m, 2) ~= 0
            stats_args{m} = names{(m + 1) / 2};
        else
            stats_args{m} = cell(n_reference_filtered, n_channels);
        end
    end
    vectors_cell_matched_i = struct(stats_args{:});
    
    % Output the results
    channel_ind = 1;
    for j = 1:n_channels
        if j == reference_index
            vectors_cell_matched_i(:, reference_index) = vectors_cell{i}{reference_index}(filter_i);
        else
            vectors_cell_matched_i(:, j) = vectors_cell{i}{j}(neighbours_in_other_series(:, channel_ind));
            channel_ind = channel_ind + 1;
        end
    end
    
    vectors_cell_matched{i} = vectors_cell_matched_i;
end

vectors_matrix = cell2mat(vectors_cell_matched);

end

