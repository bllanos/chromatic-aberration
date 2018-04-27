function [vectors_matrix] = matchByVectors(vectors_cell, vector_field, reference_index)
% MATCHBYVECTORS  Associate structures by a vector-valued field
%
% ## Syntax
% vectors_matrix = matchByVectors(vectors_cell, vector_field, reference_index)
%
% ## Description
% vectors_matrix = matchByVectors(vectors_cell, vector_field, reference_index)
%   Organize structures into an array of rows of matched structures
%
% ## Input Arguments
%
% vectors_cell -- Sets of series of structures
%   A cell vector, where each cell contains a 2D structure array with 'n'
%   columns. The cells may contain arrays with different numbers of rows.
%
%   Each cell represents a domain within which matches can be made. In
%   other words, structures in different cells cannot be matched. The
%   columns of the 2D arrays represent different data series. The purpose
%   of this function is to associate structures across data series.
%
% vector_field -- Structure field name
%   A character vector containing the name of the field by which structures
%   are to be matched. The value of this field must be a row vector of the
%   same length for all structures.
%
% reference_index -- Reference data series index
%   The column in the 2D arrays stored in `vectors_cell` which is to be
%   treated as the reference column. Matches will be made between the
%   structures in this column and the closest structures in the other
%   columns, where proximity is measured by the Euclidean distance between
%   the values of the `vector_field` field.
%
%   A reference column is required, because, for efficiency, this function
%   will not assess matches for reciprocity (i.e. whether the structures
%   associated together are all mutual nearest neighbours between columns).
%
% ## Output Arguments
%
% vectors_matrix -- Matched structures
%   A 2D structure array, created by calling 'cell2mat()' on a new version
%   of `vectors_cell`, `x`, where the elements of the array within each
%   cell, `x{i}`, are rearranged relative to those in `vectors_cell{i}`.
%
%   Each array `x{i}` is created by iterating through the elements in the
%   reference column, `vectors_cell{i}(:, reference_index)`. `x{i}(i,
%   reference_index)` is a copy of `vectors_cell{i}(i, reference_index)`,
%   whereas `x{i}(i, j)` is the closest element of `vectors_cell{i}(:, j)`
%   to `x{i}(i, j)`. Proximity is measured as the Euclidean distance
%   between the vector `x{i}(i, j).(vector_field)` and the corresponding
%   vectors of the elements of `vectors_cell{i}(:, j)`.
%
% See also cell2mat

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 27, 2018

n_cells = length(vectors_cell);
n_channels = size(vectors_cell{1}, 2);
vectors_cell_matched = cell(n_cells, 1);

for i = 1:n_cells
    vectors_reference_i = vertcat(vectors_cell{i}(:, reference_index).(vector_field));
    n_in_cell = size(vectors_reference_i, 1);
    vectors_cell_matched_i = vectors_cell{i};
    for j = 1:n_channels
        if j == reference_index
            continue;
        end
        vectors_channel_j = vertcat(vectors_cell{i}(:, j).(vector_field));
        map = zeros(n_in_cell, 1);
        for k = 1:n_in_cell
            vector_k_rep = repmat(vectors_cell{i}(k, j).(vector_field), n_in_cell, 1);
            distances = vectors_channel_j - vector_k_rep;
            distances = dot(distances, distances, 2);
            [~, map(k)] = min(distances);
        end
        vectors_cell_matched_i(:, j) = vectors_cell{i}(map, j);
    end
    vectors_cell_matched{i} = vectors_cell_matched_i;
end

vectors_matrix = cell2mat(vectors_cell_matched);

end

