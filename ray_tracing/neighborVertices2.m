function [ adj_matrix, v_adj ] = neighborVertices2(TR)
% NEIGHBORVERTICES2  Compute adjacency relationships in a triangulation
%
% ## Syntax
% adj_matrix = neighborVertices(TR)
% [ adj_matrix, v_adj ] = neighborVertices2(TR)
%
% ## Description
% adj_matrix = neighborVertices(TR)
%   Find the vertices adjacent to all vertices in the triangulation, and
%   output an adjacency matrix.
% [ adj_matrix, v_adj ] = neighborVertices2(TR)
%   Additionally returns a set of adjacency lists.
%
% ## Input Arguments
%
% TR -- Triangulation
%   An instance of MATLAB's `triangulation` class.
%
% ## Output Arguments
%
% adj_matrix -- Vertex adjacency matrix
%   `adj_matrix` is a sparse square 2D array with the same number of rows
%   as `TR.Points`. `adj_matrix(i,j)` is nonzero if there is an edge
%   bewtween the vertices `TR.Points(i, :)` and `TR.Points(j, :)` in the
%   triangulation `TR`. `adj_matrix(i,j)` is nonzero if and only if
%   `adj_matrix(j,i)` is nonzero.
%
% v_adj -- Vertex adjacency lists
%   `v_adj` is a cell vector of length n, where `v_adj{i}` is a column
%   vector containing the indices of the vertices (rows in `TR.Points`)
%   adjacent to the vertex `TR.Points(i, :)`.
%
% ## Notes
% - `neighborVertices()` is more efficient when adjacency information is
%   only needed for a small number of vertices.
% - This function presently only works for 2D triangulations
%
% See also triangulation, neighborVertices

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 12, 2017

nargoutchk(1, 2);
narginchk(1, 1);

if size(TR.Points, 2) > 2
    error('Only 2D triangulations are presently supported.')
end

n = size(TR.Points, 1);
tri = TR.ConnectivityList;
i = [
    tri(:);
    tri(:, 2);
    tri(:, 3);
    tri(:, 1)
    ];
j = [
    tri(:, 2);
    tri(:, 3);
    tri(:, 1);
    tri(:)
    ];
v = ones(length(i), 1);
adj_matrix = sparse(i, j, v, n, n);
if nargout > 1
    v_adj = cell(n, 1);
    for i = 1:n
        v_adj{i} = find(adj_matrix(:, i));
    end
end