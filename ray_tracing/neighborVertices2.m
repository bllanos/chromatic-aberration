function adj_matrix = neighborVertices2(TR)
% NEIGHBORVERTICES2  Compute adjacency relationships in a triangulation
%
% ## Syntax
% adj_matrix = neighborVertices(TR)
%
% ## Description
% adj_matrix = neighborVertices(TR)
%   Find the vertices adjacent to all vertices in the triangulation.
%
% ## Input Arguments
%
% TR -- Triangulation
%   An instance of MATLAB's `triangulation` class.
%
% ## Output Arguments
%
% adj_matrix -- Vertex adjacency matrix
%   `adj_matrix` is a square 2D logical array with the same number of rows
%   as `TR.Points`. `adj_matrix(i,j)` is `true` if there is an edge
%   bewtween the vertices `TR.Points(i, :)` and `TR.Points(j, :)` in the
%   triangulation `TR`. `adj_matrix(i,j)` is `true` if and only if
%   `adj_matrix(j,i)` is `true`.
%
% ## Notes
% - `neighborVertices()` is more efficient when adjacency information is
%   only needed for a small number of vertices.
%
% See also triangulation, neighborVertices

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 12, 2017

nargoutchk(1, 1);
narginchk(1, 1);

tri = TR.ConnectivityList;
n = size(TR.Points, 1);
d = size(TR.ConnectivityList, 2);
sz = [n n];
adj_matrix = false(sz);
for i = 1:d
    if i == d
        j = 1;
    else
        j = i + 1;
    end
    adj_matrix_ind = sub2ind(sz, tri(:, i), tri(:, j));
    adj_matrix(adj_matrix_ind) = true;
end
% Make edges undirected
adj_matrix = (adj_matrix | adj_matrix.');
end