function v_adj = neighborVertices(TR, vi)
% NEIGHBORVERTICES  Compute adjacency relationships in a triangulation
%
% ## Syntax
% v_adj = neighborVertices(TR, vi)
%
% ## Description
% v_adj = neighborVertices(TR, vi)
%   Find the vertices adjacent to the given vertices in the triangulation.
%
% ## Input Arguments
%
% TR -- Triangulation
%   An instance of MATLAB's `triangulation` class.
%
% vi -- Vertex indices
%   The indices of the vertices (indices into `TR.ConnectivityList`) whose
%   neighbouring vertices are to be determined. A vector of length n.
%
% ## Output Arguments
%
% v_adj -- Adjacent vertex indices
%   The indices of the vertices (indices into `TR.ConnectivityList`) of the
%   vertices which are connected with the vertices specified by `vi`.
%   `v_adj` is a cell vector of length n, where `v_adj{i}` is a column
%   vector containing the indices of the vertices adjacent to the vertex
%   with index `vi(i)`.
%
% See also triangulation

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 28, 2017

nargoutchk(1, 1);
narginchk(2, 2);

ti = vertexAttachments(TR,vi);
n_vi = length(vi);
v_adj = cell(n_vi, 1);
tri = TR.ConnectivityList;
for k = 1:n_vi
    v_adj_k = tri(ti{k}, :);
    v_adj_k = unique(v_adj_k(:));
    v_adj{k} = v_adj_k(v_adj_k ~= vi(k));
end
end