function [t, lp] = closestPointToLine(p, l_start, l_end)
% CLOSESTPOINTTOLINE Nearest points on a line defined by two points
%
% ## Syntax
% t = closestPointToLine(p, l_start, l_end)
% [t, lp] = closestPointToLine(p, l_start, l_end)
%
% ## Description
% t = closestPointToLine(p, l_start, l_end)
%   Returns the multiples of the increment between the two points on the
%   line, corresponding to the closest points on the line to the points in
%   `p`.
%
% lp = closestPointToLine(p, l_start, l_end)
%   Additionally returns the points on the line closest to the points in
%   `p`.
%
% ## Input Arguments
%
% p -- Points
%   An n x k array of points, where 'n' is the number of points, and 'k' is
%   the dimensionality of each point.
%
% l_start -- Line start point
%   A 1 x k vector containing the starting point of the line
%
% l_end -- Line end point
%   A 1 x k vector containing the end point of the line
%
% ## Output Arguments
%
% t -- Interpolation factor
%   `t` is the scalar such that `lp(i, :) = l_start + t(i) * (l_end -
%   l_start)`. `t` is a vector of length n.
%
% lp -- Points
%   An n x k array of points, where 'n' is the number of points, and 'k' is
%   the dimensionality of each point. All of the points are on the line
%   defined by `l_start` and `l_end`.
%
% ## References
% - https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
%
% See also pointsOnLinesByCoordinates

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 29, 2018

nargoutchk(1, 2);
narginchk(3, 3);

direction = l_end - l_start;
direction_length = sqrt(dot(direction, direction));

n = size(p, 1);
l_start_rep = repmat(l_start, n, 1);
direction_rep = repmat(direction, n, 1);
start_to_p = p - l_start_rep;
t = dot(start_to_p, direction_rep, 2) / direction_length;

if nargout > 1
    lp = l_start_rep + (t * direction_rep) / direction_length;
end

end
