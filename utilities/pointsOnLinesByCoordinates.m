function [ p ] = pointsOnLinesByCoordinates(start, endpoint, coord, coordIndex)
% POINTSONLINESBYCOORDINATES Find points on lines between two endpoints by one coordinate
%
% ## Syntax
% p = pointsOnLinesByCoordinates(start, endpoint, coord, coordIndex)
%
% ## Description
% p = pointsOnLinesByCoordinates(start, endpoint, coord, coordIndex)
%   Returns the full coordinates of points on the lines defined by the two
%   endpoints given one known coordinate of each point.
%
% ## Input Arguments
%
% start -- Line control point
%   An n x k array containing one of the two points defining each line. 'n'
%   is the number of lines, and 'k' is the dimensionality of the points.
%
% endpoint -- Line control point
%   An n x k array containing the other points defining the lines.
%
% coord -- Coordinate values
%   An n x 1 array containing values of coordinates that the output points
%   must satisfy.
%
% coordIndex -- Coordinate indices
%   An n x 1 array clarifying `coord` by indicating to which dimensions of
%   point coordinates the elements of `coord` correspond.
%
% ## Output Arguments
%
% p -- Points on lines
%   An n x k array of points. `p(i, :)` is a point on the line defined by
%   `start(i, :)` and `endpoint(i, :)` having `coord(i)` as the value of
%   its `coordIndex(i)`-th coordinate.
%
%   If such a point `p(i, :)` on the line is impossible, then `p(i, :)` is
%   equal to `nan(1, k)`.
%
% ## Tips
% - The results of this function are unchanged if the values of `start` and
%   `endpoint` are exchanged.

% Bernard Llanos
% Spring 2016 research assistantship supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 10, 2016

nargoutchk(1, 1);

coordIndexLinear = sub2ind(size(start), (1:size(start, 1))', coordIndex);
direction = endpoint - start;
t = (coord - start(coordIndexLinear)) ./ direction(coordIndexLinear);
t(isnan(t)) = 0; % Assuming 0 / 0 = NaN is the only reason that elements might be NaN
p = start + repmat(t, 1, size(start, 2)) .* direction;
% Correct numerical error by avoiding recalculation
p(coordIndexLinear) = coord;
p(isinf(t), :) = NaN; % Points not on the line

end

