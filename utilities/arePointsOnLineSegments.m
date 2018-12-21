function [ out ] = arePointsOnLineSegments(p, start, endpoint, tolerance)
% AREPOINTSONLINESEGMENTS Check if points are between the endpoints of line segments
%
% ## Syntax
% out = arePointsOnLineSegments(p, start, endpoint, tolerance)
%
% ## Description
% out = arePointsOnLineSegments(p, start, endpoint, tolerance)
%   Returns true if `p` are on the line segments between `start` and
%   `endpoint`, to within a distance of `tolerance`.
%
% ## Input Arguments
%
% p -- Test point
%   An n x k array of query points, where 'n' is the number of points, and
%   'k' is the dimensionality of each point.
%
% start -- Line segment endpoint
%   An n x k array containing one of the two endpoints defining each line
%   segment.
%
% endpoint -- Line segment endpoint
%   An n x k array containing the other endpoints defining the line
%   segments.
%
% tolerance -- Error threshold
%   A scalar threshold for the perpendicular distances from `p` to the
%   lines which, if exceeded, leads to the conclusion that `p` is not on
%   the lines.
%
% ## Output Arguments
%
% out -- Test results
%   An n x 1 logical array where `out(i)` is true if and only if `p(i, :)`
%   is on the line segment between `start(i, :)` and `endpoint(i, :)`,
%   given the value of `tolerance`.
%
%   Test results are false for points that are outside of the 0 to 1
%   interpolation region between `start` and `endpoint`, as determined by
%   testing calculated interpolation factors (for any coordinates).
%   `tolerance` is not applied to this portion of the test - The boundary
%   values are set to exactly 0 and 1.
%
% ## Tips
% - The results of this function are unchanged if the values of `start` and
%   `endpoint` are exchanged.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 10, 2016

nargoutchk(1, 1);

direction = endpoint - start;
start_to_p = p - start;
t = start_to_p ./ direction;
out = all(isfinite(t) & t >= 0 & t <= 1, 2);

% See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
dimensionality = size(direction, 2);
unit_direction = direction ./ repmat(sqrt(dot(direction, direction, 2)), 1, dimensionality);
distances = start_to_p - (repmat(dot(start_to_p, unit_direction, 2), 1, dimensionality) .* unit_direction);
distances = sqrt(dot(distances, distances, 2));

% A faster test would be to calculate if `t` is the same to within a
% tolerance over all coordinates, but this test has no absolute geometric
% meaning.
out = out & all(distances < tolerance, 2);

end
