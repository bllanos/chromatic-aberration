function [ intersections, box_plane_coord_indices, box_filter, line_filter ] = lineBoxIntersections( start, endpoint, xlim, ylim, zlim, varargin )
% LINEBOXINTERSECTIONS  Find intersection points of a line with an axis-aligned rectangular prism
%
% ## Syntax
% intersections = lineBoxIntersections( start, endpoint, xlim, ylim, zlim )
% [ intersections, box_plane_coord_indices ] = lineBoxIntersections( start, endpoint, xlim, ylim, zlim )
% [ intersections, box_plane_coord_indices, box_filter ] = lineBoxIntersections( start, endpoint, xlim, ylim, zlim )
% [ intersections, box_plane_coord_indices, box_filter, line_filter ] = lineBoxIntersections(
%   start, endpoint, xlim, ylim, zlim, tolerance
% )
%
% ## Description
% intersections = lineBoxIntersections( start, endpoint, xlim, ylim, zlim )
%   Returns the coordinates of the points of intersection of the line
%   segment with each of the 6 planes defining the axis-aligned rectangular
%   prism.
%
% [ intersections, box_plane_coord_indices ] = lineBoxIntersections( start, endpoint, xlim, ylim, zlim )
%   Additionally returns a map used for determining which intersection
%   points correspond to which planes.
%
% [ intersections, box_plane_coord_indices, box_filter ] = lineBoxIntersections( start, endpoint, xlim, ylim, zlim )
%   Additionally returns a logical array identifying intersection points
%   which are within the faces of the prism.
%
% [ intersections, box_plane_coord_indices, box_filter, line_filter ] = lineBoxIntersections(
%   start, endpoint, xlim, ylim, zlim, tolerance
% )
%   Additionally returns a logical array identifying intersection points
%   which are between the endpoints of the line segment.
%
% ## Input Arguments
%
% start -- Line segment endpoint
%   A 1 x 3 array containing one of the two endpoints defining the line
%   segment.
%
% endpoint -- Line segment endpoint
%   An 1 x 3 array containing the other endpoint defining the line segment.
%
% xlim -- Box limits in the x-dimension
%   A two-element vector specifying the start and end of the box in the
%   x-direction.
%
% ylim -- Box limits in the y-dimension
%   A two-element vector specifying the start and end of the box in the
%   y-direction.
%
% zlim -- Box limits in the z-dimension
%   A two-element vector specifying the start and end of the box in the
%   z-direction.
%
% tolerance -- Error threshold
%   A scalar threshold for the perpendicular distances from the
%   intersection points to the line defined by the endpoints which, if
%   exceeded, leads to the conclusion that a given intersection point is
%   not on the line segment. (This argument is needed to compensate for
%   numerical error, and is passed to `arePointsOnLineSegments`.)
%
% ## Output Arguments
%
% intersections -- Intersection points
%   An 6 x 3 array, whose columns store the x, y and z-coordinates of the
%   intersection points of the line defined by the points `start` and
%   `endpoint` with the planes coinciding with the faces of the rectangular
%   prism.
%
%   `intersection(i, :)` will be equal to `inf(1, 3)` if there is no such
%   intersection point (as determined by `pointsOnLinesByCoordinates`).
%
% box_plane_coord_indices -- Intersection points to plane mapping
%   A 6 x 1 array, where the `box_plane_coord_indices(i)` is the index of
%   the coordinate that is fixed along the plane corresponding to the
%   intersection point `intersections(i, :)`. (1 for x, 2 for y and 3 for
%   z.) The specific plane can be determined by looking at the value of
%   `intersections(i, box_plane_coord_indices(i))`.
%
%   For example, if `box_plane_coord_indices(i) == 3` and
%   `intersections(i, 3) == zlim(1)`, then the intersection point is on the
%   plane `z = zlim(1)`.
%
% box_filter -- Filter to find intersection points on prism faces
%   A 6 x 1 logical array, where `box_filter(i)` is true if
%   `intersections(i, :)` is within the face of the prism that is
%   coincident on the corresponding plane.
%
% line_filter -- Filter to find intersection points on the line segment
%   A 6 x 1 logical array, where `line_filter(i)` is true if
%   `intersections(i, :)` is between the two endpoints of the line segment
%   as opposed to somewhere else along the line (or not on the line at all).
%
% See also arePointsOnLineSegments, pointsOnLinesByCoordinates

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 10, 2016

nargoutchk(1, 4);
if nargout < 4
    narginchk(5, 5);
else
    narginchk(6, 6);
end

box_plane_coords = [ xlim'; ylim'; zlim' ];
box_plane_coord_indices = [ 1; 1; 2; 2; 3; 3 ];
n_planes = length(box_plane_coords);

cs = repmat(start, n_planes, 1);
ps = repmat(endpoint, n_planes, 1);
intersections = pointsOnLinesByCoordinates(cs, ps, box_plane_coords, box_plane_coord_indices);

if nargout > 2
    box_filter =...
        intersections(:, 1) >= xlim(1) & intersections(:, 1) <= xlim(2) &...
        intersections(:, 2) >= ylim(1) & intersections(:, 2) <= ylim(2) &...
        intersections(:, 3) >= zlim(1) & intersections(:, 3) <= zlim(2);
end

if nargout > 3
    tolerance = varargin{1};
    line_filter = arePointsOnLineSegments(intersections, cs, ps, tolerance);
end

end
