% UTILITIES
% Version 2.0.0 18-Jul-2019
%
% Helper functions, such as for data input/output, text processing, and geometry
% calculations.
%
% Geometry
%   arePointsOnLineSegments    - Check if points are between the endpoints of line segments
%   closestPointToLine         - Nearest points on a line defined by two points
%   lineBoxIntersections       - Find intersection points of a line with an axis-aligned rectangular prism
%   normalizePointsPCA         - Normalization of point coordinates using Principal Components Analysis
%   plotPCA                    - Plot a point cloud and its PCA approximation
%   pointsOnLinesByCoordinates - Find points on lines between two endpoints by one coordinate
%
% Data exploration
%   listFiles                  - List files matching a wildcard
%   listFilesRecursive         - List files matching a regular expression in a directory tree
%   TestMATFileChanges         - Test script to compare the '.mat' files in two directories
%
% Data input and output
%   loadImage                  - Load an image from an image file or a data file
%   loadVariables              - Load variables from a '.mat' file
%   saveImages                 - Save images to image files and data files
%   ShrinkMATFiles             - Convert large double precision variables to single precision variables in '.mat' files
%
% Data manipulation
%   clipAndRemap               - Clip values and convert to an integer type
%   mergeStructs               - Merge two structures
%   trimCommon                 - Extract unique parts of filepaths
