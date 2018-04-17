function [newPoints, T] = normalizePointsPCA(points)
%NORMALIZEPOINTSPCA Normalization of point coordinates using Principal Components Analysis
%
% Input
%
% points - An n x (d + 1) array representing the homogenous coordinates of
%          points, assumed to be scaled such that `points(:, end) = 1`.
%
% Output
%
% newPoints - Normalized point coordinates corresponding to `points`. All
%             'd' dimensions of the input points will be normalized such
%             that the centroid of the points is at the origin, and the
%             variances of the points in each principal component direction
%             is 1.
%
% T - The transform matrix which normalized the points
%     `newPoints = (T * points.').'`
%
% Created for: CMPUT 615 Lab 3: Geometry II and SFM
% Bernard Llanos
% February 17, 2017

d = size(points, 2) - 1; % Dimensionality
[coeff,score,~,~,~,mu] = pca(points(:, 1:d));
Translation = [
    eye(d) mu.';
    zeros(1, d) 1
    ];
muDist = var(score);
% Scaling factors so that variances in principal component
% directions are 1.
s = sqrt(muDist) .^ (-1);
T = [
    diag(s) zeros(d, 1);
    zeros(1, d) 1
    ] / [
    coeff zeros(d, 1);
    zeros(1, d) 1
    ] / Translation;
newPoints = (T * points.').';
end

