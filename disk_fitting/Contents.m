% DISK_FITTING
% Version 2.0.0 18-Jul-2019
%
% Disk keypoint detection for dispersion calibration (approximating dispersion
% as lateral chromatic aberration)
%
% Dispersion model generation scripts
%   RAWDiskDispersion - Chromatic aberration calibration from captured images
%
% Disk keypoint localization
%   ellipseModel      - Create functions for evaluating a parametric ellipse
%   findAndFitDisks   - Fit ellipses to blobs in an image
%   plotEllipse       - Plot an ellipse
%   refineDisk        - Improve an ellipse fit to an image blob
