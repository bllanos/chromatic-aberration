function fg_out = plotEllipse(ellipse_to_world, varargin)
% PLOTELLIPSE  Plot an ellipse
%
% ## Syntax
% plotEllipse(ellipse_to_world [, fg])
% fg = plotEllipse(ellipse_to_world)
%
% ## Description
% plotEllipse(ellipse_to_world [, fg])
%   Plot an ellipse in a new or existing figure
%
% fg = plotEllipse(ellipse_to_world)
%   Plot an ellipse in a new figure, and return a handle to the new figure
%
% ## Input Arguments
%
% ellipse_to_world -- Ellipse coordinate transformation
%   A matrix transforming 2D points, 'p', on the unit disk, to points on
%   the ellipse. Points on the unit circle will be mapped to points on the
%   boundary of the ellipse. Specifically, a point (x1, y1) and its image
%   (x2, y2) on the ellipse are related by the equation:
%
%     [x2 y2 1] = (ellipse_to_world * [x1 y1 1].').'
%
% fg -- Figure handle
%   A handle to the figure to update.
%
%   If not passed, a new figure is created.
%
% ## Output Arguments
%
% fg -- Figure handle
%   A handle to the new figure created.
%
% See also ellipseModel

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 20, 2018

narginchk(1, 2);

if isempty(varargin)
    nargoutchk(0, 1);
    fg_out = figure;
else
    nargoutchk(0, 0);
    fg_out = varargin{1};
    figure(fg_out);
end

n_samples = 20;
angle = linspace(-pi, pi, n_samples).';
p = [cos(angle), sin(angle), ones(n_samples, 1)];
p_ellipse = (ellipse_to_world * p.').';

hold on
plot(p_ellipse(:, 1), p_ellipse(:, 2), 'b-');
plot(ellipse_to_world(1, end), ellipse_to_world(2, end), 'bo');
hold off

end

