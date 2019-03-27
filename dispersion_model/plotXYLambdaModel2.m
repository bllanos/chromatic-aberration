function plotXYLambdaModel2(spatial_domain, varargin)
% PLOTXYLAMBDAMODEL2  Plot a model of dispersion in two or three variables
%
% ## Usage
%
% 'plotXYLambdaModel2()' is a version of 'plotXYLambdaModel()' that produces
% visualizations of a model of dispersion without comparison to the measurements
% used to fit the model. This function is also intended for generating
% print-worthy figures, and so makes a number of figure style adjustments.
%
% ## Syntax
% plotXYLambdaModel2(...
%     spatial_domain, lambda, lambda_ref, n_lambda_plot, dispersionfun...
% )
% plotXYLambdaModel2(...
%     spatial_domain, channel_ref, dispersionfun...
% )
%
% ## Description
% plotXYLambdaModel2(...
%     spatial_domain, lambda, lambda_ref, n_lambda_plot, dispersionfun...
% )
%   Creates figures for graphical inspection of a model of dispersion
%
% plotXYLambdaModel2(...
%     spatial_domain, channel_ref, dispersionfun...
% )
%   Identical operation, but for dispersion measured between colour channels
%   instead of wavelength bands
%
% ## Input Arguments
%
% spatial_domain -- Domain for the first two independent variables (x,y)
%   A four element vector giving the (x, y) coordinates of the top left and
%   bottom right corners of the rectangular region in space in which to sample
%   the model of dispersion. i.e. `spatial_domain` is of the form `[x_1, y_1,
%   x_2, y_2]`. This function will determine whether the y-axis needs to be
%   reversed (for image axis conventions) based on the y-order of the two corner
%   points.
%
% lambda -- Data for the third independent variable (wavelength)
%   A vector of values of the third independent variable at which the model of
%   dispersion can be evaluated.
%
% lambda_ref -- Reference wavelength
%   Dispersion vectors are calculated relative to image positions for this
%   wavelength. This input argument is only used for plot annotation.
%
% n_lambda_plot -- Wavelength sample count
%   The number of plots to generate (each for a different wavelength).
%   Wavelengths will be selected from 'lambda' using a procedure which favours
%   even spacing in the spectrum.
%
% channel_ref -- Reference colour channel
%   Dispersion vectors are calculated relative to image positions for the
%   colour channel of this index in the series of colour channels. This
%   input argument is only used for plot annotation.
%
% dispersionfun -- Model of dispersion
%   The 'dispersionfun' output argument of 'makeDispersionfun()', which
%   evaluates models of dispersion in X, Y, and lambda, or which evaluates
%   a separate set of models in X and Y for each colour channel.
%
% See also xylambdaPolyfit, xylambdaSplinefit, plotXYLambdaModel

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 26, 2019

nargoutchk(0, 0);
narginchk(3, 5);
if nargin == 4
    error('Unexpected number of input arguments.')
end

channel_mode = (nargin == 3);
if channel_mode
    index_ref = varargin{1};
    dispersionfun = varargin{2};
    n_plot = 3; % Number of RGB channels
    lambda = 1:n_plot;
else
    lambda = varargin{1};
    index_ref = varargin{2};
    n_plot = varargin{3};
    dispersionfun = varargin{4};
end

if channel_mode
    lambda_samples = lambda;
else
    lambda_samples = linspace(lambda(1), lambda(end), n_plot);
    lambda_samples_ind = zeros(n_plot, 1);
    for k = 1:n_plot
        [~, lambda_samples_ind(k)] = min(abs(lambda - lambda_samples(k)));
        lambda_samples(k) = lambda(lambda_samples_ind(k));
    end
end

max_x = max(spatial_domain([1, 3]));
max_y = max(spatial_domain([2, 4]));
min_x = min(spatial_domain([1, 3]));
min_y = min(spatial_domain([2, 4]));
domain_center = [mean(spatial_domain([1, 3])), mean(spatial_domain([2, 4]))];

% Individual plots for a sample of the wavelengths, or for all colour
% channels
xy_sampling_surface = [200 200];
xy_sampling_vectors = [10 10];
for k = 1:n_plot    
    % Evaluate the model

    
    [x_grid_s, y_grid_s] = meshgrid(linspace(min_x, max_x, xy_sampling_surface(1)), linspace(min_y, max_y, xy_sampling_surface(2)));
    dataset_grid_s = [x_grid_s(:), y_grid_s(:), repmat(lambda_samples(k), numel(x_grid_s), 1)];
    disparity_grid_s = dispersionfun(dataset_grid_s);
    disparity_grid_s_mag = sqrt(dot(disparity_grid_s, disparity_grid_s, 2));
    disparity_grid_s_sign = sign(dot(...
        disparity_grid_s, dataset_grid_s(:, 1:2) - repmat(domain_center, numel(x_grid_s), 1), 2 ...
    ));
    disparity_grid_s_mag = disparity_grid_s_mag .* disparity_grid_s_sign;
    disparity_grid_s_mag_matrix = reshape(disparity_grid_s_mag, size(x_grid_s));
    
    [x_grid_v, y_grid_v] = meshgrid(linspace(min_x, max_x, xy_sampling_vectors(1)), linspace(min_y, max_y, xy_sampling_vectors(2)));
    dataset_grid_v = [x_grid_v(:), y_grid_v(:), repmat(lambda_samples(k), numel(x_grid_v), 1)];
    disparity_grid_v = dispersionfun(dataset_grid_v);
    disparity_grid_v_mag = sqrt(dot(disparity_grid_v, disparity_grid_v, 2));
    disparity_grid_v_sign = sign(dot(...
        disparity_grid_v, dataset_grid_v(:, 1:2) - repmat(domain_center, numel(x_grid_v), 1), 2 ...
    ));
    disparity_grid_v_mag = disparity_grid_v_mag .* disparity_grid_v_sign;
        
    disparity_grid_s_theta = atan2d(disparity_grid_s(:, 2), disparity_grid_s(:, 1));
    disparity_grid_s_theta_matrix = reshape(disparity_grid_s_theta, size(x_grid_s));
    
    % Visualization of the disparity magnitude
    figure;
    hold on
    s = surf(x_grid_s, y_grid_s, disparity_grid_s_mag_matrix);
    set(s, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    c = colorbar;
    c.Label.String = 'Signed dispersion magnitude';
    quiver3(...
        dataset_grid_v(:, 1), dataset_grid_v(:, 2), disparity_grid_v_mag,...
        disparity_grid_v(:, 1), disparity_grid_v(:, 2), zeros(size(disparity_grid_v_mag)),...
        'k', 'LineWidth', 0.5 ...
    );
    hold off
    xlabel('X')
    ylabel('Y')
    axis equal
    if spatial_domain(4) > spatial_domain(2)
        axis ij;
    end
    if channel_mode
        zlabel(sprintf('Magnitude of dispersion wrt Channel %d', index_ref))
    else
        zlabel(sprintf('Magnitude of dispersion wrt \\lambda = %g', index_ref))
    end
    if channel_mode
        title(sprintf('Evaluation of the model of dispersion for Channel %d', lambda_samples(k)));
    else
        title(sprintf('Evaluation of the model of dispersion for \\lambda = %g', lambda_samples(k)));
    end
        
    % Visualization of the disparity direction
    figure;
    hold on
    s = pcolor(x_grid_s, y_grid_s, disparity_grid_s_theta_matrix);
    set(s, 'EdgeColor', 'none');
    c = colorbar;
    c.Label.String = 'Dispersion direction [degrees]';
    quiver(...
        dataset_grid_v(:, 1), dataset_grid_v(:, 2),...
        disparity_grid_v(:, 1), disparity_grid_v(:, 2),...
        'k', 'LineWidth', 0.5 ...
    );
    hold off
    xlabel('Image x-coordinate')
    ylabel('Image y-coordinate')
    axis image
    if spatial_domain(4) > spatial_domain(2)
        axis ij;
    end
    if channel_mode
        title(sprintf('Evaluation of the model of dispersion for Channel %d', lambda_samples(k)));
    else
        title(sprintf('Evaluation of the model of dispersion for \\lambda = %g', lambda_samples(k)));
    end
end
    
end