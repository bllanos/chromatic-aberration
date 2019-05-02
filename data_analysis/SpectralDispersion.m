%% Compare spectral dispersion models
% Visualize dispersion as a function of wavelength at different image
% locations
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### Model of dispersions
%
% A list of '.mat' files containing several variables, which are the output
% of 'RAWDiskDispersion.m', for example. The following variables are
% required in each file:
% - 'dispersion_data': A model of chromatic aberration, modeling the warping
%   from the reference wavelength band to the other wavelength bands.
%   `dispersion_data` can be converted to a function form using
%   `dispersionfun = makeDispersionfun(dispersion_data)`.
% - 'model_from_reference': A parameter of the above scripts, which
%   determines the frame of reference for the model of chromatic
%   aberration. It must be set to `true`.
% - 'bands': A vector containing the wavelengths or colour channel indices
%   at which the dispersion model was originally fit. Only the value of
%   'bands' from the first file is used.
%
% The following two additional variables are optional. If they are present,
% they will be used to convert between the coordinate system in which the
% model of chromatic aberration was constructed, and the image coordinate
% system. The first variable, 'model_space' is a structure with same form
% as the `model_space` input argument of 'modelSpaceTransform()'. The
% second variable, `fill`, can be omitted, in which case it defaults to
% `false`. `fill` corresponds to the `fill` input argument of
% 'modelSpaceTransform()'. Refer to the documentation of
% 'modelSpaceTransform.m' for details.
%
% ## Output
%
% ### Graphical output
%
% Figures showing the dispersion vector components as a function of
% wavelength for different positions in the image plane
%
% ### Data file output
%
% A '.mat' file containing the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 26, 2019

% List of parameters to save with results
parameters_list = {
    'image_sampling',...
    'lambda',...
    'n_sampled_pixels',...
    'plot_bound',...
    'forward_dispersion_model_filenames',...
    'names'...
};

%% Input data and parameters

% Image dimensions (number of rows, number of columns)
image_sampling = [2048, 2448];

% Spectral sampling
lambda = linspace(350, 750, 250).';

% Number of pixels to sample along a half-diagonal starting from the top
% left of the image to the image centre
n_sampled_pixels = 10;

% Bound on dispersion to use for plot axes limits
plot_bound = 6;

% Models of dispersion
forward_dispersion_model_filenames = {...
    'C:\Users\GraphicsLab\Documents\llanos\Data\20190208_ComputarLens\dispersion\spectral\RAWDiskDispersionResults_spectral_polynomial_fromReference.mat';
    'C:\Users\GraphicsLab\Documents\llanos\Data\20190208_ComputarLens\dispersion\spectral\RAWDiskDispersionResults_spectral_spline_fromReference.mat';
    'C:\Users\GraphicsLab\Documents\llanos\Results\dispersion\spectral\polynomial_newCV\RAWDiskDispersionResults_spectral_polynomial_fromReference.mat'
};

% Names for the models of dispersion
names = {
    'Polynomial, spatial deg. 4, spectral deg. 6 (old cross-validation)';
    'Spline, smoothing parameter 1.1884';
    'Polynomial, spatial deg. 3, spectral deg. 3 (new cross-validation)'
};

% Output directory for all figures and saved parameters
output_directory = 'C:\Users\GraphicsLab\Documents\llanos\Results\dispersion_plots';

%% Load dispersion models

n_models = length(forward_dispersion_model_filenames);

dispersionfuns = cell(n_models, 1);
for i = 1:n_models
    [...
        dispersion_data, bands_dispersionfun, transform_data...
    ] = loadDispersionModel(forward_dispersion_model_filenames{i}, true, true);
    if i == 1
        bands = reshape(bands_dispersionfun, [], 1);
        n_bands = length(bands);
    end

    [~, T_roi] = modelSpaceTransform(...
        image_sampling, transform_data.model_space, transform_data.fill, false...
    );
    dispersionfuns{i} = makeDispersionfun(dispersion_data, T_roi);
end

%% Plotting

sampled_pixels = [
    linspace(0.5, floor(image_sampling(2) / 2) + 0.5, n_sampled_pixels).',...
    linspace(0.5, floor(image_sampling(1) / 2) + 0.5, n_sampled_pixels).'
];

legend_str = [
    {sprintf('%s at control bands', names{1})};
    names
];

x_limits = sqrt(plot_bound ^2 / 2);
x_limits = [-x_limits, x_limits];
y_limits = x_limits;

for k = 1:n_sampled_pixels
    px = sampled_pixels(k, :);
    
    fg = figure;
    hold on
    px_lambda = [repmat(px, n_bands, 1), bands];
    d = dispersionfuns{1}(px_lambda);
    plot3(d(:, 1), d(:, 2), px_lambda(:, 3), 'ko');
    
    px_lambda = [repmat(px, size(lambda, 1), 1), lambda];
    
    for i = 1:n_models
        d = dispersionfuns{i}(px_lambda);
        plot3(d(:, 1), d(:, 2), px_lambda(:, 3), 'LineWidth', 2);
    end
    hold off
    
    grid on
    legend(legend_str{:});
    title(sprintf('Spectral dispersion at pixel (%g, %g)', px(1), px(2)));
    xlabel('Dispersion x-component');
    ylabel('Dispersion y-component');
    zlabel('\lambda');
    xlim(x_limits);
    ylim(y_limits);
    
    savefig(...
        fg, fullfile(output_directory, sprintf(...
                'dispersionAtX%dY%d.fig', floor(px(1)), floor(px(2))...
        )), 'compact'...
    );
    close(fg)
end

%% Save parameters to a file
save_variables_list = [ parameters_list, {} ];
save_data_filename = fullfile(output_directory, 'SpectralDispersionParams.mat');
save(save_data_filename, save_variables_list{:});