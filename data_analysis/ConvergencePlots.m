%% ADMM convergence analysis
% Plot variables from the 'saveIterations*' files output by
% 'baek2017Algorithm2LowMemory()'.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 25, 2019

%% Input data and parameters

% Image estimation parameters and results '.mat' file, such as output by
% 'CorrectByHyperspectralADMM.m'
parameters_filename = 'CorrectByHyperspectralADMM.mat';

% Wildcard for 'ls()' to find the 'saveIterations*.mat' files to load.
% Only the last few files will be loaded, corresponding to the steps in
% multi-stage image estimation.
input_data_wildcard = 'saveIterations*.mat';

%% Initialization

load(parameters_filename);
if iscell(bands)
    n_steps = length(bands);
else
    n_steps = 1;
    bands = {bands};
end

enabled_weights = solvePatchesSpectralOptions.reg_options.enabled;
n_weights = length(enabled_weights);
if length(solvePatchesColorOptions.reg_options.enabled) ~= n_weights ||...
        any(solvePatchesColorOptions.reg_options.enabled ~= enabled_weights)
    error('It is not clear what the set of enabled regularization terms is.');
end
nonneg = solvePatchesSpectralOptions.admm_options.nonneg;
if solvePatchesColorOptions.admm_options.nonneg ~= nonneg
    error('It is not clear whether a non-negativity constraint is enabled.');
end
n_z = sum(enabled_weights) + nonneg;

legend_str_px = cell(n_steps, 1);

n_bands_max = length(bands{end});
plot_colors_bands = jet(n_bands_max);
%plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
%plot_styles = {'--', ':', '-.'};

%% Process the data files

data_filenames = listFiles(input_data_wildcard);
data_filenames = data_filenames((end - n_steps + 1):end);

iter_offset = 0;
legend_str_z = {};

fg_penalty_changes = figure;
title('Iterations where penalty parameters changed');
xlabel('Cumulative ADMM iteration number');
            
for s = 1:n_steps
    load(data_filenames{s});
    bands_s = bands{s};
    n_bands_s = length(bands_s);
    iter_cum = (1:iter) + iter_offset;
    
    if s == 1
        fg_px = figure;
        title('Center pixel value');
        xlabel('Cumulative ADMM iteration number');
        ylabel('Intensity');
    else
        figure(fg_px);
    end
    legend_str_px{s} = cell(n_bands_s, 1);
    hold on
    for c = 1:n_bands_s
        plot(...
            iter_cum, center_px_iter(1:iter, c),...
            'Color', plot_colors_bands(c, :),...
            'LineWidth', 2,...
            'Marker', 'none',...
            'LineStyle', '-'...
            );
        legend_str_px{s}{c} = sprintf('\\lambda = %g [nm]', bands_s(c));
    end
    hold off
    
    if s == 1
        fg_pcg_iter = figure;
        title('Number of conjugate gradients iterations');
        xlabel('Cumulative ADMM iteration number');
        ylabel('Number of conjugate gradients iterations');
    else
        figure(fg_pcg_iter);
    end
    hold on
    plot(...
        iter_cum, pcg_n_iter(1:iter),...
        'Color', [0, 0, 0],...
        'LineWidth', 2,...
        'Marker', 'none',...
        'LineStyle', '-'...
        );
    hold off
    
    if s == 1
        fg_pcg_relres = figure;
        title('Conjugate gradients relative residual');
        xlabel('Cumulative ADMM iteration number');
        ylabel('Conjugate gradients relative residual');
    else
        figure(fg_pcg_relres);
    end
    hold on
    plot(...
        iter_cum, pcg_relres_iter(1:iter),...
        'Color', [0, 0, 0],...
        'LineWidth', 2,...
        'Marker', 'none',...
        'LineStyle', '-'...
        );
    hold off
    
    z_ind_plot = 0;
    plot_colors_z = jet(n_z);
    for z_ind = 1:size(R_norm_iter, 2)
        if z_ind > n_weights || enabled_weights(z_ind)
            z_ind_plot = z_ind_plot + 1;
        end
        if ~any(R_norm_iter(:, z_ind))
            continue
        end
        
        if s == 1 || ~ishghandle(fg_R(z_ind_plot))
            fg_R(z_ind_plot) = figure; %#ok<SAGROW>
            title(sprintf('ADMM primal convergence for z_{%d}', z_ind_plot));
            xlabel('Cumulative ADMM iteration number');
        else
            figure(fg_R(z_ind_plot));
        end
        hold on
        plot(...
            iter_cum, R_norm_iter(1:iter, z_ind),...
            'Color', [1, 0, 0],...
            'LineWidth', 2,...
            'Marker', 'none',...
            'LineStyle', '-'...
            );
        plot(...
            iter_cum, epsilon_pri_iter(1:iter, z_ind),...
            'Color', [0, 1, 0],...
            'LineWidth', 2,...
            'Marker', 'none',...
            'LineStyle', ':'...
        );
        hold off
        
        if s == 1 || ~ishghandle(fg_S(z_ind_plot))
            fg_S(z_ind_plot) = figure; %#ok<SAGROW>
            title(sprintf('ADMM dual convergence for z_{%d}', z_ind_plot));
            xlabel('Cumulative ADMM iteration number');
        else
            figure(fg_S(z_ind_plot));
        end
        hold on
        plot(...
            iter_cum, S_norm_iter(1:iter, z_ind),...
            'Color', [1, 0, 0],...
            'LineWidth', 2,...
            'Marker', 'none',...
            'LineStyle', '-'...
            );
        plot(...
            iter_cum, epsilon_dual_iter(1:iter, z_ind),...
            'Color', [0, 1, 0],...
            'LineWidth', 2,...
            'Marker', 'none',...
            'LineStyle', ':'...
        );
        hold off
        
        figure(fg_penalty_changes);
        legend_str_z{end + 1} = sprintf('Change in \\rho_{%d}', z_ind_plot);
        hold on
        plot(...
            iter_cum,...
            changed_penalty_parameters(1:iter, z_ind) + 3 * (z_ind_plot - 1),...
            'Color', plot_colors_z(z_ind_plot, :),...
            'LineWidth', 2,...
            'Marker', 'none',...
            'LineStyle', '-'...
            );
        hold off
    end
    
    iter_offset = iter_offset + iter;
end

% Add legends to plots
figure(fg_px);
% Also add a zero line, to make visual inspection easier
hold on
plot(...
    [1 iter_offset], [0, 0],...
    'Color', [0, 0, 0],...
    'LineWidth', 1,...
    'Marker', 'none',...
    'LineStyle', '-'...
);
hold off
legend_str_px = vertcat(legend_str_px{:});
legend(legend_str_px{:});

z_ind_plot = 0;
for z_ind = 1:size(R_norm_iter, 2)
    if z_ind > n_weights || enabled_weights(z_ind)
        z_ind_plot = z_ind_plot + 1;
    end
    if ~any(R_norm_iter(:, z_ind))
        continue
    end
    
    figure(fg_R(z_ind_plot));
    legend(sprintf('norm(r_{%d})', z_ind_plot), sprintf('\\epsilon_{%d, pri}', z_ind_plot));
    
    figure(fg_S(z_ind_plot));
    legend(sprintf('norm(s_{%d})', z_ind_plot), sprintf('\\epsilon_{%d, dual}', z_ind_plot));
end

figure(fg_penalty_changes);
legend(legend_str_z{:});