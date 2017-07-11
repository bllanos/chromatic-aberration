function [ disparity_spline, disparity_raw ] = radialChromaticAberration(...
    stats, x_fields, reference_wavelength_index, z, reference_z, varargin...
)
% RADIALCHROMATICABERRATION  Model chromatic aberration
%
% ## Syntax
% disparity_spline = radialChromaticAberration(...
%     stats, x_fields, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb, verbose]...
% )
% [ disparity_spline, disparity ] = radialChromaticAberration(...
%     stats, x_fields, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb, verbose]...
% )
%
% ## Description
% disparity_spline = radialChromaticAberration(...
%     stats, x_fields, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb, verbose]...
% )
%   Returns a spline model of chromatic aberration as a function of
%   distance from the origin, and scene depth.
%
% [ disparity_spline, disparity ] = radialChromaticAberration(...
%     stats, x_fields, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb, verbose]...
% )
%   Additionally returns the disparity values used to build the spline
%   model.
%
% ## Input Arguments
%
% stats -- Point spread function statistics
%   Statistics of point spread functions produced by a set of scene
%   features, for each wavelength, and for each depth. `stats(i, k,
%   j).(name)` is the value of the 'name' statistic corresponding to the
%   i-th scene feature, emitting light the k-th wavelength, and positioned
%   at the j-th depth.
%
%   `stats.(name)` can be a scalars or vectors, as long as all elements are
%   the same size. The lengths of differences between vectors at the
%   reference wavelength, and vectors at other wavelengths, will be used as
%   the response variable values for spline fitting.
%
%   Consequently, this function can generate trends in a variety of
%   quantities with respect to scene depth, such as image positions
%   (2-element vectors), and point spread function blur radii (scalars).
%
%   `stats` might be an array of the `stats` output argument of
%   'analyzePSF()', for example.
%
% x_fields -- Spline fitting input field names
%   A scalar structure with the same fields as `stats`, where each value is
%   a two-element cell vector of strings. `x_fields.(name)` contains
%   the fieldname in `stats` of the values to be used as one of the input
%   variables for spline fitting for the values in `stats.(name)`.
%
%   The splines output in `disparity_spline.(name)` model disparities
%   between the values in `stats.(name)` for different wavelengths. The
%   splines are thin-plate splines in a 2D parameter space, where one
%   parameter is depth, and the other parameter is the magnitudes
%   (distances from the origin) of the values in
%   `stats(:, reference_wavelength_index, :).(x_fields.(name))`.
%
% reference_wavelength_index -- Reference wavelength
%   The index into the second dimension of `stats` representing the
%   reference wavelength. Chromatic aberration will be measured between
%   statistics produced with other wavelengths and statistics produced with
%   this wavelength.
%
% z -- Scene depths
%   The z-positions of the scene elements producing the values in `stats`.
%   `z` is a vector with a length equal to `size(stats, 3)`.
%
%   `z` can be some transformation of the actual z-positions. More
%   generally, `z` is some variable describing positions, measured along
%   the optical axis, which is to be assessed as a predictor of chromatic
%   aberration.
%
% reference_z -- Reference depth
%   The z-position corresponding to a depth value of zero; An offset which
%   will be applied to the values in `z`.
%
% wavelengths -- Wavelengths corresponding to image measurements
%   The wavelengths of light corresponding to the elements of `stats`. A
%   row vector of length `size(stats, 2)`, where `wavelengths(k)` is the wavelength
%   used to generate the values in `stats(:, k, :)`. This
%   parameter is used for figure legends only, not for calculations.
%
%   Either all of `wavelengths`, `wavelengths_to_rgb`, and `verbose` must
%   be passed, or none of them must be passed.
%
% wavelengths_to_rgb -- Colour map for wavelengths
%   RGB colours to be used when plotting points representing values for the
%   different wavelengths. The k-th row of this `size(stats, 2)` x 3
%   matrix represents the RGB colour corresponding to the k-th wavelength,
%   `wavelengths(k)`.
%
%   Either all of `wavelengths`, `wavelengths_to_rgb`, and `verbose` must
%   be passed, or none of them must be passed.
%
% verbose -- Debugging and visualization controls
%   If recognized fields of `verbose` are true, corresponding graphical
%   output will be generated for debugging purposes.
%
%   In addition to fields which are simple Boolean values, there is a
%   `filter` field, whose value is a structure. If `filter.(name)` is true,
%   then visualization output will be generated for the corresponding data
%   in `stats.(name)`. If `filter.(name)` is false, or if `name` is not a
%   field of `filter`, then the data in `stats.(name)` will be excluded
%   from all graphical output.
%
%   Either all of `wavelengths`, `wavelengths_to_rgb`, and `verbose` must
%   be passed, or none of them must be passed.
%
% ## Output Arguments
%
% disparity_spline -- Thin-plate spline models of chromatic aberration
%   A set of thin-plate smoothing splines modeling chromatic aberration as
%   a function of `z`, and `r`. `r` is the distance of a value in `stats(:,
%   reference_wavelength_index, :).(x_fields.(name))` from the
%   origin.
%
%   `disparity_spline.(name)` is a cell vector of length `size(stats, 2)`,
%   where the k-th cell contains a thin-plate smoothing spline describing
%   the aberration between the k-th wavelength and the reference
%   wavelength, for the statistic 'name'.
%   `disparity_spline(reference_wavelength_index).(name)` is an empty cell.
%
%   `aberration = fnval(disparity_spline(k).(name),[r; z])` evaluates the
%   spline model at the distances from the origin in the row vector `r`,
%   and the associated depths in the row vector `z`. `aberration` is a row
%   vector of distances between point spread function statistics of type
%   'name'.
%
%   If there are too few data points to estimate spline models, the
%   corresponding elements of `disparity_spline.(name)` are empty cells.
%
% disparity_raw -- Input disparity values
%   The disparity values, calculated from the values in `stats` and `z`,
%   used to construct `disparity_spline`. `disparity_raw` has the same
%   dimensions as `stats`; `disparity_raw(i, k, j)` is the displacement
%   vector from `stats(i, reference_wavelength_index, j).(name)` to
%   `stats(i, k, j).(name)`. Therefore, this disparity vector is measured
%   for the i-th scene feature, emitting light at the k-th wavelength, and
%   positioned at the j-th depth.
%
% See also doubleSphericalLensPSF, analyzePSF, tpaps

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 27, 2017

nargoutchk(1, 2);
narginchk(5, 8);

sz = size(stats);
n_points = sz(1);
n_wavelengths = sz(2);
n_depths = sz(3);

if ~isempty(varargin)
    if length(varargin) ~= 3
        error('Unexpected number of input arguments. All of `wavelengths`, `wavelengths_to_rgb`, and `verbose` should be passed, or none should be passed.');
    else
        wavelengths = varargin{1};
        wavelengths_to_rgb = varargin{2};
        n_points_all_depths = n_points * n_depths;
        zeros_plot = zeros(n_points_all_depths, 1);
        verbose = varargin{3};
        verbose_filter = verbose.filter;
        display_raw_values = verbose.display_raw_values;
        display_raw_disparity = verbose.display_raw_disparity;
        display_disparity_splines = verbose.display_disparity_splines;
    end
else
    verbose_filter = struct();
    display_raw_values = false;
    display_raw_disparity = false;
    display_disparity_splines = false;
end

names = fieldnames(stats);
n_names = length(names);

% Preliminary processing
stats_cell = struct2cell(stats);
for i = 1:n_names
    stats_cell_i = squeeze(stats_cell(i, :, :, :));
    stats_cell_i = reshape(stats_cell_i, n_points, 1, n_wavelengths, n_depths);
    stats_mat_i = cell2mat(stats_cell_i);
    stats_mat.(names{i}) = stats_mat_i;
end

z_adjusted = z - reference_z;
z_adjusted = repelem(z_adjusted, n_points);
if size(z_adjusted, 1) > size(z_adjusted, 2)
    z_adjusted = z_adjusted.';
end

for i = 1:n_names
    name_i = names{i};
    name_display_i = replace(name_i, '_', '\_');
    verbose_filter_i = isfield(verbose_filter, name_i) && verbose_filter.(name_i);
    
    % Find the other input variable
    name_x = x_fields.(name_i);
    name_display_x = replace(name_x, '_', '\_');

    stats_mat_i = stats_mat.(name_i);
    stats_reference_i = repmat(...
        stats_mat_i(:, :, reference_wavelength_index, :),...
        1, 1, n_wavelengths, 1 ...
    );
    disparity_raw_i = stats_mat_i - stats_reference_i;
    
    stats_mat_x = stats_mat.(name_x);
    
    dimensionality_i = size(disparity_raw_i, 2);
    dimensionality_x = size(stats_mat_x, 2);
    
    if verbose_filter_i && dimensionality_i <= 2 && dimensionality_x <= 2
        
        stats_mat_i_3D = reshape(permute(stats_mat_i, [1, 4, 2, 3]), [], dimensionality_i, n_wavelengths);
        stats_mat_x_3D = reshape(permute(stats_mat_x, [1, 4, 2, 3]), [], dimensionality_x, n_wavelengths);
        disparity_raw_i_3D = reshape(permute(disparity_raw_i, [1, 4, 2, 3]), [], dimensionality_i, n_wavelengths);
        
        if display_raw_values
            figure
            hold on
            legend_strings = cell(n_wavelengths, 1);
            for k = 1:n_wavelengths
                if dimensionality_x == 1
                    y_plot = zeros_plot;
                else
                    y_plot = stats_mat_x_3D(:, 2, k);
                end
                if dimensionality_i == 1
                    v_plot = zeros_plot;
                else
                    v_plot = stats_mat_i_3D(:, 2, k);
                end
                
                quiver3(...
                    stats_mat_x_3D(:, 1, reference_wavelength_index),...
                    y_plot,...
                    z_adjusted.',...
                    stats_mat_i_3D(:, 1, k),...
                    v_plot,...
                    zeros_plot,...
                    'Color', wavelengths_to_rgb(k, :), 'AutoScale', 'on'...
                    );
                legend_strings{k} = sprintf(...
                    '\\lambda = %g nm',...
                    wavelengths(k)...
                    );
            end
            legend(legend_strings);
            title(sprintf(...
                'Values of ''%s'' plotted as (autoscaled) arrows at values of ''%s''',...
                name_display_i, name_display_x...
            ))
            if dimensionality_x == 2
                xlabel(sprintf('%s_x', name_display_x));
                ylabel(sprintf('%s_y', name_display_x));
                zlabel('Depth')
            elseif dimensionality_x == 1
                xlabel(sprintf('%s', name_display_x));
                zlabel('Depth')
            end
            hold off
        end
        
        if display_raw_disparity
            figure
            hold on
            legend_strings = cell(n_wavelengths * 2 - 1, 1);
            for k = 1:n_wavelengths
                if dimensionality_x == 1
                    y_plot = zeros_plot;
                else
                    y_plot = stats_mat_x_3D(:, 2, k);
                end

                scatter3(...
                    stats_mat_x_3D(:, 1, k), y_plot,...
                    z_adjusted, [], wavelengths_to_rgb(k, :), 'o'...
                );
                legend_strings{k} = sprintf(...
                    '''%s'' for \\lambda = %g nm',...
                    name_display_x, wavelengths(k)...
                );
            end
            k_legend = 1;
            for k = 1:n_wavelengths
                if k ~= reference_wavelength_index
                    if dimensionality_x == 1
                        y_plot = zeros_plot;
                    else
                        y_plot = stats_mat_x_3D(:, 2, reference_wavelength_index);
                    end
                    if dimensionality_i == 1
                        v_plot = zeros_plot;
                    else
                        v_plot = disparity_raw_i_3D(:, 2, k);
                    end

                    quiver3(...
                        stats_mat_x_3D(:, 1, reference_wavelength_index),...
                        y_plot,...
                        z_adjusted.',...
                        disparity_raw_i_3D(:, 1, k),...
                        v_plot,...
                        zeros_plot,...
                        'Color', wavelengths_to_rgb(k, :), 'AutoScale', 'off'...
                    );
                    legend_strings{n_wavelengths + k_legend} = sprintf(...
                        'Aberration of ''%s'' for \\lambda = %g nm',...
                        name_display_i, wavelengths(k)...
                    );
                    k_legend = k_legend + 1;
                end
            end
            legend(legend_strings);
            title(sprintf(...
                'Aberrations of ''%s'' at values of ''%s''',...
                name_display_i, name_display_x...
            ))
            if dimensionality_x == 2
                xlabel(sprintf('%s_x', name_display_x));
                ylabel(sprintf('%s_y', name_display_x));
                zlabel('Depth')
            elseif dimensionality_x == 1
                xlabel(sprintf('%s', name_display_x));
                zlabel('Depth')
            end
            hold off
        end
    elseif verbose_filter_i
        warning('`radialChromaticAberration` cannot produce a visualization of statistics with more than two dimensions.');
    end

    disparity_raw_i_radial = sqrt(dot(disparity_raw_i, disparity_raw_i, 2));
    disparity_raw_i_radial_signs = sign(dot(stats_reference_i, disparity_raw_i, 2));
    disparity_raw_i_radial_signs(disparity_raw_i_radial_signs == 0) = 1;
    disparity_raw_i_radial = disparity_raw_i_radial .* disparity_raw_i_radial_signs;
    stats_reference_x_radial = sqrt(dot(...
        stats_mat_x(:, :, reference_wavelength_index, :),...
        stats_mat_x(:, :, reference_wavelength_index, :),...
        2 ...
    ));
    stats_reference_x_radial = stats_reference_x_radial(:).';

    % Prepare data for spline fitting
    spline_predictors = [ stats_reference_x_radial; z_adjusted ];
    spline_responses = permute(squeeze(disparity_raw_i_radial), [1, 3, 2]);
    spline_responses = reshape(spline_responses, 1, [], n_wavelengths);

    spline_predictors_filter = all(isfinite(spline_predictors), 1);
    spline_predictors = spline_predictors(:, spline_predictors_filter);
    spline_responses = spline_responses(:, spline_predictors_filter, :);

    % Spline fitting
    disparity_spline_i = cell(n_wavelengths, 1);
    for k = 1:n_wavelengths
        if k ~= reference_wavelength_index

            spline_responses_k = spline_responses(:, :, k);
            spline_responses_filter = isfinite(spline_responses_k);
            spline_responses_k = spline_responses_k(spline_responses_filter);
            spline_predictors_k = spline_predictors(:, spline_responses_filter);
            spline_predictors_unique = unique(spline_predictors_k(1, :));
            sufficient_data = (length(spline_predictors_unique) > 1);
            spline_predictors_unique = unique(spline_predictors_k(2, :));
            sufficient_data = sufficient_data & (length(spline_predictors_unique) > 1);
            sufficient_data = sufficient_data & (size(spline_predictors_k, 2) > 2);

            if sufficient_data
                try
                    disparity_spline_i{k} = tpaps(...
                        spline_predictors_k,...
                        spline_responses_k...
                    );
                catch ME
                    if (strcmp(ME.identifier,'SPLINES:TPAPS:collinearsites'))
                        sufficient_data = false;
                        warning('Spline fit failed. Most likely there are insufficient data points available to construct a spline model of ''%s'' as a function of ''%s'', and depth, for \\lambda at index %d.',...
                            name_i, name_x, k...
                        )
                    else
                        rethrow(ME)
                    end
                end
            else
                warning(...
                    'Insufficient data points available to construct a spline model of ''%s'' as a function of ''%s'', and depth, for \\lambda at index %d.',...
                    name_i, name_x, k...
                )
            end

            if verbose_filter_i && display_disparity_splines
                figure
                if sufficient_data
                    pts = fnplt(disparity_spline_i{k});
                    surf(pts{1}, pts{2}, pts{3}, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
                    colorbar
                    colormap summer
                    colorbar;
                    title(sprintf(...
                        'Aberration in ''%s'' for \\lambda = %g nm',...
                        name_display_i, wavelengths(k)...
                    ))
                    legend_str = {'Thin plate spline'};
                else
                    title(sprintf(...
                        'Aberration in ''%s'' for \\lambda = %g nm (Insufficient data for spline model)',...
                        name_display_i, wavelengths(k)...
                    ))
                    legend_str = cell(1, 0);
                end
                xlabel(sprintf('Magnitude of ''%s''', name_display_x));
                ylabel('Depth');
                zlabel(sprintf('Aberration in ''%s''', name_display_i));
                hold on
                plot3(...
                    spline_predictors_k(1, :), spline_predictors_k(2, :),...
                    spline_responses_k,...
                    'ko','markerfacecolor','r'...
                    )
                legend([legend_str, 'Original data points'])
                hold off
            end
        end
    end
    
    disparity_spline.(name_i) = disparity_spline_i;
    disparity_raw.(name_i) = disparity_raw_i;
end

end