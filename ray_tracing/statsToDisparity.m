function [...
    disparity_raw, disparity_raw_radial, stats_mat...
] = statsToDisparity(...
    stats, reference_wavelength_index, varargin...
)
% STATSTODISPARITY  Calculate raw chromatic aberration vectors
%
% ## Syntax
% disparity_raw = statsToDisparity(...
%     stats, reference_wavelength_index,...
%     [, z, reference_z, x_fields, wavelengths, wavelengths_to_rgb, verbose]...
% )
% [...
%     disparity_raw, disparity_raw_radial...
% ] = statsToDisparity(...
%     stats, reference_wavelength_index,...
%     [, z, reference_z, x_fields, wavelengths, wavelengths_to_rgb, verbose]...
% )
% [...
%     disparity_raw, disparity_raw_radial, stats_mat...
% ] = statsToDisparity(...
%     stats, reference_wavelength_index,...
%     [, z, reference_z, x_fields, wavelengths, wavelengths_to_rgb, verbose]...
% )
%
% ## Description
% disparity_raw = statsToDisparity(...
%     stats, reference_wavelength_index,...
%     [, z, reference_z, x_fields, wavelengths, wavelengths_to_rgb, verbose]...
% )
%   Returns disparity vectors relating point spread function statistics
%   between wavelengths, as a function of scene depth.
%
% [...
%     disparity_raw, disparity_raw_radial...
% ] = statsToDisparity(...
%     stats, reference_wavelength_index,...
%     [, z, reference_z, x_fields, wavelengths, wavelengths_to_rgb, verbose]...
% )
%   Additionally returns signed magnitudes of disparity vectors.
%
% [...
%     disparity_raw, disparity_raw_radial, stats_mat...
% ] = statsToDisparity(...
%     stats, reference_wavelength_index,...
%     [, z, reference_z, x_fields, wavelengths, wavelengths_to_rgb, verbose]...
% )
%   Additionally returns the input `stats` argument in a different layout,
%   for easier manipulation.
%
% ## Input Arguments
%
% In the following, the term "colour channel" can be substituted for
% "wavelength". This function is intended to process data for "wavebands",
% which can be individual wavelengths, colour channels, or more general
% spectral distributions.
%
% Either all of optional input arguments must be passed, or none of them
% must be passed. The optional input arguments are used only to generate
% graphical output.
%
% stats -- Point spread function statistics
%   Statistics of point spread functions produced by a set of scene
%   features, for each wavelength, and for each depth. `stats(i, k,
%   j).(name)` is the value of the 'name' statistic corresponding to the
%   i-th scene feature, emitting light of the k-th wavelength, and
%   positioned at the j-th depth.
%
%   `stats.(name)` can be a scalars or vectors, as long as all elements are
%   the same size. `stats` might be an array of the `stats` output argument
%   of 'analyzePSF()', for example.
%
% reference_wavelength_index -- Reference wavelength
%   The index into the second dimension of `stats` representing the
%   reference wavelength. Disparity vectors will be measured between
%   statistics for other wavelengths and statistics for this wavelength.
%
% z -- Scene depths
%   The z-positions of the scene elements producing the values in `stats`.
%   `z` is a vector with a length equal to `size(stats, 3)`.
%
%   `z` can be some transformation of the actual z-positions. More
%   generally, `z` is some variable describing positions, measured along
%   the optical axis.
%
% reference_z -- Reference depth
%   The z-position corresponding to a depth value of zero; An offset which
%   will be applied to the values in `z`.
%
% x_fields -- Spline fitting input field names
%   A scalar structure with the same fields as `stats`, where each value is
%   a two-element cell vector of strings. `x_fields.(name)` contains the
%   fieldname in `stats` of the independent variables to use when plotting
%   the disparity vectors calculated from `stats.(name)`.
%
% wavelengths -- Wavelengths or colours corresponding to image measurements
%   Either a vector of wavelengths of light, or a cell vector of colour
%   channel names (character vectors), corresponding to the elements of
%   `stats`. A row vector, or cell vector, of length `size(stats, 2)`,
%   where `wavelengths(k)` is associated with the values in `stats(:, k,
%   :)`.
%
% wavelengths_to_rgb -- Colour map for wavelengths
%   RGB colours to be used when plotting points representing values for the
%   different wavelengths or colour channels. The k-th row of this
%   `size(stats, 2)` x 3 matrix represents the RGB colour corresponding to
%   the k-th wavelength or colour channel, `wavelengths(k)`.
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
% ## Output Arguments
%
% disparity_raw -- Raw disparity vectors
%   Disparity vectors, calculated from the values in `stats`.
%   `disparity_raw.(name)` has the same dimensions as `stats`;
%   `disparity_raw.(name)(i, k, j)` is the displacement vector from
%   `stats(i, reference_wavelength_index, j).(name)` to `stats(i, k,
%   j).(name)`. Therefore, this disparity vector is measured for the i-th
%   scene feature, emitting light at the k-th wavelength, and positioned at
%   the j-th depth.
%
% disparity_raw_radial -- Raw disparity vector signed magnitudes
%   The magnitudes of the vectors in `disparity_raw`, with signs indicating
%   if they are aligned with the corresponding vectors in `stats(:,
%   reference_wavelength_index, :).(name)` (positive), or in the opposite
%   direction (negative).
%
% stats_mat -- Reformatted point spread function statistics
%   A scalar structure version of `stats`, where each field stores an array
%   of values. In contrast, `stats` is a structure array, where each field
%   of each element stores a single scalar/vector value. The format of
%   `stats_mat` is the same as the format of `disparity_raw` and
%   `disparity_raw_radial`.
%
% See also doubleSphericalLensPSF, analyzePSF, radialChromaticAberration,
% tpaps

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 26, 2017

nargoutchk(1, 3);
narginchk(2, 8);

sz = size(stats);
n_points = sz(1);
n_wavelengths = sz(2);
if length(sz) < 3
    n_depths = 1;
else
    n_depths = sz(3);
end

if ~isempty(varargin)
    if length(varargin) ~= 6
        error('Unexpected number of input arguments. All of the optional should be passed, or none should be passed.');
    else
        z = varargin{1};
        reference_z = varargin{2};
        x_fields = varargin{3};
        wavelengths = varargin{4};
        have_wavelength_names = iscell(wavelengths);
        wavelengths_to_rgb = varargin{5};
        n_points_all_depths = n_points * n_depths;
        zeros_plot = zeros(n_points_all_depths, 1);
        verbose = varargin{6};
        verbose_filter = verbose.filter;
        display_raw_values = verbose.display_raw_values;
        display_raw_disparity = verbose.display_raw_disparity;
        
        z_adjusted = z - reference_z;
        z_adjusted = repelem(z_adjusted, n_points);
        if size(z_adjusted, 1) > size(z_adjusted, 2)
            z_adjusted = z_adjusted.';
        end
    end
else
    verbose_filter = struct();
    display_raw_values = false;
    display_raw_disparity = false;
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

for i = 1:n_names
    name_i = names{i};
    name_display_i = replace(name_i, '_', '\_');
    verbose_filter_i = isfield(verbose_filter, name_i) && verbose_filter.(name_i);
    
    if verbose_filter_i
        % Find the other input variable
        name_x = x_fields.(name_i);
        name_display_x = replace(name_x, '_', '\_');
        stats_mat_x = stats_mat.(name_x);
        dimensionality_x = size(stats_mat_x, 2);
    end

    stats_mat_i = stats_mat.(name_i);
    stats_reference_i = repmat(...
        stats_mat_i(:, :, reference_wavelength_index, :),...
        1, 1, n_wavelengths, 1 ...
    );
    disparity_raw_i = stats_mat_i - stats_reference_i;
    disparity_raw_i_radial = sqrt(dot(disparity_raw_i, disparity_raw_i, 2));
    disparity_raw_i_radial_signs = sign(dot(stats_reference_i, disparity_raw_i, 2));
    disparity_raw_i_radial_signs(disparity_raw_i_radial_signs == 0) = 1;
    disparity_raw_i_radial = disparity_raw_i_radial .* disparity_raw_i_radial_signs;
    
    dimensionality_i = size(disparity_raw_i, 2);
    
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
                if have_wavelength_names
                    legend_strings{k} = sprintf(...
                        '%s channel',...
                        wavelengths{k}...
                        );
                else
                    legend_strings{k} = sprintf(...
                        '\\lambda = %g nm',...
                        wavelengths(k)...
                        );
                end
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
                
                if have_wavelength_names
                    legend_strings{k} = sprintf(...
                        '''%s'' for %s channel',...
                        name_display_x, wavelengths{k}...
                        );
                else
                    legend_strings{k} = sprintf(...
                        '''%s'' for \\lambda = %g nm',...
                        name_display_x, wavelengths(k)...
                        );
                end
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
                    if have_wavelength_names
                        legend_strings{n_wavelengths + k_legend} = sprintf(...
                            'Aberration of ''%s'' for %s channel',...
                            name_display_i, wavelengths{k}...
                        );
                    else
                        legend_strings{n_wavelengths + k_legend} = sprintf(...
                            'Aberration of ''%s'' for \\lambda = %g nm',...
                            name_display_i, wavelengths(k)...
                        );
                    end
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
        warning('`statsToDisparity` cannot produce a visualization of statistics with more than two dimensions.');
    end

    disparity_raw.(name_i) = disparity_raw_i;
    disparity_raw_radial.(name_i) = disparity_raw_i_radial;
end

end