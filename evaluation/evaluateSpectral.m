function [e_spectral, varargout] = evaluateSpectral(...
    I_spectral, R_spectral, lambda, spectral_weights, options...
)
% EVALUATESPECTRAL  Compare spectral images
%
% ## Syntax
% e_spectral = evaluateSpectral(...
%     I_spectral, R_spectral, lambda, spectral_weights, options...
% )
% [e_spectral, fg_spectral] = evaluate(...
%     I_spectral, R_spectral, lambda, spectral_weights, options...
% )
%
% ## Description
%
% e_spectral = evaluateSpectral(...
%     I_spectral, R_spectral, lambda, spectral_weights, options...
% )
%   Returns a structure containing the values of different evaluation
%   metrics comparing the two spectral images, and generates corresponding
%   graphical output, if requested.
%
% [e_spectral, fg_spectral] = evaluate(...
%     I_spectral, R_spectral, lambda, spectral_weights, options...
% )
%   Additionally returns a structure of figure handles from the spectral
%   image evaluation graphical output.
%
% ## Input Arguments
%
% I_spectral -- Test spectral image
%   An h x w x c2 array containing an estimated spectral image.
%
% R_spectral -- Reference spectral image
%   An h x w x c1 array containing the ideal/true spectral image.
%
% lambda -- Wavelengths
%   A vector of length 'c1' containing the wavelengths corresponding to the
%   third dimension of `R_spectral`.
%
% spectral_weights -- Spectral resampling matrix
%   A matrix of dimensions c1 x c2 defining a mapping from the spectral
%   space of `I_spectral` to the spectral space of `R_spectral`.
%
% options -- Output options for spectral images
%   A structure which controls graphical output to figures relating to
%   the spectral images, and also determines some of the output in
%   `e_spectral`. `options` has the following fields:
%   - 'metric': The preferred evaluation metric. 'metric' is a character
%     vector with the following possible values:
%     - 'rmse': Root mean squared error
%     - 'mrae': Mean relative absolute error
%     - 'gof': Goodness-of-fit
%   - 'error_map': If `true`, then this function will produce a 'metric'
%     error map figure.
%   - 'radiance': A matrix, where each row is a four-element vector
%     describing an image patch (center pixel x-coordinate, center pixel
%     y-coordinate, width, and height). Image patches widths and heights
%     must be odd integers. For each image patch, a figure will be
%     generated containing plots of spectral information from the two
%     images, averaged within the image patch.
%   - 'radiance_fg': (Optional) A vector of figure handles. The output
%     produced for 'radiance' will be added to these figures, instead of
%     being shown in new figures. It is assumed that the existing figures
%     already have plotlines for the reference image, so only the test
%     image's information will be added.
%   - 'scanlines': A matrix, where each row is a four-element vector
%     describing a line segment (start pixel x-coordinate, start pixel
%     y-coordinate, end pixel x-coordinate, and end pixel y-coordinate).
%     For each line segment, a figure will be generated, containing a plot
%     of 'metric' error along the line.
%   - 'scanlines_fg': A vector of figure handles. 'scanlines_fg' has the
%     same role with respect to 'scanlines' as 'radiance_fg' has for
%     'radiance'.
%   - 'mi_bands': A two-element vector containing the indices of the
%     spectral bands (in the spectral space of the reference image) to
%     compute mutual information between.
%   - 'bands_diff': A two-element vector containing the indices of the
%     spectral bands (in the spectral space of the reference image) to
%     create difference images between (to show dispersion between
%     wavelengths). The band with the first index will be subtracted from
%     the band with the second index. Two figures will be created, one for
%     the reference image, and one for the test image. If this field does
%     not exist, then the figures will not be created.
%   - 'plot_color': A 3-element vector containing the RGB triplet to use
%     for plotlines for the test image. The plotlines for the reference
%     image will always be black.
%   - 'plot_marker': A character vector describing the datapoint markers to
%     use for plotlines for the test image. The plotlines for the reference
%     image will not have markers. 'plot_marker' is used for the 'Marker'
%     line property.
%   - 'plot_style': A character vector describing the line style to use for
%     plotlines for the test image. The plotlines for the reference image
%     will have solid lines. 'plot_style' is used for the 'LineStyle' line
%     property.
%   - 'reference_patch': An optional field containing a four-element vector
%     describing an image patch (center pixel x-coordinate, center pixel
%     y-coordinate, width, and height). Image patches widths and heights
%     must be odd integers. If this field exists, the average spectral
%     signal will be computed within this patch, for each of `I_spectral`
%     and `R_spectral`. `I_spectral` will be divided by its average
%     spectral signal, then multiplied by the average spectral signal of
%     `R_spectral`. 'reference_patch' is intended to allow for comparisions
%     between spectral images that differ in illumination.
%
% ## Output Arguments
%
% e_spectral -- Spectral error statistics
%   A structure with the following fields:
%   - 'psnr': The peak signal-to-noise ratio between the two images. 'psnr'
%     is a structure with the following fields:
%     - 'min': The peak signal-to-noise ratio for the band with the lowest
%       peak signal-to-noise ratio.
%     - 'mean': The average peak signal-to-noise ratio over all bands
%     - 'median': The peak signal-to-noise ratio for the band with the
%       median peak signal-to-noise ratio.
%     - 'raw': A c1-element vector, containing the peak signal-to-noise
%       ratios for each band.
%   - 'ssim': The Structural Similarity Index Measure computed between the
%     two images. 'ssim' is a structure with the following fields:
%     - 'min': The SSIM value for the band with the lowest SSIM value.
%     - 'mean': The average SSIM value over all bands
%     - 'median': The SSIM value for the band with the median mean SSIM value.
%     - 'raw': A c1-element vector, containing the SSIM values for each band.
%   - 'mi_within': A two element vector, containing the mutual information
%     between two bands in the reference image, and between two bands in
%     the test image, respectively. The indices of the bands are given in
%     `options.mi_bands`.
%   - 'mi_between': A c1-element vector, containing the mutual information
%     between corresponding bands in the reference image and the test
%     image.
%   - 'rmse': The root mean square error between the spectral information
%     at pixels in the test image and pixels in the reference image. 'rmse'
%     is a structure with the following fields:
%     - 'max': The highest RMSE value in the image
%     - 'mean': The average RMSE value over all pixels
%     - 'median': The median RMSE value over all pixels
%     - 'global': The RMSE value computed over all pixels, in contrast to
%       'mean', which is the average of per-pixel RMSE values.
%   - 'mrae': The mean relative absolute error between the spectral
%     information of pixels in the test image and in the reference image.
%     'mrae' is a structure with the following fields:
%     - 'max': The highest mean relative absolute error in the image
%     - 'mean': The average mean relative absolute error value over all
%       pixels. Note that this value is the same as would be computed
%       globally over all pixels, so there is no need for a 'global' field.
%     - 'median': The median mean relative absolute error over all pixels
%   - 'gof': The goodness-of-fit coefficient between the spectral
%     information of pixels in the test image and in the reference image.
%     'gof' is a structure with the following fields:
%     - 'min': The lowest goodness-of-fit value in the image
%     - 'mean': The average goodness-of-fit value over all pixels. This is
%       the same value which should be used as a 'global' value.
%     - 'median': The median goodness-of-fit value over all pixels
%   - 'radiance': A structure vector, with the same length as the number of
%     rows in `options.radiance`. Each element contains metrics describing
%     the image patch given by the corresponding element of
%     `options.radiance`, and has the following fields:
%     - 'rmse': The root mean square error of the average spectral radiance
%       in the test patch relative to the average spectral radiance in the
%       reference patch
%     - 'mrae': The mean relative absolute error of the average spectral
%       radiance in the test patch relative to the average spectral
%       radiance in the reference patch
%     - 'gof': The goodness-of-fit coefficient of the average spectral
%       radiance in the test patch relative to the average spectral
%       radiance in the reference patch.
%
% fg_spectral -- Spectral error evaluation figures
%   A structure with the following fields, all of which store figure
%   handles:
%   - 'error_map': A figure handle corresponding to the output triggered by
%     `options.error_map`.
%   - 'radiance': A vector of figure handles corresponding to the output
%     triggered by `options.radiance`.
%   - 'patches': A figure showing the image locations of all patches
%     described by `options.radiance`, produced only if `options` does not
%     have a 'radiance_fg' field.
%   - 'scanlines': A vector corresponding to the output triggered by
%     `options.scanlines`.
%   - 'scanlines_locations': A figure showing the image locations of all
%     line segments described by `options.scanlines`, produced only if
%     `options` does not have a 'scanlines_fg' field.
%   - 'bands_diff': A two-element vector, containing figure handles for the
%     difference images between spectral bands. The first element is for
%     the reference image, whereas the second element is for the test
%     image. The spectral bands being subtracted are given by
%     `options.bands_diff`.
%
% ## Notes
% - A border of `border` pixels (a local variable defined in the code) is
%   excluded from the images when calculating image-wide statistics, to
%   ignore artifacts in image estimation arising from the image borders.
% - Figures will not be generated if the corresponding fields of
%   `options` are absent.
% - Figures are produced with titles and axis labels, but without legends.
% - Images can be input in integer or floating-point formats. For peak
%   signal-to-noise ratio calculations, the peak value will be 1.0 for
%   floating-point formats, and the maximum possible positive integer value
%   for integer formats.
%
% ## References
% - Image borders are excluded from image similarity measurements in the
%   image demosaicking literature, such as in:
%
%   Monno, Y., Kiku, D., Tanaka, M., & Okutomi, M. (2017). "Adaptive
%     residual interpolation for color and multispectral image
%     demosaicking." Sensors (Switzerland), 17(12). doi:10.3390/s17122787
%
% - The code for calculating mutual information was retrieved from MATLAB
%   Central,
%   https://www.mathworks.com/matlabcentral/fileexchange/36538-very-fast-mutual-information-betweentwo-images
%   The function 'third_party/MI_GG/MI_GG.m' was written by Generoso
%   Giangregorio, and corresponds to the following article:
%
%   M. Ceccarelli, M. di Bisceglie, C. Galdi, G. Giangregorio, S.L. Ullo,
%     "Image Registration Using Non–Linear Diffusion", IGARSS 2008.
%
% - The idea of using mutual information to evaluate image alignment is
%   mentioned in, among other articles,
%
%   Brauers, J., Schulte, B., & Aach, T. (2008). "Multispectral
%     Filter-Wheel Cameras: Geometric Distortion Model and Compensation
%     Algorithms." IEEE Transactions on Image Processing, 17(12),
%     2368-2380. doi:10.1109/TIP.2008.2006605
%
% - The goodness-of-fit coefficient is used as a spectral error measure
%   (Equation 15) in:
%
%   Nguyen R.M.H., Prasad D.K., Brown M.S. (2014) Training-Based Spectral
%     Reconstruction from a Single RGB Image. In: Fleet D., Pajdla T.,
%     Schiele B., Tuytelaars T. (eds) Computer Vision – ECCV 2014. ECCV
%     2014. Lecture Notes in Computer Science, vol 8695. Springer, Cham
%
% - Root mean squared error and mean relative absolute error are defined
%   in:
%
%   B. Arad, O. Ben-Shahar and R. Timofte, "NTIRE 2018 challenge on
%     spectral reconstruction from RGB images," in The IEEE Conference on
%     Computer Vision and Pattern Recognition (CVPR) Workshops, 2018.
%
% See also metrics, psnr, ssim, evaluateRGB, samplingWeights, plot

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2018

narginchk(5, 5);
nargoutchk(1, 2);

border = 10;

n_bands = length(lambda);
if n_bands ~= size(spectral_weights, 1)
    error('The number of bands in `lambda` must equal the size of `spectral_weights` in its first dimension.');
end
if size(spectral_weights, 1) ~= size(R_spectral, 3)
    error('The number of rows of `spectral_weights` must equal the size of `R_spectral` in its third dimension.');
end
if size(spectral_weights, 2) ~= size(I_spectral, 3)
    error('The number of columns of `spectral_weights` must equal the size of `I_spectral` in its third dimension.');
end

class_spectral = class(I_spectral);
if ~isa(R_spectral, class_spectral)
    error('The two images do not have the same datatype.')
end
if isinteger(I_spectral)
    peak_spectral = intmax(class_spectral);
elseif isfloat(I_spectral)
    peak_spectral = 1;
else
    error('The two images are of an unexpected datatype.')
end

if strcmp(options.metric, 'rmse')
    metric_ind = 1;
elseif strcmp(options.metric, 'mrae')
    metric_ind = 2;
elseif strcmp(options.metric, 'gof')
    metric_ind = 3;
else
    error('Unrecognized value, "%s", of `options.metric`.', options.metric);
end

I_spectral = channelConversion(I_spectral, spectral_weights);

image_height = size(I_spectral, 1);
image_width = size(I_spectral, 2);
if isfield(options, 'reference_patch')
    bounds = options.reference_patch;
    if mod(bounds(3), 2) == 0 || mod(bounds(4), 2) == 0
        error('Reference patch dimensions must be odd integers.');
    end
    widthDiv2 = (bounds(3) - 1) / 2;
    heightDiv2 = (bounds(4) - 1) / 2;
    rectangle = [...
        max(1, bounds(1) - widthDiv2);
        max(1, bounds(2) - heightDiv2);
        min(image_width, bounds(1) + widthDiv2);
        min(image_height, bounds(2) + heightDiv2)
    ];

    radiance_average_I = mean(mean(I_spectral(...
        rectangle(2):rectangle(4), rectangle(1):rectangle(3), :...
    ), 1), 2);
    radiance_average_R = squeeze(mean(mean(R_spectral(...
        rectangle(2):rectangle(4), rectangle(1):rectangle(3), :...
    ), 1), 2));
    radiance_ratio = radiance_average_R ./ radiance_average_I;
    I_spectral = I_spectral .* repmat(radiance_ratio, image_height, image_width, 1);
end

I_clipped = I_spectral((border + 1):(end - border), (border + 1):(end - border), :);
R_clipped = R_spectral((border + 1):(end - border), (border + 1):(end - border), :);

e_spectral.psnr = struct('raw', zeros(n_bands, 1));
e_spectral.ssim = struct('raw', zeros(n_bands, 1));
for c = 1:n_bands
    [~, ~, e_spectral.psnr.raw(c)] = metrics(I_clipped(:, :, c), R_clipped(:, :, c), 3, peak_spectral, true);
    e_spectral.ssim.raw(c) = ssim(I_clipped(:, :, c), R_clipped(:, :, c));
end
e_spectral.psnr.min = min(e_spectral.psnr.raw);
e_spectral.psnr.mean = mean(e_spectral.psnr.raw);
e_spectral.psnr.median = median(e_spectral.psnr.raw);
e_spectral.ssim.min = min(e_spectral.ssim.raw);
e_spectral.ssim.mean = mean(e_spectral.ssim.raw);
e_spectral.ssim.median = median(e_spectral.ssim.raw);

mi_class = 'uint8';
if ~isa(I_spectral, mi_class)
    peak_int = double(intmax(mi_class));
    I_clipped_int = uint8(I_clipped * peak_int / peak_spectral);
    R_clipped_int = uint8(R_clipped * peak_int / peak_spectral);
else
    I_clipped_int = I_clipped;
    R_clipped_int = R_clipped;
end

e_spectral.mi_within(1) = MI_GG(...
    R_clipped_int(:, :, options.mi_bands(1)),...
    R_clipped_int(:, :, options.mi_bands(2))...
);
e_spectral.mi_within(2) = MI_GG(...
    I_clipped_int(:, :, options.mi_bands(1)),...
    I_clipped_int(:, :, options.mi_bands(2))...
);

e_spectral.mi_between = zeros(n_bands, 1);
for c = 1:n_bands
    e_spectral.mi_between(c) = MI_GG(...
        I_clipped_int(:, :, c),...
        R_clipped_int(:, :, c)...
    );
end

[rmse, mrae, gof, global_rmse] = metrics(I_clipped, R_clipped, 3, 0, true);
e_spectral.rmse = struct(...
    'max', max(rmse(:)),...
	'mean', mean(rmse(:)),...
	'median', median(rmse(:)),...
    'global', global_rmse...
);
e_spectral.mrae = struct(...
    'max', max(mrae(:)),...
	'mean', mean(mrae(:)),...
	'median', median(mrae(:))...
);
e_spectral.gof = struct(...
    'min', min(gof(:)),...
	'mean', mean(gof(:)),...
	'median', median(gof(:))...
);

if isfield(options, 'radiance')
    n_patches = size(options.radiance, 1);
    e_spectral.radiance = struct(...
        'rmse', cell(n_patches, 1),...
        'mrae', cell(n_patches, 1),...
        'gof', cell(n_patches, 1)...
    );
    radiance_averages_I = zeros(n_patches, n_bands);
    radiance_averages_R = zeros(n_patches, n_bands);
    patch_rectangles = zeros(n_patches, 4);
    for i = 1:n_patches
        bounds = options.radiance(i, :);
        if mod(bounds(3), 2) == 0 || mod(bounds(4), 2) == 0
            error('Patch sizes must be odd integers.');
        end
        widthDiv2 = (bounds(3) - 1) / 2;
        heightDiv2 = (bounds(4) - 1) / 2;
        patch_rectangles(i, 1) = max(1, bounds(1) - widthDiv2);
        patch_rectangles(i, 2) = max(1, bounds(2) - heightDiv2);
        patch_rectangles(i, 3) = min(image_width, bounds(1) + widthDiv2);
        patch_rectangles(i, 4) = min(image_height, bounds(2) + heightDiv2);

        radiance_average_I = squeeze(mean(mean(I_spectral(...
            patch_rectangles(i, 2):patch_rectangles(i, 4),...
            patch_rectangles(i, 1):patch_rectangles(i, 3), :...
            ), 1), 2)).';

        radiance_average_R = squeeze(mean(mean(R_spectral(...
            patch_rectangles(i, 2):patch_rectangles(i, 4),...
            patch_rectangles(i, 1):patch_rectangles(i, 3), :...
            ), 1), 2)).';
        
        patch_rectangles(i, 3) = patch_rectangles(i, 3) - patch_rectangles(i, 1) + 1;
        patch_rectangles(i, 4) = patch_rectangles(i, 4) - patch_rectangles(i, 2) + 1;
        
        radiance_averages_I(i, :) = radiance_average_I;
        radiance_averages_R(i, :) = radiance_average_R;

        [...
            e_spectral.radiance(i).rmse,...
            e_spectral.radiance(i).mrae,...
            e_spectral.radiance(i).gof...
        ] = metrics(radiance_average_I, radiance_average_R, 2, 0, false);
    end
end

% Graphical output
fg_spectral = struct;

if isfield(options, 'error_map') && options.error_map
    fg_spectral.error_map = figure;
    switch metric_ind
        case 1
            map = metrics(I_spectral, R_spectral, 3, 0, false);
            title('Spectral root-mean-square error');
        case 2
            [~, map] = metrics(I_spectral, R_spectral, 3, 0, false);
            title('Spectral mean relative absolute error');
        case 3
            [~, ~, map] = metrics(I_spectral, R_spectral, 3, 0, false);
            title('Spectral goodness-of-fit');
    end
    if metric_ind == 1
        imagesc(map);
    else
        imagesc(min(max(map, 0), 1));
    end
    colorbar;
    if metric_ind ~= 1
        caxis([0 1]);
    end 
end

background = [];
if isfield(options, 'radiance')
    for i = 1:n_patches
        if isfield(options, 'radiance_fg')
            fg_spectral.radiance(i) = options.radiance_fg(i);
            figure(fg_spectral.radiance(i));
            hold on
        else
            fg_spectral.radiance(i) = figure;
            hold on
            plot(...
                lambda, radiance_averages_R(i, :),...
                'Color', 'k', 'LineWidth', 2, 'Marker', 'none'...
            );
            xlabel('Wavelength [nm]');
            ylabel('Radiance');
            title(sprintf('Spectral radiance for image patch %d', i));
        end
        plot(...
                lambda, radiance_averages_I(i, :),...
                'Color', options.plot_color, 'LineWidth', 2,...
                'Marker', options.plot_marker, 'LineStyle', options.plot_style...
            );
        hold off
    end
    
    % Generate a figure showing all patches
    if ~isfield(options, 'radiance_fg')
        fg_spectral.patches = figure;
        % Show the band with highest average intensity
        R_spectral_mean = mean(mean(R_spectral, 1), 2);
        [~, max_ind] = max(R_spectral_mean);
        background = R_spectral(:, :, max_ind);
        background = background ./ max(max(background));
        labels = cell(n_patches, 1);
        for i = 1:n_patches
            labels{i} = sprintf('%d', i);
        end
        figure_image = insertObjectAnnotation(...
            background, 'rectangle', patch_rectangles, labels,...
            'TextBoxOpacity', 0.9, 'FontSize', 12, 'LineWidth', 2,...
            'Color', jet(n_patches)...
        );
        imshow(figure_image);
        title('Image patches for spectral radiance evaluation')
    end
end

if isfield(options, 'scanlines')
    n_lines = size(options.scanlines, 1);
    
    % Find pixels along the lines
    lines_ordering_filter = options.scanlines(:, 1) > options.scanlines(:, 3);
    lines_endpoints = options.scanlines;
    lines_endpoints(lines_ordering_filter, :) = options.scanlines(lines_ordering_filter, [3, 4, 1, 2]);
    lines_pixels_indices = cell(n_lines, 1);
    line_pixels_x = cell(n_lines, 1);
    for i = 1:n_lines
        line_pixels_x{i} = (lines_endpoints(i, 1):lines_endpoints(i, 3)).';
        n_px_line = length(line_pixels_x{i});
        line_pixels_xy = round(pointsOnLinesByCoordinates(...
            repmat(lines_endpoints(i, 1:2), n_px_line, 1),...
            repmat(lines_endpoints(i, 3:4), n_px_line, 1),...
            line_pixels_x{i}, ones(n_px_line, 1)...
        ));
        line_pixels_xy(line_pixels_xy < 1) = 1;
        filter = line_pixels_xy(:, 1) > image_width;
        line_pixels_xy(filter, 1) = image_width;
        filter = line_pixels_xy(:, 2) > image_height;
        line_pixels_xy(filter, 2) = image_height;
        lines_pixels_indices{i} = sub2ind(...
            [image_height, image_width],...
            line_pixels_xy(:, 2), line_pixels_xy(:, 1)...
        );
    end
    
    % Plot spectral radiance error
    for i = 1:n_lines
        n_px_line = length(lines_pixels_indices{i});
        lines_pixels_indices_spectral =...
            repmat(lines_pixels_indices{i}, n_bands, 1) +...
            repelem((0:(n_bands - 1)).' * (image_width * image_height), n_px_line, 1);
        radiance_R = reshape(R_spectral(lines_pixels_indices_spectral), n_px_line, n_bands);

        lines_pixels_indices_spectral =...
            repmat(lines_pixels_indices{i}, size(I_spectral, 3), 1) +...
            repelem((0:(size(I_spectral, 3) - 1)).' * (image_width * image_height), n_px_line, 1);
        radiance_I = reshape(I_spectral(lines_pixels_indices_spectral), n_px_line, []);

        
        switch metric_ind
            case 1
                line_error = metrics(radiance_I, radiance_R, 2, 0, false);
            case 2
                [~, line_error] = metrics(radiance_I, radiance_R, 2, 0, false);
            case 3
                [~, ~, line_error] = metrics(radiance_I, radiance_R, 2, 0, false);
            otherwise
                error('Unrecognized value of `metric_ind`.');
        end
        
        if isfield(options, 'scanlines_fg')
            fg_spectral.scanlines(i) = options.scanlines_fg(i);
            figure(fg_spectral.scanlines(i));
        else
            fg_spectral.scanlines(i) = figure;
            xlabel('Image x-coordinate');
            switch metric_ind
                case 1
                    ylabel('Spectral root-mean-square error');
                case 2
                    ylabel('Spectral mean relative absolute error');
                case 3
                    ylabel('Spectral goodness-of-fit');
                otherwise
                    error('Unrecognized value of `metric_ind`.');
            end
            title(sprintf('Spectral radiance error along line %d', i));
        end
        hold on
        plot(...
                line_pixels_x{i}, line_error,...
                'Color', options.plot_color, 'LineWidth', 2,...
                'Marker', options.plot_marker, 'LineStyle', options.plot_style...
            );
        hold off
        if metric_ind == 3
            ylim([0, 1]);
        end
    end
    
    % Generate a figure showing all lines
    if ~isfield(options, 'scanlines_fg')
        fg_spectral.scanlines_locations = figure;
        % Show the band with highest average intensity
        if isempty(background)
            R_spectral_mean = mean(mean(R_spectral, 1), 2);
            [~, max_ind] = max(R_spectral_mean);
            background = R_spectral(:, :, max_ind);
            background = background ./ max(max(background));
        end
        labels = cell(n_lines, 1);
        colors = jet(n_lines);
        n_channels = size(colors, 2);
        figure_image = repmat(background, 1, 1, n_channels);
        for i = 1:n_lines
            labels{i} = sprintf('Line %d', i);
            n_px_line = length(lines_pixels_indices{i});
            lines_pixels_indices_rgb =...
                repmat(lines_pixels_indices{i}, n_channels, 1) +...
                repelem((0:(n_channels - 1)).' * (image_width * image_height), n_px_line, 1);
            figure_image(lines_pixels_indices_rgb) = reshape(...
                repmat(colors(i, :), n_px_line, 1), [], 1 ...
            );
        end
        figure_image = insertText(...
            figure_image, lines_endpoints(:, 1:2), 1:n_lines,...
            'FontSize', 12, 'BoxColor', colors, 'BoxOpacity', 0.9,...
            'AnchorPoint', 'LeftTop'...
        );
        imshow(figure_image);
        title('Image lines for spectral radiance evaluation')
    end
end

if isfield(options, 'bands_diff')
    diff_band = R_spectral(:, :, options.bands_diff(2)) - R_spectral(:, :, options.bands_diff(1));
    fg_spectral.bands_diff = figure;
    imagesc(diff_band);
    colorbar;
    title(sprintf('Difference image between band %g and band %g in the reference image',...
        lambda(options.bands_diff(1)), lambda(options.bands_diff(2))));

    diff_band = I_spectral(:, :, options.bands_diff(2)) - I_spectral(:, :, options.bands_diff(1));
    fg_spectral.bands_diff(2) = figure;
    imagesc(diff_band);
    colorbar;
    title(sprintf('Difference image between band %g and band %g in the test image',...
        lambda(options.bands_diff(1)), lambda(options.bands_diff(2))));
end

if nargout > 1
    varargout{1} = fg_spectral;
end
    
end
