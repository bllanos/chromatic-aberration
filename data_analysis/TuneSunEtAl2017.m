%% Tuning the method of Sun et al. 2017
%
% A script testing different parameters for the method of Sun et al. 2017 so as
% to minimize colour fringing at image edges, but maximize colour faithfulness.
%
% The input image is an image of a ColorChecker chart, corrected for vignetting,
% along with a label image identifying the different patches as well as the
% frame around the patches. Ground truth colours will be taken to be the average
% colours in eroded versions of the patches. Chromatic aberration correction is
% evaluated in narrow regions abutting the edges of the patches (and not
% extending outside of the edges): The error between pixels and their ground
% truth colours is summed over all pixels in the evaluation region.
%
% An output image will be saved for each set of parameters tested.
%
% ## References
%
%   Sun, T., Peng, Y., & Heidrich, W. (2017). "Revisiting cross-channel
%   information transfer for chromatic aberration correction." In 2017 IEEE
%   International Conference on Computer Vision (ICCV) (pp. 3268–3276).
%   doi:10.1109/ICCV.2017.352
%
% The above method uses the following blind deconvolution method as a
% preprocessing step:
%
%   Krishnan, D., Tay, T. & Fergus, R. (2011). "Blind deconvolution using a
%   normalized sparsity measure." In IEEE Conference on Computer Vision and
%   Pattern Recognition (CVPR) (pp. 233–240).

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 8, 2019


%% Input data and parameters

% ## Parameters for Krishnan et al. 2011

% Kernel sizes must be odd integers. I will set the kernel size based on the
% estimated amount of dispersion.
kernel_sz = 9;

min_lambda = logspace(-5, 5, 30);

% Window in which to estimate the PSF: (y1, x1, y2, x2) of the top left and
% bottom right corners. (Set it to an empty array to use the entire image.)
krishnan_opts.kernel_est_win = [];

% ## Parameters for Sun et al. 2017

% PSF estimation window size, as a fraction of the image's largest dimension
psf_sz = logspace(log10(0.05), log10(0.2), 10);

% CCT implementation window size, as a fraction of the image's largest dimension
win_sz = logspace(log10(0.0016), log10(0.2), 20);

% Other parameters for Sun et al. 2017
alpha = 0.3;
beta = 0.3;
n_iter = 3;

% Index of the reference colour channel (Green)
reference_channel_index = 2;

% ## Testing parameters

% Wildcard for 'ls()' to find the label image (an image file, not a '.mat' file)
label_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/exposure_blending/d2_colorChecker30cm_unfiltered_background0_patches1to24_frame25.png';

% Wildcard for 'ls()' to find the input image to process.
% '.mat' or image files can be loaded
input_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/colorChecker_preprocessed/d2_colorChecker30cm_unfiltered_vc.mat';
%input_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/exposure_blending/d2_colorChecker30cm_unfiltered.mat';
input_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Width of the evaluation region in pixels
eval_width = 12;

bayer_pattern = 'gbrg'; % Colour-filter pattern

% Quantiles used for clipping to produce nice output images in image-format
% files
quantiles = [0.01, 0.99];

% Split the image in half to reduce memory consumption
split_image = true;

% Output directory
output_directory = '/home/llanos/Downloads/sun2017_parameter_matrix';

%% Prepare the image

label_filename = listFiles(label_wildcard);
I_label = imread(label_filename{1});

input_filename = listFiles(input_wildcard);
I_raw = loadImage(input_filename{1}, input_variable_name);

channel_mask = bayerMask(size(I_raw, 1), size(I_raw, 2), bayer_pattern);

% Crop to the colorChecker chart
n_patches = 24;
[row, col] = find(I_label == (n_patches + 1));
roi = [min(row), max(row), min(col), max(col)];
% Avoid changing the Bayer pattern
roi(mod(roi, 2) == 0) = roi(mod(roi, 2) == 0) - 1;
I_label = I_label(roi(1):roi(2), roi(3):roi(4));
I_raw = I_raw(roi(1):roi(2), roi(3):roi(4));
channel_mask = channel_mask(roi(1):roi(2), roi(3):roi(4), :);
image_sampling = [size(I_raw, 1), size(I_raw, 2)];

[~, I_out_filename] = fileparts(input_filename{1});
I_out_remapped = clipAndRemap(I_raw, 'uint8', 'quantiles', quantiles);
saveImages(...
    'image', output_directory, I_out_filename,...
    I_out_remapped, '_cropped', []...
);

% Find ground truth colours
n_channels_rgb = 3;
se_clean = strel('square', 3);
se_eval = strel('disk', eval_width);
eval_mask = false(image_sampling);
I_gt = zeros([image_sampling, n_channels_rgb]);
for pc = 1:n_patches
    patch_mask = imopen(I_label == pc, se_clean);
    patch_average_mask = imerode(patch_mask, se_eval); % This also cleans up some noise in the mask
    patch_eval_mask = logical(patch_mask - patch_average_mask);
    eval_mask = eval_mask | patch_eval_mask;
    
    for c = 1:n_channels_rgb
        I_gt_c = I_gt(:, :, c);
        I_gt_c(patch_eval_mask) = mean(I_raw(patch_average_mask & channel_mask(:, :, c)));
        I_gt(:, :, c) = I_gt_c;
    end
end
eval_mask = repmat(eval_mask, 1, 1, n_channels_rgb);

I_out_remapped = clipAndRemap(I_gt, 'uint8', 'quantiles', quantiles);
saveImages(...
    'image', output_directory, I_out_filename,...
    I_out_remapped, '_evaluationColors', []...
);

I_gt_columnar = reshape(I_gt(eval_mask), [], n_channels_rgb);

% Produce a full-colour input image
I_rgb = demosaic_ARI(...
    repmat(I_raw, 1, 1, n_channels_rgb), bayer_pattern...
);

I_out_remapped = clipAndRemap(I_rgb, 'uint8', 'quantiles', quantiles);
saveImages(...
    'image', output_directory, I_out_filename,...
    I_out_remapped, '_ariDemosaiced', []...
);

%% Test all combinations of parameters for Krishnan et al. 2011

n_kernel_sz = length(kernel_sz);
n_min_lambda = length(min_lambda);
error_krishnan = zeros(n_kernel_sz, n_min_lambda);

krishnan_opts.prescale = 1;
krishnan_opts.k_reg_wt = 1;
krishnan_opts.gamma_correct = 1;
krishnan_opts.k_thresh = 0.0;
krishnan_opts.kernel_init = 3;
krishnan_opts.delta = 0.001;
krishnan_opts.x_in_iter = 2; 
krishnan_opts.x_out_iter = 2;
krishnan_opts.xk_iter = 21;
krishnan_opts.nb_lambda = 3000;
krishnan_opts.nb_alpha = 1.0;
krishnan_opts.use_ycbcr = 1;

krishnan_opts.blur = I_rgb(:, :, reference_channel_index);

min_error = Inf;

for i = 1:n_kernel_sz
    kernel_sz_i = kernel_sz(i);
    krishnan_opts.kernel_size = kernel_sz_i;
    
    for j = 1:n_min_lambda
        min_lambda_j = min_lambda(j);
        
        krishnan_opts.min_lambda = min_lambda_j;
        
        [yorig, I_deblur_ij, kernel, opts] = ms_blind_deconv([], krishnan_opts);
        close all
        
        [~, mrae] = metrics(...
            reshape(I_deblur_ij(eval_mask(:, :, reference_channel_index)), [], 1),...
            I_gt_columnar(:, reference_channel_index),...
            2, [], true...
        );
        error_krishnan(i, j) = mean(mrae);
        
        if min_error > error_krishnan(i, j)
            I_deblur = I_deblur_ij;
            kernel_sz_opt = kernel_sz_i;
            min_lambda_opt = min_lambda_j;
            min_error = error_krishnan(i, j);
        end
        
        I_out_remapped = clipAndRemap(I_deblur_ij, 'uint8', 'quantiles', quantiles);
        saveImages(...
            'image', output_directory, I_out_filename,...
            I_out_remapped, sprintf('_kernel%g_minLambda%g', kernel_sz(i), min_lambda(j)), []...
        );
    end
end

[kernel_sz_grid, min_lambda_grid] = ndgrid(kernel_sz, min_lambda);

fg = figure;
hold on
if n_kernel_sz > 1 && n_min_lambda > 1
    surf(kernel_sz_grid, min_lambda_grid, error_krishnan);
    xlabel('Kernel size');
    ylabel('min \lambda');
    zlabel('MRAE');
elseif n_min_lambda > 1
    plot(...
        min_lambda_grid, error_krishnan, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--'...
    );
    xlabel('min \lambda');
    ylabel('MRAE');
else
    plot(...
        kernel_sz_grid, error_krishnan, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--'...
    );
    xlabel('Kernel size');
    ylabel('MRAE');
end
title(sprintf(...
    'Krishnan et al. 2011 performance on a %d x %d image',...
    image_sampling(1), image_sampling(2)...
));
hold off
savefig(...
    fg,...
    fullfile(output_directory, [I_out_filename '_error_krishnan.fig']), 'compact'...
);
close(fg);

%% Test all combinations of parameters for Sun et al. 2017 on the best image from Krishnan et al. 2011

n_win_sz = length(win_sz);
n_psf_sz = length(psf_sz);
error_sun = zeros(n_win_sz, n_psf_sz);
mi_rg = zeros(n_win_sz, n_psf_sz);
mi_gb = zeros(n_win_sz, n_psf_sz);
mi_rb = zeros(n_win_sz, n_psf_sz);

min_error = Inf;

max_image_size = max(image_sampling);
for i = 1:n_win_sz
    win_sz_i = ceil(max_image_size * win_sz(i));
    
    for j = 1:n_psf_sz
        psf_sz_j = ceil(max_image_size * psf_sz(j));
        
        I_out = zeros([image_sampling, n_channels_rgb]);
        I_out(:, :, reference_channel_index) = I_deblur;
        
        for c = 1:n_channels_rgb
            if c ~= reference_channel_index
                if split_image
                    [I_out(1:floor(end / 2), :, c), psf] = ref_deblur(...
                        I_deblur(1:floor(end / 2), :), I_rgb(1:floor(end / 2), :, c),...
                        psf_sz_j, win_sz_i, alpha, beta, n_iter...
                        );
                    [I_out((floor(end/2) + 1):end, :, c), psf] = ref_deblur(...
                        I_deblur((floor(end/2) + 1):end, :), I_rgb((floor(end/2) + 1):end, :, c),...
                        psf_sz_j, win_sz_i, alpha, beta, n_iter...
                        );
                else
                    [I_out(:, :, c), psf] = ref_deblur(...
                        I_deblur, I_rgb(:, :, c),...
                        psf_sz_j, win_sz_i, alpha, beta, n_iter...
                        );
                end
            end
        end
        
        [~, mrae] = metrics(...
            reshape(I_out(eval_mask), [], n_channels_rgb), I_gt_columnar,...
            2, [], true);
        error_sun(i, j) = mean(mrae);
        
        if min_error > error_sun(i, j)
            win_sz_opt = win_sz_i;
            psf_sz_opt = psf_sz_j;
            min_error = error_sun(i, j);
        end
        
        I_out_remapped = clipAndRemap(I_out, 'uint8', 'quantiles', quantiles);
        
        mi_rg(i, j) = MI_GG(...
            I_out_remapped(:, :, 1),...
            I_out_remapped(:, :, 2)...
        );
        mi_gb(i, j) = MI_GG(...
            I_out_remapped(:, :, 2),...
            I_out_remapped(:, :, 3)...
        );
        mi_rb(i, j) = MI_GG(...
            I_out_remapped(:, :, 1),...
            I_out_remapped(:, :, 3)...
        );

        saveImages(...
            'image', output_directory, I_out_filename,...
            I_out_remapped, sprintf('_win%g_psf%g', win_sz_i, psf_sz_j), []...
        );
    end
end

[win_sz_grid, psf_sz_grid] = ndgrid(win_sz, psf_sz);

fg = figure;
hold on
if n_win_sz > 1 && n_psf_sz > 1
    surf(win_sz_grid, psf_sz_grid, error_sun);
    xlabel('CCT window size');
    ylabel('PSF window size');
    zlabel('MRAE');
elseif n_psf_sz > 1
    plot(...
        psf_sz_grid, error_sun, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--'...
    );
    xlabel('PSF window size');
    ylabel('MRAE');
else
    plot(...
        win_sz_grid, error_sun, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--'...
    );
    xlabel('CCT window size');
    ylabel('MRAE');
end
title(sprintf(...
    'Sun et al. 2017 performance on a %d x %d image',...
    image_sampling(1), image_sampling(2)...
));
hold off
savefig(...
    fg,...
    fullfile(output_directory, [I_out_filename '_error_sun.fig']), 'compact'...
);
close(fg);

fg = figure;
hold on
if n_win_sz > 1 && n_psf_sz > 1
    surf(win_sz_grid, psf_sz_grid, mi_rg, 'FaceAlpha', 0.5, 'FaceColor', [1, 1, 0]);
    surf(win_sz_grid, psf_sz_grid, mi_gb, 'FaceAlpha', 0.5, 'FaceColor', [0, 1, 1]);
    surf(win_sz_grid, psf_sz_grid, mi_rb, 'FaceAlpha', 0.5, 'FaceColor', [1, 0, 1]);
    xlabel('CCT window size');
    ylabel('PSF window size');
    zlabel('Mutual Information');
elseif n_psf_sz > 1
    plot(...
        psf_sz_grid, mi_rg, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--', 'Color', 'y'...
    );
    plot(...
        psf_sz_grid, mi_gb, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--', 'Color', 'c'...
    );
    plot(...
        psf_sz_grid, mi_rb, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--', 'Color', 'm'...
    );
    xlabel('PSF window size');
    ylabel('Mutual Information');
else
    plot(...
        win_sz_grid, mi_rg, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--', 'Color', 'y'...
    );
    plot(...
        win_sz_grid, mi_gb, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--', 'Color', 'c'...
    );
    plot(...
        win_sz_grid, mi_rb, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--', 'Color', 'm'...
    );
    xlabel('CCT window size');
    ylabel('Mutual Information');
end
title('Sun et al. 2017 mutual information');
legend('Red-Green', 'Green-Blue', 'Red-Blue');
hold off
savefig(...
    fg,...
    fullfile(output_directory, [I_out_filename '_mi_sun.fig']), 'compact'...
);
close(fg);