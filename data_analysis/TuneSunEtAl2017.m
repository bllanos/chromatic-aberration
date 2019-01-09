%% Visualization of the differences in dispersion within colour channels
%
% A script testing different window sizes for the method of Sun et al. 2017.
%
% ## References
%
%   Sun, T., Peng, Y., & Heidrich, W. (2017). "Revisiting cross-channel
%   information transfer for chromatic aberration correction." In 2017 IEEE
%   International Conference on Computer Vision (ICCV) (pp. 3268–3276).
%   doi:10.1109/ICCV.2017.352

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 8, 2019


%% Input data and parameters

% PSF estimation window size, as a fraction of the image's largest dimension
psf_sz = logspace(log10(0.05), log10(0.5), 4);

% CCT implementation window size, as a fraction of the image's largest dimension
win_sz = logspace(log10(0.0039), log10(0.2), 20);

% Other parameters for Sun et al. 2017
alpha = 0.3;
beta = 0.3;
n_iter = 3;

% Index of the reference colour channel (Green)
reference_channel_index = 2;

% Wildcard for 'ls()' to find the true image (reference image).
% '.mat' or image files can be loaded
true_image_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180726_Demosaicking_Kodak/PNG_Richard W Franzen/kodim08.png';
true_image_variable_name = []; % Used only when loading '.mat' files

% Wildcard for 'ls()' to find the input image to process
% '.mat' or image files can be loaded
input_image_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_Kodak_demosaicing/evaluation_per_image/08_ARI.tif';
input_image_variable_name = []; % Used only when loading '.mat' files

% Names for the image
image_name = 'kodim08'; % Name used in figures. Must not contain spaces.
image_filename = image_name; % Name used in filenames. Must not contain spaces.

% Image evaluation options used with 'evaluateRGB()'
evaluation_options = struct; %struct('error_map', true);

% Output directory
output_directory = '/home/llanos/Downloads/sun2017_parameter_matrix';

%% Test all combinations of parameters

dp.evaluation.global_rgb = struct;
dp.evaluation.custom_rgb.(image_name) = evaluation_options;

fullfilename_params = fullfile(output_directory, image_filename);

true_image_filename_load = listFiles(true_image_wildcard);
I_gt = im2double(loadImage(true_image_filename_load{1}, true_image_variable_name));
image_sampling_3 = size(I_gt);

input_image_filename = listFiles(input_image_wildcard);
I_in = im2double(loadImage(input_image_filename{1}, input_image_variable_name));
if any(size(I_in) ~= image_sampling_3)
    error('The reference image must have the same dimensions as the input image.');
end

e_rgb_table = [];
n_win_sz = length(win_sz);
n_psf_sz = length(psf_sz);
error = zeros(n_win_sz, n_psf_sz);

max_image_size = max(image_sampling_3(1:2));
for i = 1:n_win_sz
    win_sz_i = ceil(max_image_size * win_sz(i));
    
    for j = 1:n_psf_sz
        psf_sz_j = ceil(max_image_size * psf_sz(j));
        
        I_out = zeros(image_sampling_3);
        I_out(:, :, reference_channel_index) = I_in(:, :, reference_channel_index);
        
        for c = 1:image_sampling_3(3)
            if c ~= reference_channel_index
                [I_out(:, :, c), psf] = ref_deblur(...
                    I_in(:, :, reference_channel_index), I_in(:, :, c),...
                    psf_sz_j, win_sz_i, alpha, beta, n_iter...
                );
            end
        end
        
        fullfilename_params_ij = sprintf('%s_%g_%g', fullfilename_params, win_sz(i), psf_sz(j));
        filename_params_ij = sprintf('%s_%g_%g', image_filename, win_sz(i), psf_sz(j));
        e_rgb_table_current = evaluateAndSaveRGB(...
            I_out, I_gt, dp, image_name, image_name,...
            fullfilename_params_ij...
        );
        error(i, j) = e_rgb_table_current{1, {'CPSNR'}};
        if ~isempty(e_rgb_table)
            e_rgb_table = union(e_rgb_table_current, e_rgb_table);
        else
            e_rgb_table = e_rgb_table_current;
        end
        saveImages(...
            output_directory, filename_params_ij,...
            I_out, '_rgb', 'I_rgb'...
        );
    end
end

writetable(...
    e_rgb_table,...
    fullfile(output_directory, [image_filename, '_evaluateRGB.csv'])...
);

[win_sz_grid, psf_sz_grid] = ndgrid(win_sz, psf_sz);

fg = figure;
hold on
if n_win_sz > 1 && n_psf_sz > 1
    surf(win_sz_grid, psf_sz_grid, error);
    xlabel('CCT window size');
    ylabel('PSF window size');
    zlabel('CPSNR');
elseif n_psf_sz > 1
    plot(...
        psf_sz_grid, error, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--'...
    );
    xlabel('PSF window size');
    ylabel('CPSNR');
else
    plot(...
        win_sz_grid, error, 'LineWidth', 2,...
        'Marker', 'o', 'LineStyle', '--'...
    );
    xlabel('CCT window size');
    ylabel('CPSNR');
end
title(sprintf(...
    'Performance on a %d x %d image',...
    image_sampling_3(1), image_sampling_3(2)...
));
hold off
savefig(...
    fg,...
    [fullfilename_params 'error.fig'], 'compact'...
);
close(fg);