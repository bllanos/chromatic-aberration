%% Ad-hoc evaluation of a set of images
%
% Compare images with a reference image

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 25, 2018

%% Input data and parameters

% Wildcard for 'ls()' to find the true image (reference image).
% '.mat' or image files can be loaded
true_image_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/OthersCode/2017_Choi_et_al_highQualityHyperspectralReconstructionUsingASpectralPrior_ACMTransGraphics/inputs/synthetic/KAIST/scene01_oneNinth.mat';
true_image_variable_name = 'img_hs'; % Used only when loading '.mat' files

% Wildcards for 'ls()' to find the comparison images. Use one wildcard per
% image.
% '.mat' or image files can be loaded
images_wildcards = {
    '/home/llanos/GoogleDrive/ThesisResearch/OthersCode/2017_Choi_et_al_highQualityHyperspectralReconstructionUsingASpectralPrior_ACMTransGraphics/outputs/recon_synthetic/scene01_oneNinth_recon.mat'...
    };
% Corresponding variable names to use when loading '.mat' files. Elements
% corresponding to files other than '.mat' files can be empty.
images_variable_names = { 'x_recon' };

% Names for the images
true_image_name = 'GT'; % Name used in figures. Must not contain spaces.
true_image_filename = 'scene01_oneNinth'; % Name used in filenames. Must not contain spaces.
images_names = { 'Choi et al. 2017' };

% `true` for RGB image evaluation, `false` for multispectral or
% hyperspectral image evaluation.
is_rgb = false;

if ~is_rgb
    % Path and filename of a '.mat' file containing the wavelengths
    % or colour channel indices corresponding to the true image
    bands_filename = '/home/llanos/GoogleDrive/ThesisResearch/OthersCode/2017_Choi_et_al_highQualityHyperspectralReconstructionUsingASpectralPrior_ACMTransGraphics/inputs/synthetic/KAIST/scene01_oneNinth.mat';
    bands_variable = 'wvls2b'; % Variable name in the above file

    % Paths and filenames of '.mat' files containing the matrices for
    % converting the images to the colour space of the true image. An
    % identity mapping will be used for empty cells
    spectral_weights_filenames = {[]};
    spectral_weights_variables = {[]}; % Variable names in the above files
end

% Evaluation options, in the format of the 'options' input argument of
% 'evaluateSpectral()' or 'evaluateRGB()' as appropriate. Fields for
% storing figure handles should not be included. The fields 'plot_*' should
% also be omitted.
if is_rgb
    evaluation_options = struct('error_map', true);
else
    evaluation_options = struct(...
        'metric', 'mrae',...
        'error_map', true,...
        'mi_bands', [4, 24],...
        'bands_diff', [4, 24],...
        'radiance', [...
            218, 179, 31, 31;
            321, 179, 31, 31
        ]...
    );
end

% Output directory
output_directory = '/home/llanos/Downloads';

%% Perform the evaluations

n_images = length(images_wildcards);
if n_images ~= length(images_names)
    error('There must be as many other images as there are names for those images.');
end
if ~is_rgb && n_images ~= length(spectral_weights_filenames)
    error('There must be as many other images as there are colour space conversion files.');
end

if ~is_rgb
    bands = loadVariables(bands_filename, bands_variable);
end

if is_rgb
    dp.evaluation.global_rgb = struct;
    dp.evaluation.custom_rgb.(true_image_name) = evaluation_options;
else
    dp.evaluation.global_spectral = struct;
    dp.evaluation.custom_spectral.(true_image_name) = evaluation_options;
    dp.evaluation.custom_spectral.(true_image_filename) = evaluation_options;
    evaluation_plot_colors = jet(n_images);
    evaluation_plot_markers = {'v', 'o', '+', '*', '<', '.', 'x', 's', 'd', '^', 'p', 'h', '>'};
    evaluation_plot_styles = {'--', ':', '-.'};
end
name_params = fullfile(output_directory, true_image_filename);

true_image_filename_load = listFiles(true_image_wildcard);
I_gt = loadImage(true_image_filename_load{1}, true_image_variable_name);

e_rgb_table = [];
e_spectral_table = [];
fg_spectral = struct;
for i = 1:n_images
    image_filename = listFiles(images_wildcards{i});
    I = loadImage(image_filename{1}, images_variable_names{i});
    if is_rgb
        e_rgb_table_current = evaluateAndSaveRGB(...
            I, I_gt, dp, true_image_name, images_names{i},...
            name_params...
        );
        if ~isempty(e_rgb_table)
            e_rgb_table = union(e_rgb_table_current, e_rgb_table);
        else
            e_rgb_table = e_rgb_table_current;
        end
    else
        if isempty(spectral_weights_filenames{i})
            spectral_weights = eye(length(bands));
        else
            spectral_weights = loadVariables(spectral_weights_filenames{i}, spectral_weights_variables{i});
        end
        
        dp.evaluation.global_spectral.plot_color =...
            evaluation_plot_colors(i, :);
        dp.evaluation.global_spectral.plot_marker =...
            evaluation_plot_markers{...
            mod(i - 1, length(evaluation_plot_markers)) + 1 ...
            };
        dp.evaluation.global_spectral.plot_style =...
            evaluation_plot_styles{...
            mod(i - 1, length(evaluation_plot_styles)) + 1 ...
            };
        
        [e_spectral_table_current, fg_spectral] = evaluateAndSaveSpectral(...
            I, I_gt, bands, spectral_weights,...
            dp, true_image_name, images_names{i},...
            name_params, fg_spectral...
            );
        if ~isempty(e_spectral_table)
            e_spectral_table = union(e_spectral_table_current, e_spectral_table);
        else
            e_spectral_table = e_spectral_table_current;
        end
    end
end

%% Output files common to all images

if ~isempty(e_rgb_table)
    writetable(...
        e_rgb_table,...
        fullfile(output_directory, [true_image_filename, '_evaluateRGB.csv'])...
    );
end
if ~isempty(e_spectral_table)
    writetable(...
        e_spectral_table,...
        fullfile(output_directory, [true_image_filename, '_evaluateSpectral.csv'])...
    );
    % Also save completed figures
    evaluateAndSaveSpectral(output_directory, dp, true_image_filename, images_names, fg_spectral);
end