%% Create montages of image patches
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Description
%
% Outputs montages for visual evaluation of small patches in images. All input
% image files ('.mat' files are not expected) should be of the same numeric
% class.
%
% ## Notes
% - This script assumes that all image patches will lie within the image borders
%   (even with padding around the patches).
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
% File created March 11, 2019

% Presently, this script just saves some key information, not a complete set of
% parameters
parameters_list = {};

%% Input data and parameters

patch_size = [101, 101]; % Odd integers: Number of columns, number of rows

% Images to evaluate, as filename prefixes and evaluation patches
images = struct(...
    'd2_book_', struct(...
        'flip', true,...
        'patches', [... % Patches are listed as (x, y) image indices of the centre pixel
            723, 1676; % Smallest doll
            766, 610; % Pattern on book cover
            2131, 1455; % Painting inside glass bottle
            1961, 930; % Text reflected in CD
        ]...
    ),...
    'd2_glass_', struct(...
        'flip', true,...
        'patches', [...
            449, 883; % Dot seen through the crystal ball
            1028, 830; % Dots in shadow of crystal ball
            2083, 557; % Dots seen through star trophy
        ]...
    ),...
    'd2_ship_', struct(...
        'flip', true,...
        'patches', [...
            1662, 533; % Rock, noise pattern, and dolphin on upper right
            496, 1500; % Post in a sail
            1650, 946; % Stern of ship just below dolphin
        ]...
    )...
);

% Algorithms to evaluate, given as filename middle portions
% Each row in the inner cell arrays will become a row in the montage
% Algorithms of the form '${ALG}_krishnan' will be handled by loading the result
% of algorithm 'ALG', followed by blind deconvolution using the method of
% Krishnan et al. 2011, and then chromatic aberration correction using the
% method of Sun et al. 2017.
deconv_suffix = '_krishnan';
algorithm_groups = {...
    {...
        'bilinear_rgb', 'bilinear_cct_rgb', 'bilinear_rgb_krishnan', 'bilinear_channelWarp_rgb';    
        'MATLABdemosaic_rgb', 'MATLABdemosaic_cct_rgb', 'MATLABdemosaic_rgb_krishnan', 'MATLABdemosaic_channelWarp_rgb';
        'ARI_rgb', 'ARI_cct_rgb', 'ARI_rgb_krishnan', 'ARI_channelWarp_rgb';
        'bands6_L2L2NonNeg_DMfw_latent', 'bands6_L1L1NonNeg_DMfw_latent', 'bands6_L2NonNeg_DMfw_latent', 'bands6_L1NonNeg_DMfw_latent'
    };
    {...
        'RGB_L2L2NonNeg_DMfw_ignoreDispersion_rgb', 'RGB_L1L1NonNeg_DMfw_ignoreDispersion_rgb', 'RGB_L2NonNeg_DMfw_ignoreDispersion_rgb', 'RGB_L1NonNeg_DMfw_ignoreDispersion_rgb';
        'RGB_L2L2NonNeg_DMfw_rgb', 'RGB_L1L1NonNeg_DMfw_rgb', 'RGB_L2NonNeg_DMfw_rgb', 'RGB_L1NonNeg_DMfw_rgb';
        'bands6_L2L2NonNeg_DMfw_latent', 'bands6_L1L1NonNeg_DMfw_latent', 'bands6_L2NonNeg_DMfw_latent', 'bands6_L1NonNeg_DMfw_latent';
        'bands6_L2L2NonNeg_DMfw_ignoreDispersion_latent', 'bands6_L1L1NonNeg_DMfw_ignoreDispersion_latent', 'bands6_L2NonNeg_DMfw_ignoreDispersion_latent', 'bands6_L1NonNeg_DMfw_ignoreDispersion_latent'...
    }...
};

% Padding to add for the methods of Krishnan et al. 2011 and Sun et al. 2017.,
% to reduce image border artifacts
padding = 16;
    
% Filename input postfix
input_postfix = '_sRGB.tif';

n_channels_rgb = 3;

% Directory containing the input images
input_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_sRGB';

% Refer to the MATLAB documentation on "Figure Properties"
% Output figure paper size, in inches
output_size_page = [10, 10];
% Output figure width, in inches
output_width = 7;
% Horizontal and vertical offsets from the lower left corner of the page, in
% inches
output_margin = [0.25, 0.25];

% Output directory
output_directory = '/home/llanos/Downloads/montage';

% ## Parameters for Krishnan et al. 2011 and Sun et al. 2017
run('SetFixedParameters.m')

%% Preprocessing

image_prefixes = fieldnames(images);
n_images = length(image_prefixes);
n_groups = length(algorithm_groups);

% Find the class of the image data
I_filename = fullfile(input_directory, sprintf('%s%s%s', image_prefixes{1}, algorithm_groups{1}{1}, input_postfix));
I = imread(I_filename);
class_images = class(I);

half_width = floor((patch_size - 1) / 2);
output_postfix = '.pdf';
font_size = max(12, floor(0.05 * max(patch_size)));
text_offset = [font_size * 2, font_size];

%% Process the images

for i = 1:n_images
    prefix = image_prefixes{i};
    flip = images.(prefix).flip;
    patches = images.(prefix).patches;
    
    for g = 1:n_groups
        group = algorithm_groups{g};        
        n_rows = size(group, 1);
        n_cols = size(group, 2);
        n_algorithms = numel(group);
        labels = cell(n_algorithms, 1);
        algorithm_filenames = cell(n_algorithms, 1);
        
        I_montage = zeros([n_rows * patch_size(2), n_cols * patch_size(1), n_channels_rgb], class_images);
        
        for pc = 1:size(patches, 1)
            patch = patches(pc, :);
            roi_image = [
                patch(2) - half_width(2), patch(2) + half_width(2),...
                patch(1) - half_width(1), patch(1) + half_width(1)
            ];
            algorithm_index = 1;
            
            for col = 1:n_cols
                for row = 1:n_rows
                    algorithm = group{algorithm_index};
                    use_deconv = (strfind(algorithm, deconv_suffix) == (length(algorithm) - length(deconv_suffix) + 1));
                    if use_deconv
                        algorithm_filenames{algorithm_index} = algorithm(1:(end - length(deconv_suffix)));
                    else
                        algorithm_filenames{algorithm_index} = algorithm;
                    end
                    algorithm_filenames{algorithm_index} = fullfile(...
                        input_directory,...
                        sprintf('%s%s%s', prefix, algorithm_filenames{algorithm_index}, input_postfix)...
                    );
                    I = imread(algorithm_filenames{algorithm_index});

                    if ~isa(I, class_images)
                        error('The image "%s" has an unexpected datatype.', algorithm_filenames{algorithm_index})
                    end
                    if size(I, 3) ~= n_channels_rgb
                        error('The image "%s" does not have %d channels.', n_channels_rgb)
                    end
                    
                    roi_montage = [
                        ((row - 1) * patch_size(2)) + 1, row * patch_size(2),...
                        ((col - 1) * patch_size(1)) + 1, col * patch_size(1)...
                    ];

                    if use_deconv
                        krishnan_opts = krishnan2011Options;
                        sub_image_input = im2double(I(...
                            (roi_image(1) - padding):(roi_image(2) + padding),...
                            (roi_image(3) - padding):(roi_image(4) + padding),...
                            :...
                        ));
                        krishnan_opts.blur = sub_image_input(:, :, sun2017Options.reference_channel_index);
                        psf = cell(n_channels_rgb, 1);
                        [...
                            ~, I_deblur, psf{sun2017Options.reference_channel_index}, ~...
                        ] = ms_blind_deconv([], krishnan_opts);
                        close all % 'ms_blind_deconv()' opens figures
                        
                        image_sampling_patch = size(I_deblur);
                        sub_image = zeros([image_sampling_patch, n_channels_rgb]);
                        sub_image(:, :, sun2017Options.reference_channel_index) = I_deblur;
                        
                        psf_sz = min(image_sampling_patch);
                        
                        for c = 1:n_channels_rgb
                            if c ~= sun2017Options.reference_channel_index
                                [sub_image(:, :, c), psf{c}] = ref_deblur(...
                                    I_deblur, sub_image_input(:, :, c),...
                                    psf_sz, sun2017Options.win_sz, sun2017Options.alpha, sun2017Options.beta, sun2017Options.iter...
                                    );
                            end
                        end

                        % Save data for later examination, if desired
                        save(...
                            fullfile(output_directory, sprintf('%spatch%d_%s_data.mat', prefix, pc, algorithm)),...
                            'psf_sz', 'sun2017Options', 'psf', 'krishnan_opts', 'padding', 'roi_image', 'sub_image'...
                        );
                    
                        if isinteger(I_montage)
                            peak_rgb = double(intmax(class_images));
                            sub_image = sub_image * peak_rgb;
                        end
                        sub_image = cast(sub_image, class_images);
                        sub_image = sub_image((padding + 1):(end - padding), (padding + 1):(end - padding), :);
                    else
                        sub_image = I(roi_image(1):roi_image(2), roi_image(3):roi_image(4), :);
                    end
                    
                    if flip
                        sub_image = rot90(sub_image, 2);
                    end
                    I_montage(roi_montage(1):roi_montage(2), roi_montage(3):roi_montage(4), :) = sub_image;
                    
                    algorithm_index = algorithm_index + 1;
                end
            end
            
            output_filename_prefix = fullfile(output_directory, sprintf('%spatch%d_group%d', prefix, pc, g));
            
            % Output
            fg = figure;
            imshow(I_montage);
            set(fg, 'PaperUnits', 'inches');
            set(fg, 'PaperSize', output_size_page);
            output_height = output_width * size(I_montage, 1) / size(I_montage, 2);
            set(fg, 'PaperPosition', [output_margin(1) output_margin(2) output_width, output_height]);
            print(...
                fg, sprintf('%s%s', output_filename_prefix, output_postfix),...
                '-dpdf', '-r300'...
            );
        
            % Annotate and output
            for algorithm_index = 1:n_algorithms
                labels{algorithm_index} = sprintf('(%s)', char(double('a') + (algorithm_index - 1)));
            end
            labels = reshape(reshape(labels, n_rows, n_cols).', [], 1);
            algorithm_index = 1;
            for col = 1:n_cols
                for row = 1:n_rows
                    x = col * patch_size(1) - text_offset(1);
                    y = row * patch_size(2) - text_offset(2);
                    text(x, y, labels{algorithm_index}, 'FontSize', font_size);
                    algorithm_index = algorithm_index + 1;
                end
            end
            print(...
                fg, sprintf('%s_annotated%s', output_filename_prefix, output_postfix),...
                '-dpdf', '-r300'...
            );
            close(fg);
            
            % Output a legend
            labels_table = reshape(reshape(labels, n_rows, n_cols).', [], 1);
            group_table = reshape(group.', [], 1);
            algorithm_filenames_table = reshape(reshape(algorithm_filenames, n_rows, n_cols).', [], 1);
            t = table(labels_table, group_table, algorithm_filenames_table, 'VariableNames', {'Label', 'Algorithm', 'SourceFile'});
            writetable(...
                t, sprintf('%s_legend.csv', output_filename_prefix)...
            );
        end
    end
end