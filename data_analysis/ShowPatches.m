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
inset_fraction = 0.18; % Dimensions of region to extract for the inset, as a fraction of the patch size
inset_magnification = 2; % Integer factor by which to magnify the region extracted for the inset

% Images to evaluate, as filename prefixes and evaluation patches
%
% Patches are listed as (x, y) image indices of the centre pixel, followed by
% the (x, y) image indices of the centre pixel for the inset. The inset must lie
% within the borders of the patch.
images = struct(...
    'd1_colorChecker30cm_', struct(...
        'flip', false,...
        'patches', [...
            1963, 1526, 1981, 1503; % 'mm' sign
            286, 1471, 290, 1471 % 'X-rite'
        ]...
    ),...
    'd2_book_', struct(...
        'flip', true,...
        'patches', [...
            723, 1676, 726, 1688; % Smallest doll
            766, 610, 792, 634; % Pattern on book cover
            2131, 1455, 2094, 1436; % Painting inside glass bottle
            1789, 1366, 1750, 1379; % Box for glass bottle
            437, 196, 437, 196 % Edge between two colours on book cover
        ]...
    ),...
    'd2_glass_', struct(...
        'flip', true,...
        'patches', [...
            927, 1071, 912, 1070; % Dots in shadow of crystal ball
            2361, 743, 2353, 730 % Dot beside star trophy
        ]...
    ),...
    'd2_ship_', struct(...
        'flip', true,...
        'patches', [...
            633, 465, 640, 465; % Specular highlight on dolphin's head
            1662, 533, 1677, 521; % Rock, noise pattern, and dolphin on upper right
            2363, 1359, 2363, 1359 % Noise pattern
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
        'bands6_L2NonNeg_DMfw_latent', 'bands8_L2NonNeg_DMfw_latent';
        'bands6_L1NonNeg_DMfw_latent', 'bands8_L1NonNeg_DMfw_latent'
    };
    {...
        'bilinear_rgb', 'bilinear_cct_rgb', 'bilinear_rgb_krishnan', 'bilinear_channelWarp_rgb', 'bands8_Lap2NonNeg_DMfw_latent';
        'ARI_rgb', 'ARI_cct_rgb', 'ARI_rgb_krishnan', 'ARI_channelWarp_rgb', 'bands8_L1NonNeg_DMfw_latent'
    };
    {...
        'RGB_L1SpatialLap2NonNeg_DMfw_rgb', 'RGB_Lap1NonNeg_DMfw_rgb', 'RGB_Lap2NonNeg_DMfw_rgb', 'RGB_L2NonNeg_DMfw_rgb', 'RGB_L1NonNeg_DMfw_rgb';
        'bands8_L1SpatialLap2NonNeg_DMfw_latent', 'bands8_Lap1NonNeg_DMfw_latent', 'bands8_Lap2NonNeg_DMfw_latent', 'bands8_L2NonNeg_DMfw_latent', 'bands8_L1NonNeg_DMfw_latent'
    };
    {...
        'RGB_L1NonNeg_DMfw_rgb_unwarped', 'RGB_Lap2NonNeg_DMfw_rgb_unwarped',  'bands8_Lap2NonNeg_DMfw_latent_unwarped', 'bands8_L1NonNeg_DMfw_latent_unwarped', ;
        'RGB_L1NonNeg_DMfw_rgb', 'RGB_Lap2NonNeg_DMfw_rgb', 'bands8_Lap2NonNeg_DMfw_latent', 'bands8_L1NonNeg_DMfw_latent';
    }
};

% Padding to add for the methods of Krishnan et al. 2011 and Sun et al. 2017.,
% to reduce image border artifacts
padding = 16;

% Filename input postfix
input_postfix = '.tif';

n_channels_rgb = 3;

% Directory containing the input images
input_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190610_ThesisResults_CopiedImages/images_MHomogCorrectedRGB';

% Directory containing the linear, uncorrected images, for patch-specific
% processing with Krishnan et al. 2011 and Sun et al. 2017.
input_directory_uncorrected = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/run_on_dataset_allEstimatedImages_MATFiles';
input_variable_uncorrected = 'I_rgb'; % Used only when loading '.mat' files 
input_postfix_uncorrected = '.mat';

% Path and filename of a '.mat' file containing the conversion matrix for
% raw colour channels to XYZ
xyz_weights_filename = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190421_ComputarLens_revisedAlgorithms/CalibrateColorCorrectionData_unfiltered.mat';
xyz_weights_variable = 'M_homog'; % Variable name in the above file
% Whitepoint to use for XYZ to sRGB conversion
whitepoint = [1, 1, 1];

% Refer to the MATLAB documentation on "Figure Properties"
% Output figure paper size, in inches
output_size_page = [15, 15];
% Output figure width, in inches
output_width = 9.5;
% Horizontal and vertical offsets from the lower left corner of the page, in
% inches
output_margin = [0.25, 0.25];

% Output directory
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190610_ThesisResults_CopiedImages/montage';

% ## Parameters for Krishnan et al. 2011 and Sun et al. 2017
run('SetFixedParameters.m')

%% Preprocessing

image_prefixes = fieldnames(images);
n_images = length(image_prefixes);
n_groups = length(algorithm_groups);

xyz_weights = loadVariables(xyz_weights_filename, xyz_weights_variable);

% Find the class of the image data
I_filename = fullfile(input_directory, sprintf('%s%s%s', image_prefixes{1}, algorithm_groups{1}{1}, input_postfix));
I = imread(I_filename);
class_images = class(I);

inset_size = floor(inset_fraction * patch_size);
inset_size(mod(inset_size, 2) == 0) = inset_size(mod(inset_size, 2) == 0) + 1;
inset_size_output = inset_magnification * inset_size;
roi_inset_output = [
    patch_size(2) - inset_size_output(2) + 1, patch_size(2),...
    1, inset_size_output(1)
];
half_width = floor((patch_size - 1) / 2);
inset_half_width = floor((inset_size - 1) / 2);
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
            patch = patches(pc, 1:2);
            roi_image = [
                patch(2) - half_width(2), patch(2) + half_width(2),...
                patch(1) - half_width(1), patch(1) + half_width(1)
            ];
            inset = patches(pc, 3:4) - roi_image([3, 1]) + 1;
            roi_inset = [
                inset(2) - inset_half_width(2), inset(2) + inset_half_width(2),...
                inset(1) - inset_half_width(1), inset(1) + inset_half_width(1)
            ];
            algorithm_index = 1;
            
            for col = 1:n_cols
                for row = 1:n_rows
                    algorithm = group{algorithm_index};
                    use_deconv = (strfind(algorithm, deconv_suffix) == (length(algorithm) - length(deconv_suffix) + 1));
                    if use_deconv
                        algorithm_filenames{algorithm_index} = algorithm(1:(end - length(deconv_suffix)));
                        algorithm_filenames{algorithm_index} = fullfile(...
                            input_directory_uncorrected,...
                            sprintf('%s%s%s', prefix, algorithm_filenames{algorithm_index}, input_postfix_uncorrected)...
                        );
                        I = loadImage(algorithm_filenames{algorithm_index}, input_variable_uncorrected);
                    else
                        algorithm_filenames{algorithm_index} = algorithm;
                        algorithm_filenames{algorithm_index} = fullfile(...
                            input_directory,...
                            sprintf('%s%s%s', prefix, algorithm_filenames{algorithm_index}, input_postfix)...
                        );
                        I = imread(algorithm_filenames{algorithm_index});
                        if ~isa(I, class_images)
                            error('The image "%s" has an unexpected datatype.', algorithm_filenames{algorithm_index})
                        end
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
                        sub_image_input = I(...
                            (roi_image(1) - padding):(roi_image(2) + padding),...
                            (roi_image(3) - padding):(roi_image(4) + padding),...
                            :...
                        );
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
                    
                        % Colour conversion
                        sub_image = channelConversion(sub_image, xyz_weights);
                        sub_image = xyz2rgb(...
                            sub_image, 'ColorSpace', 'srgb',...
                            'WhitePoint', whitepoint,...
                            'OutputType', class_images...
                        );
                        sub_image = sub_image((padding + 1):(end - padding), (padding + 1):(end - padding), :);
                    else
                        sub_image = I(roi_image(1):roi_image(2), roi_image(3):roi_image(4), :);
                    end
                    
                    % Prepare inset
                    sub_image_inset = sub_image(roi_inset(1):roi_inset(2), roi_inset(3):roi_inset(4), :);
                    sub_image_inset = imresize(...
                        sub_image_inset, [inset_size_output(2) inset_size_output(1)],...
                        'nearest'...
                    );
                    % Add a small border
                    sub_image_inset(1, :, :) = 0;
                    sub_image_inset(end, :, :) = 0;
                    sub_image_inset(:, 1, :) = 0;
                    sub_image_inset(:, end, :) = 0;
                    
                    if flip
                        sub_image = rot90(sub_image, 2);
                        sub_image_inset = rot90(sub_image_inset, 2);
                    end
                    sub_image(roi_inset_output(1):roi_inset_output(2), roi_inset_output(3):roi_inset_output(4), :) = sub_image_inset;
                    

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
            labels = reshape(reshape(labels, n_cols, n_rows).', [], 1);
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