%% Visual inspection of the relationships between colour channels in images
%
% Produce colour clouds from images, and find their PCA components
%
% See also relativeSensitivity

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 21, 2018

% List of parameters to save with results
parameters_list = {
        'range',...
        'filter_gradient',...
        'bayer_pattern'...
    };

%% Input data and parameters

% Wildcard for 'ls()' to find the images to process. All images are
% expected to be in one directory.
images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181130_LightBox/preprocessed/blended_averaged/pos1_1mmDots_400nm.mat';
images_variable_name = 'I_raw'; % Used only when loading '.mat' files

% Range of pixel values to consider
range = [0, 0.95];

% Filter to low-gradient regions of the images
filter_gradient = true;

% Colour-filter pattern
%
% Single-channel images are assumed to be colour-filter array images, whereas
% three-channel images are assumed to be full-colour or demosaiced images.
bayer_pattern = 'gbrg';

% ## Output directory
output_directory = '/home/llanos/Downloads';

%% Find and process the images

image_filenames = listFiles(images_wildcard);
n_images = length(image_filenames);

n_channels_rgb = 3;
reference_index = 2; % The Green channel is the reference channel
channel_names = {'Red', 'Green', 'Blue'};
label_str = cell(1, n_channels_rgb);
for c = 1:n_channels_rgb
    label_str{c} = sprintf('%s pixel values', channel_names{c});
end
legend_str = {'Pixel values', 'First component', 'Second component', 'Third component'};

components = zeros(n_channels_rgb, n_channels_rgb, n_images);
means = zeros(n_channels_rgb, 1, n_images);
for i = 1:n_images
    I = loadImage(image_filenames{i}, images_variable_name);
    if size(I, 3) == 1
        I_rgb = bayerDownsample(I, bayer_pattern);
        is_raw = true;
    elseif size(I, 3) == n_channels_rgb
        I_rgb = I;
        is_raw = false;
    else
        error('Unexpected number of channels in image "%s".', image_filenames{i});
    end
    
    pixels = reshape(I_rgb, [], n_channels_rgb);
    
    % Filter by intensity
    filter = all(pixels > range(1) & pixels < range(2), 2);
    
    % Filter by gradient
    if filter_gradient
        if is_raw
            channels_filter = false(n_channels_rgb, 1);
            channels_filter(reference_index) = true;
            I_green = bilinearDemosaic(I, bayer_pattern, channels_filter);
        else
            I_green = I(:, :, reference_index);
        end
        G = imresize(imgradient(I_green), 0.5, 'bilinear');
        level = graythresh(G);
        gradient_mask = ~imbinarize(G, level);
        filter = filter & reshape(gradient_mask, [], 1);
    end
    
    pixels = pixels(filter, :);
    
    [components(:, :, i), ~, ~, ~, ~, means(:, :, i)] = pca(pixels);
    [~, name] = fileparts(image_filenames{i});
    title_str = sprintf('Colour cloud for image ''%s''', name);
    fg = plotPCA(pixels, components(:, :, i), means(:, :, i), label_str, legend_str, title_str);
    savefig(...
        fg,...
        fullfile(output_directory, [name '_colorCloud.fig']), 'compact'...
    );
    close(fg);
end

%% Save parameters, intermediate variables, and output data
save_variables_list = [ parameters_list, {...
        'image_filenames',...
        'components',...
        'means',...
    } ];
save_data_filename = fullfile(output_directory, 'ColorCloudData.mat');
save(save_data_filename, save_variables_list{:});