%% Measurement of chromatic aberration
% Measure chromatic aberration across the image plane, using probability
% distribution divergence formulations.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% Refer to the first code section below.
%
% ## Output
%
% TODO

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 4, 2017

%% Input data and parameters

% Directory containing the input images
in_directory = 'C:\Users\llanos\Google Drive\ThesisResearch\Data and Results\20170808_OpticalTableMatrix\averaged';

% Wildcard for image filenames
wildcard = 'd44_a1.4_far_disksWhite.tif';

% Colour-filter pattern used to interpret raw images
align = 'gbrg';

% Chromatic aberration measurement techniques
f_div = {
    @censusDivergence
    };
f_div_names = {'Census transform'};

% Window definition for chromatic aberration measurement
% Refer to the documentation of `divergenceMap()` for details
p = 0.01;
if p == 0
    half_width = 5;
    half_height = 5;
    std_weights = 2 * [half_width half_height] + 1;
else
    std_weights = 3;
end

% Which pairs of channels to measure chromatic aberration between?
% The first number in each pair (row) is the "reference" channel
color_channel_pairs = [
    2, 1;
    2, 3;
    1, 3
    ];
color_channel_names = {'Red', 'Green', 'Blue'};

% ## Debugging Flags
show_sampled_image = true;
show_chromatic_aberration = true;

%% Process input images
names = ls(fullfile(in_directory, wildcard));
n_img = size(names, 1);
n_f = length(f_div);
n_pairs = size(color_channel_pairs, 1);

if p == 0
    window_name = sprintf('(%d x %d window)', std_weights(1), std_weights(2));
else
    window_name = sprintf('(\\sigma = %g Gaussian window, p = %g)', std_weights, p);
end

for i = 1:n_img
    I_raw = imread(fullfile(in_directory, names(i, :)));
    I_downsampled = bayerDownsample( I_raw, align );
    if show_sampled_image
        figure; imshow(I_downsampled / max(max(max(I_downsampled))));
        title('Original image, downsampled and approximately demosaiced')
    end
    for k = 1:n_f
        f_div_k = f_div{k};
        for c = 1:n_pairs
            c_1 = color_channel_pairs(c, 1);
            c_2 = color_channel_pairs(c, 2);
            I_1 = I_downsampled(:, :, c_1);
            I_2 = I_downsampled(:, :, c_2);
            div_map = divergenceMap( I_1, I_2, f_div_k, std_weights, p );
            if show_chromatic_aberration
                figure;
                imagesc(div_map)
                colorbar
                xlabel('X');
                ylabel('Y');
                colormap jet
                c_bar = colorbar;
                c_bar.Label.String = 'Chromatic aberration measurement';
                title(sprintf(...
                    '%s measurements on %s relative to %s %s',...
                    f_div_names{k},...
                    color_channel_names{c_2}, color_channel_names{c_1},...
                    window_name...
                    ));
            end
        end
    end
end