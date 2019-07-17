%% Demo data images preprocessing script
%
% Script used to reduce the volume of data stored in the repository.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 17, 2019

%% Parameters

directory = fullfile('.', 'demo_data', 'captured_images');

crop_regex = '^d1_.+$';
scale_regex = '^d2_.+$';

image_size = [2048, 2448];
image_center = (image_size / 2);
top_left_roi = floor(image_center * 2/3 + image_size * 1/3);
% Preserve Bayer pattern
top_left_roi(mod(top_left_roi, 2) == 0) = top_left_roi(mod(top_left_roi, 2) == 0) + 1;
bottom_right_roi = ceil(image_center * 1/3 + image_size * 2/3);
bottom_right_roi(mod(bottom_right_roi, 2) == 1) = bottom_right_roi(mod(bottom_right_roi, 2) == 1) + 1;
crop_roi = [top_left_roi(1), bottom_right_roi(1), top_left_roi(2), bottom_right_roi(2)];

scale = 1/6;
smaller_size = floor(scale * image_size);
% Preserve Bayer pattern
smaller_size(mod(smaller_size, 2) == 1) = smaller_size(mod(smaller_size, 2) == 1) + 1;

% Colour-filter pattern
bayer_pattern = 'gbrg';

%% Processing

filepaths = listFilesRecursive(crop_regex, directory);

for f = 1:length(filepaths)
    fprintf('File "%s"\n', filepaths{f});
    I = imread(filepaths{f});
    I = I(crop_roi(1):crop_roi(2), crop_roi(3):crop_roi(4));
    I = im2uint8(I);
    imwrite(I, filepaths{f});
end

filepaths = listFilesRecursive(scale_regex, directory);

for f = 1:length(filepaths)
    fprintf('File "%s"\n', filepaths{f});
    I = im2double(imread(filepaths{f}));
    I = bilinearDemosaic(I, bayer_pattern);
    I = imresize(I, smaller_size, 'bilinear');
    I = mosaic(I, bayer_pattern);
    I = im2uint8(I);
    imwrite(I, filepaths{f});
end