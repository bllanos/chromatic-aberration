%% Convert output data from the modified version of Choi et al. 2017's code
%
% Output radiance and colour images that can be compared to the original
% images converted for input to the modified code from Choi et al. 2017.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## References
% - Choi, I., Jeon, D. S., Gutierrez, D., & Kim, M. H. (2017).
%   "High-Quality Hyperspectral Reconstruction Using a Spectral Prior." ACM
%   Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 218:1-13.
%   10.1145/3130800.3130810

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 28, 2018

% List of parameters to save with results
parameters_list = {
    'conversion_filename'...
};

%% Input data and parameters

% Wildcard for 'ls()' to find the input *reflectance* images
% '.mat' or image files can be loaded
input_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190130_KAIST_crop/ChoiEtAl2017/*recon.mat';

% Conversion information created by 'Choi2017Input.m'
conversion_filename = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/preprocessed/Choi2017Input.mat';

% Output directory for all images and saved parameters
output_directory = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190130_KAIST_crop/ChoiEtAl2017';

%% Load conversion data

variables_required = { 'spectral_weights_output', 'color_weights' };
load(conversion_filename, variables_required{:});
if ~all(ismember(variables_required, who))
    error('One or more of the required colour/spectral space conversion variables is not loaded.')
end

%% Output images

image_filenames = listFiles(input_images_wildcard);
n_images = length(image_filenames);

for i = 1:n_images
    [I, name] = loadImage(image_filenames{i}, 'x_recon');
    I = double(I);
    
    I_radiance = channelConversion(I, spectral_weights_output, 3);
    I_rgb = channelConversion(I, color_weights, 3);
    
    name_postfix = [name, '_choiOutConverted'];

    saveImages(...
        output_directory, name_postfix,...
        I_radiance, '_latent', 'I_latent',...
        I_rgb, '_rgb', 'I_rgb'...
    );
end

%% Save parameters and common data to a file
save_variables_list = [ parameters_list, {...
        'image_filenames'...
    } ];
save_data_filename = fullfile(output_directory, 'Choi2017Output.mat');
save(save_data_filename, save_variables_list{:});