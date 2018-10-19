function [dataset_params] = describeDataset(name)
% DESCRIBEDATASET  Create a structure describing a dataset for evaluation
%
% ## Syntax
% dataset_params = describeDataset(name)
%
% ## Description
% dataset_params = describeDataset(name)
%   Returns a structure describing the dataset
%
% ## Input Arguments
%
% name -- Dataset name
%   A character vector containing the name of a dataset. `name` must be one
%   of the recognized names listed below.
%
% ## Output Arguments
%
% dataset_params -- Dataset description
%   A structure describing how to load the dataset, which image estimation
%   algorithms to run on the dataset, and how to evaluate the results of
%   the image estimation algorithms. `dataset_params` has the following
%   fields, all of which are either empty, or store character vectors:
%   - 'raw_images_wildcard': A wildcard for 'ls()' to find the input RAW
%     images of the dataset. If empty, the RAW images are to be generated
%     from RGB images.
%   - 'raw_images_variable': The variable name to use for loading RAW
%     images from '.mat' files, where applicable.
%   - 'rgb_images_wildcard': A wildcard for 'ls()' to find the true RGB
%     images of the dataset. If empty, the RGB images are to be generated
%     from spectral images.
%   - 'rgb_images_variable': The variable name to use for loading RGB
%     images from '.mat' files, where applicable.
%   - 'spectral_images_wildcard': A wildcard for 'ls()' to find the true
%     spectral images of the dataset. If empty, the dataset lacks spectral
%     images.
%   - 'spectral_images_variable': The variable name to use for loading
%     spectral images from '.mat' files, where applicable.
%   - 'spectral_reflectances': A Boolean field, required only if the
%     dataset has spectral images. If 'spectral_reflectances' is `true`,
%     then the spectral images are reflectance images. Otherwise, the
%     spectral images are radiance images. Reflectance images need to be
%     preprocessed by multiplying them by the spectral power distribution
%     of an illuminant.
%   - 'dispersion_rgb_forward': The filename and path of the model of
%     dispersion of the RGB channels relative to the reference channel,
%     stored in a '.mat' file. The dispersion model is a function of
%     coordinates in the reference channel (usually the Green channel). If
%     empty, any correction of lateral chromatic aberration in the RGB
%     domain will be based on priors only, not on calibration data.
%   - 'dispersion_rgb_reverse': Similar to 'dispersion_rgb_forward', but
%     the dispersion model is a function of coordinates in the given
%     channel, not the reference channel.
%   - 'dispersion_spectral_reverse': The filename and path of the model of
%     dispersion of wavelength bands relative to the reference band, stored
%     in a '.mat' file. The dispersion model is a function of coordinates
%     in the given band. If empty, any correction of lateral chromatic
%     aberration in the spectral domain will be based on priors only, not
%     on calibration data.
%   - 'color_map': The filename and path of the colour space conversion
%     '.mat' data file. This data file is used to convert spectral
%     information to RGB (or other input colour space), and has the form
%     described in the documentation of 'CorrectByHyperspectralADMM.m'.
%   - 'wavelengths': The filename and path of a '.mat' file containing a
%     variable, `bands`, storing the wavelengths which determine the
%     spectral sampling of spectral images. This field is required only if
%     the dataset has spectral images.
%   - 'evaluation': Evaluation parameters, a structure with the following
%     fields:
%     - 'global_rgb': A structure of the form of the `options` input
%       argument of 'evaluateRGB()', describing default evaluation options
%       for all RGB images.
%     - 'custom_rgb': A structure, where the value of each field is of the
%       form of the `options` input argument of 'evaluateRGB()', and
%       describes custom evaluation options for an RGB image. The
%       fieldnames of 'custom_rgb' are the filenames (excluding file
%       extensions and common suffixes) of the RGB image files that will be
%       loaded, or of the spectral images that will be converted to RGB
%       images. Fields not present in the structures should default to the
%       values given by 'global_rgb'.
%     - 'global_spectral': A structure of the form of the `options` input
%       argument of 'evaluateSpectral()', describing default evaluation
%       options for all spectral images. Fields for storing figure handles
%       should not be included. The fields 'plot_*' should also be
%       omitted. 'global_spectral' is needed only for datasets with
%       spectral images.
%     - 'custom_spectral': A structure, where the value of each field is
%       of the form of the `options` input argument of
%       'evaluateSpectral()', and describes custom evaluation options for a
%       spectral image. The fieldnames of 'custom_spectral' are the
%       filenames (excluding file extensions and common suffixes) of the
%       spectral image files that will be loaded. Fields not present in the
%       structures should default to the values given by 'global_spectral'.
%       Fields for storing figure handles should not be included. The
%       fields 'plot_*' should also be omitted. 'custom_spectral' is needed
%       only for datasets with spectral images.
%
% ## Recognized (high-quality) datasets
% - 'kodak': The Kodak Lossless True Color Image Suite dataset, often
%   used to evaluate demosaicking algorithms.
%   - Source: http://r0k.us/graphics/kodak/
%   - Maintainer: Richard W Franzen
% - 'kaist': The KAIST Dataset of Hyperspectral Reflectance Images
%   - Source: http://vclab.kaist.ac.kr/siggraphasia2017p1/kaistdataset.html
%   - Reference:
%     Choi, I., Jeon, D. S., Nam, G., Gutierrez, D., & Kim, M. H. (2017).
%       "High-Quality Hyperspectral Reconstruction Using a Spectral Prior."
%       ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6),
%       218:1–13. doi:10.1145/3130800.3130810
%
% There are some datasets defined in the code used for preliminary testing
% only, such as:
% - 'kaist-crop': The color checker chart from Scene 30 of the KAIST
%   dataset.
%
% ## Notes
% - In order to compare results across images, the 'mi_bands' field of the
%   `options` input argument of 'evaluateSpectral()' should be set globally
%   for the dataset. In other words, it should not be a field of any of the
%   values of `dataset_params.evaluation.custom_spectral`. This will affect
%   how 'mergeSpectralTables()' aggregates evaluation results across
%   images.
%
% See also evaluateRGB, evaluateSpectral, mergeSpectralTables

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 27, 2018

nargoutchk(1, 1);
narginchk(1, 1);

if strcmp(name, 'kodak')
    dataset_params.raw_images_wildcard = [];
    dataset_params.raw_images_variable = [];
    dataset_params.rgb_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180726_Demosaicking_Kodak/PNG_Richard W Franzen/*.png';
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = [];
    dataset_params.spectral_images_variable = [];
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.color_map = [];
    dataset_params.wavelengths = [];
    dataset_params.evaluation = struct(...
        'global_rgb', struct,...
        'custom_rgb', struct(...
            'kodim01', struct('error_map', true))...
        );
elseif strcmp(name, 'kaist')
    dataset_params.raw_images_wildcard = [];
    dataset_params.raw_images_variable = [];
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/*.exr';
    dataset_params.spectral_images_variable = [];
    dataset_params.spectral_reflectances = true;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/SonyColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/wavelengths.mat';
elseif strcmp(name, '20180817_TestSpectralDataset')
    dataset_params.raw_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/*raw.mat';
    dataset_params.raw_images_variable = 'I_raw';
    dataset_params.rgb_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/*3.mat';
    dataset_params.rgb_images_variable = 'I_3';
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/*hyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/BimaterialImagesData.mat';
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/NikonD5100ColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/BimaterialImagesData.mat';
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'error_map', true,...
            'mi_bands', [1, 23],...
            'bands_diff', [1, 23]...
            ),...
        'custom_spectral', struct(...
            'lacelike_0016', struct(...
                'radiance', [181, 75, 11, 11; 195, 397, 11, 11],...
                'scanlines', [124, 214, 273, 380; 266, 371, 50, 435]...
            ),...
            'wrinkled_0143_grey', struct(...
                'radiance', [129, 121, 15, 15; 155, 317, 5, 5],...
                'scanlines', [66, 172, 270, 186]...
            )...
        )...
    );
elseif strcmp(name, 'kaist-crop')
    dataset_params.raw_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181019_KAIST_ForPDFAResearchDayPoster/dataset/*cropped_raw.mat';
    dataset_params.raw_images_variable = 'I_raw';
    dataset_params.rgb_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181019_KAIST_ForPDFAResearchDayPoster/dataset/*cropped_rgb.mat';
    dataset_params.rgb_images_variable = 'I_rgb';
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181019_KAIST_ForPDFAResearchDayPoster/dataset/*cropped_latent.mat';
    dataset_params.spectral_images_variable = 'I_latent';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181019_KAIST_ForPDFAResearchDayPoster/dataset/CIE1931ColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181019_KAIST_ForPDFAResearchDayPoster/dataset/wavelengths.mat';
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'error_map', true,...
            'mi_bands', [4, 24],... % 450 nm and 650 nm
            'bands_diff', [4, 24]...
            ),...
        'custom_spectral', struct(...
            'scene30_cropped_latent', struct(...
                'radiance', [ % Colour checker patches
                    117, 112, 51, 51;
                    225, 117, 51, 51;
                    331, 115, 51, 51;
                    435, 116, 51, 51;
                    542, 123, 51, 51;
                    649, 121, 51, 51;
                    112, 217, 51, 51;
                    220, 221, 51, 51;
                    328, 224, 51, 51;
                    434, 227, 51, 51;
                    540, 228, 51, 51;
                    648, 229, 51, 51;
                    109, 329, 51, 51;
                    217, 330, 51, 51;
                    325, 329, 51, 51;
                    431, 330, 51, 51;
                    540, 332, 51, 51;
                    645, 334, 51, 51;
                    111, 430, 51, 51;
                    217, 435, 51, 51;
                    323, 435, 51, 51;
                    431, 436, 51, 51;
                    541, 435, 51, 51;
                    643, 434, 51, 51
                ]...
            )...
        )...
    );
else
    error('Unrecognized dataset name.');
end

end

