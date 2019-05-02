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
%
%   - 'raw_images_wildcard': A wildcard for 'ls()' to find the input RAW
%     images of the dataset. If empty, the RAW images are to be generated
%     from RGB images.
%   - 'raw_images_variable': The variable name to use for loading RAW
%     images from '.mat' files, where applicable.
%
%   - 'rgb_images_wildcard': A wildcard for 'ls()' to find the true RGB
%     images of the dataset. If empty, the RGB images are to be generated
%     from spectral images.
%   - 'rgb_images_variable': The variable name to use for loading RGB
%     images from '.mat' files, where applicable.
%
%   - 'spectral_images_wildcard': A wildcard for 'ls()' to find the true
%     spectral images of the dataset. If empty, the dataset lacks spectral
%     images.
%   - 'spectral_images_variable': The variable name to use for loading
%     spectral images from '.mat' files, where applicable.
%
%   - 'spectral_reflectances': A Boolean field, required only if the
%     dataset has spectral images. If 'spectral_reflectances' is `true`,
%     then the spectral images are reflectance images. Otherwise, the
%     spectral images are radiance images. Reflectance images need to be
%     preprocessed by multiplying them by the spectral power distribution
%     of an illuminant.
%
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
%   - 'is_aberrated': A logical scalar indicating if the images in the
%     dataset are affected by dispersion. If `true`, the models of
%     dispersion will not be used when generating "ground truth" RGB and/or
%     RAW images. Moreover, the number of evaluations for dispersion-aware
%     image estimation algorithms will be doubled in cases where the
%     dataset includes models of dispersion: The dataset images will still
%     be used to evaluate images estimated with dispersion correction, even
%     though the comparison is not strictly valid. The dataset images will
%     also be used to evaluated the versions of the estimated images with
%     dispersion re-applied.
%
%   - 'color_map': The filename and path of the colour space conversion
%     '.mat' data file. This data file is used to convert spectral
%     information to RGB (or other input colour space), and has the form
%     described in the documentation of 'CorrectByHyperspectralADMM.m'.
%   - 'fix_bands': An optional field containing a logical scalar. If this
%     field is present and `true`, then spectral image estimation must
%     estimate all spectral bands loaded from the data file referred to by
%     'color_map', rather than estimating a different set of spectral
%     bands, determined by the parameters in 'SetFixedParameters.m'. This
%     setting is only used when 'color_map' describes spectral to colour
%     conversion, not conversion between colour channels.
%
%   - 'wavelengths': The filename and path of a '.mat' file containing a
%     variable, `bands`, storing the wavelengths which determine the
%     spectral sampling of spectral images. This field is required only if
%     the dataset has spectral images.
%
%   - 'patch_size': A two-element vector containing the height and width,
%     respectively, of the image patches to use for patch-wise image estimation.
%     Patch dimensions should be even integers to simplify handling of colour
%     filter array images.
%   - 'padding': A scalar containing the width of the border around each image
%     patch. During patch-wise image estimation, the actual patch size will be
%     `patch_size + (2 * padding)` (except as cropped to the image borders). The
%     padding around each patch will be discarded when combining the results
%     from each patch into the final image. Larger values of 'padding' result in
%     increasing overlap between the areas used to estimate each image patch.
%     'padding' should be an even integer to simplify handling of colour filter
%     array images.
%
%   - 'params_patches': A structure, where the value of each field is a matrix
%     of image patches to be used for selecting regularization weights (or other
%     parameters for image estimation algorithms). The columns of each matrix
%     contain the center pixel x-coordinates, and center pixel y-coordinates,
%     respectively, of the image patches. The fieldnames of 'params_patches' are
%     the filenames (excluding file extensions and common suffixes) of the image
%     files in the dataset. The same patches should be used for both colour and
%     spectral versions of images. Therefore, if applicable, both versions of
%     each image must have filenames that map to the same field of
%     'params_patches'.
%
%     This field is optional and, if it is present, it need not have fields for
%     every image. Random patches should be selected when there are no patches
%     provided for an image by the dataset description.
%
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
%       should not be included. The fields 'plot_*', 'radiance',
%       'scanlines', and 'reference_patch' should also be omitted.
%       'global_spectral' is needed only for datasets with spectral images.
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
%   - 'choi_rgb_wildcard': A wildcard for 'ls()' to find RGB images
%     estimated by the method of Choi et al. 2017. This field is optional.
%   - 'choi_spectral_wildcard': A wildcard for 'ls()' to find spectral
%     images estimated by the method of Choi et al. 2017. This field is
%     optional.
%
% ## Recognized (high-quality) datasets
% - 'kodak': The Kodak Lossless True Color Image Suite dataset, often
%   used to evaluate demosaicking algorithms.
%   - Source: http://r0k.us/graphics/kodak/
%   - Maintainer: Richard W Franzen
% - 'kaist-crop': The KAIST Dataset of Hyperspectral Reflectance Images (cropped
%   versions)
%   - Source: http://vclab.kaist.ac.kr/siggraphasia2017p1/kaistdataset.html
%   - Reference:
%     Choi, I., Jeon, D. S., Nam, G., Gutierrez, D., & Kim, M. H. (2017).
%       "High-Quality Hyperspectral Reconstruction Using a Spectral Prior."
%       ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6),
%       218:1–13. doi:10.1145/3130800.3130810
%
% There are some datasets defined in the code used for preliminary testing
% only, such as:
% - 'kaist-crop-2': The color checker chart from Scene 30 of the KAIST
%   dataset.
% - 'choi-test': A dataset used to test the conversion from reflectance
%   images to radiance images, and to test modifications to the method of
%   Choi et al. 2017 (cited above).
%
% ## Notes
% - In order to compare results across images, the 'mi_bands' field of the
%   `options` input argument of 'evaluateSpectral()' should be set globally
%   for the dataset. In other words, it should not be a field of any of the
%   values of `dataset_params.evaluation.custom_spectral`. This will affect
%   how 'mergeSpectralTables()' aggregates evaluation results across
%   images.
%
% ## References
% - Choi, I., Jeon, D. S., Gutierrez, D., & Kim, M. H. (2017).
%   "High-Quality Hyperspectral Reconstruction Using a Spectral Prior." ACM
%   Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6), 218:1-13.
%   10.1145/3130800.3130810
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
    dataset_params.is_aberrated = true;
    dataset_params.color_map = [];
    dataset_params.wavelengths = [];
    dataset_params.patch_size = [30 30];
    dataset_params.padding = 10;
    dataset_params.evaluation = struct(...
        'global_rgb', struct,...
        'custom_rgb', struct(...
            'kodim01', struct('error_map', true))...
        );
    
elseif strcmp(name, 'kaist-crop')
    dataset_params.raw_images_wildcard = [];
    dataset_params.raw_images_variable = [];
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/cropped/*_reflectance.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = true;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.is_aberrated = true;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/NikonD5100ColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/wavelengths.mat';
    dataset_params.patch_size = [128 128];
    dataset_params.padding = 16;
    % One patch is on the ColorChecker chart, to the top right of the white
    % patch, and the other is in another location with strong edges.
    dataset_params.params_patches = struct(...
        'four', [2231, 1400],... %[576, 1736; 2231, 1400],...
        'ten', [1998, 614],... %[587, 992; 1998, 614],...
        'twentyOne', [1583, 1080],... %[304, 1457; 1583, 1080],...
        'twentySix', [2174, 708]... %[182, 1241; 2174, 708]...
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [1, 31],... % 420 nm and 720 nm
            'bands_diff', [1, 31]...
        ),...
        'custom_spectral', struct(...
            'four', struct(... %'reference_patch', [350, 1737, 101, 101],...
                'radiance', [
                    1937, 704, 15, 15;
                    2764, 456, 15, 15;
                    1979, 155, 15, 15;
                    1154, 398, 15, 15;
                    1144, 1400, 15, 15;
                    1898, 1755, 15, 15;
                    2703, 1470, 15, 15
                ],...
                'scanlines', [1941, 648, 1890, 1817]...
            ),...
            'ten', struct(... %'reference_patch', [234, 976, 101, 101],...
                'radiance', [ % ColorChecker
                    539, 724, 15, 15;
                    644, 726, 15, 15;
                    748, 727, 15, 15;
                    854, 727, 15, 15;
                    963, 727, 15, 15;
                    1069, 728, 15, 15;
                    532, 833, 15, 15;
                    641, 833, 15, 15;
                    751, 834, 15, 15;
                    859, 834, 15, 15;
                    963, 834, 15, 15;
                    1072, 835, 15, 15;
                    528, 941, 15, 15;
                    642, 941, 15, 15;
                    752, 941, 15, 15;
                    857, 941, 15, 15;
                    965, 941, 15, 15;
                    1070, 941, 15, 15;
                    532, 1045, 15, 15;
                    640, 1046, 15, 15;
                    749, 1046, 15, 15;
                    855, 1046, 15, 15;
                    963, 1046, 15, 15;
                    1069, 1047, 15, 15
                ],...
                'scanlines', [ % Disc
                    1630, 186, 1817, 1003;
                    1912, 84, 2308, 795
                ]...
            ),...
            'twentyOne', struct(... %'reference_patch', [669, 330, 101, 101],...
                'radiance', [
                    1910, 397, 15, 15; % Red on helmet
                    1998, 490, 15, 15; % Gold on mask
                    2591, 382, 15, 15; % Blue hat
                    2573, 703, 15, 15; % Blue eyes
                    2459, 865, 15, 15; % Yellow beak
                    2770, 1226, 15, 15; % White body
                    2673, 1497, 11, 11 % Green "grass"
                ],...
                'scanlines', [ % Rows of ColorChecker chart
                    117, 987, 1250, 1008;
                    110, 1168, 1249, 1189;
                    110, 1360, 1244, 1372;
                    104, 1541, 1240, 1557
                ]...
            ),...
            'twentySix', struct(... %'reference_patch', [339, 629, 101, 101],...
                'radiance', [
                    1103, 1089, 15, 15; % Left leaf
                    2123, 251, 15, 15; % Top leaf
                    2103, 552, 15, 15; % Centre of top blossom
                    2281, 540, 9, 9; % Petal of top blossom
                    1457, 631, 9, 9; % Petal of top left blossom
                    1304, 892, 9, 9; % Petal of left blossom
                    1593, 1002, 9, 9; % Petal of centre blossom
                    2090, 883, 9, 9 % Petal of right blossom
                ],...
                'scanlines', [1225, 1306, 2130, 229]... % Path through bouquet
            )...
        )...
    );
    dataset_params.choi_rgb_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190130_KAIST_crop/ChoiEtAl2017/*_recon_choiOutConverted_rgb.mat';
    dataset_params.choi_spectral_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190130_KAIST_crop/ChoiEtAl2017/*_recon_choiOutConverted_latent.mat';
    
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
    dataset_params.is_aberrated = false;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/NikonD5100ColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20180817_TestSpectralDataset/dataset/BimaterialImagesData.mat';
    dataset_params.patch_size = [64 64];
    dataset_params.padding = 16;
    dataset_params.params_patches = struct(...
        'lacelike_0016', [181, 75; 195, 397],...
        'wrinkled_0143_grey', [129, 121]...
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
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

elseif strcmp(name, 'kaist-crop-2')
    dataset_params.raw_images_wildcard = [];
    dataset_params.raw_images_variable = [];
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181019_KAIST_ForPDFAResearchDayPoster/dataset/*cropped_latent.mat';
    dataset_params.spectral_images_variable = 'I_latent';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.is_aberrated = true;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/NikonD5100ColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_highQualityHyperspectralReconstructionUsingASpectralPrior_LCTFSystem/wavelengths.mat';
    dataset_params.patch_size = [256 256];
    dataset_params.padding = 16;
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
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

elseif strcmp(name, 'choi-test')
    dataset_params.raw_images_wildcard = [];
    dataset_params.raw_images_variable = [];
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181127_TestingChoiEtAl2017/Original/colorIDs_hyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.is_aberrated = true;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181127_TestingChoiEtAl2017/NikonD5100ColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181127_TestingChoiEtAl2017/Original/BimaterialImagesData.mat';
    dataset_params.patch_size = [64 64];
    dataset_params.padding = 16;
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [4, 24],...
            'bands_diff', [4, 24]...
            ),...
        'custom_spectral', struct(...
            'colorIDs_hyper', struct(...
                'radiance', [ % Colour checker patches
                10    10    11    11
                30    10    11    11
                50    10    11    11
                70    10    11    11
                90    10    11    11
                110    10    11    11
                10    30    11    11
                30    30    11    11
                50    30    11    11
                70    30    11    11
                90    30    11    11
                110    30    11    11
                10    50    11    11
                30    50    11    11
                50    50    11    11
                70    50    11    11
                90    50    11    11
                110    50    11    11
                10    70    11    11
                30    70    11    11
                50    70    11    11
                70    70    11    11
                90    70    11    11
                110    70    11    11
                ]...
            )...
        )...
    );
    dataset_params.choi_rgb_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181127_TestingChoiEtAl2017/ChoiEtAl2017_OutputConverted/recon_choiOutConverted_rgb.mat';
    dataset_params.choi_spectral_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20181127_TestingChoiEtAl2017/ChoiEtAl2017_OutputConverted/recon_choiOutConverted_latent.mat';
    
elseif strcmp(name, '20190107_DiskPattern_rawFromSpectral')
    dataset_params.raw_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/disks47cm_raw.mat';
    dataset_params.raw_images_variable = 'I_raw';
    dataset_params.rgb_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/disks47cm_d3.mat';
    dataset_params.rgb_images_variable = 'I_3';
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/disks47cm_dHyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion/RAWDiskDispersionResults_RGB_spline_fromReference.mat';
    dataset_params.dispersion_rgb_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion/RAWDiskDispersionResults_RGB_spline_fromNonReference.mat';
    dataset_params.dispersion_spectral_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion/RAWDiskDispersionResults_spectral_spline_fromNonReference.mat';
    dataset_params.is_aberrated = true;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/sensor.mat';
    dataset_params.fix_bands = true;
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/sensor.mat';
    dataset_params.patch_size = [256 256];
    dataset_params.padding = 16;
    dataset_params.params_patches = struct(...
        'disks47cm_dHyper', [2274, 237; 339, 304; 350 1833; 2274, 1896; 1198, 1007]...
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [1, 7],...
            'bands_diff', [1, 7]...
            ),...
        'custom_spectral', struct(...
            'disks47cm_dHyper', struct(...
                'radiance', [
                    2280, 220, 17, 17;
                    233, 304, 17, 17;
                    350 1833, 17, 17;
                    2322, 1905, 17, 17;
                    1198, 1007, 17, 17;
                    2245, 254, 17, 17;
                    260, 330, 17, 17;
                    373, 1805, 17, 17;
                    2290, 1871, 17, 17;
                    1169, 1036, 17, 17
                ],...
                'scanlines', [233, 304, 2322, 1905; 350, 1833, 2280, 220]...
            )...
        )...
    );

elseif strcmp(name, '20190107_DiskPattern_rawCaptured')
    dataset_params.raw_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/preprocessed_images/exposure_blended/disks47cm_unfiltered.mat';
    dataset_params.raw_images_variable = 'I_raw';
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/disks47cm_dHyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion/RAWDiskDispersionResults_RGB_spline_fromReference.mat';
    dataset_params.dispersion_rgb_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion/RAWDiskDispersionResults_RGB_spline_fromNonReference.mat';
    dataset_params.dispersion_spectral_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/dispersion/RAWDiskDispersionResults_spectral_spline_fromNonReference.mat';
    dataset_params.is_aberrated = true;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/SonyColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190107_DiskPattern_real/channel_scaling/sensor.mat';
    dataset_params.patch_size = [256 256];
    dataset_params.padding = 16;
    dataset_params.params_patches = struct(...
        'disks47cm_dHyper', [2274, 237; 339, 304; 350 1833; 2274, 1896; 1198, 1007]...
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [1, 7],...
            'bands_diff', [1, 7]...
            ),...
        'custom_spectral', struct(...
            'disks47cm_dHyper', struct(...
                'reference_patch', [1224, 1024, 17, 17],...
                'radiance', [
                    2280, 220, 17, 17;
                    233, 304, 17, 17;
                    350 1833, 17, 17;
                    2322, 1905, 17, 17;
                    1198, 1007, 17, 17;
                    2245, 254, 17, 17;
                    260, 330, 17, 17;
                    373, 1805, 17, 17;
                    2290, 1871, 17, 17;
                    1169, 1036, 17, 17
                ],...
                'scanlines', [233, 304, 2322, 1905; 350, 1833, 2280, 220]...
            )...
        )...
    );

elseif strcmp(name, '20190208_ComputarLens_rawCaptured_ignoreDispersion')
    dataset_params.raw_images_wildcard = 'C:\Users\GraphicsLab\Documents\llanos\Results\20190208_ComputarLens\dataset\exposure_blending\*_unfiltered.mat';
    dataset_params.raw_images_variable = 'I_raw';
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = 'C:\Users\GraphicsLab\Documents\llanos\Results\20190208_ComputarLens\dataset\channel_scaling\*_dHyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = [];
    dataset_params.dispersion_rgb_reverse = [];
    dataset_params.dispersion_spectral_reverse = [];
    dataset_params.is_aberrated = true;
    dataset_params.color_map = 'C:\Users\GraphicsLab\Documents\llanos\Results\20190208_ComputarLens\dataset\SonyColorMapData.mat';
    dataset_params.wavelengths = 'C:\Users\GraphicsLab\Documents\llanos\Results\20190208_ComputarLens\dataset\channel_scaling\sensor.mat';
    dataset_params.patch_size = [128 128];
    dataset_params.padding = 16;
    dataset_params.params_patches = struct(...
        'd1_colorChecker30cm', [146, 141],... % "colorchecker CLASSIC" text
        'd1_disks32cmV1', [2164, 228],... % Dot near top right
        'd1_disks32cmV2', [2176, 225],... % Dot near top right
        'd2_book', [1180, 205],... % "Wiley"
        'd2_colorChecker30cm', [336, 1470],... % "X" symbol beside "X-rite"
        'd2_glass', [449, 883],... % Dot seen through the crystal ball
        'd2_ship', [2234, 1132]... % Noise pattern
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [1, 7],...
            'bands_diff', [1, 7]...
            ),...
        'custom_spectral', struct(...
            'd1_colorChecker30cm', struct(... % The white square of the ColorChecker is overexposed, unfortunately
                'reference_patch', [695, 1303, 15, 15],... % Square right of the white ColorChecker square
                'radiance', [ % ColorChecker patches
                    296, 344, 15, 15;
                    641, 353, 15, 15;
                    980, 362, 15, 15;
                    1337, 371, 15, 15;
                    1697, 356, 15, 15;
                    2045, 365, 15, 15;
                    332, 688, 15, 15;
                    664, 691, 15, 15;
                    989, 700, 15, 15;
                    1341, 696, 15, 15;
                    1679, 706, 15, 15;
                    2027, 718, 15, 15;
                    371, 1003, 15, 15;
                    678, 1007, 15, 15;
                    1004, 1015, 15, 15;
                    1326, 1020, 15, 15;
                    1660, 1020, 15, 15;
                    1994, 1021, 15, 15;
                    395, 1308, 15, 15;
                    695, 1303, 15, 15; % Also the reference patch!
                    1014, 1316, 15, 15;
                    1325, 1320, 15, 15;
                    1641, 1327, 15, 15;
                    1970, 1330, 15, 15
                ],...
                'scanlines', [
                    296, 344, 2045, 365; % ColorChecker rows
                    332, 688, 2027, 718;
                    371, 1003, 1994, 1021;
                    395, 1308, 1970, 1330;
                    1620, 1510, 1980, 1515 % Across the 'mm' ruler gradations
                ]...
            ),...
            'd1_disks32cmV1', struct(...
                'reference_patch', [1223, 989, 15, 15],... % White in the image centre
                'radiance', [
                    1223, 989, 15, 15; % Also the reference patch!
                    1185, 1024, 15, 15; % Dot near image centre
                    2286, 112, 15, 15; % White near top right
                    2327, 73, 15, 15; % Dot closer to the corner
                    1817, 490, 15, 15; % White 1/2 way to top right
                    1851, 528, 15, 15; % Dot nearby
                    144, 175, 15, 15; % White near top left
                    180, 212, 15, 15; % Dot closer to the centre
                    577, 474, 15, 15; % White 1/2 way to top left
                    616, 511, 15, 15; % Dot closer to the centre
                    94, 1945, 15, 15; % White near bottom left
                    124, 1915, 15, 15; % Dot closer to the centre
                    603, 1517, 15, 15; % White 1/2 way to bottom left
                    635, 1487, 15, 15; % Dot closer to the centre
                    2286, 1944, 15, 15; % White near bottom right
                    2318, 1979, 15, 15; % Dot closer to the corner
                    1838, 1548, 15, 15; % White 1/2 way to bottom right
                    1802, 1517, 15, 15 % Dot closer to the centre
                ],...
                'scanlines', [
                    26, 56, 2406, 73; % Top left to top right
                    26, 56, 125, 1912; % Top left to bottom left
                    125, 1912, 2252, 1972; % Bottom left to bottom right
                    2252, 1972, 2406, 73; % Bottom right to top right
                    26, 56, 2252, 1972; % Top left to bottom right
                    125, 1912, 2406, 73 % Bottom left to top right
                ]...
            ),...
            'd1_disks32cmV2', struct(...
                'reference_patch', [1229, 1049, 15, 15],... % White in the image centre
                'radiance', [
                    1229, 1049, 15, 15; % Also the reference patch!
                    1193, 1016, 15, 15; % Dot near image centre
                    2293, 188, 15, 15; % White near top right
                    2259, 147, 15, 15; % Dot closer to the corner
                    1821, 488, 15, 15; % White 1/2 way to top right
                    1858, 526, 15, 15; % Dot nearby
                    152, 171, 15, 15; % White near top left
                    111, 131, 15, 15; % Dot closer to the corner
                    590, 548, 15, 15; % White 1/2 way to top left
                    628, 582, 15, 15; % Dot closer to the centre
                    208, 1893, 15, 15; % White near bottom left
                    237, 1865, 15, 15; % Dot closer to the centre
                    627, 1492, 15, 15; % White 1/2 way to bottom left
                    596, 1521, 15, 15; % Dot closer to the corner
                    2329, 1910, 15, 15; % White near bottom right
                    2361, 1942, 15, 15; % Dot closer to the corner
                    1830, 1529, 15, 15; % White 1/2 way to bottom right
                    1860, 1559, 15, 15 % Dot closer to the corner
                ],...
                'scanlines', [
                    32, 51, 2426, 63; % Top left to top right
                    32, 51, 180, 1920; % Top left to bottom left
                    180, 1920, 2222, 1994; % Bottom left to bottom right
                    2222, 1994, 2426, 63; % Bottom right to top right
                    32, 51, 2222, 1994; % Top left to bottom right
                    180, 1920, 2426, 63 % Bottom left to top right
                ]...
            ),...
            'd2_book', struct(...
                'reference_patch', [928, 940, 9, 9],... % White text on book
                'radiance', [
                    911, 107, 15, 15; % Orange on book near top of image
                    398, 1108, 15, 15; % Blue around book title
                    608, 1471, 15, 15; % Black above "Second Edition"
                    1118, 1763, 15, 15; % Red on largest doll
                    944, 1688, 7, 7; % Purple on second-largest doll
                    2216, 1150, 75, 75 % Large patch of wood
                ],...
                'scanlines', [
                    914, 1034, 946, 1034; % Going across the 'I' in "Image"
                ]...
            ),...
            'd2_colorChecker30cm', struct(...
                'reference_patch', [749, 1306, 15, 15],... % Square right of the white ColorChecker square
                'radiance', [ % ColorChecker patches
                    362, 338, 15, 15;
                    703, 342, 15, 15;
                    1043, 347, 15, 15;
                    1396, 345, 15, 15;
                    1751, 350, 15, 15;
                    2108, 359, 15, 15;
                    383, 676, 15, 15;
                    716, 688, 15, 15;
                    1046, 687, 15, 15;
                    1398, 694, 15, 15;
                    1734, 696, 15, 15;
                    2072, 703, 15, 15;
                    419, 994, 15, 15;
                    735, 1004, 15, 15;
                    1058, 1005, 15, 15;
                    1387, 1017, 15, 15;
                    1719, 1019, 15, 15;
                    2052, 1023, 15, 15;
                    446, 1297, 15, 15;
                    749, 1306, 15, 15; % Also the reference patch!
                    1061, 1315, 15, 15;
                    1375, 1321, 15, 15;
                    1704, 1324, 15, 15;
                    2033, 1336, 15, 15
                ],...
                'scanlines', [ % ColorChecker rows
                    362, 338, 2108, 359;
                    383, 676, 2072, 703;
                    419, 994, 2052, 1023;
                    446, 1297, 2033, 1336;
                    1673, 1512, 2034, 1515 % Across the 'mm' ruler gradations
                ]...
            ),...
            'd2_glass', struct(...
                'reference_patch', [1243, 1042, 15, 15],... % White patch of background around image centre
                'radiance', [
                    488, 170, 15, 15; % Grey substrate
                    1574, 1510, 15, 15; % White patch beside star
                    1724, 1454, 15, 15; % White patch seen through star
                    1786, 1585, 15, 15; % White patch seen through star, in deep shadow of star
                    1857, 1585, 15, 15; % White patch seen through star, in light shadow of star
                    2039, 1581, 15, 15; % White patch in light shadow of star
                    2208, 1585, 15, 15; % White patch in deep shadow of star
                ],...
                'scanlines', [
                    1281, 1234, 1398, 1328; % Line through a rainbow passing through a dot
                    1902, 1582, 1942, 1655; % Line through a rainbow seen through the star
                ]...
            ),...
            'd2_ship', struct(...
                'reference_patch', [1888, 1577, 15, 15],... % White background of noise pattern
                'radiance', [
                    1818, 1869, 15, 15; % White background seen through laminate
                    1810, 1985, 15, 15; % White background below laminate
                ],...
                'scanlines', [
                    963, 1013, 1394, 1235; % Cut through dorsal fin of dolphin
                    1312, 909, 1581, 324 % Cut through dolphin and rock above it in the image
                ]...
            )...
        )...
    );

elseif strcmp(name, '20190208_ComputarLens_rawCaptured_dispersion')
    dataset_params.raw_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/exposure_blending/*_unfiltered.mat';
    dataset_params.raw_images_variable = 'I_raw';
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/channel_scaling/*_dHyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/rgb/full_image/RAWDiskDispersionResults_RGB_polynomial_fromReference.mat';
    dataset_params.dispersion_rgb_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/rgb/full_image/RAWDiskDispersionResults_RGB_polynomial_fromNonReference.mat';
    dataset_params.dispersion_spectral_reverse = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dispersion/spectral/full_image/RAWDiskDispersionResults_spectral_polynomial_fromNonReference.mat';
    dataset_params.is_aberrated = true;
    dataset_params.color_map = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/SonyColorMapData.mat';
    dataset_params.wavelengths = '/home/llanos/GoogleDrive/ThesisResearch/Results/20190208_ComputarLens/dataset/channel_scaling/sensor.mat';
    dataset_params.patch_size = [128 128];
    dataset_params.padding = 16;
    dataset_params.params_patches = struct(...
        'd1_colorChecker30cm', [146, 141],... % "colorchecker CLASSIC" text
        'd1_disks32cmV1', [2164, 228],... % Dot near top right
        'd1_disks32cmV2', [2176, 225],... % Dot near top right
        'd2_book', [1180, 205],... % "Wiley"
        'd2_colorChecker30cm', [336, 1470],... % "X" symbol beside "X-rite"
        'd2_glass', [449, 883],... % Dot seen through the crystal ball
        'd2_ship', [2234, 1132]... % Noise pattern
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [1, 7],...
            'bands_diff', [1, 7]...
            ),...
        'custom_spectral', struct(...
            'd1_colorChecker30cm', struct(... % The white square of the ColorChecker is overexposed, unfortunately
                'reference_patch', [695, 1303, 15, 15],... % Square right of the white ColorChecker square
                'radiance', [ % ColorChecker patches
                    296, 344, 15, 15;
                    641, 353, 15, 15;
                    980, 362, 15, 15;
                    1337, 371, 15, 15;
                    1697, 356, 15, 15;
                    2045, 365, 15, 15;
                    332, 688, 15, 15;
                    664, 691, 15, 15;
                    989, 700, 15, 15;
                    1341, 696, 15, 15;
                    1679, 706, 15, 15;
                    2027, 718, 15, 15;
                    371, 1003, 15, 15;
                    678, 1007, 15, 15;
                    1004, 1015, 15, 15;
                    1326, 1020, 15, 15;
                    1660, 1020, 15, 15;
                    1994, 1021, 15, 15;
                    395, 1308, 15, 15;
                    695, 1303, 15, 15; % Also the reference patch!
                    1014, 1316, 15, 15;
                    1325, 1320, 15, 15;
                    1641, 1327, 15, 15;
                    1970, 1330, 15, 15
                ],...
                'scanlines', [
                    296, 344, 2045, 365; % ColorChecker rows
                    332, 688, 2027, 718;
                    371, 1003, 1994, 1021;
                    395, 1308, 1970, 1330;
                    1620, 1510, 1980, 1515 % Across the 'mm' ruler gradations
                ]...
            ),...
            'd1_disks32cmV1', struct(...
                'reference_patch', [1223, 989, 15, 15],... % White in the image centre
                'radiance', [
                    1223, 989, 15, 15; % Also the reference patch!
                    1185, 1024, 15, 15; % Dot near image centre
                    2286, 112, 15, 15; % White near top right
                    2327, 73, 15, 15; % Dot closer to the corner
                    1817, 490, 15, 15; % White 1/2 way to top right
                    1851, 528, 15, 15; % Dot nearby
                    144, 175, 15, 15; % White near top left
                    180, 212, 15, 15; % Dot closer to the centre
                    577, 474, 15, 15; % White 1/2 way to top left
                    616, 511, 15, 15; % Dot closer to the centre
                    94, 1945, 15, 15; % White near bottom left
                    124, 1915, 15, 15; % Dot closer to the centre
                    603, 1517, 15, 15; % White 1/2 way to bottom left
                    635, 1487, 15, 15; % Dot closer to the centre
                    2286, 1944, 15, 15; % White near bottom right
                    2318, 1979, 15, 15; % Dot closer to the corner
                    1838, 1548, 15, 15; % White 1/2 way to bottom right
                    1802, 1517, 15, 15 % Dot closer to the centre
                ],...
                'scanlines', [
                    26, 56, 2406, 73; % Top left to top right
                    26, 56, 125, 1912; % Top left to bottom left
                    125, 1912, 2252, 1972; % Bottom left to bottom right
                    2252, 1972, 2406, 73; % Bottom right to top right
                    26, 56, 2252, 1972; % Top left to bottom right
                    125, 1912, 2406, 73 % Bottom left to top right
                ]...
            ),...
            'd1_disks32cmV2', struct(...
                'reference_patch', [1229, 1049, 15, 15],... % White in the image centre
                'radiance', [
                    1229, 1049, 15, 15; % Also the reference patch!
                    1193, 1016, 15, 15; % Dot near image centre
                    2293, 188, 15, 15; % White near top right
                    2259, 147, 15, 15; % Dot closer to the corner
                    1821, 488, 15, 15; % White 1/2 way to top right
                    1858, 526, 15, 15; % Dot nearby
                    152, 171, 15, 15; % White near top left
                    111, 131, 15, 15; % Dot closer to the corner
                    590, 548, 15, 15; % White 1/2 way to top left
                    628, 582, 15, 15; % Dot closer to the centre
                    208, 1893, 15, 15; % White near bottom left
                    237, 1865, 15, 15; % Dot closer to the centre
                    627, 1492, 15, 15; % White 1/2 way to bottom left
                    596, 1521, 15, 15; % Dot closer to the corner
                    2329, 1910, 15, 15; % White near bottom right
                    2361, 1942, 15, 15; % Dot closer to the corner
                    1830, 1529, 15, 15; % White 1/2 way to bottom right
                    1860, 1559, 15, 15 % Dot closer to the corner
                ],...
                'scanlines', [
                    32, 51, 2426, 63; % Top left to top right
                    32, 51, 180, 1920; % Top left to bottom left
                    180, 1920, 2222, 1994; % Bottom left to bottom right
                    2222, 1994, 2426, 63; % Bottom right to top right
                    32, 51, 2222, 1994; % Top left to bottom right
                    180, 1920, 2426, 63 % Bottom left to top right
                ]...
            ),...
            'd2_book', struct(...
                'reference_patch', [928, 940, 9, 9],... % White text on book
                'radiance', [
                    911, 107, 15, 15; % Orange on book near top of image
                    398, 1108, 15, 15; % Blue around book title
                    608, 1471, 15, 15; % Black above "Second Edition"
                    1118, 1763, 15, 15; % Red on largest doll
                    944, 1688, 7, 7; % Purple on second-largest doll
                    2216, 1150, 75, 75 % Large patch of wood
                ],...
                'scanlines', [
                    914, 1034, 946, 1034; % Going across the 'I' in "Image"
                ]...
            ),...
            'd2_colorChecker30cm', struct(...
                'reference_patch', [749, 1306, 15, 15],... % Square right of the white ColorChecker square
                'radiance', [ % ColorChecker patches
                    362, 338, 15, 15;
                    703, 342, 15, 15;
                    1043, 347, 15, 15;
                    1396, 345, 15, 15;
                    1751, 350, 15, 15;
                    2108, 359, 15, 15;
                    383, 676, 15, 15;
                    716, 688, 15, 15;
                    1046, 687, 15, 15;
                    1398, 694, 15, 15;
                    1734, 696, 15, 15;
                    2072, 703, 15, 15;
                    419, 994, 15, 15;
                    735, 1004, 15, 15;
                    1058, 1005, 15, 15;
                    1387, 1017, 15, 15;
                    1719, 1019, 15, 15;
                    2052, 1023, 15, 15;
                    446, 1297, 15, 15;
                    749, 1306, 15, 15; % Also the reference patch!
                    1061, 1315, 15, 15;
                    1375, 1321, 15, 15;
                    1704, 1324, 15, 15;
                    2033, 1336, 15, 15
                ],...
                'scanlines', [ % ColorChecker rows
                    362, 338, 2108, 359;
                    383, 676, 2072, 703;
                    419, 994, 2052, 1023;
                    446, 1297, 2033, 1336;
                    1673, 1512, 2034, 1515 % Across the 'mm' ruler gradations
                ]...
            ),...
            'd2_glass', struct(...
                'reference_patch', [1243, 1042, 15, 15],... % White patch of background around image centre
                'radiance', [
                    488, 170, 15, 15; % Grey substrate
                    1574, 1510, 15, 15; % White patch beside star
                    1724, 1454, 15, 15; % White patch seen through star
                    1786, 1585, 15, 15; % White patch seen through star, in deep shadow of star
                    1857, 1585, 15, 15; % White patch seen through star, in light shadow of star
                    2039, 1581, 15, 15; % White patch in light shadow of star
                    2208, 1585, 15, 15; % White patch in deep shadow of star
                ],...
                'scanlines', [
                    1281, 1234, 1398, 1328; % Line through a rainbow passing through a dot
                    1902, 1582, 1942, 1655; % Line through a rainbow seen through the star
                ]...
            ),...
            'd2_ship', struct(...
                'reference_patch', [1888, 1577, 15, 15],... % White background of noise pattern
                'radiance', [
                    1818, 1869, 15, 15; % White background seen through laminate
                    1810, 1985, 15, 15; % White background below laminate
                ],...
                'scanlines', [
                    963, 1013, 1394, 1235; % Cut through dorsal fin of dolphin
                    1312, 909, 1581, 324 % Cut through dolphin and rock above it in the image
                ]...
            )...
        )...
    );

elseif strcmp(name, '20190421_ComputarLens_dHyper_dispersion')
    dataset_params.raw_images_wildcard = [];
    dataset_params.raw_images_variable = [];
    dataset_params.rgb_images_wildcard = [];
    dataset_params.rgb_images_variable = [];
    dataset_params.spectral_images_wildcard = 'C:\Users\GraphicsLab\Documents\llanos\Results\channel_scaling\*_dHyper.mat';
    dataset_params.spectral_images_variable = 'I_hyper';
    dataset_params.spectral_reflectances = false;
    dataset_params.dispersion_rgb_forward = 'C:\Users\GraphicsLab\Documents\llanos\Data\20190208_ComputarLens\dispersion\rgb\full_image\RAWDiskDispersionResults_RGB_polynomial_fromReference.mat';
    dataset_params.dispersion_rgb_reverse = 'C:\Users\GraphicsLab\Documents\llanos\Data\20190208_ComputarLens\dispersion\rgb\full_image\RAWDiskDispersionResults_RGB_polynomial_fromNonReference.mat';
    dataset_params.dispersion_spectral_reverse = 'C:\Users\GraphicsLab\Documents\llanos\Results\dispersion\spectral\polynomial_newCV\RAWDiskDispersionResults_spectral_polynomial_fromNonReference.mat';
    dataset_params.is_aberrated = true;
    dataset_params.color_map = 'C:\Users\GraphicsLab\Documents\llanos\Results\channel_scaling\sensor.mat';
    dataset_params.wavelengths = 'C:\Users\GraphicsLab\Documents\llanos\Results\channel_scaling\sensor.mat';
    dataset_params.patch_size = [128 128];
    dataset_params.padding = 16;
    dataset_params.params_patches = struct(...
        'd1_colorChecker30cm', [146, 141],... % "colorchecker CLASSIC" text
        'd1_disks32cmV1', [2164, 228],... % Dot near top right
        'd1_disks32cmV2', [2176, 225],... % Dot near top right
        'd2_book', [1180, 205],... % "Wiley"
        'd2_colorChecker30cm', [336, 1470],... % "X" symbol beside "X-rite"
        'd2_glass', [449, 883],... % Dot seen through the crystal ball
        'd2_ship', [2234, 1132]... % Noise pattern
    );
    dataset_params.evaluation = struct(...
        'global_rgb', struct('error_map', true),...
        'custom_rgb', struct,...
        'global_spectral', struct(...
            'metric', 'mrae',...
            'error_map', true,...
            'mi_bands', [1, 7],...
            'bands_diff', [1, 7]...
            ),...
        'custom_spectral', struct(...
            'd1_colorChecker30cm', struct(... % The white square of the ColorChecker is overexposed, unfortunately
                'radiance', [ % ColorChecker patches
                    296, 344, 15, 15;
                    641, 353, 15, 15;
                    980, 362, 15, 15;
                    1337, 371, 15, 15;
                    1697, 356, 15, 15;
                    2045, 365, 15, 15;
                    332, 688, 15, 15;
                    664, 691, 15, 15;
                    989, 700, 15, 15;
                    1341, 696, 15, 15;
                    1679, 706, 15, 15;
                    2027, 718, 15, 15;
                    371, 1003, 15, 15;
                    678, 1007, 15, 15;
                    1004, 1015, 15, 15;
                    1326, 1020, 15, 15;
                    1660, 1020, 15, 15;
                    1994, 1021, 15, 15;
                    395, 1308, 15, 15;
                    695, 1303, 15, 15; % Also the reference patch!
                    1014, 1316, 15, 15;
                    1325, 1320, 15, 15;
                    1641, 1327, 15, 15;
                    1970, 1330, 15, 15
                ],...
                'scanlines', [
                    296, 344, 2045, 365; % ColorChecker rows
                    332, 688, 2027, 718;
                    371, 1003, 1994, 1021;
                    395, 1308, 1970, 1330;
                    1620, 1510, 1980, 1515 % Across the 'mm' ruler gradations
                ]...
            ),...
            'd1_disks32cmV1', struct(...
                'radiance', [
                    1223, 989, 15, 15; % Also the reference patch!
                    1185, 1024, 15, 15; % Dot near image centre
                    2286, 112, 15, 15; % White near top right
                    2327, 73, 15, 15; % Dot closer to the corner
                    1817, 490, 15, 15; % White 1/2 way to top right
                    1851, 528, 15, 15; % Dot nearby
                    144, 175, 15, 15; % White near top left
                    180, 212, 15, 15; % Dot closer to the centre
                    577, 474, 15, 15; % White 1/2 way to top left
                    616, 511, 15, 15; % Dot closer to the centre
                    94, 1945, 15, 15; % White near bottom left
                    124, 1915, 15, 15; % Dot closer to the centre
                    603, 1517, 15, 15; % White 1/2 way to bottom left
                    635, 1487, 15, 15; % Dot closer to the centre
                    2286, 1944, 15, 15; % White near bottom right
                    2318, 1979, 15, 15; % Dot closer to the corner
                    1838, 1548, 15, 15; % White 1/2 way to bottom right
                    1802, 1517, 15, 15 % Dot closer to the centre
                ],...
                'scanlines', [
                    26, 56, 2406, 73; % Top left to top right
                    26, 56, 125, 1912; % Top left to bottom left
                    125, 1912, 2252, 1972; % Bottom left to bottom right
                    2252, 1972, 2406, 73; % Bottom right to top right
                    26, 56, 2252, 1972; % Top left to bottom right
                    125, 1912, 2406, 73 % Bottom left to top right
                ]...
            ),...
            'd1_disks32cmV2', struct(...
                'radiance', [
                    1229, 1049, 15, 15; % Also the reference patch!
                    1193, 1016, 15, 15; % Dot near image centre
                    2293, 188, 15, 15; % White near top right
                    2259, 147, 15, 15; % Dot closer to the corner
                    1821, 488, 15, 15; % White 1/2 way to top right
                    1858, 526, 15, 15; % Dot nearby
                    152, 171, 15, 15; % White near top left
                    111, 131, 15, 15; % Dot closer to the corner
                    590, 548, 15, 15; % White 1/2 way to top left
                    628, 582, 15, 15; % Dot closer to the centre
                    208, 1893, 15, 15; % White near bottom left
                    237, 1865, 15, 15; % Dot closer to the centre
                    627, 1492, 15, 15; % White 1/2 way to bottom left
                    596, 1521, 15, 15; % Dot closer to the corner
                    2329, 1910, 15, 15; % White near bottom right
                    2361, 1942, 15, 15; % Dot closer to the corner
                    1830, 1529, 15, 15; % White 1/2 way to bottom right
                    1860, 1559, 15, 15 % Dot closer to the corner
                ],...
                'scanlines', [
                    32, 51, 2426, 63; % Top left to top right
                    32, 51, 180, 1920; % Top left to bottom left
                    180, 1920, 2222, 1994; % Bottom left to bottom right
                    2222, 1994, 2426, 63; % Bottom right to top right
                    32, 51, 2222, 1994; % Top left to bottom right
                    180, 1920, 2426, 63 % Bottom left to top right
                ]...
            ),...
            'd2_book', struct(...
                'radiance', [
                    911, 107, 15, 15; % Orange on book near top of image
                    398, 1108, 15, 15; % Blue around book title
                    608, 1471, 15, 15; % Black above "Second Edition"
                    1118, 1763, 15, 15; % Red on largest doll
                    944, 1688, 7, 7; % Purple on second-largest doll
                    2216, 1150, 75, 75 % Large patch of wood
                ],...
                'scanlines', [
                    914, 1034, 946, 1034; % Going across the 'I' in "Image"
                ]...
            ),...
            'd2_colorChecker30cm', struct(...
                'radiance', [ % ColorChecker patches
                    362, 338, 15, 15;
                    703, 342, 15, 15;
                    1043, 347, 15, 15;
                    1396, 345, 15, 15;
                    1751, 350, 15, 15;
                    2108, 359, 15, 15;
                    383, 676, 15, 15;
                    716, 688, 15, 15;
                    1046, 687, 15, 15;
                    1398, 694, 15, 15;
                    1734, 696, 15, 15;
                    2072, 703, 15, 15;
                    419, 994, 15, 15;
                    735, 1004, 15, 15;
                    1058, 1005, 15, 15;
                    1387, 1017, 15, 15;
                    1719, 1019, 15, 15;
                    2052, 1023, 15, 15;
                    446, 1297, 15, 15;
                    749, 1306, 15, 15; % Also the reference patch!
                    1061, 1315, 15, 15;
                    1375, 1321, 15, 15;
                    1704, 1324, 15, 15;
                    2033, 1336, 15, 15
                ],...
                'scanlines', [ % ColorChecker rows
                    362, 338, 2108, 359;
                    383, 676, 2072, 703;
                    419, 994, 2052, 1023;
                    446, 1297, 2033, 1336;
                    1673, 1512, 2034, 1515 % Across the 'mm' ruler gradations
                ]...
            ),...
            'd2_glass', struct(...
                'radiance', [
                    488, 170, 15, 15; % Grey substrate
                    1574, 1510, 15, 15; % White patch beside star
                    1724, 1454, 15, 15; % White patch seen through star
                    1786, 1585, 15, 15; % White patch seen through star, in deep shadow of star
                    1857, 1585, 15, 15; % White patch seen through star, in light shadow of star
                    2039, 1581, 15, 15; % White patch in light shadow of star
                    2208, 1585, 15, 15; % White patch in deep shadow of star
                ],...
                'scanlines', [
                    1281, 1234, 1398, 1328; % Line through a rainbow passing through a dot
                    1902, 1582, 1942, 1655; % Line through a rainbow seen through the star
                ]...
            ),...
            'd2_ship', struct(...
                'radiance', [
                    1818, 1869, 15, 15; % White background seen through laminate
                    1810, 1985, 15, 15; % White background below laminate
                ],...
                'scanlines', [
                    963, 1013, 1394, 1235; % Cut through dorsal fin of dolphin
                    1312, 909, 1581, 324 % Cut through dolphin and rock above it in the image
                ]...
            )...
        )...
    );

else
    error('Unrecognized dataset name.');
end

end

