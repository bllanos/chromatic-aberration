%% Evaluate demosaicking, spectral reconstruction, and chromatic aberration correction
% For an image of an X-Rite ColorChecker CLASSIC colour chart.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% ### ColorChecker image
% 
% A raw colour filter array image of a ColorChecker chart must be provided. The
% image is expected to have been preprocessed, such as using
% 'PreprocessRAWImages.m', so that it does not need to be linearized after being
% loaded. For an image format file, the image will be loaded with the Image
% Processing Toolbox 'imread()' function. For a '.mat' file, the variable to be
% loaded must be provided in the script parameters.
%
% Along with the captured image, two single-channel label images are required.
% The first is used for colour calibration, vignetting correction, and
% evaluation of colour-based chromatic aberration correction algorithms. For the
% 24 patches of the chart, the corresponding pixels in the label image must have
% values equal to the patch indices. The uniformly-coloured portions of the
% frame surrounding the patches, used to calibrate vignetting, must be labelled
% 25. Lastly, the rest of the image should have a different label (e.g. zero).
%
% The second label image is used for evaluation of spectral-based chromatic
% aberration correction algorithms. It should have all 24 patches labelled with
% their indices, but need not have other pixels labelled. The boundaries of the
% patches in the first label image should be accurate in the reference channel
% of colour-based models of chromatic aberration (e.g. the Green channel). In
% contrast, the boundaries of the patches in the second label image should be
% accurate in the spectral band used as the reference band for spectral models
% of chromatic aberration.
%
% ### Spectral reflectances and conversion to colour
%
% #### Spectral reflectances
% A '.csv' file containing a header row, a first column for wavelength
% values, and remaining columns for relative spectral reflectances of each
% patch in the ColorChecker.
%
% #### CIE tristimulus functions
% A '.mat' file containing a variable 'xyzbar', which can be used as the
% `C` input argument of 'cieSpectralToColor()'.
%
% ### Estimated images
%
% Several filepath wildcards in the script parameters below are used to locate
% the reconstructed RGB and spectral images to be evaluated. These images should
% have been generated using the raw image of the ColorChecker chart as input.
% '.mat' or image format files can be loaded. In the output files describing the
% evaluation results, the names of the algorithms being evaluated will be the
% unique portions of the image filenames.
%
% ## Output
%
% ### Data file output
%
% #### Intermediate data and parameters
% A '.mat' file containing the following variables, as appropriate:
% - 'bands_estimated': A vector containing the wavelengths of the spectral
%   bands used in the estimated images.
% - 'bands_measured': A vector containing the wavelengths of the spectral
%   bands used by the reference reflectance data.
% - 'spectral_weights': A matrix for converting pixels in the spectral
%   space of the estimated hyperspectral images to the spectral space of
%   the reference reflectance data.
% - 'algorithms': A structure containing the names of the algorithms which were
%   evaluated, extracted from the names of the input estimated images.
%   'algorithms' has the following fields:
%   - 'spectral': A cell vector of character vectors containing the names of the
%     algorithms which produced the estimated spectral images.
%   - 'color': A cell vector of character vectors containing the names of the
%     algorithms which produced the estimated colour images.
% - 'algorithm_filenames': A structure of the same form as 'algorithms'
%   containing the names and paths of the input estimated images.
% - 'input_image_filename': A character vector containing the name and path of the
%   reference raw image of the ColorChecker.
% - 'input_label_filenames': A two-element cell vector of character vectors
%   containing the names and paths of the first and second label images for the
%   ColorChecker.
% 
% Additionally, the file contains the values of all parameters listed in
% `parameters_list`.
%
% The file is saved as 'EvaluateColorCheckerData.mat'.
%
% #### Evaluation results
%
% For each estimated ColorChecker patch, RGB error metrics and spectral error
% metrics, are output in the form of CSV files. Each CSV file contains results
% for all (RGB or spectral) estimated images. The estimated images will be
% listed in the files under the unique portions of their filenames. RGB error
% metrics are saved as '*_evaluateRGB.csv', whereas the spectral error metrics
% are saved as '*_evaluateSpectral.csv', where '*' is the number of the patch.
%
% Error metrics are also aggregated across images, and saved as
% 'EvaluateColorChecker_evaluateRGB.csv' and
% 'EvaluateColorChecker_evaluateSpectral.csv'.
%
% ## Detailed description of processing
%
% Ideal XYZ colours corresponding to the ColorChecker patches are created by
% integrating over the products of the measured spectral reflectances and the
% CIE tristimulus functions. These colours could be converted to RGB, if
% desired, using a whitepoint of [1, 1, 1] (for an equal energy radiator).
%
% A model of vignetting in the image of the ColorChecker is created from the
% variation in intensity of the ColorChecker's frame. The model is applied to
% the measured spectral reflectances, and ideal XYZ colours, to obtain "ground
% truth" images.  A 3 x 3 colour correction matrix is then fit between the XYZ
% colour image and the raw image of the ColorChecker, and output for later use
% in converting the estimated images to the sRGB colour space, for example.
%
% For RGB image evaluate a ground truth image is created from the raw image of
% the ColorChecker as follows:
% - The image is corrected for vignetting using the above model of vignetting.
% - For each ColorChecker patch, the averages of the Red, Green, and Blue pixels
%   in a square around the centroid of the patch are computed.
% - These average values are then replicated to every pixel in the patch, and
%   the vignetting model is applied to them to simulate an image without
%   spectral dispersion at the edges of the ColorChecker patches.
%
% For each ColorChecker patch, an independent evaluation is performed for all
% reconstructed images as follows:
% - As the illuminant is unknown, the average spectral signal in the
%   reconstructed image is calculated within a reference patch. The
%   reconstructed image is then multiplied by the ratio of the corresponding
%   average signal in the true image to this average signal.
% - The centroid of the patch is located, and a square is cropped from both the
%   true and reconstructed spectral images. 'evaluateAndSaveSpectral()' is
%   called on these patches to evaluate spectral reconstruction of the
%   homogenous area.
% - The border region of the patch is located, and its pixels are reshaped into
%   1D images from the true and reconstructed spectral images.
%   'evaluateAndSaveSpectral()' is called on these patches to evaluate spectral
%   reconstruction of the border region.
% - The above two steps are repeated for the RGB version of the reconstructed
%   image, using 'evaluateAndSaveRGB()'. In this case, the reference image is
%   created from the raw image of the ColorChecker as described above, and so
%   there is no need to rescale one image's channels to match the other's.
%
% The input raw image of the ColorChecker is also evaluated, by filling in `NaN`
% values for the colour channel values which were not sensed because of the
% colour filter array.
%
% The evaluation results for all patches are merged into summaries.
%
% Note that the extraction of patch centre and border regions uses the
% appropriate label image for the reconstructed image being evaluated.
%
% ## References
%
% The colour correction matrix is calculated based on the description in:
%
%   Karaimer, Hakki C., & Brown, Michael S. (2018). "Improving Color
%   Reproduction Accuracy on Cameras," In IEEE Conference on Computer Vision and
%   Pattern Recognition (CVPR) (pp. 6440–6449).

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created March 7, 2019

% List of parameters to save with results
parameters_list = {};

%% Input data and parameters


%% Save parameters and additional data to a file
save_variables_list = [ parameters_list, {...
} ];
save_data_filename = fullfile(output_directory, ['EvaluateColorCheckerData.mat']);
save(save_data_filename, save_variables_list{:});