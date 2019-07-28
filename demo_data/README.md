# Demo/Walkthrough Instructions
Bernard Llanos

Supervised by Dr. Y.H. Yang

University of Alberta, Department of Computing Science

File created July 17, 2019

## Notes
- All MATLAB scripts are expected to be executed while MATLAB's current
  folder is the root folder of the codebase (the parent folder of the folder
  containing this file).

## High-dynamic range image synthesis
Create high-dynamic range images from raw images captured using different
exposure times.

1. Create a sub-directory called 'averaged_images' in this directory.
2. Run [../raw_image_preprocessing/PreprocessRAWImages.m](../raw_image_preprocessing/PreprocessRAWImages.m)
   - A series of figures will open to show correlations between pixel intensities
     recorded under different exposures. The figures will allow you to assess
     whether the image sensor has a linear response to exposure time.
   - Images output to './averaged_images/' are dark frame-subtracted versions of
     the input images from [./captured_images/images_filtered_light/disks/](./captured_images/images_filtered_light/disks/).
   - Exposure-blended versions of the images from './averaged_images/' are saved
     in [./hdr_averaged_images/](./hdr_averaged_images/). The '.mat' files contain the high-quality,
     versions of the images, whereas the '.tif' files are the human-viewable
     versions of the images. '_raw01.tif' images are scaled versions of the
     exposure-blended images, whereas '_rgb01.tif' images have additionally
     been demosaiced. The colours of the demosaiced images are uncorrected,
     and are not necessarily meaningful, because the images were captured under
     narrowband-filtered light.

## Multispectral image synthesis
Create multispectral images by stacking images captured under light filtered
by different optical bandpass filters. Simultaneously calibrate the relative
spectral sensitivity of the RGB camera's colour channels at the wavelengths
of light corresponding to each filter.

1. Run [../raw_image_preprocessing/RAWBandImagesToDataset.m](../raw_image_preprocessing/RAWBandImagesToDataset.m)
   - Figures will open to show correlations between colour channels under each
     of the different filtered illuminations. An ideal image sensor would have
     linear relationships between intensities in different colour channels
     under narrowband-filtered illumination.
   - The last figures to open show the calibrated relative spectral sensitivity
     functions of the camera.
   - Images are output to [./multispectral_images/](./multispectral_images/):
     - '.mat' files are the high-quality versions of the images
     - Images with filenames ending in '3_01.tif' are versions of synthetic
       colour images for human viewing.
     - Images with filenames ending in 'nm.tif' are versions of the individual
       bands in the multispectral images for human viewing.
     - The types of images are described in the documentation comments of
       [../raw_image_preprocessing/RAWBandImagesToDataset.m](../raw_image_preprocessing/RAWBandImagesToDataset.m)

## Lateral chromatic aberration calibration
Calibrate models of lateral chromatic aberration to use for chromatic
aberration correction or simulation.

### Modelling dispersion from disk keypoints
1. It should not be necessary for the demo images provided, but to assist
   with disk localization, the image can be corrected for intensity non-uniformity,
   and a mask can be provided to eliminate regions without disk keypoints.
   Sample data has been provided for both tasks:
   - [./hdr_averaged_images/d1_disks32cmV2_maskVignetting.png](./hdr_averaged_images/d1_disks32cmV2_maskVignetting.png) is a mask used
     for intensity non-uniformity (vignetting) correction. White marks pixels
     which should all have the same intensity. Note that the black regions extend
     outside the borders of the disks, to prevent blurring of the disks from
     affecting vignetting calibration.
   - [./hdr_averaged_images/d1_disks32cmV2_maskDisks.png](./hdr_averaged_images/d1_disks32cmV2_maskDisks.png) is a mask used for
     disk localization. White marks pixels in regions which contain disks.
     It is normally not necessary to make such an accurate mask. Note that the
     white regions extend outside the borders of the disks.
   - Both of the masks are single-channel images.
2. To create models of spectral dispersion, run [../disk_fitting/RAWDiskDispersion.m](../disk_fitting/RAWDiskDispersion.m)
   - There will be warnings about rank deficient matrices, because the
     cross-validation procedure is currently configured to test polynomial
     degrees (in [../dispersion_model/xylambdaPolyfit.m](../dispersion_model/xylambdaPolyfit.m)) that are much higher
     than necessary for the sample data.
   - Figures will be generated showing the models of dispersion produced
     for image warping to either apply, or correct aberration. (Therefore,
     it will seem like two identical sets of figures have been produced.)
3. In preparation for creating models of warping between colour channels, edit
   [../disk_fitting/RAWDiskDispersion.m](../disk_fitting/RAWDiskDispersion.m):
   - Replace `rgb_mode = false;` with `rgb_mode = true;`.
   - Replace `input_images_wildcard = './demo_data/hdr_averaged_images/*disks*nm.mat';`
     with `input_images_wildcard = './demo_data/multispectral_images/d1_disks32cmV2_raw.mat';`.
4. Run [../disk_fitting/RAWDiskDispersion.m](../disk_fitting/RAWDiskDispersion.m) to create models of warping
   between colour channels. The output will be similar to the output in Step 2.
   Note that the vignetting and disk mask images will not be used, because they
   do not have filepaths that will be found based on the input image. Model
   generation should succeed regardless.

#### Notes
- To see the results of disk detection, set more of the verbosity options in
  [../disk_fitting/RAWDiskDispersion.m](../disk_fitting/RAWDiskDispersion.m) to `true`.
- If disk detection is extremely difficult, one can manually isolate disks using
  an accurate segmentation mask, and set `findAndFitDisks_options.mask_as_threshold`
  in [../disk_fitting/RAWDiskDispersion.m](../disk_fitting/RAWDiskDispersion.m) to `true`.
  [../disk_fitting/RAWDiskDispersion.m](../disk_fitting/RAWDiskDispersion.m) can then use the mask as the initial
  segmentation of disks from the background.
- Output files are saved in [./dispersion_models/disk_fitting/](./dispersion_models/disk_fitting/).

### Modelling dispersion from image registration
1. To create models of spectral dispersion, run [../dispersion_model/RegistrationDispersion.m](../dispersion_model/RegistrationDispersion.m)
   - There will be warnings about the image registration failing to converge.
     These warnings are to be expected for the input image, as it contains large
     textureless regions.
   - There will be warnings about rank deficient matrices, because the
     cross-validation procedure is currently configured to test polynomial
     degrees (in [../dispersion_model/xylambdaPolyfit.m](../dispersion_model/xylambdaPolyfit.m)) that are much higher
     than necessary for the sample data.
   - Figures will be generated showing the models of dispersion produced
     for image warping to either apply, or correct aberration. (Therefore,
     it will seem like two identical sets of figures have been produced.)
2. In preparation for creating models of warping between colour channels, edit
   [../dispersion_model/RegistrationDispersion.m](../dispersion_model/RegistrationDispersion.m):
   - Replace `rgb_mode = false;` with `rgb_mode = true;`.
   - Replace `input_images_wildcard = './demo_data/multispectral_images/*colorChecker30cm*_dHyper.mat';`
     with `input_images_wildcard = './demo_data/multispectral_images/*colorChecker30cm*_d3.mat';`.
   - Replace `input_images_variable_name = 'I_hyper';`
     with `input_images_variable_name = 'I_3';`.
3. Run [../dispersion_model/RegistrationDispersion.m](../dispersion_model/RegistrationDispersion.m) to create models of warping
   between colour channels. The output will be similar to the output in Step 1.

#### Notes
- The dispersion models generated by this demo will be different from the ones
  generated during the previous demo. The previous demo operated on a
  cropped image, containing only a portion of the field of view. The current demo
  is run on a downscaled image that captures the entire field of view.
- Output files are saved in [./dispersion_models/registration/](./dispersion_models/registration/).

## Colour correction

1. [../data_analysis/CalibrateColorCorrection.m](../data_analysis/CalibrateColorCorrection.m) can produce colour correction
   models, but presently requires the CIE tristimulus functions (`xyzbar_filename`),
   which are not included in this repository. The output of the script has been
   saved to [./color_correction/CalibrateColorCorrectionData.mat]().
2. Run [../data_analysis/ColorCorrection.m](../data_analysis/ColorCorrection.m) to generate white-balanced and
   gamma-corrected images in [./color_correction/](./color_correction/) from the images in
   [./multispectral_images/](./multispectral_images/). The output images have names ending in '_wb.tif'.
3. If the third-party colour correction code mentioned in the top-level README
   file is available, replace `use_chromadapt = true;`
   with `use_chromadapt = false;` in [../data_analysis/ColorCorrection.m](../data_analysis/ColorCorrection.m),
   and then re-run the script. The colours in the output images (having names
   ending in '_M_homog.tif') should be more saturated than those generated by
   white-balancing.

### Notes
- The corrected versions of the multispectral images in this demo are identical
  to the corrected versions of the RGB images, because the RGB images were
  simulated from the multispectral images, and because the colour correction
  is performed in RGB space.

## Image reconstruction

### Chromatic aberration correction by image warping
Approximately correct lateral chromatic aberration by warping colour channels.

1. Run [../aberration_correction/CorrectByWarping.m](../aberration_correction/CorrectByWarping.m).
2. The output colour image (in [./aberration_correction/](./aberration_correction/)) is in the raw
   colour space of the camera. In order to correct its colours, modify
   [../data_analysis/ColorCorrection.m](../data_analysis/ColorCorrection.m) as follows:
   - Replace `spectral_wildcard = './demo_data/multispectral_images/*dHyper.mat';`
     with `spectral_wildcard = [];`.
   - Replace `color_wildcard = './demo_data/multispectral_images/*d3.mat';`
     with `color_wildcard = './demo_data/aberration_correction/d1_disks32cmV2_d3_unwarped.mat';`.
   - Replace `color_variable_name = 'I_3';` with
     `color_variable_name = 'I_unwarped';`.
3. Run [../data_analysis/ColorCorrection.m](../data_analysis/ColorCorrection.m). Compare the output file,
   [./color_correction/d1_disks32cmV2_d3_unwarped_wb.tif](./color_correction/d1_disks32cmV2_d3_unwarped_wb.tif)
   with [./color_correction/d2_colorChecker30cm_d3_wb.tif](./color_correction/d2_colorChecker30cm_d3_wb.tif).
   The first image should have reduced colour fringes around the disks
   relative to the second.

### Chromatic aberration correction by spectral reconstruction
Correct lateral chromatic aberration in the spectral domain by simultaneously
recovering a multispectral image corresponding to the input image.

1. Run [../aberration_correction/CorrectByHyperspectralADMM.m](../aberration_correction/CorrectByHyperspectralADMM.m). Note that the
   script will take time to run (about 5-10 minutes).
2. As done above, prepare to produce a colour-corrected result. Modify
   [../data_analysis/ColorCorrection.m](../data_analysis/ColorCorrection.m) as follows:
   - Set `spectral_wildcard` to `'./demo_data/aberration_correction/*latent.mat'`.
   - Replace `spectral_variable_name = 'I_hyper';` with
     `spectral_variable_name = 'I_latent';`.
   - Replace `color_wildcard = './demo_data/multispectral_images/*d3.mat';`
     with `color_wildcard = [];`.
3. Run [../data_analysis/ColorCorrection.m](../data_analysis/ColorCorrection.m). Compare the output file,
   [./color_correction/d1_disks32cmV2_raw_patch64x64_pad8_weightsTarget99And99_latent_wb.tif](./color_correction/d1_disks32cmV2_raw_patch64x64_pad8_weightsTarget99And99_latent_wb.tif)
   with [./color_correction/d1_disks32cmV2_d3_unwarped_wb.tif](./color_correction/d1_disks32cmV2_d3_unwarped_wb.tif).
   The two images will both have reduced colour fringes around the disks.
   The first image will have some false colours from regularization.
   The severity of the false colours depends on the camera spectral sensitivity,
   and on the spectral characteristics of the scene, as described in
   the thesis associated with this codebase. (See Figure 5.17, Section 5.5.2,
   and Section 6.3 in the thesis.) In natural images, as opposed to images
   synthesized from images taken under filtered light, results may improve.