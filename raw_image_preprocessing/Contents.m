% RAW_IMAGE_PREPROCESSING
% Version 2.0.0 18-Jul-2019
%
% Dataset generation from raw colour-filter array images, and utility functions
% for colour-filter array image processing.
%
% Image dataset generation
%   blendExposures         - Combine images taken under different exposure settings
%   darkSubtract           - Read raw Bayer pattern images and perform averaging and dark frame subtraction
%   findAndGroupImages     - List images and optionally group them by spectral band
%   PreprocessRAWImages    - Real image dataset preparation
%   RAWBandImagesToDataset - Assemble images for spectral bands into spectral and colour images
%   relativeSensitivity    - Fit scaling factors between colour channels
%
% Colour-filter array images
%   bayerDownsample        - Create subresolution image from colour-filter array data
%   bayerMask              - Create a logical array to index colour channels in a raw image
%   bilinearDemosaic       - Demosaic an image by bilinear interpolation
%   imreadRAW              - Read a raw Bayer pattern image
%   mosaic                 - Mosaic colour images to produce colour-filter array images
%   offsetBayerPattern     - Determine the Bayer pattern for an offset image
%
% Non-high-dynamic range image generation (unused)
%   AverageRAWImages       - Average RAW images
%   dirreadRAW             - Read and write raw Bayer pattern images in a directory