% SAMPLING
% Version 2.0.0 18-Jul-2019
%
% Sampling theory utility functions and conversion between sampling and colour
% spaces.
%
% Sampling rate calculation
%   bandlimit           - Find the approximate bandlimits of discrete signals
%   findSampling        - Find an optimal sampled representation for spectral data
%   TestSamplingWeights - Test script for the 'findSampling()' and 'bandlimit()'
%
% Resampling
%   resampleArrays      - Interpolate array data to match array dimensions
%   resamplingWeights   - Create a resampling matrix operator
%
% Colour space conversion
%   channelConversion   - Convert an image to a different colour/spectral space
%   colorWeights        - Create a colour conversion matrix operator
%
% Filter kernels for resampling signals
%   delta               - Delta distribution interpolation kernel (i.e. no interpolation)
%   gaussian            - Gaussian interpolation kernel
%   triangle            - Triangle interpolation kernel (linear interpolation)
%
% Other
%   integrationWeights  - Calculate weights for numerical integration
%   scaleSignals        - Find scaling factors to make a set of signals more uniform
%   TestScaleSignals    - Test script for the 'scaleSignals()' function
