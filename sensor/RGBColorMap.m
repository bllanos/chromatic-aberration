%% Identity RGB colour conversion
% Create a mapping from the RGB colour space to itself, to be used when
% estimating latent RGB images using 'CorrectByHyperspectralADMM.m', for
% example.
%
% ## Usage
% Run the script.
%
% ## Input
%
% None
%
% ## Output
%
% ### Colour conversion data
%
% A '.mat' file containing the following variables:
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th colour channel (Red, Green, or Blue) to light of the j-th
%   colour channel (Red, Green, or Blue).
% - 'bands': A vector containing the colour channel indices corresponding
%   to the columns of `sensor_map`.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 30, 2018

bands = (1:3).';
sensor_map = eye(length(bands));

save_variables_list = {...
        'bands',...
        'sensor_map'...
    };
uisave(save_variables_list,'RGBColorMapData');