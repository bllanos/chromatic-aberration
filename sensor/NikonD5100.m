%% Nikon D5100 camera
% Save the Nikon D5100 spectral response curves in a format used by other
% scripts.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% The Nikon D5100 spectral sensitivity data were retrieved from
% https://spectralestimation.wordpress.com/data/, and correspond to the
% following publication:
%
%   Darrodi, M. M., Finlayson, G., Goodman, T., & Mackiewicz, M. (2015).
%     "Reference data set for camera spectral sensitivity estimation."
%     Journal of the Optical Society of America A, 32(3), 381-391.
%     doi:10.1364/JOSAA.32.000381
%
% ## Output
%
% ### Graphical output
% - A plot of the spectral response functions.
%
% ### Sensor quantum efficiency data
%
% A '.mat' file containing the following variables:
% - 'channel_mode': A Boolean value set to `false` to indicate that the
%   data in `sensor_map` represents spectral sensitivities, not colour
%   channel mappings
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th spectral response function (Red, Green, Blue) to light of
%   the j-th wavelength.
% - 'bands': A vector containing the wavelengths corresponding to the
%   columns of `sensor_map`.
%
% Additionally, the file contains the values of all parameters in the first
% section of the script below, for reference. (Specifically, those listed
% in `parameters_list`, which should be updated if the set of parameters is
% changed.)

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 3, 2018

% List of parameters to save with results
parameters_list = {
        'data_source',...
        'channel_mode'...
    };

%% Input data and parameters

data_source = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180802_referenceDataSetForCameraSpectralSensitivityEstimation/nikon.csv';
% Data represents spectral sensitivities, not colour channel mappings
channel_mode = false;

%% Load the data

qe_table = readtable(data_source);
bands = qe_table{:, 1};
sensor_map = qe_table{:, 2:end};

%% Visualization

figure;
hold on
plot(bands, sensor_map(:, 1), 'r');
plot(bands, sensor_map(:, 2), 'g');
plot(bands, sensor_map(:, 3), 'b');
hold off
xlabel('Wavelength [nm]')
ylabel('Response')
title('Nikon D5100 spectral sensitivity functions')
legend('Red', 'Green', 'Blue')

%% Save to a file
sensor_map = sensor_map.';
save_variables_list = [ parameters_list, {...
        'sensor_map',...
        'bands'...
    } ];
uisave(save_variables_list,'NikonD5100ColorMapData');