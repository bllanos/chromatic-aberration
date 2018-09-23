%% CIE 1931 camera
% Save the CIE 1931 colour matching functions in a camera response function
% format, to serve as a reference camera.
%
% ## Usage
% Modify the parameters, the first code section below, then run.
%
% ## Input
%
% The CIE 1931 Standard (2°) Observer spectral tristimulus functions were
% retrieved from Table 1 of the ASTM E308 standard and saved as
% '/home/llanos/GoogleDrive/ThesisResearch/Data and Results/20180614_ASTM_E308/Table1_CIE1931_2DegStandardObserver.csv'.
%
% ## Output
%
% ### Graphical output
% - A plot of the tristimulus functions.
%
% ### Sensor quantum efficiency data
%
% A '.mat' file containing the following variables:
% - 'channel_mode': A Boolean value set to `false` to indicate that the
%   data in `sensor_map` represents spectral sensitivities, not colour
%   channel mappings
% - 'sensor_map': A 2D array, where `sensor_map(i, j)` is the sensitivity
%   of the i-th tristimulus function (x-bar, y-bar or z-bar) to light of
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
% File created June 15, 2018

% List of parameters to save with results
parameters_list = {
        'data_source',...
        'channel_mode'...
    };

%% Input data and parameters

data_source = '/home/llanos/GoogleDrive/ThesisResearch/Data/20180614_ASTM_E308/Table1_CIE1931_2DegStandardObserver.csv';
% Data represents spectral sensitivities, not colour channel mappings
channel_mode = false;

%% Load the data

xyzbar_table = readtable(data_source);
bands = xyzbar_table{:, 1};
sensor_map = xyzbar_table{:, 2:end};

%% Visualization

figure;
hold on
plot(bands, sensor_map(:, 1), 'r');
plot(bands, sensor_map(:, 2), 'g');
plot(bands, sensor_map(:, 3), 'b');
hold off
xlabel('Wavelength [nm]')
ylabel('Response')
title('CIE 1931 colour matching functions')
legend('x-bar', 'y-bar', 'z-bar')

%% Save to a file
sensor_map = sensor_map.';
save_variables_list = [ parameters_list, {...
        'sensor_map',...
        'bands'...
    } ];
uisave(save_variables_list,'CIE1931ColorMapData');