%% Test script for the 'sellmeierDispersion()' function
% 'sellmeierDispersion()' should be capable of reproducing the plot
% at https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 26, 2018

lambda = linspace(200, 1200, 1000);
% Constants for SCHOTT N-BK7 glass retrieved from the SCHOTT glass
% datasheet provided at
% https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT
constants.B_1 = 1.03961212;
constants.B_2 = 0.231792344;
constants.B_3 = 1.01046945;
constants.C_1 = 0.00600069867;
constants.C_2 = 0.0200179144;
constants.C_3 = 103.560653;

n = sellmeierDispersion(lambda, constants);

figure;
plot(lambda, n, 'k');
xlabel('Wavelength [nm]')
ylabel('Refractive index')
title('SCHOTT N-BK7 refractive index')