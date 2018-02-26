%% Test script for the 'sonyQuantumEfficiency()' function
% 'sonyQuantumEfficiency()' should be capable of reproducing the plot
% '../FL3_GE_50S5C_quantumEfficiencyData.png'

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 26, 2018

lambda = linspace(200, 1200, 1000);
qe = sonyQuantumEfficiency(lambda);
qe_percent = qe * 100;

figure;
hold on
plot(lambda, qe_percent(:, 1), 'r');
plot(lambda, qe_percent(:, 2), 'g');
plot(lambda, qe_percent(:, 3), 'b');
hold off
xlabel('Wavelength [nm]')
ylabel('Quantum efficiency [%]')
title('Sony ICX655, 2/3" colour quantum efficiency')
legend('Red', 'Green', 'Blue')