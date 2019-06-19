%% Test script for the 'scaleSignals()' function

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 26, 2018

%% Run qualitative test

bands = (200:1200).';
signals = sonyQuantumEfficiency(bands).';

[scale, signals_scaled] = scaleSignals(signals, 0.05, [1, 2, 1]);

%% Visualization

figure;
hold on
plot(bands, signals(1, :), 'r');
plot(bands, signals_scaled(1, :), 'Color', 'r', 'LineStyle', ':');
plot(bands, signals(2, :), 'g');
plot(bands, signals_scaled(2, :), 'Color', 'g', 'LineStyle', ':');
plot(bands, signals(3, :), 'b');
plot(bands, signals_scaled(3, :), 'Color', 'b', 'LineStyle', ':');
hold off
xlabel('Wavelength [nm]')
ylabel('Quantum efficiency [%]')
title('Sony ICX655, 2/3" colour quantum efficiency')
legend('Red', 'Red - Scaled', 'Green', 'Green - Scaled', 'Blue', 'Blue - Scaled')