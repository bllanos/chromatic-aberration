function [qe] = sonyQuantumEfficiency(lambda)
% SONYQUANTUMEFFICIENCY Query Sony ICX655, 2/3" colour quantum efficiency data
%
% ## Syntax
% qe = sonyQuantumEfficiency(lambda)
%
% ## Description
% qe = sonyQuantumEfficiency(lambda)
%   Returns Red, Green, and Blue quantum efficiencies for the given
%   wavelengths
%
% ## Input Arguments
%
% lambda -- Wavelengths
%   A vector containing wavelengths of light, measured in nanometres
%
% ## Output Arguments
%
% qe -- Sony ICX655, 2/3" colour quantum efficiencies
%   A 3-column matrix with the same number of rows as there are elements in
%   `lambda`. `qe(i, :)` contains the Red, Green, and Blue colour channels'
%   quantum efficiencies for the wavelength `lambda(i)`. The data is drawn
%   from the plot of quantum efficiency curves for the Sony ICX655, 2/3"
%   sensor saved in 'FL3_GE_50S5C_quantumEfficiencyData.png'.
%
% ## Data Source
%
% FLIR. (Jan. 27, 2017). FLIR FLEA3 GigE Vision Imaging Performance
%   Specification. version 1.1, [Online]. Available:
%   https://www.ptgrey.com/support/downloads/10109 (visited on 05/08/2017).
%
% ## Notes
% - Out of concern for possible copyright implications, the image
%   'FL3_GE_50S5C_quantumEfficiencyData.png' is not part of the repository.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 23, 2018

nargoutchk(1, 1);
narginchk(1, 1);

I = imread('FL3_GE_50S5C_quantumEfficiencyData.png');

% Pixel coordinates of the plot axes
x_axis_lim = [99, 931];
y_axis_lim = [429, 73];

% Measurement unit coordinates of the plot axes
x_axis_lim_lambda = [300, 1100];
y_axis_lim_qe = [0, 0.6];

% Start and end x-coordinates of segments along which the curves are
% occluded
r_segment_x = [
    186, 189;
    207, 210;
    370, 376;
    397, 399;
    542, 931
    ];
g_segment_x = [
    99, 185;
    306, 309;
    531, 931
    ];
b_segment_x = zeros(0, 2);
segments_x = {r_segment_x, g_segment_x, b_segment_x};
    
% Start and end y-coordinates of segments along which the curves are
% occluded
r_segment_y = [
    422, 412;
    395, 399;
    428, 425;
    277, 265;
    425, 429
    ];
g_segment_y = [
    423, 421;
    272, 261;
    426, 429
    ];
b_segment_y = zeros(0, 2);
segments_y = {r_segment_y, g_segment_y, b_segment_y};

threshold = 5;

x = (lambda - x_axis_lim_lambda(1)) / (x_axis_lim_lambda(2) - x_axis_lim_lambda(1));
x = round(x_axis_lim(1) + x * (x_axis_lim(2) - x_axis_lim(1)));

n_points = length(lambda);
n_channels = 3;
qe = zeros(n_points, n_channels);
for i = 1:n_points
    x_i = x(i);
    if x_i >= x_axis_lim(1) && x_i <= x_axis_lim(2)
        pixels = I(y_axis_lim(2):y_axis_lim(1), x_i, :);
        pixels = reshape(pixels, [], n_channels);
        for c = 1:n_channels
            segments_x_c = segments_x{c};
            done = false;
            for j = 1:size(segments_x_c, 1)
                if x_i >= segments_x_c(j, 1) && x_i <= segments_x_c(j, 2)
                    t_x = (x_i - segments_x_c(j, 1)) / (segments_x_c(j, 2) - segments_x_c(j, 1));
                    y_i = (1 - t_x) * segments_y{c}(j, 1) + t_x * segments_y{c}(j, 2);
                    done = true;
                    break;
                end
            end
            if ~done
                y_filter = [
                    pixels(:, c) - pixels(:, mod(c, n_channels) + 1),...
                    pixels(:, c) - pixels(:, mod(c + 1, n_channels) + 1) 
                    ];
                y_filter = all(y_filter >= threshold, 2);
                y_i = find(y_filter, 1, 'last');
                if isempty(y_i)
                    error('No suitable y-value found for \\lambda = %g', lambda(i));
                end
                y_i = y_axis_lim(2) + y_i - 1; 
            end
            t_y = (y_i - y_axis_lim(1)) / (y_axis_lim(2) - y_axis_lim(1));
            qe(i, c) = (1 - t_y) * y_axis_lim_qe(1) + t_y * y_axis_lim_qe(2);
        end
    end
end

end

