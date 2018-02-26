function [n] = sellmeierDispersion(lambda, constants)
% SELLMEIERDISPERSION Calculate indices of refraction from the Sellmeier dispersion formula
%
% ## Syntax
% n = sellmeierDispersion(lambda, constants)
%
% ## Description
% n = sellmeierDispersion(lambda, constants)
%   Returns indices of refraction for the given wavelengths
%
% ## Input Arguments
%
% lambda -- Wavelengths
%   A vector containing wavelengths of light, measured in nanometres
%
% constants -- Constants in the Sellmeier dispersion formula
%   A structure containing the values of the constants in the Sellmeier
%   dispersion formula, which are properties of the material. `constants`
%   has the following fields, defined at
%   https://en.wikipedia.org/wiki/Sellmeier_equation
%   - B_1
%   - B_2
%   - B_3
%   - C_1
%   - C_2
%   - C_3
%
%   The units of the 'C' constants are square micrometres.
%
% ## Output Arguments
%
% n -- Refractive indices
%   Indices of refraction calculated from the Sellmeier dispersion formula
%   given the material constants in `constants`. `n(i)` is the index of
%   refraction for light of wavelength `lambda(i)`.
%
% ## References
% - MATLAB code for the Sellmeier dispersion formula at
%   RefractiveIndex.INFO:
%   https://refractiveindex.info/?shelf=glass&book=BK7&page=SCHOTT

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 26, 2018

nargoutchk(1, 1);
narginchk(2, 2);

lambda_sq = (lambda / 1000) .^ 2; % Convert to micrometres
n = sqrt(...
    1 +...
    constants.B_1 ./ (1 - (constants.C_1 ./ lambda_sq)) +...
    constants.B_2 ./ (1 - (constants.C_2 ./ lambda_sq)) +...
    constants.B_3 ./ (1 - (constants.C_3 ./ lambda_sq))...
    );
end