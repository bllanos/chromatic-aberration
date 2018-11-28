function [lambda_Rad, Rad, varargout] = reflectanceToRadiance(lambda_L, L, lambda_Ref, Ref, varargin)
% REFLECTANCETORADIANCE Convert spectral reflectances to spectral radiances
%
% ## Syntax
% [lambda_Rad, Rad] = reflectanceToRadiance(...
%     lambda_L, L, lambda_Ref, Ref...
% )
% [lambda_Rad, Rad, Rad_normalized] = reflectanceToRadiance(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, C_ind, int_method]...
% )
% [...
%   lambda_Rad, Rad, Rad_normalized, lambda_C_resampled, C_resampled...
% ] = reflectanceToRadiance(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, C_ind, int_method]...
% )
%
% ## Description
% [lambda_Rad, Rad] = reflectanceToRadiance(...
%     lambda_L, L, lambda_Ref, Ref...
% )
%   Returns radiances for the given reflectances seen under the given
%   illuminant
%
% [lambda_Rad, Rad, Rad_normalized] = reflectanceToRadiance(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, C_ind, int_method]...
% )
%   Additionally returns radiances normalized by the sensor's response to
%   the illuminant.
%
% [...
%   lambda_Rad, Rad, Rad_normalized, lambda_C_resampled, C_resampled...
% ] = reflectanceToRadiance(...
%     lambda_L, L, lambda_Ref, Ref, lambda_C, C [, C_ind, int_method]...
% )
%   Additionally returns the version of the sensor's response functions
%   which can be used for colour calculation, along with the wavelengths at
%   which it was resampled. Five output arguments must be requested.
%
% ## Input Arguments
%
% lambda_L -- Illuminant wavelengths
%   A vector containing wavelengths of light at which the spectral power
%   distribution of the illuminant has been sampled. Units are nanometres.
%
% L -- Illuminant
%   A vector the same length as `lambda_L` containing the spectral power
%   distribution of the illuminant. `L(i)` is the spectral power at the
%   wavelength `lambda_L(i)`.
%
% lambda_Ref -- Reflectance wavelengths
%   A vector containing wavelengths of light at which the spectral
%   reflectances have been sampled. Units are nanometres.
%
% Ref -- Spectral reflectances
%   A matrix of size length(lambda_Ref) x n, where 'n' is the number of
%   different spectral reflectances. The rows of `Ref` correspond to the
%   wavelengths in `lambda_Ref`. `Ref(i, j)` is the reflectance of the j-th
%   sample at the i-th wavelength. Reflectances should be fractions in the
%   range [0, 1].
%
% lambda_C -- Reference wavelength values
%   A vector of wavelengths at which the 'C(lambda)' distributions were
%   sampled. If `C` are the CIE 1931 colour matching functions, an
%   even sampling of 5 nm is required, spanning the range from 360 nm
%   (inclusive) to 780 nm (inclusive), in order for this function to
%   compute normalized radiances according to the ASTM E308 standard.
%
% C -- Sensor response functions
%   A matrix, where `C(i, j)` contains the value of the j-th sensor
%   response function evaluated at the wavelength `lambda_C(i)`. `C(:,
%   C_ind)` is used to normalize radiances by the sensor response to the
%   illuminant, such as when computing colours according to the ASTM E308
%   standard. The sensor response functions are assumed to be zero outside
%   of the range of `lambda_C`.
%
% C_ind -- Sensor response function index
%   The column index of the sensor response function in `C` to use for
%   normalizing radiances. For instance, if `C` are the CIE 1931 colour
%   matching functions, `C_ind` should be 2.
%
%   Defaults to 1 if not passed.
%
% int_method -- Numerical integration method
%   The numerical integration method to use when computing the radiance
%   normalization constant. `int_method` is passed to
%   `integrationWeights()` as its `method` input argument.
%
%   Defaults to 'rect' if not passed.
%
% ## Output Arguments
%
% lambda_Rad -- Radiance wavelength values
%   A vector with a length equal to the size of the first dimension of
%   `Rad`, containing the wavelengths at which the radiances in `Rad` and
%   `Rad_normalized` have been sampled. The input reflectances and
%   illuminant spectral power distribution are resampled according to this
%   new wavelength array, so that they can be multiplied together to
%   produce the output radiances.
%
% Rad -- Radiances
%   An n-column matrix, where the columns correspond to the columns of
%   `Ref`. The i-th column is the radiance of the i-th reflectance,
%   computed by taking the product of the reflectance with the illuminant
%   spectral power distribution.
%
% Rad_normalized -- Normalized radiances
%   A version of `Rad` which has been normalized by the response of the
%   `C` function to the illuminant spectral power distribution.
%   Normalized radiances are required when computing the corresponding
%   colours.
%
% lambda_C_resampled -- Resampled reference wavelength values
%   A vector containing the wavelength values in `lambda_Rad` that are
%   inside the interval represented by `lambda_C`.
%
% C_resampled -- Resampled sensor response function
%   A version of `C` which has been resampled according to the
%   `lambda_C_resampled` series of wavelengths. The first and last values
%   of `C_resampled` have been adjusted to compensate for any truncation of
%   the domain of wavelengths `lambda_C` represented by `lambda_Rad`.
%   `C_resampled` should be used instead of `C` when calculating CIE
%   tristimulus responses to the radiances in `Rad_normalized`, according
%   to the ASTM E308 standard.
%
% ## References
% - ASTM E308-17 Standard Practice for Computing the Colors of Objects by
%   Using the CIE System, ASTM International, West Conshohocken, PA, 2017,
%   https://doi.org/10.1520/E0308-17
%
% See also ciedIlluminant, reflectanceToColor, cieSpectralToColor,
% resampleArrays, integrationWeights

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 15, 2018

    % Add the weights for wavelengths missing from the full range used in
    % the ASTM E308 standard. See Section 7.3.2.2 of the standard.
    function C_resampled = sumEnds(lambda, lambda_resampled, C, C_resampled)
        spacing = mean(diff(lambda_resampled));
        spacing_reference = 5;
        if any(diff(lambda) ~= spacing_reference) || lambda(1) ~= 360 || lambda(end) ~= 780
            warning('`lambda_C` does not have the form required for this function to compute normalized radiances according to the ASTM E308 standard.')
        end
        ind = find(lambda < lambda_resampled(1));
        if ~isempty(ind)
            if length(ind) > 1
                lower_missing_weights = sum(spacing_reference .* C(ind(1:(end - 1)), :), 1);
            else
                lower_missing_weights = zeros(1, size(C, 2));
            end
            lower_missing_weights = lower_missing_weights +...
                mod(lambda_resampled(1), spacing_reference) * C(ind(end), :);
            C_resampled(1, :) = C_resampled(1, :) + lower_missing_weights ./ spacing;
        end
        ind = find(lambda > lambda_resampled(end));
        if ~isempty(ind)
            if length(ind) > 1
                upper_missing_weights = sum(spacing_reference .* C(ind(2:end), :), 1);
            else
                upper_missing_weights = zeros(1, size(C, 2));
            end
            upper_missing_weights = upper_missing_weights +...
                (spacing_reference - mod(lambda_resampled(end), spacing_reference)) * C(ind(1), :);
            C_resampled(end, :) = C_resampled(end, :) + upper_missing_weights ./ spacing;
        end
    end

if isempty(varargin)
    nargoutchk(2, 2);
    narginchk(4, 4);
    normalize = false;
elseif length(varargin) == 2 || length(varargin) == 3 || length(varargin) == 4
    nargoutchk(3, 5);
    if nargout == 4
        error('Either three or five output arguments must be requested.');
    end
    narginchk(6, 8);
    normalize = true;
    lambda_C = varargin{1};
    C = varargin{2};
    int_method = 'rect';
    if length(varargin) > 2
        C_ind = varargin{3};
        
        if length(varargin) > 3
            int_method = varargin{4};
        end
    else
        C_ind = 1;
    end
else
    error('Either zero, two, three, or four optional input arguments must be provided, not %d.', length(varargin));
end

[...
    L_resampled, Ref_resampled,...
    lambda_Rad...
] = resampleArrays(...
    lambda_L, L, lambda_Ref, Ref,...
    'spline'...
    );
Rad = Ref_resampled .* repmat(L_resampled, 1, size(Ref_resampled, 2));

% Normalize radiances
if normalize
    ybar = C(:, C_ind);
    [...
        L_resampled2, ybar_resampled,...
        lambda_resampled2...
    ] = resampleArrays(...
        lambda_L, L, lambda_C, ybar,...
        'spline'...
        );
    weights = integrationWeights(lambda_resampled2, int_method);
    ybar_resampled = sumEnds(lambda_C, lambda_resampled2, ybar, ybar_resampled);

    % Compute the normalization constant
    N = dot(L_resampled2 .* ybar_resampled, weights);

    Rad_normalized = Rad ./ N;
    
    varargout{1} = Rad_normalized;
    
    if nargout > 3
        [C_resampled, lambda_C_resampled] = resampleArrays(...
            lambda_C, C, lambda_Rad,...
            'spline'...
        );
        varargout{2} = lambda_C_resampled;
        C_resampled = sumEnds(lambda_C, lambda_Rad, C, C_resampled);
        varargout{3} = C_resampled;
    end
end

end

