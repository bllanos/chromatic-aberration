function I = clipAndRemap(I, class_name, method, limits)
% CLIPANDREMAP  Clip values and convert to an integer type
%
% ## Syntax
% I = clipAndRemap(I, class_name, method, limits)
%
% ## Description
%
% I = clipAndRemap(I, class_name, method, limits)
%   Returns the array as an array of the specified type, with values clipped and
%   rescaled according to the given method.
%
% ## Input Arguments
%
% I -- Data
%   An array containing data to convert from one type to another
%
% class_name -- New type
%   A character vector containing the name of the new class for the data. The
%   class is expected to be an integer class.
%
% method -- Clipping method
%   A character vector, which must be one of the following options:
%   - 'quantiles': `limits` is interpreted as a pair of quantiles. The data
%     will be clipped to range between the lower and upper quantiles given
%     in `limits` before being rescaled.
%   - 'values': `limits` is interpreted as a pair of values. The data will be
%     clipped to the range between the lower and upper values provided in
%     `limits` before being rescaled.
%
% limits -- Clipping arguments
%   A two element vector giving lower and upper limits, respectively, on the
%   data to be enforced before rescaling. `limits` is interpreted as described
%   in the documentation of `method` above.
%
% ## Output Arguments
%
% I -- Output data
%   A version of the input argument `I` cast to class `class_name`. Before
%   casting, `I` is first clipped to the limits given in `limits`, then rescaled
%   so that the lower and upper limits are mapped to zero, and to the upper
%   limit of the type `class_name`, respectively.
%
% See also cast

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 16, 2019

narginchk(4, 4);
nargoutchk(1, 1);

if strcmp(method, 'quantiles')
    limit_values = quantile(I, limits, 'all');
elseif strcmp(method, 'values')
    limit_values = limits;
else
    error('Unrecognized value of `method`, "%s".', method);
end
I(I < limit_values(1)) = limit_values(1);
I(I > limit_values(2)) = limit_values(2);

min_value = 0;
max_value = intmax(class_name);
I = ((I - limit_values(1)) .* ((max_value - min_value) ./ (limit_values(2) - limit_values(1)))) + min_value;
I = cast(I, class_name);

end