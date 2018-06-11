function [Y1_resampled, varargout] = resampleArrays(x1, Y1, x2, varargin)
% RESAMPLEARRAYS Interpolate array data to match array dimensions
%
% ## Syntax
% [Y_new, x_new] = resampleArrays(x1, Y1, x2 [, method, extrapolation])
% [Y1, Y2, x_new] = resampleArrays(x1, Y1, x2, Y2 [, method, extrapolation])
%
% ## Description
% [Y_new, x_new] = resampleArrays(x1, Y1, x2 [, method, extrapolation])
%   Converts the `Y1` array from the sampling of `x1` to that of `x_new`.
%   One or two output arguments can be requested.
%
% [Y1, Y2, x_new] = resampleArrays(x1, Y1, x2, Y2 [, method, extrapolation])
%   Resamples either `Y1` or `Y2` such that both arrays have the finest
%   sampling. Two or three output arguments can be requested.
%
% ## Input Arguments
%
% x1 -- Sampling points for the first array
%   A vector containing values representing the sample locations at which
%   slices of `Y1` were obtained. `x1(i)` corresponds to the data in
%   `Y1(:,...,:,i,:,...,:)`. The first dimension in `Y1` with a size
%   matching the length of `x1` is assumed to correspond to the values in
%   `x1`.
%
% Y1 -- First array
%   An N-d array corresponding to the sample values in `x1`.
%
% x2 -- Second set of sampling points
%   In the syntax where `Y2` is passed as a non-empty array, `x2` is to
%   `Y2` as `x1` is to `Y1`. Otherwise, `x2` is the new set of sampling
%   locations at which the data in `Y1` must be resampled.
%
% Y2 -- Second array
%   An N-d array corresponding to the sample values in `x2`.
%
% method -- Interpolation method
%   The `method` input argument for the MATLAB 'interp1()' function, to use
%   when resampling array data.
%
% extrapolation -- Extrapolation settings
%   The `extrapolation` input argument for the MATLAB 'interp1()' function,
%   to use when resampling array data.
%
%   If `extrapolation` is not passed, or is 'none', no extrapolation will
%   be performed, and resampling will be bounded to the sample points in
%   the interval `[max(min(x1), min(x2)), min(max(x1), max(x2))]`.
%
% ## Output Arguments
%
% Y_new -- Resampled version of Y1
%   If extrapolation is active (see above), `Y_new` is an N-d array similar
%   to the input argument `Y1`, but resampled so that the array has a size
%   of `length(x2)` along the dimension corresponding to `x1`. `Y_new` is
%   produced by resampling each line of data in `Y1` using the MATLAB
%   'interp1()' function. Lines are taken along the dimension corresponding
%   to `x1`. Each line is resampled as follows:
%     ```
%     Y_new(j_1, j_2, ..., j_(d-1), :, j_(d + 1), ..., j_n) = interp1(...
%                  x1,...
%                  Y1(j_1, j_2, ..., j_(d-1), :, j_(d + 1), ..., j_n)...
%                  x2,...
%                  method, extrapolation...
%                );
%     ```
%   where the set of `j` values index a possible line of data in `Y1`, and
%   `d` is the dimension corresponding to `x1`.
%
%   If extrapolation is disabled, the locations used for resampling are the
%   values in `x_new`.
%
% x_new -- Resampling locations
%   If extrapolation is enabled, `x_new` is equal to `x1` or `x2`,
%   whichever is longer. Otherwise, `x_new` consists of the values from
%   `x1` or in `x2` which are in the interval `[max(min(x1), min(x2)),
%   min(max(x1), max(x2))]`, whichever set of values (the set from `x1` or
%   the set from `x2`) is more numerous.
%
% Y1 -- Resampled version of Y1
%   In the syntax where two data arrays are provided, and extrapolation is
%   enabled, `Y1` is equal to the input argument `Y1` if `length(x1)` is
%   greater than or equal to `length(x2)`. Otherwise, again with
%   extrapolation enabled, `Y1` is resampled so that the dimension
%   corresponding to `x1` now has a size equal to the length of `x2`. The
%   resampling is performed as described in the documentation of `Y_new`
%   above.
%
%   If extrapolation is disabled, `Y1` is resampled according to the new
%   sampling locations `x_new`. 
%
% Y2 -- Resampled version of Y2
%   A version of the input argument `Y2`, produced in the same way as
%   described in the documentation of the output argument `Y1` above.
%
% ## Caution
% Ensure that the dimension of `Y1` corresponding to `x1` is the first
% dimension of `Y1` with a size equal to the length of `x1` (and likewise
% for `Y2` and `x2`), otherwise resampling will occur along the wrong
% dimension.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 10, 2018

% Parse input arguments
narginchk(3, 6);
n_varargs = length(varargin);
Y2 = [];
method = 'spline';
method_ind = 0;
extrapolation = 'none';
if n_varargs > 0
    if isStringScalar(varargin{1}) || ischar(varargin{1})
        method_ind = 1;
        method = varargin{1};
        nargoutchk(1, 2);
    else
        Y2 = varargin{1};
        nargoutchk(2, 3);
    end
    
    if n_varargs > 1
        if isStringScalar(varargin{2}) || ischar(varargin{2})
            if method_ind
                extrapolation = varargin{2};
            else
                method = varargin{2};
            end
        elseif method_ind && isscalar(varargin{2}) && isnumeric(varargin{2})
            extrapolation = varargin{2};
        else
            error('Expected `method` as a second optional argument');
        end
        
        if n_varargs > 2
            extrapolation = varargin{3};
        end
    end
end

do_extrapolation = true;
if isStringScalar(extrapolation) || ischar(extrapolation)
    if ~strcmp(extrapolation, 'extrap') && ~strcmp(extrapolation, 'none')
        error('`extrapolation` must be ''extrap'' if it is a character vector or a string scalar.');
    end
    do_extrapolation = ~strcmp(extrapolation, 'none');
elseif ~(isscalar(varargin{2}) && isnumeric(varargin{2}))
    error('`extrapolation` must be a numeric scalar if it is neither a character vector nor a string scalar.');
end

Y2_passed = ~isempty(Y2);

n_x1 = length(x1);
n_x2 = length(x2);
Y1_resampled = Y1;
Y2_resampled = Y2;
x_new = x1;
if (n_x1 ~= n_x2) || ~all(x1(:) == x2(:))
    if do_extrapolation
        resample_Y = [(n_x2 > n_x1) (Y2_passed & n_x2 < n_x1)];
        if resample_Y(1)
            x_new = x2;
        end
    else
        interval = [max(min(x1), min(x2)), min(max(x1), max(x2))];
        x1_new = x1(x1 >= interval(1) & x1 <= interval(2));
        x2_new = x2(x2 >= interval(1) & x2 <= interval(2));
        n_x1_new = length(x1_new);
        n_x2_new = length(x2_new);
        use_x2 = n_x2_new > n_x1_new;
        if use_x2
            x_new = x2_new;
        else
            x_new = x1_new;
        end
        resample_Y = [
            (use_x2 || (n_x1_new ~= n_x1) || ~all(x1(:) == x1_new(:))),...
            (Y2_passed && (~use_x2 || (n_x2_new ~= n_x2) || ~all(x2(:) == x2_new(:))))
        ];
    end
    n_x_new = length(x_new);
    
    for i = 1:2
        if ~resample_Y(i)
            continue;
        end
        if i == 1
            Z = Y1;
            x_src = x1;
            n_xsrc = n_x1;
        else
            Z = Y2;
            x_src = x2;
            n_xsrc = n_x2;
        end

        sz = size(Z);
        dim_Z = find(sz == n_xsrc);
        sz_new = sz;
        sz_new(dim_Z) = n_x_new;
        permute_vector = [dim_Z, 1:(dim_Z - 1), (dim_Z + 1):ndims(Z)];
        Z_unwrapped = permute(Z, permute_vector);
        Z_resampled = zeros(prod(sz_new), 1);
        n_groups = numel(Z) / n_xsrc;
        offset_old = 1;
        offset_new = 1;
        for j = 1:n_groups
            if do_extrapolation
                Z_resampled(offset_new:(offset_new + n_x_new - 1)) = interp1(...
                        x_src,...
                        Z_unwrapped(offset_old:(offset_old + n_xsrc - 1)),...
                        x_new,...
                        method, extrapolation...
                    );
            else
                Z_resampled(offset_new:(offset_new + n_x_new - 1)) = interp1(...
                        x_src,...
                        Z_unwrapped(offset_old:(offset_old + n_xsrc - 1)),...
                        x_new,...
                        method...
                    );
            end
           offset_old = offset_old + n_xsrc;
           offset_new = offset_new + n_x_new;
        end
        reshape_vector = size(Z_unwrapped);
        Z_resampled = reshape(Z_resampled, [n_x_new, reshape_vector(2:end)]);
        Z_resampled = ipermute(Z_resampled, permute_vector);

        if i == 1
            Y1_resampled = Z_resampled;
        else
            Y2_resampled = Z_resampled;
        end
    end
end

if nargout > 1
    if Y2_passed
        varargout{1} = Y2_resampled;
    else
        varargout{1} = x_new;
    end
    if nargout > 2
        varargout{2} = x_new;
    end
end
end

