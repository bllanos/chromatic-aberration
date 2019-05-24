function varargout = loadVariables(filename, varnames, empty, none)
% LOADVARIABLES Load variables from a '.mat' file
%
% ## Syntax
% [val1, val2, ... ] = loadVariables(filename, varnames [, empty, none])
%
% ## Description
% [val1, val2, ... ] = loadVariables(filename, varnames [, empty, none])
%   Returns the variables loaded from the file.
%
% ## Input Arguments
%
% filename -- Filename
%   A character vector or string containing the filename and path of the '.mat'
%   file, including the file extension.
%
% varnames -- Variable names
%   The name of the single variable to be loaded from the file (a character
%   vector or string scalar), or a cell vector or string vector of variable
%   names.
%
% empty -- Allow empty variables
%   A flag indicating whether an empty value of the variable is allowed
%   (`true`), or whether an error should be raised if the variable has an empty
%   value (`false`). Defaults to `false` if not passed. If there are multiple
%   variables to load, `empty` can optionally be a logical vector of
%   corresponding per-variable flags.
%
% none -- Allow missing variables
%   A flag indicating whether the variable can be missing from the file
%   (`true`), or whether an error should be raised if the variable is missing
%   (`false`). Defaults to `false` if not passed. If `none` is `true`, and the
%   variable is not found in the file, `val` will be empty (`[]`). If there are
%   multiple variables to load, `empty` can optionally be a logical vector of
%   corresponding per-variable values.
%
% ## Output Arguments
%
% val1, val2, ... -- Variable values
%   The results of loading variables named in `varnames` from the '.mat' file.
%
% ## Notes
% - An error will be thrown if any variables are allowed to be absent from the
%   '.mat' file, but are expected to be non-empty.
%
% See also load, loadImage

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 23, 2019


narginchk(2, 4);

if ischar(varnames)
    varnames = {varnames};
elseif ~(iscell(varnames) || isstring(varnames))
    error('`varnames` must be a character vector, string scalar, string array, or cell vector.');
end
n_vars = length(varnames);
nargoutchk(n_vars, n_vars);

if nargin < 3
    empty = false(n_vars, 1);
elseif length(empty) < n_vars
    if ~isscalar(empty)
        error('`empty` must have the same length as `varnames`, or must be a scalar.');
    end
    empty = repmat(empty, n_vars, 1);
end

if nargin < 4
    none = false(n_vars, 1);
elseif length(none) < n_vars
    if ~isscalar(none)
        error('`none` must have the same length as `varnames`, or must be a scalar.');
    end
    none = repmat(none, n_vars, 1);
end

if any(none & ~empty)
    error('Variables cannot be allowed to be absent, yet expected to be non-empty.');
end

file = matfile(filename, 'Writable', false);
varargout = cell(n_vars, 1);
for i = 1:n_vars
    try
        val = file.(varnames{i});
    catch ME
        switch ME.identifier
            case 'MATLAB:MatFile:VariableNotInFile'
                if none(i)
                    val = [];
                else
                    rethrow(ME)
                end
            otherwise
                rethrow(ME)
        end
    end
    if isempty(val) && ~empty(i)
        error('Empty value of `%s` loaded from "%s".', varnames{i}, filename);
    end
    varargout{i} = val;
end

end