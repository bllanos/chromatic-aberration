function [ names ] = listFiles(wildcard)
% LISTFILES  List files matching a wildcard
%
% ## Syntax
% names = listFiles(wildcard)
%
% ## Description
% names = listFiles(wildcard)
%   Return file names matching the wildcard as a cell array.
%
% ## Input Arguments
%
% wildcard -- Input filename wildcard
%   A wildcard expression for `ls()`, containing the path, if the desired
%   files are located in a different directory from the current directory.
%   (i.e. The wildcard is used as-is - This function does not add any
%   contextual information.)
%
% ## Output Arguments
%
% names -- Filepaths and names
%   A cell vector, where each element is a character containing the name
%   and path of a file in the output of `ls()`.
%
% ## Notes
% - This function treats Windows systems as a special case, and assumes
%   that `ls()` has Linux-specific behaviour on all other systems.
%
% See also ls

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created January 23, 2018

nargoutchk(1,1)
narginchk(1,1)

if ispc
    filepath = fileparts(wildcard);
    raw_names = ls(wildcard);
    names = cell(size(raw_names, 1), 1);
    for i = 1:size(raw_names, 1)
        names{i} = fullfile(filepath, strtrim(raw_names(i, :)));
    end
else
    names = strtrim(strsplit(ls(wildcard), {'  ','\f','\n','\r','\t','\v'}));
    names = names(1:(end - 1)); % There is always a terminating newline
end

end


