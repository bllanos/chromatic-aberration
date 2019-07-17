function [ files ] = listFilesRecursive(regex, directory)
% LISTFILESRECURSIVE  List files matching a regular expression in a directory tree
%
% ## Syntax
% names = listFilesRecursive(regex)
% names = listFilesRecursive(regex, directory)
%
% ## Description
% names = listFilesRecursive(regex)
%   Returns file names matching the regular expression as a cell array. The
%   current directory and all of its subdirectories are recursively
%   searched.
%
% names = listFilesRecursive(regex, directory)
%   Returns file names matching the regular expression as a cell array. The
%   given directory and all of its subdirectories are recursively searched.
%
% ## Input Arguments
%
% regex -- Input filename regular expression
%   A regular expression to be matched with filenames (and not with their paths).
%
% directory -- Starting directory path
%   A character vector or string scalar containing the path of the root
%   directory for the recursive search for files with names matching
%   `regex`.
%
% ## Output Arguments
%
% names -- Filepaths and names
%   A cell column vector, where each element is a character containing the
%   name and path of a file with a name matching `regex`. Directories are
%   not included in the output, even if their names match `regex`.
%
% ## Notes
% - Note that the input argument of 'listFiles()', in contrast, is not a regular
%   expression, but an expression of the form handled by the
%   platform-specific 'ls' function.
% - Regular expression matching is case-sensitive.
%
% See also listFiles, dir, ls

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 17, 2019

nargoutchk(1,1)
narginchk(1,2)

% Explore directories by depth-first search
if nargin > 1
    subdirectory_paths = {directory};
else
    subdirectory_paths = {'.'};
end
subdirectory_stack = 0;
subdirectory_listings = { dir(subdirectory_paths{1}) };
if isempty(subdirectory_listings{1})
    error('"%s" is not a valid directory.', subdirectory_paths{1});
end

files = {};
done = false;
while ~done
    subdirectory_stack(end) = subdirectory_stack(end) + 1;
    current_name = subdirectory_listings{end}(subdirectory_stack(end)).name;
    if subdirectory_listings{end}(subdirectory_stack(end)).isdir 
        if ~any(strcmp(current_name, {'.', '..'}))
            % Push the subdirectory
            subdirectory_paths{end + 1} = fullfile(subdirectory_paths{end}, current_name);
            subdirectory_stack(end + 1) = 0;
            subdirectory_listings{end + 1} = dir(subdirectory_paths{end});
        end
    elseif ~isempty(regexp(current_name, regex, 'once'))
        files{end + 1, 1} = fullfile(subdirectory_paths{end}, current_name);
    end
    while subdirectory_stack(end) == length(subdirectory_listings{end})
        % Pop the subdirectory
        subdirectory_paths = subdirectory_paths(1:(end - 1));
        subdirectory_stack = subdirectory_stack(1:(end - 1));
        subdirectory_listings = subdirectory_listings(1:(end - 1));
        if isempty(subdirectory_stack)
            done = true;
            break
        end
    end
end

end