function s_out = trimCommon(s_in, varargin)
% TRIMCOMMON  Extract unique parts of filepaths
%
% ## Syntax
% s_out = trimCommon(s_in [, ends])
%
% ## Description
% s_out = trimCommon(s_in [, ends])
%   Returns a cell vector of the unique portions of the input filepaths.
%
% ## Input Arguments
%
% s_in -- Filenames
%   A cell vector of character vectors representing filenames and paths.
%
% ends -- Ends to trim
%   A two-element logical vector, where `ends(1)` indicates whether the
%   common start of a filename (excluding the directory path) should be
%   trimmed, and `ends(2)` indicates whether the common ending of a
%   filename (excluding the extension) should be trimmed. Regardless of the
%   value of `ends`, common paths and common extensions will be trimmed.
%
%   Defaults to `[true, true]` if not passed.
%
% ## Output Arguments
%
% s_out -- Trimmed character vectors
%   A cell vector the same dimensions as `s_in`, containing trimmed
%   versions of the character vectors in `s_in`. Starting and ending
%   substrings common to all elements of `s_in` are absent from the
%   elements of `s_out`, depending on the value of `ends` (see above).
%
%   If the input cell vector has length 1, the function will attempt to
%   strip a file extension and path from its first element, returning just
%   the filename without the extension.
%
% See also fileparts

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 20, 2018

narginchk(1, 2);
nargoutchk(1, 1);

ends = true(2, 1);
if ~isempty(varargin)
    ends = varargin{1};
end

if length(s_in) == 1
    s_out = cell(1, 1);
    [~, s_out{1}] = fileparts(s_in{1});
    return;
end

start_ind = length(s_in{1});
end_ind = start_ind;

c_1_start = s_in{1};
c_1_end = s_in{1};
for i = 2:length(s_in)
    c_2 = s_in{i};
    len_2 = length(c_2);
    if start_ind
        if start_ind > len_2
            start_ind = len_2;
            c_1_start = c_1_start(1:start_ind);
        end
        test = (c_2(1:start_ind) ~= c_1_start);
        find_result = find(test, 1);
        if ~isempty(find_result)
            start_ind = find_result - 1;
        end
        c_1_start = c_2(1:start_ind);
    end
    if end_ind
        if end_ind > len_2
            end_ind = len_2;
            c_1_end = c_1_end((end - end_ind + 1):end);
        end
        test = (c_2((end - end_ind + 1):end) ~= c_1_end);
        find_result = find(test, 1, 'last');
        if ~isempty(find_result)
            end_ind = end_ind - find_result;
        end
        c_1_end = c_2((end - end_ind + 1):end);
    end
end
start_ind = start_ind + 1;

[path, s_out{1}, ext] = fileparts(s_in{1});
if ~ends(1)
    if isempty(path)
        start_ind = 1;
    else
        % Add 2 to account for the path separator
        start_ind = min(length(path) + 2, start_ind);
    end
end
if ~ends(2)
    if isempty(ext)
        end_ind = 0;
    else
        % `ext` includes the dot
        end_ind = min(length(ext), end_ind);
    end
end

s_out = cell(size(s_in));
for i = 1:length(s_in)
    s_out{i} = s_in{i}((start_ind):(end - end_ind));
end
    
end
