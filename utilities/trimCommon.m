function s_out = trimCommon(s_in)
% TRIMCOMMON  Extract unique parts of character vectors
%
% ## Syntax
% s_out = trimCommon(s_in)
%
% ## Description
% s_out = trimCommon(s_in)
%   Returns a cell vector of the unique portions of the input character
%   vectors.
%
% ## Input Arguments
%
% s_in -- Character vectors
%   A cell vector of character vectors.
%
% ## Output Arguments
%
% s_out -- Trimmed character vectors
%   A cell vector the same dimensions as `s_in`, containing trimmed
%   versions of the character vectors in `s_in`. Starting and ending
%   substrings common to all elements of `s_in` are absent from the
%   elements of `s_out`.
%
% See also fileparts

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 20, 2018

narginchk(1, 1);
nargoutchk(1, 1);

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

s_out = cell(size(s_in));
for i = 1:length(s_in)
    s_out{i} = s_in{i}((start_ind + 1):(end - end_ind));
end
    
end
