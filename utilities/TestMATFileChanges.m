%% Test script to compare the '.mat' files in two directories
% This script checks whether there are any '.mat' files found in one
% directory which are not in the other directory (by matching filenames).
% It then iterates through the '.mat' files in one directory. For each
% '.mat' file, it checks if the values of the numeric variables it contains
% match the values of the variables of the same names in the corresponding
% '.mat' file in the other directory. It also checks whether the two '.mat'
% files contain the same lists of variables.
%
% This script is useful for testing if the output of a script has changed,
% following some modifications to the codebase.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 13, 2018

%% parameters
clc;
clear;

dir_1 = '/home/llanos/Downloads/old';
dir_2 = '/home/llanos/Downloads/new';

%% Test for identical filenames

mat_ext = '.mat';
mat_wildcard = ['*' mat_ext];

files_1 = listFiles(fullfile(dir_1, mat_wildcard));
files_2 = listFiles(fullfile(dir_2, mat_wildcard));

if length(files_1) ~= length(files_2)
    error('The directories contain different numbers of %s files.', mat_ext)
end

filenames_1 = cell(length(files_1), 1);
for f = 1:length(files_1)
    [~, filename_1, ~] = fileparts(files_1{f});
    [~, filename_2, ~] = fileparts(files_2{f});
    if ~strcmp(filename_1, filename_2)
        error('The file %s is not found in the second directory.', files_1{f})
    end
    filenames_1{f} = filename_1;
end

%% Test for identical file contents

for f = 1:length(files_1)
    fprintf('Looking at the versions of file "%s" across the two directories...\n', filenames_1{f});
    m1 = matfile(files_1{f}, 'Writable', false);
    m2 = matfile(files_2{f}, 'Writable', false);
    vars1 = whos(m1);
    vars2 = whos(m2);
    
    if length(vars1) ~= length(vars2)
        error('The versions contain different numbers of variables.')
    end

    for v = 1:length(vars1)
        if ~strcmp(vars1(v).name, vars2(v).name)
            error('The variable %s is not found in the second version.', vars1(v).name)
        end
    end
    
    for v = 1:length(vars1)
        var1 = vars1(v);
        var2 = vars2(v);
        name = var1.name;
        sz1 = var1.size;
        sz2 = var2.size;
        if any(sz1 ~= sz2)
            error('The variable %s has different sizes.', name)
        end
        if ~strcmp(var1.class, var2.class)
            error('The variable %s has different classes: %s and %s.', name, var1.class, var2.class)
        end
        if any(strcmp(var1.class, {'single', 'double', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'}))
            val1 = m1.(name);
            val2 = m2.(name);
            if any(val1(:) ~= val2(:))
                error('The numeric variable %s has different values.', name);
            end
        end
    end     
end