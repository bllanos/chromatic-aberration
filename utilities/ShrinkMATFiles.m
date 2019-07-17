%% Convert large double precision variables to single precision variables in '.mat' files
% This script recursively processes the subdirectories of the given
% directory. For each '.mat' file found, it loads class 'double' variables
% larger than 1 MiB, converts them to class 'single', and then updates the
% '.mat' file with the new variable values.
%
% This script is used to reduce the volume of data stored.

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 17, 2019

%% Parameters

directory = '.';

%% Processing

mat_regex = '^.+\.mat$';
min_size = 2^20;
filepaths = listFilesRecursive(mat_regex, directory);

for f = 1:length(filepaths)
    fprintf('File "%s"\n', filepaths{f});
    file = matfile(filepaths{f}, 'Writable', true);
    vars = whos(file);
    
    for v = 1:length(vars)
        var = vars(v);
        name = var.name;
        bytes1 = var.bytes;
        if strcmp(var.class, 'double') && ~var.sparse && bytes1 > min_size
            file.(name) = single(file.(name));
            fprintf('\tUpdated variable "%s"\n', name);
        end
    end     
end