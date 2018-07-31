function [I, name, ext] = loadImage(filename, varargin)
% LOADIMAGE Load an image from an image file or a data file
%
% ## Syntax
% I = loadImage(filename [, variable_name])
% [I, name] = loadImage(filename [, variable_name])
% [I, name, ext] = loadImage(filename [, variable_name])
%
% ## Description
% I = loadImage(filename [, variable_name])
%   Returns the image loaded from the file
% [I, name] = loadImage(filename [, variable_name])
%   Additionally returns the name of the image file, excluding its
%   extension
% [I, name, ext] = loadImage(filename [, variable_name])
%   Additionally returns the file extension of the image file
%
% ## Input Arguments
%
% filename -- Dispersion model filename
%   A character vector containing the filename and path of the '.mat' or
%   image format file containing the image.
%
% variable_name -- Image variable name
%   The name of the variable to be loaded from the file, should `filename`
%   refer to a '.mat' file. (Otherwise `variable_name` can be empty.)
%
% ## Output Arguments
%
% I -- Image
%   The result of loading the variable 'variable_name' from the '.mat' file
%   referred to by `filename`, or the result of calling 'imread()' on the
%   image file referred to by `filename`.
%
% name -- Filename excluding path and extension
%   The name of the file referred to by `filename`.
%
% ext -- File format
%   The extension of the file referred to by `filename`.
%
% See also saveImages, imread, load

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2018

nargoutchk(1, 3);
narginchk(1, 2);

mat_ext = '.mat';

[~, name, ext] = fileparts(filename);
if strcmp(ext, mat_ext)
    if isempty(varargin)
        error('A variable name must be given to load input images from %s files.', mat_ext);
    else
        variable_name = varargin{1};
    end
    load(filename, variable_name);
    if exist(variable_name, 'var')
        I = eval(variable_name);
    else
        error(...
            'The input image variable %s was not loaded from %s.',...
            variable_name, filename...
            );
    end
else
    I = imread(filename);
end

end

