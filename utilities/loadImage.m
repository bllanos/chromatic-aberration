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
%   image format file containing the image, including the file extension.
%
% variable_name -- Image variable name
%   The name of the variable to be loaded from the file, should `filename`
%   refer to a '.mat' file. (Otherwise `variable_name` can be empty.)
%
% ## Output Arguments
%
% I -- Image
%   The result of loading, depending on the type of file referred to by `filename`,
%   - The variable 'variable_name' from a '.mat' file
%     - An error will be thrown if `variable_name` is not a numeric
%       variable.
%     - Integer variables will be converted to double precision using
%       'im2double()'.
%     - 'single' variables will be converted to double precision using
%       'double()'.
%   - The result of calling 'imread()' on an image file.
%     Images loaded with 'imread()' are converted to double precision using
%     'im2double()'.
%   - The non-RGB channels of an OpenEXR format hyperspectral image
%     ('.exr'), in the order in which they are given in the list of
%     channels returned by 'exrinfo()'.
%
% name -- Filename excluding path and extension
%   The name of the file referred to by `filename`.
%
% ext -- File format
%   The extension of the file referred to by `filename`.
%
% ## References
% - The 'exr*()' functions are from a GitHub repository of MATLAB bindings
%   for the OpenEXR file format: https://github.com/KAIST-VCLAB/openexr-matlab
%   The repository is connected with the following article:
%
%   Choi, I., Jeon, D. S., Nam, G., Gutierrez, D., & Kim, M. H. (2017).
%     "High-Quality Hyperspectral Reconstruction Using a Spectral Prior."
%     ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2017), 36(6),
%     218:1–13. doi:10.1145/3130800.3130810
%
% See also saveImages, imread, load, im2double

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2018

nargoutchk(1, 3);
narginchk(1, 2);

mat_ext = '.mat';
exr_ext = '.exr';

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
        if ~isnumeric(I)
            error(...
                'The input image variable `%s` in %s does not have a numeric class.',...
                variable_name, filename...
            );
        elseif isinteger(I)
            I = im2double(I);
        elseif isa(I, 'single')
            I = double(I);
        end
    else
        error(...
            'The input image variable %s was not loaded from %s.',...
            variable_name, filename...
            );
    end
elseif strcmp(ext, exr_ext)
    disp('Loading an OpenEXR image...');
    info = exrinfo(filename);
    channel_names = info.channels;
    channel_names_filter = ~(...
        strcmp('R',channel_names) |...
        strcmp('G',channel_names) |...
        strcmp('B',channel_names)...
    );
    channel_names = channel_names(channel_names_filter);
    I_map = exrreadchannels(filename, channel_names{:});
    I = zeros([info.size, length(channel_names)]);
    for i = 1:length(channel_names)
        I(:, :, i) = I_map(channel_names{i});
    end
    fprintf('\tDone\n');
else
    I = im2double(imread(filename));
end

end

