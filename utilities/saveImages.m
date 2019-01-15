function filepaths = saveImages(varargin)
% SAVEIMAGES Save images to image files and data files
%
% ## Syntax
% saveImages(...
%   directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
% saveImages(...
%   'data', directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
% saveImages(...
%   'image', directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
% filepaths = saveImages(____)
%
% ## Description
% saveImages(...
%   directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
%   Saves the images to '.mat' files and image files.
%
% saveImages(...
%   'data', directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
%   Saves the images to '.mat' files only.
%
% saveImages(...
%   'image', directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
%   Saves the images to image files only, except that images with more than
%   three channels are saved to '.mat' files only.
%
% filepaths = saveImages(____)
%   Returns the paths and names of the saved files.
%
% ## Input Arguments
%
% directory -- Output directory
%   A character vector containing the path of the directory in which to
%   save the files. `directory` cannot be `'data'` or `'image'`.
%
% name -- Partial filename
%   The portion common to all output file names (does not include a file
%   extension).
%
% I1, I2, ... -- Images
%   Images with arbitrary numbers of channels, to be saved to files.
%
% I1_postfix, I2_postfix, ... -- Filename postfixes
%   The portion unique to the file names of the corresponding images, to be
%   affixed to the ends of file names (does not include a file extension).
%
% I1_variable, I2_variable, ... -- Image variable name
%   The variable names to use when saving the corresponding images to
%   '.mat' files.
%
% ## Output
%
% Up to two files will be generated for each image having one or three
% channels, depending on what kinds of output files were requested: One
% '.mat' file, where `IN` is stored under the variable name `IN_variable`,
% and one image file. The filenames of these files will be postfixed with
% 'IN_postfix'. Images having other numbers of channels will only be saved
% to '.mat' files, and will be saved to '.mat' files even if the first
% input argument is 'image'.
%
% filepaths -- File paths and names
%   A cell array of character vectors containing the full filepaths and
%   names of the output files. The elements of `filepaths` will be ordered
%   first by image, and then by output file type ('.mat' file and/or image
%   file).
%
% See also loadImage, imwrite, save

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 31, 2018

nargoutchk(0, 1);
if isempty(varargin)
    error('No input arguments passed.')
elseif ~isStringScalar(varargin{1}) && ~ischar(varargin{1})
    error('The first input argument must be a character vector or a string scalar.');
end
data_only = strcmp(varargin{1}, 'data');
images_only = strcmp(varargin{1}, 'image');
data_and_images = ~data_only && ~images_only;
if data_and_images
    offset = 2;
else
    offset = 3;
end
if (nargin - offset) < 3
    error('At least one image, filename postfix, and variable name triplet must be provided.');
end
if ~isStringScalar(varargin{offset - 1}) && ~ischar(varargin{offset - 1})
    error('The `directory` input argument must be a character vector or a string scalar.');
end
directory = varargin{offset - 1};
if ~isStringScalar(varargin{offset}) && ~ischar(varargin{offset})
    error('The `name` input argument must be a character vector or a string scalar.');
end
name = varargin{offset};
if mod(nargin - offset, 3) ~= 0
    error('Three arguments must be passed for each image to be saved.');
end

img_ext = '.tif';
mat_ext = '.mat';
image_arguments = varargin((offset + 1):end);
n_images = round(length(image_arguments) / 3);

output = struct;
filepaths = cell(n_images * 2, 1);
used_filepaths = false(n_images * 2, 1);
for i = 1:n_images
    filename = fullfile(directory, [name image_arguments{i * 3 - 1}]);
    I = image_arguments{i * 3 - 2};
    sz = size(I);
    if length(sz) < 2
        error('Image %d is one-dimensional.', i)
    elseif length(sz) > 3
        error('Image %d is four (or more) dimensional.', i)
    end
    can_write_image = (length(sz) == 2 || sz(3) == 1 || sz(3) == 3);
    if can_write_image && ~data_only
        filepaths{i * 2} = [filename img_ext];
        imwrite(I, filepaths{i * 2});
        used_filepaths(i * 2) = true;
    end
    if ~can_write_image || ~images_only
        I1_variable = image_arguments{i * 3};
        output.(I1_variable) = I;
        filepaths{(i * 2) - 1} = [filename mat_ext];
        save(filepaths{(i * 2) - 1}, '-struct', 'output', I1_variable);
        used_filepaths(i * 2 - 1) = true;
    end
end
filepaths = filepaths(used_filepaths);

end

