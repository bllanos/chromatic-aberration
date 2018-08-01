function saveImages(directory, name, varargin)
% SAVEIMAGES Save images to image files and data files
%
% ## Syntax
% saveImages(...
%   directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
%
% ## Description
% saveImages(...
%   directory, name [, I1, I1_postfix, I1_variable, I2, I2_postfix,
%   I2_variable, ...] ...
% )
%   Saves the images to '.mat' files and image files.
%
% ## Input Arguments
%
% directory -- Output directory
%   A character vector containing the path of the directory in which to
%   save the files.
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
% Two files will be generated for each image having one or three channels:
% one '.mat' file, where `IN` is stored under the variable name
% `IN_variable`, and one image file. The filenames of these files will
% be postfixed with 'IN_postfix'. Images having other numbers of channels
% will only be saved to '.mat' files.
%
% See also loadImage, imwrite, save

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 31, 2018

nargoutchk(0, 0);
if nargin < 5
    error('At least one image must be provided for writing to files.');
elseif mod(nargin - 2, 3) ~= 0
    error('Three arguments must be passed for each image to be saved.');
end

img_ext = '.tif';
mat_ext = '.mat';
n_images = round((nargin - 2) / 3);

output = struct;
for i = 1:n_images
    filename = fullfile(directory, [name varargin{(i-1) * 3 + 2}]);
    I = varargin{(i-1) * 3 + 1};
    sz = size(I);
    if length(sz) < 2
        error('Image %d is one-dimensional.', i)
    elseif length(sz) > 3
        error('Image %d is four-dimensional.', i)
    end
    if length(sz) == 2 || sz(3) == 1 || sz(3) == 3
        imwrite(I, [filename img_ext]);
    end
    I1_variable = varargin{(i-1) * 3 + 3};
    output.(I1_variable) = I;
    save([filename mat_ext], '-struct', 'output', I1_variable);
end

end

