function [ output_files ] = dirreadRAW(...
    in_directory, out_directory, ext, wildcard, regex, ops, varargin...
)
% DIRREADRAW  Read and write raw Bayer pattern images in a directory
%
% ## Syntax
% dirreadRAW(...
%   in_directory, out_directory, ext, wildcard, regex, ops [, varargin]...
% )
% output_files = dirreadRAW(...
%   in_directory, out_directory, ext, wildcard, regex, ops [, varargin]...
% )
%
% ## Description
% dirreadRAW(...
%   in_directory, out_directory, ext, wildcard, regex, ops, varargin...
% )
%   Load and process images, then save them to the output directory
% output_files = dirreadRAW(...
%   in_directory, out_directory, ext, wildcard, regex, ops, varargin...
% )
%   Additionally returns the filepaths and names of the output images
%
% ## Input Arguments
%
% in_directory -- Input directory
%   A character vector containing the path of the directory containing
%   image files to load.
%
% out_directory -- Output directory
%   A character vector containing the path of the directory where output
%   image files are to be saved.
%
% ext -- Output filename extension
%   The filename extension for output images, which `imwrite()` will use to
%   determine the output file format. `ext` does not include the separating
%   dot ('.').
%
% wildcard -- Input filename wildcard
%   A wildcard expression for `ls()`. `wildcard` determines which files in
%   the input directory will be processed.
%
% regex -- Deduplicating regular expression
%   A regular expression which will be removed from the filenames of the
%   input image files. The resulting filenames will be grouped, and a
%   single output image will be generated for each unique filename. The
%   output image will be the average of the processed versions of the
%   grouped files. To disable grouping, pass an empty value for `regex`,
%   such as `''`, or `[]`.
%
%   Note: Filename extensions are removed prior to the regular expression
%   operation.
%
% ops -- Desired processing
%   A structure describing the processing to be performed on the image.
%   `ops` is passed directly to `imreadRAW()`.
%
% varargin -- Additional arguments for `imreadRAW()`
%   The remaining input arguments of `imreadRAW()`, after `ops`.
%
% ## Output Arguments
%
% output_files -- Output image filenames
%   A cell vector of character vectors. Each element contains the full
%   filename and path of an output image saved to the output directory.
%
% ## Notes
% - Only one output image will be requested during each call to
%   `imreadRAW()`. Therefore, if `ops.demosaic` is `true`, this function
%   will generate full colour images. Otherwise, it will generate raw
%   colour filter array data.
% - Note that image averaging induced by `regex` will occur after any
%   requested linearization.
%
% See also ls, imreadRAW, imwrite

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2017

nargoutchk(0,1)
% `narginchk` is done by `imreadRAW`

% Find all filenames
names = ls(fullfile(in_directory, wildcard));
n = size(names, 1);
cleaned_names = cell(n, 1);
for i = 1:n
    cleaned_names{i} = strtrim(names(i, :));
end
cleaned_names = string(cleaned_names);

% Remove filename extensions
processed_names = cleaned_names;
for i = 1:n
    [~, processed_names{i}, ~] = fileparts(processed_names{i});
end

% Find unique filename roots
if isempty(regex)
    unique_names = processed_names;
    unique_name_indices = 1:n;
else
    processed_names = regexprep(processed_names, regex, '');
    [unique_names, ~, unique_name_indices] = unique(processed_names);
end
m = length(unique_names);

% Process images in batches by name root
output_files = cell(m, 1);
for i = 1:m
    names_i = cleaned_names(unique_name_indices == i);
    I = imreadRAW( fullfile(in_directory, char(names_i(1))), ops, varargin{:} );
    n_m = length(names_i);
    for j = 2:n_m
        I = I + imreadRAW( fullfile(in_directory, char(names_i(j))), ops, varargin{:} );
    end
    I = I ./ n_m;
    output_files{i} = fullfile(out_directory, [char(unique_names(i)) '.' ext]);
    imwrite(I, output_files{i});
end

end