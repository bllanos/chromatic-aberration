function [ output_files ] = darkSubtract(dir, var_name, wildcards, regex)
% DARKSUBTRACT  Read raw Bayer pattern images and perform averaging and dark frame subtraction
%
% ## Syntax
% output_files = darkSubtract(dir, var_name, wildcards, regex)
%
% ## Description
% output_files = darkSubtract(dir, var_name, wildcards, regex)
%   Load and process images, then save them to the output directory and
%   return the names and paths of the output files. (No output arguments
%   can also be requested.)
%
% ## Input Arguments
%
% dir -- Output directories
%   A structure with the following fields:
%   - 'out_averaged': A character vector containing the path of the
%     directory in which to store averaged dark corrected images. If this
%     field is absent, images will not be averaged, although dark frame
%     images will still be averaged before being subtracted from the other
%     images.
%   - 'out_single': A character vector containing the path of the
%     directory in which to store non-averaged dark corrected images. If
%     this field is absent, these images will not be output.
%
% var_name -- Output variable name
%   A character vector containing the variable name of the variable used to
%   store output images in '.mat' files. (Images are output as '.mat'
%   files, not as image files.)
%
% wildcards -- Input image filename and path wildcards
%   A structure with the following fields:
%   - 'in':  A character vector containing a wildcard expression for
%     `ls()`, used to locate the images that are to be processed.
%   - 'dark':  A character vector containing a wildcard expression for
%     `ls()`, used to locate dark frame image files to load. If this field
%     is absent, dark frame subtraction will be disabled.
%
% regex -- Regular expressions
%   A structure with the following fields:
%   - 'dedup': A regular expression which will be removed from the
%     filenames of the input image files (including dark frame images). The
%     resulting filenames will be grouped into groups of identical
%     filenames. Dark frame images will be averaged by group before being
%     subtracted. If 'dir.out_averaged' exists, an average (dark corrected)
%     image will be generated for each group of non-dark frame images.
%   - 'dark_match': A regular expression which will be used to match dark
%     frame images with the other images. The portions of filenames
%     matching this regular expression must be identical between the dark
%     frame and the other image in order for the dark frame to be
%     subtracted from the other image. Note that matching is done after
%     'dedup' has been used to simplify filenames.
%
%   Filename paths and  extensions are removed prior to regular expression
%   operations.
%
% ## Output Arguments
%
% output_files -- Output image filenames
%   A structure with the following fields:
%   - 'out_averaged': A cell vector of character vectors containing the
%     names and paths of the averaged dark corrected images. This field is
%     present only if 'dir.out_averaged' exists.
%   - 'out_single': A cell vector of character vectors containing the
%     names and paths of the non-averaged dark corrected images. This field
%     is present only if 'dir.out_single' exists.
%
% ## Notes
% - Image averaging and dark frame subtraction will occur after any
%   linearization of the RAW images.
%
% See also listFiles, imreadRAW, dirreadRAW

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 3, 2018

nargoutchk(0, 1)
narginchk(4, 4)

do_average = isfield(dir, 'out_averaged');
do_output_single = isfield(dir, 'out_single');
if ~do_average && ~do_output_single
    error('No output images requested.');
end
do_subtract = isfield(wildcards, 'dark');

% Find all filenames
other_paths = listFiles(wildcards.in);
n_other = length(other_paths);
all_paths = {other_paths};
if do_subtract
    dark_paths = listFiles(wildcards.dark);
    all_paths = [all_paths, {dark_paths}];
end
n_path_types = length(all_paths);

% Remove filename paths and extensions
names = cell(n_path_types, 1);
for i = 1:n_path_types
    names{i} = cell(size(all_paths{i}));
    for j = 1:length(names{i})
        [~, names{i}{j}, ~] = fileparts(all_paths{i}{j});
    end
end

% Find unique filename roots
other_grouping = regexprep(names{1}, regex.dedup, '');
[other_grouping, ~, other_grouping_indices] = unique(other_grouping);
n_other_grouping = length(other_grouping);
if do_subtract
    dark_grouping = regexprep(names{2}, regex.dedup, '');
    [dark_grouping, ~, dark_grouping_indices] = unique(dark_grouping);
    n_dark_grouping = length(dark_grouping);
end

% Match dark frames with other images
if do_subtract
    other_tokens = regexp(other_grouping, regex.dark_match, 'tokens', 'forceCellOutput');
    other_substr = cell(n_other_grouping, 1);
    for i = 1:n_other_grouping
        other_substr{i} = cat(2, other_tokens{i}{1}{:});
    end
    dark_tokens = regexp(dark_grouping, regex.dark_match, 'tokens', 'forceCellOutput');
    dark_substr = cell(n_dark_grouping, 1);
    for i = 1:n_dark_grouping
        dark_substr{i} = cat(2, dark_tokens{i}{1}{:});
    end
    group_pairs = zeros(n_other_grouping, 1);
    for i = 1:n_other_grouping
        pairs_i = strcmp(other_substr{i}, dark_substr);
        if sum(pairs_i) ~= 1
            error('Expected one dark frame group for the group "%s", not %d.',...
                other_grouping{i}, sum(pairs_i));
        end
        group_pairs(i) = find(pairs_i);
    end
end

% Options for RAW image loading
ops.linearize = true;
ops.demosaic = false;
ops.convertColor = false;
ops.wb = false;

% Process images in batches by name root
output_files = struct;
if do_output_single
    output_files.out_single = cell(n_other, 1);
end
if do_average
    output_files.out_averaged = cell(n_other_grouping, 1);
end
for i = 1:n_other_grouping
    % Calculate the averaged dark frame corresponding to this group
    if do_subtract
        dark_names_i = dark_paths(dark_grouping_indices == group_pairs(i));
        I_dark = imreadRAW(dark_names_i{1}, ops);
        n_i = length(dark_names_i);
        for j = 2:n_i
            I_dark = I_dark + imreadRAW(dark_names_i{j}, ops);
        end
        I_dark = I_dark ./ n_i;
    end
    
    ind_i = find(other_grouping_indices == i);
    n_i = length(ind_i);
    
    I = imreadRAW(other_paths{ind_i(1)}, ops);
    if do_output_single
        if do_subtract
            I_out = I - I_dark;
        else
            I_out = I;
        end
        output_files.out_single{ind_i(1)} = saveImages(...
            'data', dir.out_single, names{1}{ind_i(1)}, I_out, '', var_name...
        );
    end
    I_sum = I;
    for j = 2:n_i
        I = imreadRAW(other_paths{ind_i(j)}, ops);
        if do_output_single
            if do_subtract
                I_out = I - I_dark;
            else
                I_out = I;
            end
            output_files.out_single{ind_i(j)} = saveImages(...
                'data', dir.out_single, names{1}{ind_i(j)}, I_out, '', var_name...
            );
        end
        if do_average
            I_sum = I_sum + I;
        end
    end
    if do_average
        I = I_sum ./ n_m;
        if do_subtract
            I_out = I - I_dark;
        else
            I_out = I;
        end
        output_files.out_averaged{i} = saveImages(...
            'data', dir.out_averaged, other_grouping{i}, I_out, '', var_name...
        );
    end
end

end