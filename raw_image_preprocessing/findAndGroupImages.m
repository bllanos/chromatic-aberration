function [...
    grouped_filenames, path, group_names, bands...
] = findAndGroupImages(wildcard, varargin)
% FINDANDGROUPIMAGES  List images and optionally group them by spectral band
%
% ## Syntax
% [...
%   grouped_filenames, path, group_names...
% ] = findAndGroupImages(wildcard)
% [...
%   grouped_filenames, path, group_names, bands...
% ] = findAndGroupImages(wildcard, regex)
%
% ## Description
%
% [...
%   grouped_filenames, path, group_names...
% ] = findAndGroupImages(wildcard)
%   Returns one to three output arguments: The files corresponding to the
%   wildcard, the directory containing the files, and the names of the
%   files (stripped of paths and extensions).
%
% [...
%   grouped_filenames, path, group_names, bands...
% ] = findAndGroupImages(wildcard, regex)
%   Returns one to four output arguments: The files corresponding to the
%   wildcard, in groups, the directory containing the files, the names
%   of the groups, and the spectral bands discovered in the filenames.
%
% ## Input Arguments
%
% wildcard -- Filename and path wildcards
%   A character vector containing a wildcard expression for `ls()`, used to
%   locate the files that are to be processed. The wildcard is expected to
%   locate files in only one directory.
%
% regex -- Regular expression
%   An character vector containing a regular expression matching the
%   portion of filenames containing spectral band information. `regex` must
%   contain a capturing group to extract the spectral band. For example,
%   `'_(\d+)nm'` will match `'_600nm'`, and will capture the token `'600'`.
%
%   It is expected that only a single substring in each filename (without
%   the path and file extension) matches `regex`.
%
% ## Output Arguments
%
% grouped_filenames -- Filenames and paths
%   A cell vector of cell vectors of character vectors containing the full
%   filenames and paths of the files discovered using `wildcard`. Each cell
%   contains a group of filepaths, where a group is a set of filepaths
%   which become identical when the portion matching `regex` is removed.
%   All groups are expected to contain the same number of files, and all
%   groups are expected to contain the same unique portions matching
%   `regex`. In other words, there should be one image for each spectral
%   band in each scene.
%
%   If `regex` is not passed, `grouped_filenames` is a cell vector of
%   single cells containing character vectors. In other words, files are
%   not grouped together.
%
% path -- Directory
%   A character vector containing the path to the directory containing the
%   files referred to by `grouped_filenames`.
%
% group_names -- Partial filenames
%   A cell vector with the same length as `grouped_filenames`. Each cell
%   contains a character vector obtained by stripping the path and file
%   extension from a filepath in the corresponding cell of
%   `grouped_filenames`, and then removing the portion of the filename
%   matching `regex`. In other words, `group_names` contains the partial
%   filenames which are common to each group of filenames.
%
% bands -- Spectral bands
%   A vector with the same length as the contents of each cell of
%   `grouped_filenames`, containing the spectral bands captured from
%   filenames using `regex`. `bands` is a sorted vector of wavelengths.
%
% See also darkSubtract, listFiles

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 12, 2018

narginchk(1, 2);

do_group = false;
if ~isempty(varargin)
    if isempty(varargin{1})
        error('If `regex` is passed, it must not be empty.');
    end
    do_group = true;
    regex = varargin{1};
end
if do_group
    nargoutchk(1, 4);
else
    nargoutchk(1, 3);
end

filenames = listFiles(wildcard);
n_images = length(filenames);

% Remove filename paths and extensions
names = cell(n_images, 1);
for i = 1:n_images
    [path, names{i}, ~] = fileparts(filenames{i});
end

if do_group
    % Find wavelength information
    bands = zeros(n_images, 1);
    tokens = regexp(names, regex, 'tokens', 'forceCellOutput');
    for i = 1:n_images
        if isempty(tokens{i}) || isempty(tokens{i}{1}{1})
            error('No wavelength information found in filename "%s".', filenames{i});
        end
        bands(i) = str2double(tokens{i}{1}{1});
        if ~isfinite(bands(i))
            error('Error converting wavelength information, "%s", from filename "%s", to a number.',...
                tokens{i}{1}{1}, filenames{i});
        end
    end
    [bands, ~, bands_filenames_map] = unique(bands);
    n_bands = length(bands);

    % Group images according to scene
    names_stripped = regexprep(names, regex, '');
    [group_names, ~, group_names_indices] = unique(names_stripped);
    n_groups = length(group_names);
    
    % Make sure there are images for every band in every scene, and sort
    % images by group, then by band
    grouped_filenames = cell(n_groups, 1);
    for g = 1:n_groups
        group_filter = (group_names_indices == g);
        grouped_filenames{g} = filenames(group_filter);
        bands_in_group = bands_filenames_map(group_filter);
        if (length(bands_in_group) ~= n_bands) || any(sort(bands_in_group) ~= (1:n_bands).')
            error('Not all spectral bands are represented exactly once in scene "%s".',...
                group_names{g});
        end
        [~, sorting_map] = sort(bands_in_group);
        grouped_filenames{g} = grouped_filenames{g}(sorting_map);
    end
else
    grouped_filenames = cell(n_images, 1);
    for i = 1:n_images
        grouped_filenames{i} = filenames(i);
    end
    group_names = names;
end

end