function td = mergeDispersionModelROI(varargin)
% MERGEDISPERSIONMODELROI Find the intersection of dispersion model domains
%
% ## Syntax
% td = mergeDispersionModelROI(td_1 [td_2, td_3, ...])
%
% ## Description
% td = mergeDispersionModelROI(td_1 [td_2, td_3, ...])
%   Returns dispersion model coordinate system conversion information with
%   a domain that is the intersection of the domains of the input
%   coordinate system conversion information structures.
%
% ## Input Arguments
%
% td_1, td_2, td_3, ... -- Spatial coordinate conversion information structures
%   Structures of the form of the `transform_data` output argument of
%   'loadDispersionModel()'. All structures must differ only in their
%   `model_space.corners` fields.
%
% ## Output Arguments
%
% td -- Spatial coordinate conversion information
%   A structure which is a copy of `td_1`, if the intersection of the
%   regions defined by the `model_space.corners` fields of the input
%   structures is empty. Otherwise, `td.model_space.corners` defines the
%   region corresponding to the intersection of the regions defined by the
%   input structures.
%
% See also loadDispersionModel, makeDispersionForImage, modelSpaceTransform

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created December 14, 2018

nargoutchk(1, 1);
if nargin < 1
    error('At least one dispersion model coordinate conversion information structure must be passed.');
end

td = varargin{1};
if strcmp(td.model_space.system, 'image')
    is_geometric = false;
elseif strcmp(td.model_space.system, 'geometric')
    is_geometric = true;
else
    error('Unrecognized value of the `model_space.system` field in the first input structure.')
end
poly_current = polyshape(...
    td.model_space.corners([1, 1, 2, 2]),...
    td.model_space.corners([3, 4, 4, 3])...
);

for i = 2:length(varargin)
    td_i = varargin{i};
    if td.fill && td_i.fill
        % Compatible - The region of interest is the entire image domain
    elseif ~td.fill && ~td_i.fill
        if ~strcmp(td.model_space.system, td_i.model_space.system)
            error('All structures must have the same `model_space.system` strings.');
        end
        if is_geometric && (td.model_space.pixel_size ~= td_i.model_space.pixel_size)
            error('The structures do not have the same pixel sizes.')
        end
        poly_i = polyshape(...
            td_i.model_space.corners([1, 1, 2, 2]),...
            td_i.model_space.corners([3, 4, 4, 3])...
        );
        poly_current = intersect(poly_current, poly_i);
    else
        error('All models of dispersion must have the same `fill` fields.');
    end
end

if area(poly_current)
    vertices = poly_current.Vertices;
    if is_geometric
        td.model_space.corners = [
            min(vertices(:, 1)), max(vertices(:, 2));
            max(vertices(:, 1)), min(vertices(:, 2))
        ];
    else
        td.model_space.corners = [
            min(vertices(:, 1)), min(vertices(:, 2));
            max(vertices(:, 1)), max(vertices(:, 2))
        ];
    end
end

end