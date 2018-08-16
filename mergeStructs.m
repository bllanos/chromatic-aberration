function sc = mergeStructs(sa, sb, intersect, overwrite)
% MERGESTRUCTS  Merge two structures
%
% ## Syntax
% sc = mergeStructs(sa, sb, intersect, overwrite)
%
% ## Description
% sc = mergeStructs(sa, sb, intersect, overwrite)
%   Returns a structure with fields from both input structures
%
% ## Input Arguments
%
% sa -- First structure
%   A structure
%
% sb -- Second structure
%   A structure
%
% intersect -- Take intersection or union of fields
%   If `true`, the output structure `sc` will have only those fields which
%   are common to `sa` and `sb`. Otherwise, `sc` will have all of the
%   fields of `sa` and all of the fields of `sb`.
%
% overwrite -- Overwrite common fields
%   If `true`, the values of the fields of `sc` that are present on both
%   `sa` and `sb` will be the same as the values of those fields on `sb`.
%   Otherwise, they will take their values from `sa`.
%
% ## Output Arguments
%
% sc -- Merged structure
%   A structure whose set of fields is the intersection or union of the
%   fields in `sa` and the fields in `sb` (depending on the value of
%   `intersect`). `overwrite` determines the values of the fields in cases
%   where both input structures have values for the same fields. Fields for
%   which only one of the input structure has a value obtain their values
%   from that input structure.
%
% See also struct

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 15, 2018

narginchk(4, 4);
nargoutchk(1, 1);

fields_a = fieldnames(sa);
fields_b = fieldnames(sb);
if intersect
    fields_c = intersect(fields_a, fields_b);
else
    fields_c = union(fields_a, fields_b);
end
sc = struct;
for f = 1:length(fields_c)
    field = fields_c{f};
    if overwrite && isfield(sb, field)
        sc.(field) = sb.(field);
    elseif isfield(sa, field)
        sc.(field) = sa.(field);
    elseif isfield(sb, field)
        sc.(field) = sb.(field);
    end
end
    
end
