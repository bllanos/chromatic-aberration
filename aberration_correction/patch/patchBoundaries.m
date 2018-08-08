function [ patch_lim, trim ] = patchBoundaries(...
    image_sampling, patch_size, padding, corner...
)
% PATCHBOUNDARIES  Find the inner and outer boundaries of an image patch
%
% ## Syntax
% [ patch_lim, trim ] = patchBoundaries(...
%   image_sampling, patch_size, padding, corner...
% )
%
% ## Description
% [ patch_lim, trim ] = patchBoundaries(...
%   image_sampling, patch_size, padding, corner...
% )
%   Returns the indices bounding the patch in the image, and the indices
%   within the patch of its inner region
%
% ## Input Arguments
%
% image_sampling -- Image pixel dimensions
%   A two-element vector containing the image height and width in pixels.
%
% patch_size -- Patch size
%   A two-element vector containing the height and width, respectively, of
%   the image patch (a rectangular patch). `patch_size` does not include
%   padding, but is the dimensions of the patch's inner region.
%
% padding -- Patch border width
%   A scalar containing the pixel width of the border surrounding the
%   patch's inner region (on all sides, if not constrained by the image
%   boundaries).
%
% corner -- Patch corner location
%   A two-element vector containing the row and column indices,
%   respectively, of the top-left corner of the patch in the image.
%
% ## Output Arguments
%
% The output arguments describe how the patch with the given size and
% surrounding border is clipped by the image boundaries.
%
% patch_lim -- Patch bounding indices
%   A 2 x 2 array, with the following elements:
%   - `patch_lim(1, 1)` is the row index of the top left corner of the
%     patch in the image.
%   - `patch_lim(2, 1)` is the row index of the bottom right corner of the
%     patch in the image.
%   - `patch_lim(1, 2)` is the column index of the top left corner of the
%     patch in the image.
%   - `patch_lim(2, 2)` is the column index of the bottom right corner of
%     the patch in the image.
%
% trim -- Patch inner region bounding indices
%   A 2 x 2 array, with the following elements, describing the region of
%   the patch inside the boundary determined by `padding` and by the image
%   boundaries:
%   - `trim(1, 1)` is the row index in the patch of the top left corner of
%     the patch's inner region.
%   - `trim(2, 1)` is the row index in the patch of the bottom right corner
%     of the patch's inner region.
%   - `trim(1, 2)` is the column index in the patch of the top left corner
%     of the patch's inner region.
%   - `trim(2, 2)` is the column index in the patch of the bottom right
%     corner of the patch's inner region.
%
% See also solvePatches, solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 7, 2018

nargoutchk(2, 2);
narginchk(4, 4);

patch_lim = [
    corner(1) - padding, corner(2) - padding;
    corner(1) + patch_size(1) + padding - 1, corner(2) + patch_size(2) + padding - 1
    ];
trim = [padding + 1, padding + 1];
if patch_lim(1, 1) < 1
    trim(1, 1) = trim(1, 1) + patch_lim(1, 1) - 1;
    patch_lim(1, 1) = 1;
end
if patch_lim(1, 2) < 1
    trim(1, 2) = trim(1, 2) + patch_lim(1, 2) - 1;
    patch_lim(1, 2) = 1;
end
trim = [trim; trim + patch_size - 1];
if patch_lim(2, 1) > image_sampling(1)
    trim(2, 1) = trim(2, 1) + min(0, image_sampling(1) - patch_lim(2, 1) + padding);
    patch_lim(2, 1) = image_sampling(1);
end
if patch_lim(2, 2) > image_sampling(2)
    trim(2, 2) = trim(2, 2) + min(0, image_sampling(2) - patch_lim(2, 2) + padding);
    patch_lim(2, 2) = image_sampling(2);
end

end