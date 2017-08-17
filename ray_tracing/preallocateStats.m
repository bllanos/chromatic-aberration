function [ empty_stats ] = preallocateStats( sz )
% PREALLOCATESTATS  Return an empty PSF statistics structure array
%
% ## Syntax
% empty_stats = preallocateStats( sz )
%
% ## Description
% empty_stats = preallocateStats( sz )
%   Returns a structure array, for preallocation.
%
% ## Input Arguments
%
% sz -- Structure array dimensions
%   A vector containing the dimensions of the `empty_stats` output array,
%   such that `all(size(empty_stats) == sz)` is true.
%
% ## Output Arguments
%
% empty_stats -- Empty statistics structure array
%   A structure array with empty field values, but with the same fields as
%   the `stats` output argument of 'analyzePSF()' or 'analyzePSFImage()'.
%   `empty_stats` is useful for preallocating a structure array, to be
%   filled with calls to 'analyzePSF()' or 'analyzePSFImage()'.
%
% See also analyzePSF, analyzePSFImage

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 17, 2017

nargoutchk(1, 1);
narginchk(1, 1);

empty_stats = struct(...
    'mean_position', cell(sz),...
    'mean_value', cell(sz),...
    'max_position', cell(sz),...
    'max_value', cell(sz),...
    'radius', cell(sz)...
    );

end

