function I_out = channelConversion(I_in, sensitivity, varargin)
% CHANNELCONVERSION  Convert an image to a different colour/spectral space
%
% ## Syntax
% I_out = channelConversion(I_in, sensitivity [, dim])
%
% ## Description
% I_out = channelConversion(I_in, sensitivity [, dim])
%   Returns the image in a different colour or spectral space.
%
% ## Input Arguments
%
% I_in -- Input image
%   A 2D or 3D array containing an image. See below in the documentation of
%   `dim` for valid formats of 2D arrays.
%
% sensitivity -- Colour or spectral space conversion matrix
%   A 2D array, where `sensitivity(i, j)` is the sensitivity of the i-th
%   colour channel or spectral band of `I_out` to the j-th input colour
%   channel or spectral band of `I_in`. `sensitivity` is a matrix mapping
%   colours in `I_in` to colours in `I_out`. `sensitivity` must account for
%   any numerical intergration that is part of colour conversion.
%
% dim -- Spectral dimension
%   The dimension corresponding to spectral information in the input image.
%   `dim` must be passed to disambiguate two possible cases wherein `I_in`
%   has two dimensions:
%   - `I_in` may be a column vector, where channels/bands are stored in
%     sequence, as in the case where `I_in` represents the unwrapping of a
%     3D array. Values in `I_in` are ordered first by pixel index, then by
%     channel/band. In this case, `dim` should be `1`.
%   - `I_in` may be a 2D array, where rows represent pixels, and columns
%     represent channels/bands. In this case, `dim` should
%     be `2`.
%   - `I_in` may be a single channel/band image represented as a 2D array,
%     of pixel values, in which case `dim` should be `3`.
%
%   Defaults to `3` if not passed, which is also the value corresponding to
%   `I_in` having three dimensions, as in the case of a conventional
%   multi-channel/band image.
%
% ## Output Arguments
%
% I_out -- Converted image
%   The equivalent of the input image in the new spectral or colour space,
%   generated using the colour space conversion data in `sensitivity`. The
%   size of `I_out` will differ from the size of `I_in` if `sensitivity` is
%   not square, but the layout of `I_out` will have the same interpretation
%   as that of `I_in`.
%
% See also channelConversionMatrix, imageFormation

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 20, 2018

nargoutchk(1, 1);
narginchk(2, 3);

if all(ndims(I_in) ~= [2 3])
    error('`I_in` must be a 2D or 3D array.');
end

if isempty(varargin)
    dim = 3;
else
    dim = varargin{1};
end
if all(dim ~= 1:3)
    error('`dim` can only be an integer from `1` to `3`.');
end
if dim ~= 3 && ndims(I_in) == 3
    error('For a 3D array `I_in`, `dim` is expected to be `3`.');
end

n_channels_in = size(sensitivity, 2);
if dim == 1
    if size(I_in, 2) ~= 1
        error('`I_in` must be a column vector if `dim` is `1`.');
    end
    n_px = size(I_in, 1) / n_channels_in;
    if n_px ~= round(n_px)
        error('`I_in` does not have the right number of values to be compatible with `sensitivity`.');
    end
    n_px = round(n_px);
elseif size(I_in, dim) ~= n_channels_in
    error('The number of columns of `sensitivity` must equal the size of `I_in` in its %d-th dimension (`dim`).', dim);
end

switch dim
    case 1
        I_out = reshape(...
            (sensitivity * (reshape(I_in, n_px, n_channels_in).')).',...
            [], 1 ...
        );
    case 2
        I_out = (sensitivity * (I_in.')).';
    case 3
        I_out = reshape(...
            (sensitivity * (reshape(I_in, [], n_channels_in).')).',...
            size(I_in, 1), size(I_in, 2), size(sensitivity, 1)...
        );
    otherwise
        error('Unrecognized value of `dim`.');
end
    
end