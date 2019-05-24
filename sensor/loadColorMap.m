function [sensor_map, channel_mode, bands] = loadColorMap(filename, channel_mode_expected)
% LOADCOLORMAP Load colour space conversion data
%
% ## Syntax
% sensor_map = loadColorMap(filename [, channel_mode])
% [sensor_map, channel_mode] = loadColorMap(filename [, channel_mode])
% [sensor_map, channel_mode, bands] = loadColorMap(filename [, channel_mode])
%
% ## Description
% sensor_map = loadColorMap(filename [, channel_mode])
%   Returns the `sensor_map` variable loaded from the file.
%
% [sensor_map, channel_mode] = loadColorMap(filename [, channel_mode])
%   Additionally returns the `channel_mode` variable loaded from the file.
%
% [sensor_map, channel_mode, bands] = loadColorMap(filename [, channel_mode])
%   Additionally returns the `bands` variable loaded from the file.
%
% ## Input Arguments
%
% filename -- Filename
%   A character vector containing the filename and path of the '.mat' file from
%   which to load the data. `filename` must include the file extension.
%
% channel_mode -- Expected value of `channel_mode`
%   A logical scalar indicating the expected value of `channel_mode` in the
%   file. An error will be thrown if the value of `channel_mode` in the file
%   does not match this value. If the expected value is empty, or `channel_mode`
%   is not passed as an input argument, then the value of `channel_mode` in the
%   file will not be validated.
%
% ## Output Arguments
%
% sensor_map -- Channel sensitivity data
%   A 2D array, where `sensor_map(i, j)` is the sensitivity of the i-th colour
%   channel or spectral band of an output colour space to the j-th colour
%   channel or spectral band of an input colour space. For example, `sensor_map`
%   is a matrix mapping discretized spectral power distributions to RGB colours.
%
% channel_mode -- Type of sensitivity data
%   A Boolean value indicating whether the input colour space is a set of
%   colour channels (true) or a set of spectral bands (false).
%
% bands -- Input colour space
%   A vector containing the wavelengths or colour channel indices corresponding
%   to the second dimension of 'sensor_map'. The file need not contain this
%   variable if `bands` is not requested as an output argument.
%
% See also loadVariables

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created May 23, 2019

nargoutchk(1, 3);
narginchk(1, 2);

if nargin < 2
    channel_mode_expected = [];
end

if ~isempty(channel_mode_expected) || nargout > 1
    channel_mode = loadVariables(filename, 'channel_mode');
    if ~isempty(channel_mode_expected)
        if channel_mode_expected ~= channel_mode
            error(...
                '`channel_mode` in "%s" is expected to be `logical(%d)`.',...
                filename, channel_mode_expected...
            );
        end
    end
end

if nargout > 2
    [sensor_map, bands] = loadVariables(filename, {'sensor_map', 'bands'});
else
    sensor_map = loadVariables(filename, 'sensor_map');
end

end