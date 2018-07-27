function [dispersion_data, bands, transform_data] = loadDispersionModel(filename, mfr_expected, varargin)
% LOADDISPERSIONMODEL Load dispersion information from a file
%
% ## Syntax
% dispersion_data = loadDispersionModel(filename, mfr_expected [, bands_req])
% [dispersion_data, bands] = loadDispersionModel(filename, mfr_expected [, bands_req])
% [dispersion_data, bands, transform_data] = loadDispersionModel(filename, mfr_expected [, bands_req])
%
% ## Description
% dispersion_data = loadDispersionModel(filename, mfr_expected [, bands_req])
%   Returns the dispersion model from the file
%
% [dispersion_data, bands] = loadDispersionModel(filename, mfr_expected [, bands_req])
%   Additionally returns the wavelengths or colour channels to be used for
%   evaluating the dispersion model.
%
% [dispersion_data, bands, transform_data] = loadDispersionModel(filename, mfr_expected [, bands_req])
%   Additionally returns the data needed to apply the dispersion model in
%   an image's coordinate space.
%
% ## Input Arguments
%
% filename -- Dispersion model filename
%   A character vector containing the filename and path of the '.mat' file
%   containing the dispersion model.
%
% mfr_expected -- Expected frame of reference
%   A Boolean indicating whether the dispersion model is expected to be in
%   the frame of reference of the reference wavelength/colour channel
%   (`true`), or in the frame of reference of the wavelength/colour channel
%   for which the dispersion is being calculated (`false`). This function
%   will check if the dispersion model is in the expected frame of
%   reference, and will throw an error if the expectation is violated.
%
% bands_req -- `bands` variable requirement
%   A Boolean scalar indicating whether (`true`) or not (`false`) there
%   must be a non-empty `bands` variable in the '.mat' file. Defaults to
%   `true` if not passed.
%
% ## Output Arguments
%
% dispersion_data -- Dispersion model
%   A model of dispersion, modelling the warping from the reference
%   wavelength band or colour channel to the other wavelength bands or
%   colour channels. `dispersion_data` can be converted to a function form
%   using `dispersionfun = makeDispersionfun(dispersion_data)`.
%
% bands -- Wavelength bands or colour channels
%   A vector containing the wavelengths or colour channel indices to use
%   when evaluating the functional form of the dispersion model. `bands` is
%   empty if it is not found in the '.mat' file, and if `bands_req` is
%   `false`. Otherwise, if a non-empty value of `bands` is not found in the
%   '.mat' file, an error is thrown.
%
% transform_data -- Spatial coordinate conversion information
%   If the dispersion model file contains a `model_space` variable,
%   `model_space` is packaged with some other fields to form a structure
%   with the following fields:
%   - 'model_space': A structure with same form as the `model_space` input
%     argument of 'modelSpaceTransform()'. It describes the coordinate
%     system, and domain, of the dispersion model.
%   - 'fill': The `fill` input argument of 'modelSpaceTransform()', which
%     defaults to `false` if it is not present in the dispersion model
%     file. `fill` indicates whether the dispersion model should be
%     resized, so that its domain fills the image.
%
%   `transform_data` is empty (`[]`) if the dispersion model file does not
%   contain a `model_space` variable.
%
% See also makeDispersionfun, modelSpaceTransform

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 26, 2018

nargoutchk(1, 3);
narginchk(2, 3);

bands_req = true;
if ~isempty(varargin)
    bands_req = varargin{1};
end

dispersion_data = [];
model_from_reference = [];
bands = [];
model_space = [];
fill = [];

model_variables = { 'dispersion_data', 'model_from_reference', 'bands' };
variables_transform = { 'model_space', 'fill' };
load(filename, model_variables{:}, variables_transform{:});
if isempty(dispersion_data) || isempty(model_from_reference)
    error('One or more of the dispersion model variables is not loaded.')
end
if bands_req && isempty(bands)
    error('The `bands` variable is empty or is not loaded.')
end
if mfr_expected ~= model_from_reference
    error('Dispersion model is in the wrong frame of reference.')
end

if ~isempty(model_space)
    transform_data.model_space = model_space;
    if ~isempty(fill)
        transform_data.fill = fill;
    else
        transform_data.fill = false;
    end
else
    transform_data = [];
end

end

