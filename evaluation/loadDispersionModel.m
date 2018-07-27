function [dispersion_data, bands, transform_data] = loadDispersionModel(filename, mfr_expected)
% LOADDISPERSIONMODEL Load dispersion information from a file
%
% ## Syntax
% dispersion_data = loadDispersionModel(dispersion_model_filename, mfr_expected)
% [dispersion_data, bands] = loadDispersionModel(dispersion_model_filename, mfr_expected)
% [dispersion_data, bands, transform_data] = loadDispersionModel(dispersion_model_filename, mfr_expected)
%
% ## Description
% dispersion_data = loadDispersionModel(dispersion_model_filename, mfr_expected)
%   Returns the dispersion model from the file
%
% [dispersion_data, bands] = loadDispersionModel(dispersion_model_filename, mfr_expected)
%   Additionally returns the wavelengths or colour channels to be used for
%   evaluating the dispersion model.
%
% [dispersion_data, bands, transform_data] = loadDispersionModel(dispersion_model_filename, mfr_expected)
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
%   when evaluating the functional form of the dispersion model.
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
narginchk(2, 2);

dispersion_data = [];
model_from_reference = [];
bands = [];
model_space = [];
fill = [];

model_variables_required = { 'dispersion_data', 'model_from_reference', 'bands' };
model_variables_transform = { 'model_space', 'fill' };
load(filename, model_variables_required{:}, model_variables_transform{:});
if any([isempty(dispersion_data), isempty(model_from_reference), isempty(bands)])
    error('One or more of the dispersion model variables is not loaded.')
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

