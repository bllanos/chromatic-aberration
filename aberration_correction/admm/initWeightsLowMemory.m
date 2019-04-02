function out = initWeightsLowMemory(I_in, dispersion_matrix, len_I)
% INITWEIGHTSLOWMEMORY  Allocate memory for 'weightsLowMemory()'
%
% ## Syntax
% out = initWeightsLowMemory(I_in, dispersion_matrix, len_I)
%
% ## Description
% out = initWeightsLowMemory(I_in, dispersion_matrix, len_I)
%   Returns a structure containing arrays to be used by 'weightsLowMemory()'
%
% ## Input Arguments
%
% I_in -- True image structure
%   Refer to the documentation of weightsLowMemory.m. `I_in` can be empty
%   (`[]`).
%
% dispersion_matrix -- Model of dispersion
%   `dispersion_matrix` can be empty (`[]`), if there is no model of dispersion.
%   Otherwise, `dispersion_matrix` must be a matrix for warping `I`, the
%   estimated latent image, to the space of `I_in.I`, which is affected by
%   dispersion. `dispersion_matrix` must also convert `I` to the colour space of
%   `I_in.I`. If `I_in.I` is not affected by dispersion, then
%   `dispersion_matrix` should be empty.
%
% len_I -- Image size
%   The number of values in the image being estimated.
%
% ## Output Arguments
%
% out -- Preallocated arrays and intermediate data
%   The `in` input/output argument of 'weightsLowMemory()'. Refer to the
%   documentation of weightsLowMemory.m.
%
% See also weightsLowMemory

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created October 29, 2018

nargoutchk(1, 1);
narginchk(3, 3);

out = struct;
if ~isempty(I_in)
    if size(I_in.spectral_weights, 1) ~= size(I_in.I, 3)
        error('The number of rows of `I_in.spectral_weights` must equal the size of `I_in.I` in its third dimension.');
    end
    image_sampling = [size(I_in.I, 1), size(I_in.I, 2)];
    if isempty(dispersion_matrix)
        out.Omega_Phi = channelConversionMatrix(image_sampling, I_in.spectral_weights);
    else
        out.Omega_Phi = dispersion_matrix;
    end
    out.I_est = zeros(numel(I_in.I), 1);
end
out.I_init = zeros(len_I, 1);

end