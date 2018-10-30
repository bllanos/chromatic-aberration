function out = initWeightsLowMemory(I_in)
% INITWEIGHTSLOWMEMORY  Allocate memory for 'weightsLowMemory()'
%
% ## Syntax
% out = initWeightsLowMemory(I_in)
%
% ## Description
% out = initWeightsLowMemory(I_in)
%   Returns a structure containing arrays to be used by 'weightsLowMemory()'
%
% ## Input Arguments
%
% I_in -- True image structure
%   Refer to the documentation of weightsLowMemory.m. `I_in` can be empty
%   (`[]`).
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
narginchk(1, 1);

out = struct;
if ~isempty(I_in)
    image_sampling = [size(I_in.I, 1), size(I_in.I, 2)];
    out.Omega_Phi = channelConversionMatrix(image_sampling, I_in.spectral_weights);
    out.I_est = zeros(numel(I_in.I), 1);
end

end