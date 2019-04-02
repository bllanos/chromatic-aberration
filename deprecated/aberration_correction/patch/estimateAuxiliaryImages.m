% A helper function for 'solvePatchesAligned()'
%
% See also solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 6, 2018

function [I_3D_out, images] = estimateAuxiliaryImages(...
        I_3D, dispersion_matrix_patch, padding_filter, I_size,...
        sensitivity, align, n_images...
    )

n_bands = size(I_3D, 3);
padding_filter_spectral = repmat(padding_filter, 1, 1, n_bands);
I_3D_out = reshape(I_3D(padding_filter_spectral), [I_size, n_bands]);
if n_images > 0
    images = cell(1, n_images);
    n_channels_rgb = 3;

    I = I_3D(:);
    I_rgb = channelConversion(I, sensitivity, 1);
    padding_filter_rgb = repmat(padding_filter, 1, 1, n_channels_rgb);
    images{1} = reshape(I_rgb(padding_filter_rgb), [I_size, n_channels_rgb]);

    if n_images > 1
        if ~isempty(dispersion_matrix_patch)
            I_warped = dispersion_matrix_patch * I;
        else
            I_warped = I;
        end
        J_full = channelConversion(I_warped, sensitivity, 1);
        J_full = reshape(J_full(padding_filter_rgb), [I_size, n_channels_rgb]);
        images{2} = J_full;

        if n_images > 2
            images{3} = mosaic(J_full, align);
            
            if n_images > 3
                images{4} = reshape(I_warped(padding_filter_spectral), [I_size, n_bands]);
            end
        end
    end
else
    images = [];
end
end