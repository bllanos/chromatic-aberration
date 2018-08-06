% A helper function for 'solvePatchesAligned()'
%
% See also solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 6, 2018

function [I_3D_out, images] = estimateAuxiliaryImages(...
        I_3D, dispersion_matrix_patch,...
        padding_filter, I_size,...
        sensitivity, lambda, int_method,...
        align, target_patch, n_images...
    )

n_bands = length(lambda);
I_3D_out = reshape(I_3D(repmat(padding_filter, 1, 1, n_bands)), [I_size, n_bands]);
if n_images > 0
    images = cell(1, n_images);
    do_integration = ~(...
        isempty(int_method) || strcmp(int_method, 'none')...
    );
    n_channels_rgb = 3;

    I = I_3D(:);
    image_sampling_patch = [size(I_3D, 1), size(I_3D, 2)];
    if do_integration
        Omega = channelConversionMatrix(...
            image_sampling_patch, sensitivity, lambda, int_method...
            );
    else
        Omega = channelConversionMatrix(image_sampling_patch, sensitivity);
    end
    I_rgb = Omega * I;
    padding_filter_rgb = repmat(padding_filter, 1, 1, n_channels_rgb);
    images{1} = reshape(I_rgb(padding_filter_rgb), [I_size, n_channels_rgb]);

    if n_images > 1
        if ~isempty(dispersion_matrix_patch)
            J_full = Omega * dispersion_matrix_patch * I;
        else
            J_full = Omega * I;
        end
        J_full = reshape(J_full(padding_filter_rgb), [I_size, n_channels_rgb]);
        images{2} = J_full;

        if n_images > 2
            images{3} = mosaic(J_full, offsetBayerPattern(target_patch, align));
        end
    end
else
    images = [];
end
end