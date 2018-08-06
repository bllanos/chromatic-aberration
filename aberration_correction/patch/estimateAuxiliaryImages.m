% A helper function for 'solvePatchesAligned()'
%
% See also solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 6, 2018

function images = estimateAuxiliaryImages(...
        I_3D, dispersion_matrix_patch,...
        sensitivity, lambda, int_method,...
        align, target_patch, n_images...
    )
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
    images{1} = reshape(...
        I_rgb,...
        image_sampling_patch(1), image_sampling_patch(2), n_channels_rgb...
    );

    if n_images > 1
        if has_dispersion
            J_full = Omega * dispersion_matrix_patch * I;
        else
            J_full = Omega * I;
        end
        images{2} = reshape(...
            J_full,...
            image_sampling_patch(1), image_sampling_patch(2), n_channels_rgb...
        );

        if n_images > 2
            M = mosaicMatrix(...
                image_sampling_patch,...
                offsetBayerPattern(target_patch, align)...
                );
            J_est = M * J_full;
            images{3} = reshape(J_est, image_sampling_patch);
        end
    end
end
end