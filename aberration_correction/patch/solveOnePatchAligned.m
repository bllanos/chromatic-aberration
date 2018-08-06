% A helper function for 'solvePatchesAligned()'
%
% See also solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 6, 2018

function [I, dispersion_matrix_patch, varargout] = solveOnePatchAligned(...
        J, align, dispersionfun, sensitivity,...
        lambda, patch_size, padding, f, f_args, corner...
)
image_sampling = [size(J, 1), size(J, 2)];

% Find the boundaries of the patch
patch_lim = [
    corner(1) - padding, corner(2) - padding;
    corner(1) + patch_size(1) + padding - 1, corner(2) + patch_size(2) + padding - 1
    ];
trim = [padding + 1, padding + 1];
if patch_lim(1, 1) < 1
    trim(1, 1) = trim(1, 1) + patch_lim(1, 1) - 1;
    patch_lim(1, 1) = 1;
end
if patch_lim(1, 2) < 1
    trim(1, 2) = trim(1, 2) + patch_lim(1, 2) - 1;
    patch_lim(1, 2) = 1;
end
trim = [trim; trim + patch_size - 1];
if patch_lim(2, 1) > image_sampling(1)
    trim(2, 1) = trim(2, 1) + min(0, image_sampling(1) - patch_lim(2, 1) + padding);
    patch_lim(2, 1) = image_sampling(1);
end
if patch_lim(2, 2) > image_sampling(2)
    trim(2, 2) = trim(2, 2) + min(0, image_sampling(2) - patch_lim(2, 2) + padding);
    patch_lim(2, 2) = image_sampling(2);
end

% Construct arguments for the image estimation algorithm
image_sampling_f = diff(patch_lim, 1, 1) + 1;
align_f = offsetBayerPattern(patch_lim(1, :), align);
has_dispersion = ~isempty(dispersionfun);
if has_dispersion
    dispersion_f = dispersionfunToMatrix(...
        dispersionfun, lambda, image_sampling_f, image_sampling_f,...
        [0, 0, image_sampling_f(2), image_sampling_f(1)], true, flip(corner) - 1 ...
    );
else
    dispersion_f = [];
end
J_f = J(patch_lim(1, 1):patch_lim(2, 1), patch_lim(1, 2):patch_lim(2, 2), :);

% Solve for the output patch
varargout = cell(nargout - 2, 1);
[I_f, varargout{:}] = f(...
    image_sampling_f, align_f, dispersion_f, sensitivity, lambda,...
    J_f, f_args{:}...
);

% Remove padding
n_bands = length(lambda);
padding_filter_I = false([image_sampling_f, n_bands]);
padding_filter_I((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2)), :) = true;
padding_filter_I = reshape(padding_filter_I, [], 1, 1);
I = I_f(padding_filter_I);
I = reshape(I, [diff(trim, 1, 1) + 1, n_bands]);
if has_dispersion
    n_channels_J = size(J, 3);
    padding_filter_J = false([image_sampling_f, n_channels_J]);
    padding_filter_J((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2)), :) = true;
    padding_filter_J = reshape(padding_filter_J, [], 1, 1);
    dispersion_matrix_patch = dispersion_f(padding_filter_J, padding_filter_I);
else
    dispersion_matrix_patch = [];
end
end