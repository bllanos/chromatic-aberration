function [ I_3D, image_bounds, varargout ] = solvePatchesAlignedSerial(...
    J, align, dispersionfun, sensitivity,...
    lambda, options, f, f_args, varargin...
    )
% SOLVEPATCHESALIGNEDSERIAL  Run an image estimation algorithm on image patches
%
% ## Usage
%
% 'solvePatchesAlignedSerial()' is a serial version of
% 'solvePatchesAligned()', and so has even less memory footprint.
%
% The syntax of the two functions is the same; Refer to the documentation o
% of 'solvePatchesAligned()'.
%
% See also solvePatches, solvePatchesAligned

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 8, 2018

narginchk(8, 9);

verbose = true;
if verbose
    tic
end

single_patch = false;
if ~isempty(varargin)
    target_patch = varargin{1};
    single_patch = true;
end

has_dispersion = ~isempty(dispersionfun);
if has_dispersion && ~isa(dispersionfun, 'function_handle')
    error('`dispersionfun` must be a function handle.');
end

varargout = cell(1, nargout - 2);
n_auxiliary_images = min(nargout - 2, 4);
image_sampling = [size(J, 1), size(J, 2)];
image_bounds = [0, 0, image_sampling(2), image_sampling(1)];
int_method = options.int_method;

if single_patch
    [ patch_lim, trim ] = patchBoundaries(...
        image_sampling, options.patch_size, options.padding, target_patch...
    );

    % Construct arguments for the image estimation algorithm
    if isempty(align)
        align_f = [];
    else
        align_f = offsetBayerPattern(patch_lim(1, :), align);
    end
    image_sampling_f = diff(patch_lim, 1, 1) + 1;
    if has_dispersion
        dispersion_matrix_patch = dispersionfunToMatrix(...
            dispersionfun, lambda, image_sampling_f, image_sampling_f,...
            [0, 0, image_sampling_f(2), image_sampling_f(1)], true, flip(target_patch) - 1 ...
        );
    else
        dispersion_matrix_patch = [];
    end
    J_f = J(patch_lim(1, 1):patch_lim(2, 1), patch_lim(1, 2):patch_lim(2, 2), :);

    % Solve for the output patch
    [I_3D, varargout{(n_auxiliary_images + 1):end}] = f(...
        image_sampling_f, align_f, dispersion_matrix_patch, sensitivity, lambda,...
        J_f, f_args{:}...
    );

    padding_filter = false(image_sampling_f);
    padding_filter((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2))) = true;
    [I_3D, varargout(1:n_auxiliary_images)] = estimateAuxiliaryImages(...
            I_3D, dispersion_matrix_patch, padding_filter, diff(trim, 1, 1) + 1,...
            sensitivity, lambda, int_method,...
            align_f, n_auxiliary_images...
        );

else
    if verbose
        disp('Serial processing of patches...');
    end
    
    i_vector = 1:options.patch_size(1):image_sampling(1);
    n_i_vector = length(i_vector);
    j_vector = 1:options.patch_size(2):image_sampling(2);
    n_j_vector = length(j_vector);
    n_patches = n_i_vector * n_j_vector;
    patches_I = cell(n_i_vector, n_j_vector);
    patches_auxiliary = cell(n_i_vector, n_j_vector, n_auxiliary_images);
    
    for j = 1:n_j_vector
        for i = 1:n_i_vector
            corner = [i_vector(i), j_vector(j)];
            [ patch_lim, trim ] = patchBoundaries(...
                image_sampling, options.patch_size, options.padding, corner...
            );
            patch_J = J(patch_lim(1, 1):patch_lim(2, 1), patch_lim(1, 2):patch_lim(2, 2), :);
            
            if isempty(align)
                align_f = [];
            else
                align_f = offsetBayerPattern(patch_lim(1, :), align);
            end
            image_sampling_f = diff(patch_lim, 1, 1) + 1;
            if has_dispersion
                dispersion_matrix_patch = dispersionfunToMatrix(...
                    dispersionfun, lambda, image_sampling_f, image_sampling_f,...
                    [0, 0, image_sampling_f(2), image_sampling_f(1)], true,...
                    [corner(2), corner(1)] - 1 ...
                );
            else
                dispersion_matrix_patch = [];
            end

            % Solve for the output patch
            patches_I{i, j} = f(...
                image_sampling_f, align_f, dispersion_matrix_patch, sensitivity, lambda,...
                patch_J, f_args{:}...
            );

            padding_filter = false(image_sampling_f);
            padding_filter((trim(1, 1)):(trim(2, 1)), (trim(1, 2)):(trim(2, 2))) = true;
            patch_trimmed_size = diff(trim, 1, 1) + 1;
            [patches_I{i, j}, patches_auxiliary(i, j, :)] = estimateAuxiliaryImages(...
                    patches_I{i, j}, dispersion_matrix_patch, padding_filter,...
                    patch_trimmed_size,...
                    sensitivity, lambda, int_method,...
                    align_f, n_auxiliary_images...
            );
            if verbose
                fprintf('\tProcessed patch %d of %d\n', i + (j-1) * n_i_vector, n_patches);
            end
        end
    end
    
    if verbose
        fprintf('\tDone.\n');
        disp('Recombining results from patches...');
    end
    
    % Recombine patches
    I_3D = cell2mat(patches_I);
    for im = 1:n_auxiliary_images
        varargout{im} = cell2mat(patches_auxiliary(:, :, im));
    end
    
    if verbose
        fprintf('\tDone.\n');
        toc
    end
end
    
end