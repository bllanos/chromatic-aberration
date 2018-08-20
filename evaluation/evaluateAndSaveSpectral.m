function varargout = evaluateAndSaveSpectral(varargin)
% EVALUATEANDSAVESPECTRAL  Compare spectral images and save the comparison results
%
% ## Syntax
% e_spectral_table = evaluateAndSaveSpectral(...
%     I_spectral, R_spectral, lambda,...
%     dp, I_name, alg_name, name_params [, fg_spectral]...
% )
% [e_spectral_table, fg_spectral] = evaluateAndSaveSpectral(...
%     I_spectral, R_spectral, lambda,...
%     dp, I_name, alg_name, name_params [, fg_spectral]...
% )
% evaluateAndSaveSpectral(directory, dp, I_name, all_alg_names, fg_spectral)
%
% ## Description
% e_spectral_table = evaluateAndSaveSpectral(...
%     I_spectral, R_spectral, lambda,...
%     dp, I_name, alg_name, name_params [, fg_spectral]...
% )
%   Returns a table containing quantitative comparisons between the two
%   images, and saves all graphical comparisons to files.
%
% [e_spectral_table, fg_spectral] = evaluateAndSaveSpectral(...
%     I_spectral, R_spectral, lambda,...
%     dp, I_name, alg_name, name_params [, fg_spectral]...
% )
%   Additionally returns some of the graphical comparisons instead of
%   saving them to files. The comparisons figures returned in `fg_spectral`
%   are those which can be augmented with results for additional test
%   images.
%
% evaluateAndSaveSpectral(directory, dp, I_name, all_alg_names, fg_spectral)
%   Saves the graphical comparisons that can contain results for multiple
%   test images to files, and closes all of the corresponding figures. The
%   advantage of using this syntax over the first syntax is that
%   `all_alg_names` allows this function to add legends to the saved
%   figures.
%
% ## Input Arguments
%
% I_spectral -- Estimated spectral image
%   An h x w x c array containing an estimated spectral image.
%
% R_spectral -- Reference spectral image
%   An h x w x c array containing the ideal/true spectral image.
%
% lambda -- Wavelengths
%   A vector of length 'c' containing the wavelengths corresponding to the
%   third dimension of `I_spectral` and `R_spectral`.
%
% dp -- Dataset description
%   A structure output by 'describeDataset()' providing information about
%   which comparisons to generate for the image having the name `I_name`.
%
% I_name -- Image name
%   A character vector containing the name of the image being estimated,
%   used to query `dp` for special evaluations to be conducted on it.
%   `I_name` is also used as the prefix of the filename for saving figures
%   (as MATLAB figure files) that pertain to the reference image, and are
%   not specific to the estimated image.
%
% alg_name -- Algorithm name
%   A character vector describing the image estimation algorithm and its
%   parameters, used to label the output in `e_spectral_table`.
%
% all_alg_names -- All algorithm names
%   A cell vector of character vectors to be used for legends in figures
%   which contain plotlines for multiple test images. `all_alg_names` does
%   not include an entry for the reference image.
%
% name_params -- Image partial filename
%   A character vector containing the base path and filename (excluding the
%   file extension) describing the estimated image. Graphical evaluations
%   performed on the estimated image will be saved to MATLAB figure files
%   that are given filepaths starting with this string and ending with a
%   prefix describing the type of evaluation, followed by the file
%   extension.
%
% fg_spectral -- Spectral error evaluation figures
%   A structure, which is similar to the `fg_spectral` output argument of
%   'evaluateSpectral()', but fields that refer to figures that can only
%   display results for single test image-reference image pairs are
%   ignored. If passed, some of the present image comparisons augment the
%   corresponding figures in this structure, rather than being shown in new
%   figures. The `fg_spectral` input and output arguments of this function
%   can be used to create comparison figures between more than one test
%   image and the reference image.
%
%   When passed as a third input argument, in the third syntax above, the
%   figures in `fg_spectral` which can contain results from multiple test
%   images are saved to files, and then are closed.
%
% directory -- Output directory
%   The directory in which to save figures pertaining to multiple test
%   images.
%
% ## Output Arguments
%
% e_spectral_table -- Spectral error statistics
%   A table form of the `e_spectral` structure returned by
%   'evaluateSpectral()' when invoked on the estimated and reference
%   images.
%
% fg_spectral -- Spectral error evaluation figures
%   A structure, which is similar to the `fg_spectral` output argument of
%   'evaluateSpectral()', but which only contains fields that refer to
%   figures that can display results for multiple test images. The
%   `fg_spectral` input and output arguments of this function can be used
%   to create comparison figures between more than two images, by passing
%   the output of one call to this function as the input to another call.
%
% ## Side Effects
% - The figures created by 'evaluateSpectral()' are closed, other than
%   those referred to by the `fg_spectral` output argument of this
%   function. (If `fg_spectral` is not requested as an output argument, all
%   figures are closed.)
%
% See also evaluateSpectral, evaluateAndSaveRGB, describeDataset,
% writetable

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 17, 2018

short_call = (nargin == 5);
if short_call
    nargoutchk(0, 0);
    directory = varargin{1};
    dp = varargin{2};
    I_name = varargin{3};
    all_alg_names = varargin{4};
    fg_spectral = varargin{5};
    save_multi_comparisons = true;
    I_filename = fullfile(directory, I_name);
else
    narginchk(7, 8);
    nargoutchk(0, 2);
    save_multi_comparisons = (nargout < 2);
    all_alg_names = [];
    
    I_spectral = varargin{1};
    R_spectral = varargin{2};
    lambda = varargin{3};
	dp = varargin{4};
    I_name = varargin{5};
    alg_name = varargin{6};
    name_params = varargin{7};
    if length(varargin) > 7
        fg_spectral = varargin{8};
    else
        fg_spectral = struct;
    end
end

options = dp.evaluation.global_spectral;
if isfield(dp.evaluation.custom_spectral, I_name)
    options = mergeStructs(...
        dp.evaluation.global_spectral,...
        dp.evaluation.custom_spectral.(I_name), false, true...
    );
end
    
if ~short_call
    filepath = fileparts(name_params);
    I_filename = fullfile(filepath, I_name);

    % Perform the image comparisons
    if isfield(fg_spectral, 'radiance')
        options.radiance_fg = fg_spectral.radiance;
    end
    if isfield(fg_spectral, 'scanlines')
        options.scanlines_fg = fg_spectral.scanlines;
    end
    [e_spectral, fg_spectral_current] = evaluateSpectral(...
        I_spectral, R_spectral, lambda, options...
    );

    % Save figures for this algorithm to files
    if isfield(fg_spectral_current, 'error_map')
        [~, max_ind] = max(e_spectral.mse.raw);
        savefig(...
            fg_spectral_current.error_map,...
            [name_params sprintf('_diffBand%d.fig', max_ind)], 'compact'...
        );
        close(fg_spectral_current.error_map);
    end
    if isfield(fg_spectral_current, 'bands_diff')
        savefig(...
            fg_spectral_current.bands_diff(2),...
            [...
                name_params,...
                sprintf(...
                    '_diffBand%dAnd%d.fig',...
                    options.bands_diff(1),...
                    options.bands_diff(2)...
                )...
            ], 'compact'...
        );
        close(fg_spectral_current.bands_diff(2));
    end
    
    % Save figures that do not pertain to any algorithms to files
    if isfield(fg_spectral_current, 'patches')
        savefig(...
            fg_spectral_current.patches,...
            [I_filename '_evalPatchLocations.fig'], 'compact'...
        );
        close(fg_spectral_current.patches);
    end
    if isfield(fg_spectral_current, 'scanlines_locations')
        savefig(...
            fg_spectral_current.scanlines_locations,...
            [I_filename '_evalLineLocations.fig'], 'compact'...
        );
        close(fg_spectral_current.scanlines_locations);
    end
    if isfield(fg_spectral_current, 'bands_diff')
        savefig(...
            fg_spectral_current.bands_diff(1),...
            [...
                I_filename,...
                sprintf(...
                    '_diffBand%dAnd%d.fig',...
                    options.bands_diff(1),...
                    options.bands_diff(2)...
                )...
            ], 'compact'...
        );
        close(fg_spectral_current.bands_diff(1));
    end
    
    % Pass on figures for all algorithms
    if ~save_multi_comparisons
        if isfield(fg_spectral_current, 'radiance')
            fg_spectral.radiance = fg_spectral_current.radiance;
        end
        if isfield(fg_spectral_current, 'scanlines')
            fg_spectral.scanlines = fg_spectral_current.scanlines;
        end
        varargout{2} = fg_spectral;
    end

    % Produce table output
    e_spectral_formatted = struct(...
        'MSE_max', e_spectral.mse.max,...
        'MSE_mean', e_spectral.mse.mean,...
        'MSE_median', e_spectral.mse.median,...
        'PSNR_min', e_spectral.psnr.min,...
        'PSNR_mean', e_spectral.psnr.mean,...
        'PSNR_median', e_spectral.psnr.median,...
        'SSIM_min', e_spectral.ssim.min,...
        'SSIM_mean', e_spectral.ssim.mean,...
        'SSIM_median', e_spectral.ssim.median,...
        sprintf('MI_bands%dAnd%d_Reference', options.mi_bands(1), options.mi_bands(2)),...
            e_spectral.mi_within(1),...
        sprintf('MI_bands%dAnd%d', options.mi_bands(1), options.mi_bands(2)),...
            e_spectral.mi_within(2),...
        'RMSE_max', e_spectral.rmse.max,...
        'RMSE_mean', e_spectral.rmse.mean,...
        'RMSE_median', e_spectral.rmse.median,...
        'GOF_min', e_spectral.gof.min,...
        'GOF_mean', e_spectral.gof.mean,...
        'GOF_median', e_spectral.gof.median...
    );
    n_bands = length(lambda);
    for c = 1:n_bands
        e_spectral_formatted.(sprintf('MI_band%d', c)) = e_spectral.mi_between(c);
    end
    if isfield(e_spectral, 'radiance')
        for i = 1:length(e_spectral.radiance)
            e_spectral_formatted.(sprintf('Patch%d_RMSE', i)) = e_spectral.radiance(i).rmse;
            e_spectral_formatted.(sprintf('Patch%d_GOF', i)) = e_spectral.radiance(i).gof;
        end
    end
    varargout{1} = struct2table(e_spectral_formatted, 'RowNames', {alg_name});
end

% Save figures for all algorithms to files
if save_multi_comparisons
    add_legend = ~isempty(all_alg_names);
    if add_legend
        all_alg_names = [{'Ground truth'}; reshape(all_alg_names, [], 1)];
    end
    if isfield(fg_spectral, 'radiance')
        for i = 1:length(fg_spectral.radiance)
            figure(fg_spectral.radiance(i));
            if add_legend
                legend(all_alg_names{:});
            end
            savefig(...
                fg_spectral.radiance(i),...
                [...
                    I_filename,...
                    sprintf(...
                        '_evalPatchX%dY%dW%dH%d.fig',...
                        options.radiance(i, 1),...
                        options.radiance(i, 2),...
                        options.radiance(i, 3),...
                        options.radiance(i, 4)...
                    )...
                ], 'compact'...
            );
            close(fg_spectral.radiance(i));
        end
    end
    if isfield(fg_spectral, 'scanlines')
        for i = 1:length(fg_spectral.scanlines)
            figure(fg_spectral.scanlines(i).rmse);
            if add_legend
                legend(all_alg_names{2:end});
            end
            savefig(...
                fg_spectral.scanlines(i).rmse,...
                [...
                    I_filename,...
                    sprintf(...
                        '_evalLineRMSE_X%dY%dX%dY%d.fig',...
                        options.scanlines(i, 1),...
                        options.scanlines(i, 2),...
                        options.scanlines(i, 3),...
                        options.scanlines(i, 4)...
                    )...
                ], 'compact'...
            );
            close(fg_spectral.scanlines(i).rmse);
            
            figure(fg_spectral.scanlines(i).gof);
            if add_legend
                legend(all_alg_names{2:end});
            end
            savefig(...
                fg_spectral.scanlines(i).gof,...
                [...
                    I_filename,...
                    sprintf(...
                        '_evalLineGOF_X%dY%dX%dY%d.fig',...
                        options.scanlines(i, 1),...
                        options.scanlines(i, 2),...
                        options.scanlines(i, 3),...
                        options.scanlines(i, 4)...
                    )...
                ], 'compact'...
            );
            close(fg_spectral.scanlines(i).gof);
        end
    end
end

end
