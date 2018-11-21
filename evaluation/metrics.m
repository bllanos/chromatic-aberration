function [rmse, mrae, gof_or_psnr, global_rmse] = metrics(I, R, dim, peak, filter)
% METRICS  Compute image evaluation metrics
%
% ## Syntax
% rmse = metrics(I, R, dim, ~, filter)
% [rmse, mrae] = metrics(I, R, dim, ~, filter)
% [rmse, mrae, gof] = metrics(I, R, dim, 0, filter)
% [rmse, mrae, psnr] = metrics(I, R, dim, peak, filter)
% [____, global_rmse] = metrics(I, R, dim, [0 | peak], filter)
%
% ## Description
%
% rmse = metrics(I, R, dim, ~, filter)
%   Returns the per-pixel root mean squared error between the two images
%
% [rmse, mrae] = metrics(I, R, dim, ~, filter)
%   Additionally returns the mean relative absolute error between the two
%   images
%
% [rmse, mrae, gof] = metrics(I, R, dim, 0, filter)
%   Additionally returns either the per-pixel spectral goodness-of-fit,
%   between the two images.
%
% [rmse, mrae, psnr] = metrics(I, R, dim, peak, filter)
%   Additionally returns either the global peak signal-to-noise ratio
%   between the two images.
%
% [____, global_rmse] = metrics(I, R, dim, [0 | peak], filter)
%   Additionally returns the global root mean squared error.
%
% ## Input Arguments
%
% I -- Test image
%   An array containing test image data.
%
% R -- Reference image
%   An array of the same dimensions as `I` containing reference image data.
%
% dim -- Spectral dimension
%   The dimension corresponding to spectral information in the input
%   images. For example, `I` and `R` may be 2D arrays, where rows represent
%   pixels and columns represent spectral values at pixels. In this case,
%   `dim` should be `2`.
%
% peak -- Peak signal value
%   The peak value to use when computing peak signal-to-noise ratio.
%
% filter -- Filter out NaN and Inf values
%   Remove non-finite values from the per-pixel output arguments (using
%   MATLAB's 'isfinite()' function. Note that doing so destroys the spatial
%   structure of the output arguments, resulting in vectors being returned,
%   instead of arrays with dimensions based on the dimensions of `I`.
%   Per-pixel output arguments are `rmse`, `mrae`, and `gof`. `rmse` will
%   actually never contain non-finite values.
%
% ## Output Arguments
%
% rmse -- Root mean squared error
%   The root mean-squared error between `I` and `R`, evaluated per pixel
%   (in other words, computed along dimension `dim`).
%
% mrae -- Mean relative average error
%   The mean relative average error between `I` and `R`, evaluated per
%   pixel (in other words, computed along dimension `dim`).
%
% gof -- Spectral goodness-of-fit
%   The goodness-of-fit coefficients between `I` and `R`, evaluated per
%   pixel (in other words, computed along dimension `dim`). It does not
%   make sense to compute this value for monochromatic images, hence the
%   choice between computing `gof` and computing `psnr`.
%
% psnr -- Peak signal-to-noise ratio
%   The peak signal-to-noise ratio between `I` and `R`. `psnr` can be
%   computed when the size of `I` and `R` along dimension `dim` is `1`.
%   Otherwise an error will be thrown.
%
% global_rmse -- Image-wide root mean squared error
%   The root mean-squared error between `I` and `R`, evaluated by treating
%   all values in `I` and `R` equally, and pooling their squared errors
%   before taking the square root.
%
% ## References
% - The goodness-of-fit coefficient is used as a spectral error measure
%   (Equation 15) in:
%
%   Nguyen R.M.H., Prasad D.K., Brown M.S. (2014) Training-Based Spectral
%     Reconstruction from a Single RGB Image. In: Fleet D., Pajdla T.,
%     Schiele B., Tuytelaars T. (eds) Computer Vision – ECCV 2014. ECCV
%     2014. Lecture Notes in Computer Science, vol 8695. Springer, Cham
%
% - Root mean squared error and mean relative absolute error are defined
%   in:
%
%   B. Arad, O. Ben-Shahar and R. Timofte, "NTIRE 2018 challenge on
%     spectral reconstruction from RGB images," in The IEEE Conference on
%     Computer Vision and Pattern Recognition (CVPR) Workshops, 2018.
%
% See also immse, psnr, ssim, evaluateSpectral

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created November 20, 2018

narginchk(5, 5);
nargoutchk(1, 4);

diffIR = abs(I - R);
mse = dot(diffIR, diffIR, dim) ./ size(diffIR, dim);
rmse = sqrt(mse);
if nargout > 1
    mrae = mean(diffIR ./ R, dim);
    if filter
        mrae = mrae(isfinite(mrae));
    end
    if nargout > 2
        global_mse = mean(mse(:));

        if peak == 0
            gof_or_psnr = abs(dot(I, R, dim)) ./...
                sqrt(dot(I, I, dim) .* dot(R, R, dim));
            if filter
                gof_or_psnr = gof_or_psnr(isfinite(gof_or_psnr));
            end
        elseif size(I, dim) ~= 1
            error('Cannot compute PSNR for a multi-channel input.');
        else
            gof_or_psnr = 10 * log10((peak ^ 2) / global_mse);
        end

        if nargout > 3
            global_rmse = sqrt(global_mse);
        end
    end
end
end