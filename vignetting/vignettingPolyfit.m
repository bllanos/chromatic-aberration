function [ polyfun, polyfun_data ] = vignettingPolyfit(...
    I, mask, max_degree_xy, varargin...
)
% VIGNETTINGPOLYFIT  Fit a polynomial model of image intensity
%
% ## Syntax
% [ polyfun, polyfun_data ] = vignettingPolyfit(...
%     I, mask, max_degree_xy [, align, verbose]...
% )
%
% ## Description
% [ polyfun, polyfun_data ] = vignettingPolyfit(...
%     I, mask, max_degree_xy [, align, verbose]...
% )
%   Returns a polynomial model of vignetting in terms of image position.
%   One to two output arguments can be requested; Optionally returns the
%   data used to construct the function handle form of the model of
%   vignetting.
%
% ## Input Arguments
%
% I -- Image
%   A 2D array representing either a RAW colour filter array image, or a
%   single-channel image.
%
% mask -- Constant region mask
%   A 2D logical array the same size as `I`, containing a mask showing the
%   parts of the image which, in the absence of vignetting, should have the
%   same intensity (one intensity for each colour channel, if `I`
%   represents a colour filter array image).
%
% max_degree_xy -- Maximum degree
%   The maximum degree of each of the image 'x' and 'y' coordinates in the
%   set of permissible polynomial models.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern of
%   the image. For example, 'gbrg'. The model of vignetting will be fit
%   using only pixels in the Green channel.
%
%   `align` has the same form as the `sensorAlignment` input argument of
%   `demosaic()`.
%
% verbose -- Debugging and visualization controls
%   If `verbose` is `true`, graphical output will be generated for
%   debugging purposes.
%
% ## Output Arguments
%
% polyfun -- Polynomial model
%   A function which takes an input 2D array 'xy', where the two columns
%   represent image x and y-coordinates. The output of the function is a
%   column vector, 'factor', containing multiplicative factors representing
%   vignetting at the given image positions. 'factor' is an evaluation of
%   the polynomial in 'x' and 'y' for the vignetting factor. The image `I`
%   can be corrected for vignetting by dividing pixels by the corresponding
%   elements of 'factor'. If the image is a single-channel image, the
%   corrected pixels will have values close to one. If the image is a
%   colour-filter array image, the corrected Green pixels will have values
%   close to one, and pixels in other colour channels should differ by a
%   per-channel approximately constant scaling factor from adjacent Green
%   pixels.
%
% polyfun_data -- Polynomial model data
%   A structure which can be used to obtain a function with the same
%   behaviour as `polyfun` by calling `polyfun =
%   makeVignettingfun(polyfun_data)`. `polyfun_data` has the following
%   fields:
%   - T_points: The normalization transformation applied to input spatial
%     coordinates prior to evaluating the model. A 3 x 3 matrix.
%   - T_factor_inv: The inverse normalization transformation applied to
%     the output of the model to obtain vignetting factors. A 2 x 2 matrix.
%   - powers: A 3D array where the first dimension is of size one, the
%     second dimension indexes different powers in the polynomial model,
%     and the third dimension indexes independent variables (x, y). The
%     array contains the exponents applied to each independent variable in
%     the polynomial model.
%   - n_powers: The number of terms in the polynomial model
%   - coeff: A vector of coefficients of the terms in the polynomial model
%   - degree_xy: The degree of the polynomial model in each of the spatial
%     independent variables (x, y)
%
% ## Algorithm
%
% This function constructs polynomials for the vignetting factor having
% degrees from zero to `max_degree_xy` in each of 'x' and 'y'. It selects
% the degree of the final polynomial model (a bivariate polynomial) by
% cross-validation.
%
% The input data (both the independent and dependent variables) are
% centered and scaled using 'normalizePointsPCA()' prior to polynomial
% fitting.
%
% ## References
% - MATLAB Help page for the polyfit() function
% - T. Hastie and R. Tibshirani. The Elements of Statistical Learning: Data
%   Mining, Inference, and Prediction, 2nd Edition. New York: Springer,
%   2009.
% - Rudakova, V. & Monasse, P. (2014). "Precise correction of lateral
%   chromatic aberration in images" (Guanajuato). 6th Pacific-Rim Symposium
%   on Image and Video Technology, PSIVT 2013. Springer Verlag.
%   doi:10.1007/978-3-642-53842-1_2
%   - This function is inspired by their use of a general bivariate
%     polynomial to model quantities which vary over the image plane.
% - Mannan, F. & Langer, M. S. (2016). "Blur calibration for depth from
%   defocus." In J. Guerrero (Ed.), 13th Conference on Computer and Robot
%   Vision, CRV 2016 (pp. 281-288). Institute of Electrical and Electronics
%   Engineers Inc. doi:10.1109/CRV.2016.62
%   - The general approach of using a parametric model to correct vignetting is
%     inspired by this article.
%
% See also makeVignettingfun, normalizePointsPCA, xylambdaPolyfit, crossval

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 21, 2019

    function [powers, n_powers] = powersarray(deg_xy)
        exponents_xy = 0:deg_xy;
        [x_powers_grid, y_powers_grid] = ndgrid(...
            exponents_xy, exponents_xy...
        );
        powers = cat(3,...
            reshape(x_powers_grid, 1, []),...
            reshape(y_powers_grid, 1, [])...
        );
        n_powers = size(powers, 2);
    end

    function y_est = predfun(vandermonde_matrix_train, y_train, vandermonde_matrix_test)
        coeff = vandermonde_matrix_train \ y_train;
        y_est = vandermonde_matrix_test * coeff;
    end

    function mse = crossvalfun(x_train, x_test)
        x_train_3d = repmat(permute(x_train(:, 1:2), [1 3 2]), 1, n_powers, 1);
        powers_3d = repmat(powers, size(x_train, 1), 1, 1);
        vandermonde_matrix_train = prod(x_train_3d .^ powers_3d, 3);
        
        x_test_3d = repmat(permute(x_test(:, 1:2), [1 3 2]), 1, n_powers, 1);
        powers_3d = repmat(powers, size(x_test, 1), 1, 1);
        vandermonde_matrix_test = prod(x_test_3d .^ powers_3d, 3);
        
        y_est_x = predfun(vandermonde_matrix_train, x_train(:, 3), vandermonde_matrix_test);
        mse = sum((x_test(:, 3) - y_est_x) .^ 2) / (2 * length(y_est_x));
    end

nargoutchk(1, 2);
narginchk(4, 5);

% Parse input arguments
is_raw = false;
verbose = false;
if ~isempty(varargin)
    is_raw = (isStringScalar(varargin{1}) || ischar(varargin{1}));
    if is_raw
        align = varargin{1};
        if length(varargin) > 1
            verbose = varargin{2};
        end
    else
        verbose = varargin{1};
    end
end

if size(I, 3) ~= 1
    error('The image `I` is not a RAW image or a single-channel image.');
end
image_sampling = size(I);
if is_raw
    % Only use the Green channel
    mask_bayer = bayerMask(image_sampling(1), image_sampling(2), align);
    mask = mask & mask_bayer(:, :, 2);
end

% Enumerate the positions of all pixels
px_ind = find(mask);
[px_row, px_col] = ind2sub(image_sampling, px_ind);
x = px_col - 0.5; % Place coordinates at pixel centres
y = px_row - 0.5;

dataset = [x y I(px_ind)];
n_points = size(dataset, 1);

polyfun_data = struct;

powers = [];
n_powers = [];
mean_mse = zeros(max_degree_xy + 1, 1);
std_mse = zeros(max_degree_xy + 1, 1);

% Normalize the data for numerical stability
[dataset_normalized_points, polyfun_data.T_points] = normalizePointsPCA([dataset(:, 1:2), ones(n_points, 1)]);
[dataset_normalized_scale, T_factor] = normalizePointsPCA([dataset(:, end), ones(n_points, 1)]);
polyfun_data.T_factor_inv = inv(T_factor);
dataset_normalized = [
    dataset_normalized_points(:, 1:(end-1)),...
    dataset_normalized_scale(:, 1:(end-1))
    ];

% Use cross validation to find the optimal polynomial complexity
deg_xy_vector = 0:max_degree_xy;
for deg_xy = deg_xy_vector
    [powers, n_powers] = powersarray(deg_xy);
    mse = crossval(@crossvalfun, dataset_normalized); % Defaults to 10-fold cross validation
    mean_mse(deg_xy + 1)  = mean(mse);
    std_mse(deg_xy + 1)  = std(mse);
end

% "Often a 'one-standard error' rule is used with cross-validation, in
% which we choose the most parsimonious model whose error is no more than
% one standard error above the error of the best model."
%
% From T. Hastie and R. Tibshirani, 2009.
[min_mse, ind] = min(mean_mse);
min_mse_std = std_mse(ind);
choice_mse = mean_mse - min_mse_std - min_mse;
possible_choices = (choice_mse <= 0);
polyfun_data.degree_xy = min(deg_xy_vector(possible_choices));

% Visualization
if verbose
    mse_minus1std = mean_mse - min_mse_std;
    mse_plus1std = mean_mse + min_mse_std;
    figure;
    hold on
    plot(deg_xy_vector, mse_minus1std, 'r-');
    plot(deg_xy_vector, mse_plus1std, 'r-');
    plot(deg_xy_vector(~possible_choices), mean_mse(~possible_choices), 'go');
    plot(deg_xy_vector(possible_choices), mean_mse(possible_choices), 'bx');
    scatter(polyfun_data.degree_xy, mean_mse(polyfun_data.degree_xy + 1), 'filled')
    hold off
    xlabel('Polynomial degree')
    ylabel('Mean square error')
    legend(...
        '-1 std', '+1 std',...
        'Mean square error >1 std. dev. from min',...
        'Mean square error <=1 std. dev. from min',...
        'Selected model'...
    );

    title('Cross validation error estimates');
end

% Fit the final model using all data
[polyfun_data.powers, polyfun_data.n_powers] = powersarray(polyfun_data.degree_xy);
x_train_final = repmat(permute(dataset_normalized(:, 1:2), [1 3 2]), 1, polyfun_data.n_powers, 1);
powers_final_rep = repmat(polyfun_data.powers, n_points, 1, 1);
vandermonde_matrix_final = prod(x_train_final .^ powers_final_rep, 3);
polyfun_data.coeff = vandermonde_matrix_final \ dataset_normalized(:, 3);

polyfun = makeVignettingfun(polyfun_data);
end