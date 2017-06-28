function [ disparity_spline, disparity_raw ] = radialChromaticAberration(...
    X_image, reference_wavelength_index, z, reference_z, varargin...
)
% RADIALCHROMATICABERRATION  Model chromatic aberration
%
% ## Syntax
% disparity_spline = radialChromaticAberration(...
%     X_image, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb]...
% )
% [ disparity_spline, disparity ] = radialChromaticAberration(...
%     X_image, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb]...
% )
%
% ## Description
% disparity_spline = radialChromaticAberration(...
%     X_image, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb]...
% )
%   Returns a spline model of chromatic aberration as a function of
%   distance from the image centre, and scene depth.
%
% [ disparity_spline, disparity ] = radialChromaticAberration(...
%     X_image, reference_wavelength_index,...
%     z, reference_z [, wavelengths, wavelengths_to_rgb]...
% )
%   Additionally returns the disparity values used to build the spline
%   model.
%
% ## Input Arguments
%
% X_image -- Image positions
%   Image locations produced by a set of scene features, for each
%   wavelength, and for each depth. `X_image_real(i, :, k, j)` is the
%   2D position on the image plane corresponding to the i-th scene feature,
%   emitting light the k-th wavelength, and positioned at the j-th depth.
%
% reference_wavelength_index -- Reference wavelength
%   The index into the third dimension of `X_image` representing the
%   reference wavelength. Chromatic aberration will be measured between
%   images produced with other wavelengths and the image produced with
%   this wavelength.
%
% z -- Scene depths
%   The z-positions of the scene elements producing the images in
%   `X_image`. `z` is a vector with a length equal to `size(X_image, 4)`.
%
%   `z` can be some transformation of the actual z-positions. More
%   generally, `z` is some variable describing positions, measured along
%   the optical axis, which is to be assessed as a predictor of chromatic
%   aberration.
%
% reference_z -- Reference depth
%   The z-position corresponding to a depth value of zero; An offset which
%   will be applied to the values in `z`.
%
% wavelengths -- Wavelengths corresponding to image measurements
%   The wavelengths of light corresponding to the elements of `X_image`. A
%   row vector of length `size(X_image, 3)`, where `wavelengths(k)` is the wavelength
%   used to generate the image locations in `X_image(:, :, k, :)`. This
%   parameter is used for figure legends only, not for calculations.
%
%   If both `wavelengths` and `wavelengths_to_rgb` are passed, graphical
%   output will be generated.
%
% wavelengths_to_rgb -- Colour map for wavelengths
%   RGB colours to be used when plotting points representing values for the
%   different wavelengths. The k-th row of this `size(X_image, 3)` x 3
%   matrix represents the RGB colour corresponding to the k-th wavelength,
%   `wavelengths(k)`.
%
%   If both `wavelengths` and `wavelengths_to_rgb` are passed, graphical
%   output will be generated.
%
% ## Output Arguments
%
% disparity_spline -- Thin-plate spline models of chromatic aberration
%   A set of thin-plate smoothing splines modeling chromatic aberration as
%   a function of `z`, and `r`. `r` is the 2D distance in the image plane,
%   of an image point from the origin, where the image point was formed by
%   light at the reference wavelength.
%
%   `disparity_spline` is a cell vector of length `size(X_image, 3)`, where
%   the k-th cell contains a thin-plate smoothing spline describing the
%   aberration between the k-th wavelength and the reference wavelength.
%   `disparity_spline(reference_wavelength_index)` is an empty cell.
%
%   `aberration = fnval(disparity_spline{k},[r; z])` evaluates the spline
%   model at the radial positions in the row vector `r`, and the associated
%   depths in the row vector `z`. `aberration` is a row vector of distances
%   measured in the radial direction.
%
% disparity_raw -- Input disparity values
%   The disparity values calculated from the locations in `X_image` and
%   `z`, used to construct `disparity_spline`. `disparity_raw` has the same
%   dimensions as `X_image`; `disparity_raw(i, :, k, j)` is the
%   displacement on the image plane of `X_image(i, :, k, j)`, relative to
%   `X_image(i, :, reference_wavelength_index, j)`. This disparity value is
%   measured for the i-th scene feature, emitting light at the k-th
%   wavelength, and positioned at the j-th depth.
%
% See also doubleSphericalLensPSF, tpaps

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 27, 2017

nargoutchk(1, 2);
narginchk(4, 6);

if ~isempty(varargin)
    if length(varargin) ~= 2
        error('Unexpected number of input arguments. Note that both `wavelengths` and `wavelengths_to_rgb` should be passed, or neither should be passed.');
    else
        wavelengths = varargin{1};
        wavelengths_to_rgb = varargin{2};
        verbose = true;
    end
else
    verbose = false;
end

n_wavelengths = size(X_image, 3);
X_image_reference = repmat(X_image(:, :, reference_wavelength_index, :), 1, 1, n_wavelengths, 1);
disparity_raw = X_image - X_image_reference;

z_adjusted = z - reference_z;
n_points = size(X_image, 1);
z_adjusted = repelem(z_adjusted, n_points);
if size(z_adjusted, 1) > size(z_adjusted, 2)
    z_adjusted = z_adjusted.';
end

if verbose
    X_image_matrix = reshape(permute(X_image, [1, 4, 2, 3]), [], 2, n_wavelengths);
    
    figure
    hold on
    legend_strings = cell(n_wavelengths, 1);
    for k = 1:n_wavelengths
        scatter3(X_image_matrix(:, 1, k), X_image_matrix(:, 2, k), z_adjusted, [], wavelengths_to_rgb(k, :), 'o');
        legend_strings{k} = sprintf('Image locations for \\lambda = %g nm', wavelengths(k));
    end
    legend(legend_strings);
    title('Input image locations')
    xlabel('X');
    ylabel('Y');
    zlabel('Depth')
    hold off
end

disparity_raw_radial = sqrt(dot(disparity_raw, disparity_raw, 2));
X_image_reference_radial = sqrt(dot(...
    X_image(:, :, reference_wavelength_index, :),...
    X_image(:, :, reference_wavelength_index, :),...
    2 ...
));
X_image_reference_radial = X_image_reference_radial(:).';

spline_predictors = [ X_image_reference_radial; z_adjusted ];
spline_responses = permute(squeeze(disparity_raw_radial), [1, 3, 2]);
spline_responses = reshape(spline_responses, 1, [], n_wavelengths);

disparity_spline = cell(n_wavelengths, 1);
for k = 1:n_wavelengths
    if k ~= reference_wavelength_index 
        disparity_spline{k} = tpaps(...
            spline_predictors,...
            spline_responses(:, :, k)...
        );
    end
end

end