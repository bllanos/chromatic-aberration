function [ disparity_spline, disparity_raw ] = radialChromaticAberration(...
    X_image, reference_wavelength_index, z, reference_z, varargin...
)
% RADIALCHROMATICABERRATION  Model chromatic aberration
%
% ## Syntax
% disparity_spline = radialChromaticAberration(...
%     X_image, reference_wavelength_index, z, reference_z [, verbose]...
% )
% [ disparity_spline, disparity ] = radialChromaticAberration(...
%     X_image, reference_wavelength_index, z, reference_z [, verbose]...
% )
%
% ## Description
% disparity_spline = radialChromaticAberration(...
%     X_image, reference_wavelength_index, z, reference_z [, verbose]...
% )
%   Returns a spline model of chromatic aberration as a function of
%   distance from the image centre, and scene depth.
%
% [ disparity_spline, disparity ] = radialChromaticAberration(...
%     X_image, reference_wavelength_index, z, reference_z [, verbose]...
% )
%   Additionally returns the evaluation of the spline model at the input
%   image points.
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
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
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
narginchk(4, 5);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

n_wavelengths = size(X_image, 3);
X_image_reference = repmat(X_image(:, :, reference_wavelength_index, :), 1, 1, n_wavelengths, 1);
disparity_raw = X_image - X_image_reference;
disparity_raw_radial = sqrt(dot(disparity_raw, disparity_raw, 2));
X_image_reference_radial = sqrt(dot(...
    X_image(:, :, reference_wavelength_index, :),...
    X_image(:, :, reference_wavelength_index, :),...
    2 ...
));
X_image_reference_radial = X_image_reference_radial(:).';
z_adjusted = z - reference_z;
z_adjusted = repelem(z_adjusted, size(X_image, 1));
if size(z_adjusted, 1) > size(z_adjusted, 2)
    z_adjusted = z_adjusted.';
end
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