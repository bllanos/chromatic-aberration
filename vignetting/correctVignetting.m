function [ I ] = correctVignetting(I, vignettingfun)
% CORRECTVIGNETTING  Use a model of vignetting to correct an image for vignetting
%
% ## Syntax
% I = correctVignetting(I, vignettingfun)
%
% ## Description
% I = correctVignetting(I, vignettingfun)
%   Returns a version of the image corrected based on the given vignetting
%   model.
%
% ## Input Arguments
%
% I -- Image
%   A 2D or 3D array representing either an image.
%
% vignettingfun -- Vignetting model
%   A function which takes an input 2D array 'xy', where the two columns
%   represent image x and y-coordinates. The output of the function is a
%   column vector, 'factor', containing multiplicative factors representing
%   vignetting at the given image positions. An image can be corrected for
%   vignetting by dividing pixels by the corresponding elements of
%   'factor'.
%
% ## Output Arguments
%
% I -- Corrected image
%   A version of the input argument `I` created by dividing pixels by the
%   vignetting factors computed for their positions.
%
% See also vignettingPolyfit, makeVignettingfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 21, 2019

nargoutchk(1, 1);
narginchk(2, 2);

% Enumerate the positions of all pixels
image_sampling = [size(I, 1), size(I, 2)];
[X, Y] = meshgrid(1:image_sampling(2), 1:image_sampling(1));
xy = [X(:) Y(:)] - 0.5; % Place coordinates at pixel centres

% Calculate vignetting factors
factors = vignettingfun(xy);

% Output the corrected image
I = I .* repmat(1 ./ reshape(factors, image_sampling), 1, 1, size(I, 3));

end

