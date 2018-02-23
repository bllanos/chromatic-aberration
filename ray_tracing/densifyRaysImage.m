function [ I, mask ] = densifyRaysImage(...
    image_position, ray_irradiance, image_bounds, image_sampling, varargin...
)
% DENSIFYRAYSIMAGE  Generate an image from discrete samples of ray irradiance
%
% ## Syntax
% I = densifyRaysImage(...
%     image_position, ray_irradiance, image_bounds, image_sampling [, verbose]...
% )
% [ I, mask ] = densifyRaysImage(...
%     image_position, ray_irradiance, image_bounds, image_sampling [, verbose]...
% )
%
% ## Description
% I = densifyRaysImage(...
%     image_position, ray_irradiance, image_bounds, image_sampling [, verbose]...
% )
%   Returns the image approximated by summing ray irradiances within pixels
% [ I, mask ] = densifyRaysImage(...
%     image_position, ray_irradiance, image_bounds, image_sampling [, verbose]...
% )
%   Additionally returns a mask representing the nonzero portions of the
%   image
%
% ## Input Arguments
%
% image_position -- Image coordinates of rays
%   The points of intersection of the light paths traced through the lens
%   system with the image plane. `image_position` is a two-column array,
%   with the columns containing image x, and y coordinates, respectively.
%
% ray_irradiance -- Ray irradiance
%   The incident irradiance produced by individual rays, at the points of
%   intersection of the light paths with the image plane.
%
%   `ray_irradiance` is a vector, where the i-th element corresponds to the
%   i-th row of `image_position`.
%
%   Refer to the documentation of 'doubleSphericalLens.m' for more details.
%
% image_bounds -- Image domain
%   The rectangular domain of the image to be produced. `image_bounds` is a
%   vector containing the following elements:
%   1 - The x-coordinate of the bottom left corner of the image
%   2 - The y-coordinate of the bottom left corner of the image
%   3 - The width of the image (size in the x-dimension)
%   4 - The height of the image (size in the y-dimension)
%
% image_sampling -- Image resolution
%   The number of pixels at which to sample in the domain specified by
%   `image_bounds`. `image_sampling` is a two-element integer vector
%   containing the image height, and width, respectively, measured in
%   pixels.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% I -- Image
%   The irradiance pattern produced on the image plane by the rays traced
%   through the lens system. The image boundaries, and pixel resolution
%   are defined by the `image_bounds`, and `image_sampling` input
%   arguments, respectively.
%
%   Presently, pixel values are equal to the sums of the irradiances of the
%   rays intersecting their areas.
%
% mask -- Region of interest
%   A binary image describing which pixels in the image received any
%   incident rays.
%
% See also doubleSphericalLens, analyzePSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 18, 2017

nargoutchk(1, 2);
narginchk(4, 5);

if ~isempty(varargin)
    verbose = varargin{1};
else
    verbose = false;
end

% Plot the input data
if verbose
    figure
    hold on
    % Camera sensor
    surf(...
        [image_bounds(2), image_bounds(2) + image_bounds(3); image_bounds(2), image_bounds(2) + image_bounds(3)],...
        [image_bounds(1) + image_bounds(4), image_bounds(1) + image_bounds(4); image_bounds(1), image_bounds(1)],...
        [0, 0; 0, 0],...
        'EdgeColor', 'k', 'FaceAlpha', 0.4, 'FaceColor', 'g' ...
    );
    scatter(...
        image_position(:, 1),...
        image_position(:, 2),...
        [], ray_irradiance, 'filled'...
    )
    hold off
    xlabel('X');
    ylabel('Y');
    c = colorbar;
    c.Label.String = 'Exit ray irradiance';
    title('Intersection points with the image')
end

% Convert world coordinates to pixel indices
image_position_px = round([
    image_sampling(2) * ...
    (image_position(:, 1) - image_bounds(1)) / image_bounds(3),...
    image_sampling(1) * ...
    (image_bounds(2) + image_bounds(4) - image_position(:, 2)) / image_bounds(4)
    ]);
image_position_px = image_position_px(...
    image_position_px(:, 1) >= 1 & image_position_px(:, 1) <= image_sampling(2) &...
    image_position_px(:, 2) >= 1 & image_position_px(:, 2) <= image_sampling(1), :...
    );
image_position_px = sub2ind(image_sampling, image_position_px(:, 2), image_position_px(:, 1));

% Generate the image
I = zeros(image_sampling);
for i = 1:length(image_position_px)
    I(image_position_px(i)) = I(image_position_px(i)) + ray_irradiance(i);
end

if verbose
    figure
    imagesc(I);
    colorbar
    xlabel('X');
    ylabel('Y');
    zlabel('Irradiance');
    c = colorbar;
    c.Label.String = 'Irradiance';
    title('Estimated output image pixels')
end

if nargout > 1
    mask = (I > 0);
end

end

