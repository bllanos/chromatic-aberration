function [roi, T_roi] = modelSpaceTransform(image_size, model_space, varargin)
% MODELSPACETRANSFORM  Create a transformation from pixel to dispersion model coordinates
%
% ## Syntax
% roi = modelSpaceTransform(image_size, model_space [, fill, for_raw])
% [roi, T_roi] = modelSpaceTransform(image_size, model_space [, fill, for_raw])
%
% ## Description
% roi = modelSpaceTransform(image_size, model_space [, fill, for_raw])
%   Returns a region in the image corresponding to the valid range of the
%   dispersion model.
%
% [roi, T_roi] = modelSpaceTransform(image_size, model_space [, fill, for_raw])
%   Additionally returns a transformation matrix for converting image pixel
%   coordinates within the region to the coordinate space of the dispersion
%   model.
%
% ## Input Arguments
%
% image_size -- Image pixel dimensions
%   A two-element vector containing the image height and width in pixels.
%   The image is the "source image" - Coordinates in this image are to be
%   mapped to coordinates in the reference frame of the dispersion model.
%
% model_space -- Model attributes
%   A structure describing the domain of a dispersion model, having the
%   following fields:
%   - 'system': A character vector specifying whether the dispersion model
%     is constructed in an image coordinate system ('image'), or in a
%     geometrical optics coordinate system ('geometric'). The former has
%     its origin at the top left corner of the image, and has a y-axis
%     which points downwards. The latter has its origin at the image
%     centre, and has a y-axis which points upwards.
%   - 'corners': The first and second rows of this 2 x 2 matrix contain the
%     (x,y) coordinates of the top left and bottom right corners of the
%     region, respectively, in which the model is valid. The coordinates
%     are in the same frame of reference as the model.
%   - 'pixel_size': A field present if 'system' is 'geometric'. A
%     scalar giving the size of a pixel in the frame of reference of the
%     model.
%
% fill -- Flag enabling "fill" transformation construction
%   If `fill` is `true`, the output transformation, `T_roi`, will map pixel
%   coordinates in the source image to the coordinate frame of the model as
%   though the model's valid domain is exactly the source image. If `fill`
%   is `false`, no distortion of the source image to the model's valid
%   domain will occur.
%
%   Defaults to `false` if not passed.
%
% for_raw -- Flag for RAW images
%   If `for_raw` is `true`, `roi` will be created such that the colour
%   filter array pattern of the cropped image is the same as the colour
%   filter array pattern of the full image, and so such that the cropped
%   image is a valid colour filter array image (i.e. having even pixel
%   dimensions).
%
%   Defaults to `false` if not passed.
%
% ## Output Arguments
%
% roi -- Model domain in image subscripts
%   The valid domain of the model that lies within the source image, given
%   in the coordinate frame of the source image, except in pixel indices
%   (array subscripts) rather than pixel coordinates. `roi` is a 4-element
%   vector, with the following elements:
%   - The y-subscript of the top of the domain
%   - The y-subscript of the bottom of the domain
%   - The x-subscript of the left of the domain
%   - The x-subscript of the right of the domain
%
%   In other words, the model's domain corresponds to the region
%   `I(roi(1):roi(2), roi(3):roi(4))`, where `I` is the source image in
%   matrix form. The elements of `roi` are clipped to correspond to
%   positions within the source image.
%
% T_roi -- Coordinate transformation
%   A 3 x 3 transformation matrix for converting pixel coordinates to model
%   space coordinates. The homogenous vector `[x; y; 1]`, with 'x' and 'y'
%   in pixel coordinates, maps to the homogenous vector `T_roi * [x; y; 1]` in
%   the coordinate space of the dispersion model. 'x' and 'y' are relative
%   to the sub-image defined by `roi`. For instance, the point `T * [0.5,
%   0.5, 1]` is the transformation of the pixel at the top left of the
%   sub-image.
%
%   `T_roi` allows for evaluating the model of dispersion at points within
%   the model's valid domain. If `fill` is `true`, `T_roi` can be applied
%   to any point within the source image, as the source image is mapped
%   bijectively onto the model's valid domain.
%
% See also pixelsToWorldTransform, makeDispersionfun

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 1, 2018

nargoutchk(1, 2);
narginchk(2, 4);

fill = false;
for_raw = false;
if ~isempty(varargin)
    fill = varargin{1};
    if length(varargin) > 1
        for_raw = varargin{2};
    end
end

if strcmp(model_space.system, 'image')
    T = eye(3);
    roi = [...
        model_space.corners(1, 2), model_space.corners(2, 2),...
        model_space.corners(1, 1), model_space.corners(2, 1)
    ];
elseif strcmp(model_space.system, 'geometric')
    T = pixelsToWorldTransform(image_size, model_space.pixel_size);
    corners = [model_space.corners, ones(size(model_space.corners, 1), 1)];
    corners_image = (T \ corners.').';
    corners_image = corners_image(:, 1:2) ./ repmat(corners_image(:, 3), 1, 2);
    roi = [...
        corners_image(1, 2), corners_image(2, 2),...
        corners_image(1, 1), corners_image(2, 1)
    ];
else
    error('Unrecognized value of `model_space.system`.')
end

T_roi = T * [
    1, 0, roi(3);
    0, 1, roi(1);
    0, 0, 1
];

if fill
    T_roi = T_roi * [
            (roi(4) - roi(3)) / (image_size(2) - 1), 0, 0;
            0, (roi(2) - roi(1)) / (image_size(1) - 1), 0;
            0, 0, 1
        ] * [
            1, 0, -0.5;
            0, 1, -0.5;
            0, 0, 1
        ];
    roi = [1, image_size(1), 1, image_size(2)];
else
    if roi(3) < 1
        T_roi = T_roi * [
            1, 0, -roi(3);
            0, 1, 0;
            0, 0, 1
        ];
        roi(3) = 0.5;
    end
    if roi(1) < 1
        T_roi = T_roi * [
            1, 0, 0;
            0, 1, -roi(1);
            0, 0, 1
        ];
        roi(1) = 0.5;
    end
    if roi(2) > image_size(1)
        roi(2) = image_size(1) - 0.5;
    end
    if roi(4) > image_size(2)
        roi(4) = image_size(2) - 0.5;
    end
    
    roi = round(roi + 0.5); % Convert from pixel coordinates to pixel indices
    
    if for_raw
        % The following code assumes that the region of interest is at
        % least 4 pixels by 4 pixels.
        if mod(roi(3), 2) ~= 1
            T_roi = T_roi * [
                1, 0, 1;
                0, 1, 0;
                0, 0, 1
            ];
            roi(3) = roi(3) + 1;
        end
        if mod(roi(1), 2) ~= 1
            T_roi = T_roi * [
                1, 0, 0;
                0, 1, 1;
                0, 0, 1
            ];
            roi(1) = roi(1) + 1;
        end
        if mod(roi(2) - roi(1), 2) ~= 1
            roi(2) = roi(2) - 1;
        end
        if mod(roi(4) - roi(3), 2) ~= 1
            roi(4) = roi(4) - 1;
        end
    end
end

end

