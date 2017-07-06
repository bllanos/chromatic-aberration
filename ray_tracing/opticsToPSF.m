function [ psfFn ] = opticsToPSF( imageFn, U, U_prime, lens_radius, z_film )
% OPTICSTOPSF  Theoretical point spread functions from lens geometry
%
% ## Syntax
% psfFn = opticsToPSF( imageFn, U, U_prime, lens_radius, z_film )
%
% ## Description
% psfFn = opticsToPSF( imageFn, U, U_prime, lens_radius, z_film )
%   Returns a function for calculating point spread function statistics
%   corresponding to object locations.
%
% ## Input Arguments
%
% imageFn -- Imaging function
%   The `imageFn` output argument of 'opticsFromLens()'.
%
% U -- First principal plane
%   The z-coordinate of the first principal plane of the lens corresponding
%   to `imageFn`.
%
% U_prime -- Second principal plane
%   The z-coordinate of the second principal plane of the lens
%   corresponding to `imageFn`.
%
% lens_radius -- Lens radius
%   The radius of the lens (i.e. half the height of the lens when viewed
%   edge-on)
%
% z_film -- Film position
%   The z-coordinate of the image plane.
%
% ## Output Arguments
%
% psfFn -- Point spread function simulation function
%   A function of a single argument, `X_object`, which is an n x 3 array,
%   where the columns contain the x, y, and z-coordinates, respectively of
%   object points. The object points are assumed to be located at higher
%   z-coordinates than the front surface of the lens.
%
%   `psfFn` returns a n x 1 structure array. Each element is a structure of
%   the form output by 'analyzePSF()', describing the blur circle,
%   predicted by geometric optics, produced on the image plane by the given
%   object point.
%
% ## Notes
% - Image irradiances returned by `psfFn` are relative, not absolute
%   values.
% - Coordinate system:
%   - The front of the lens is spherical, with the centre of the sphere at
%     the origin.
%   - The radii of both faces of the lens are positive when the lens is
%     biconvex.
%   - The positive z-axis points towards the front of the lens, along
%     the optical axis, assuming the front of the lens is convex.
%
% ## References
% - Section 10.3 of B.K.P Horn, Robot Vision. Cambridge, Massachusetts: MIT
%   Press, 1986.
%
% See also opticsFromLens, analyzePSF

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created July 4, 2017

d_film = U_prime - z_film;
incident_normal = [ 0, 0, 1 ];
incident_position = [ 0, 0, U ];

    function stats = psfStatistics(X_object)
        n_points = size(X_object, 1);
        X_image = imageFn(X_object);
        
        % Adjust magnification to approximately account for the actual image plane
        % location
        U_prime_rep = repmat(U_prime, n_points, 1);
        magnification_correction = ...
            repmat(d_film, n_points, 1) ./ ...
            (U_prime_rep - X_image(:, 3));
        principal_point_prime = [zeros(n_points, 2), U_prime_rep];
        X_image_rays = X_image - principal_point_prime;
        mean_position = principal_point_prime + (X_image_rays .* magnification_correction);
        mean_position = mean_position(:, 1:2);
        
        % Calculate irradiance based on Section 10.3 of "Robot Vision"
        incident_position_rep = repmat(incident_position, n_points, 1);
        incident_direction = incident_position_rep - X_object;
        incident_distance_sq = dot(incident_direction, incident_direction, 2);
        incident_direction = incident_direction ./ repmat(...
            sqrt(incident_distance_sq), 1, 3 ...
        );
        incident_cosine = -dot(...
            repmat(incident_normal, n_points, 1), incident_direction, 2 ...
        );
        mean_value = incident_cosine .^ 4;
        
        % Size of blur circle
        z_film_rep = repmat(z_film, n_points, 1);
        radius = lens_radius .* abs(...
            (z_film_rep - X_image(:, 3)) ./ X_image(:, 3)...
        );
    
        mean_position = num2cell(mean_position, 2);
        mean_value = num2cell(mean_value, 2);
        stats = struct(...
            'mean_position', mean_position,...
            'mean_value', mean_value,...
            'max_position', mean_position,...
            'max_value', mean_value,...
            'radius', num2cell(radius, 2)...
        );
    end

psfFn = @psfStatistics;

end