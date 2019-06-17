function J = shotNoise(I, n)
% SHOTNOISE  Add photon shot noise to a spectral image
%
% ## Syntax
% J = shotNoise(I)
% J_n = shotNoise(I, n)
%
% ## Description
% J = shotNoise(I)
%   Returns a version of the input image with added noise.
%
% J_n = shotNoise(I, n)
%   Returns 'n' versions of the input image with added noise.
%
% ## Input Arguments
%
% I -- Input image
%   A spectral image, in the form of a 3D array, where the first two dimensions
%   represent spatial dimensions, and the third dimension represents the
%   spectral dimension.
%
% n -- Number of output images
%   The number of noisy images to generate.
%
% ## Output Arguments
%
% J -- Noisy image
%   An image produced by adding noise to image `I`. This function simulates
%   photon shot noise that scales with the square root of the spectral intensity
%   [1]. The noise in different spectral bands is independent. `J` has the same
%   dimensions as `I`.
%
% J_n -- Noisy images
%   In the syntax where `n` is passed as an input argument, `J_n` is a 4D array,
%   with a size of `n` in its fourth dimension. `J_n(:, :, :, i)` is the i-th
%   noisy image corresponding to `I`.
%
% ## Notes
% - Use 'addNoise()' to add more realistic noise to a simulated raw image. This
%   function only simulates photon shot noise, whereas 'addNoise()' is intended
%   to simulate more realistic noise, arising not only from photon shot noise,
%   but also from other sources of noise in the image sensor.
% - Photon shot noise follows a Poisson distribution with a standard deviation
%   equal to the square root of the signal [1]. As the Poisson distribution is a
%   discrete distribution, this function first discretizes the input image.
%   Presently, the discretization has a resolution of 12 bits.
% - Negative values in the input image are clipped to zero.
%
% ## References
% 1 - Martinec, E. (2008). "Noise, dynamic range and bit depth in digital
%     SLR." Retrieved from
%     http://theory.uchicago.edu/∼ejm/pix/20d/tests/noise/
%
% See also addNoise

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created June 17, 2019

nargoutchk(1, 1);
narginchk(1, 2);

I(I < 0) = 0;
% Scale the image and convert it to integer values, because 
I_max = max(I);
resolution = 2^12; % 12-bit analog-to-digital conversion
I = floor(I .* (resolution ./ I_max));

if nargin < 2
    J = poissrnd(I);
else
    J = poissrnd(repmat(I, 1, 1, 1, n));
end

J = J .* (I_max / resolution);

end
