function w = gaussian(x)
% GAUSSIAN  Gaussian interpolation kernel
%
% ## Syntax
% w = gaussian(x)
%
% ## Description
% w = gaussian(x)
%   Returns the interpolation weight for a sample at location `x`.
%
% ## Input Arguments
%
% x -- Sample location
%   An array containing the locations of samples being interpolated. The
%   location of the interpolation point is the origin, `x = 0`.
%
% ## Output Arguments
%
% w -- Interpolation weights
%   An array with the same dimensions as `x` containing interpolation
%   weights for the samples at the locations in `x`. The weights are such
%   that, when convolved with samples spaced by increments of one unit,
%   this function performs Gaussian interpolation of the samples.
%
% ## Detailed Description
%
% The Gaussian function implemented in this function is intended to have an
% approximate bandlimit of 0.5 cycles/unit. If a Gaussian function has a
% standard deviation of 'sigma' in the spatial domain, its Fourier transform is
% a Gaussian function with a standard deviation of `1 / (2 * pi * sigma)` (See,
% for example, https://en.wikipedia.org/wiki/Gaussian_filter). Its power
% spectrum is a Gaussian function with a standard deviation of `1 / (2 * sqrt(2)
% * pi * sigma)`. If I define the bandlimit as the number 'p', where 95 percent
% of the area under the curve of the power spectrum is between '-p' and 'p',
% then `p = norminv(0.975)` times the standard deviation of the power spectrum.
% The standard deviation of the Gaussian function in the spatial domain must
% therefore be `sigma = norminv(0.975) / (sqrt(2) * pi)`, since `0.5 =
% norminv(0.975) / (2 * sqrt(2) * pi * sigma)`.
%
% In practice, I found that a standard deviation of one (corresponding to an
% approximate bandlimit of 0.2188 cycles/unit) gives slightly better spectral
% reconstruction results because of improved stability at the edges of the
% visible spectrum (i.e. by setting `findSamplingOptions.interpolant = @normpdf`
% instead of `findSamplingOptions.interpolant = @gaussian` in
% 'SetFixedParameters.m').
%
% I truncate the Gaussian to 6 standard deviations, to avoid having a kernel
% with infinite support, based on the discussion at
% https://en.wikipedia.org/wiki/Scale_space_implementation#The_sampled_Gaussian_kernel
%
% See also resamplingWeights

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created February 10, 2019

narginchk(1, 1);
nargoutchk(1, 1);

std_dev = norminv(0.975) / (sqrt(2) * pi);
w = normpdf(x, 0, std_dev);
% Truncate to 6 standard deviations
w(abs(x) > 6 * std_dev) = 0;
    
end
