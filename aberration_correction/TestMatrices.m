%% Test script for matrix generation functions

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created April 5, 2019

%% Create a series of small matrices to test boundary cases

closeness_threshold = 1e-8;

n_dim = 3;
combinations = unique(nchoosek(repelem(1:6, n_dim), n_dim), 'rows');
n_perm = factorial(n_dim);
sizes = zeros(size(combinations, 1) * n_perm, n_dim);
for i = 1:size(combinations, 1)
    sizes(((i - 1) * n_perm + 1):(i * n_perm), :) = perms(combinations(i, :));
end
sizes = unique(sizes, 'rows');

for i = 1:size(sizes, 1)
    sz = sizes(i, :);
    M = rand(sz);
    M_lin = reshape(M, [], 1);
    M_padded = padarray(M, [1, 1], 'replicate');
    
    % Spatial gradient
    [gx,gy] = imgradientxyz(M_padded, 'intermediate');
    gx = gx(2:(end-1), 2:(end-1), :);
    gy = gy(2:(end-1), 2:(end-1), :);
    gdiag_1 = (M_padded(1:(end - 2), 3:end, :) - M) ./ sqrt(2);
    gdiag_2 = (M_padded(3:end, 3:end, :) - M) ./ sqrt(2);
    
    [Gxy, Gdiag] = spatialGradient(sz);
    if ~all([gx(:); gy(:)] == Gxy * M_lin)
        fprintf('spatialGradient() xy-gradient failed for size [%d, %d, %d].\n', sz(1), sz(2), sz(3));
    end
    if ~all(abs([gdiag_1(:); gdiag_2(:)] - Gdiag * M_lin) < closeness_threshold)
        fprintf('spatialGradient() diagonal gradient failed for size [%d, %d, %d].\n', sz(1), sz(2), sz(3));
    end
    
    % Mosaicing matrix
    if sz(3) == 3
        for align_cell = {'rggb', 'bggr', 'grbg', 'gbrg'}
            align = align_cell{1};
            M_raw = mosaic(M, align);
            MOp = mosaicMatrix(sz(1:2), align);
            if ~all(M_raw(:) == MOp * M_lin)
                fprintf('mosaicMatrix() failed for size [%d, %d, %d].\n', sz(1), sz(2), sz(3));
            end
        end
    end
    
    % Channel conversion
    color_map = rand(sz(end) + 1, sz(end));
    M_color = channelConversion(M, color_map, 3);
    colorMatrix = channelConversionMatrix(sz(1:2), color_map);
    if ~all(abs(M_color(:) - colorMatrix * M_lin) < closeness_threshold)
        fprintf('channelConversionMatrix() failed for size [%d, %d, %d].\n', sz(1), sz(2), sz(3));
    end
    
    % Spatial Laplacian
    laplacian = 4 * M - (...
        M_padded(1:(end - 2), 2:(end - 1), :) +...
        M_padded(3:end, 2:(end - 1), :) +...
        M_padded(2:(end - 1), 1:(end - 2), :) +...
        M_padded(2:(end - 1), 3:end, :)...
    );
    L = spatialLaplacian(sz);
    if ~all(abs(laplacian(:) - L * M_lin) < closeness_threshold)
        fprintf('spatialLaplacian() failed for size [%d, %d, %d].\n', sz(1), sz(2), sz(3));
    end
    
    % Spectral gradient
    for full_GLambda = [true, false]
        glambda_lin = reshape(diff(M, 1, 3), [], 1);
        if full_GLambda
            glambda_lin = [glambda_lin; zeros(sz(1) * sz(2), 1)];
        elseif sz(3) == 1
            continue;
        end
        GLambda = spectralGradient(sz, full_GLambda);
        if ~all(glambda_lin == GLambda * M_lin)
            fprintf('spectralGradient() failed for size [%d, %d, %d].\n', sz(1), sz(2), sz(3));
        end
    end
end