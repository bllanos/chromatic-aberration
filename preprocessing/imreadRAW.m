function [ I_bayer, I_color ] = imreadRAW( filename, ops, varargin )
% IMREADRAW  Read a raw Bayer pattern image
%
% ## Syntax
% I_bayer = imreadRAW( filename, ops [, verbose] )
% I_color = imreadRAW( filename, ops, align [, wb, verbose] )
% [ I_bayer, I_color ] = imreadRAW( filename, ops, align [, wb, verbose] )
%
% ## Description
% I_bayer = imreadRAW( filename, ops [, verbose] )
%   Returns a preprocessed colour filter array image loaded from the file
% I_color = imreadRAW( filename, ops, align [, wb, verbose] )
%   Returns a full-colour image created from the file
% [ I_bayer, I_color ] = imreadRAW( filename, ops, align [, wb, verbose] )
%   Returns both colour filter array, and full-colour images
%
% ## Input Arguments
%
% filename -- Image filename
%   A character vector containing the full filename and path of the image
%   file to load. The image is expected to be a TIFF or DNG (Digital
%   Negative) file, as both can be loaded using MATLAB's `Tiff` class.
%
% ops -- Desired processing
%   A structure describing the processing to be performed on the image.
%   `ops` has the following Boolean fields:
%   - linearize: Linearize raw sensor counts using the
%     'SubIFDs.LinearizationTable' tag of the image file, if the tag exists
%   - demosaic: Perform demosaicing using MATLAB's `demosaic` function, to
%     generate `I_color`
%   - convertColor: Convert the image to the sRGB colour space, if the
%     image file contains the tag 'ColorMatrix2'
%   - wb: Perform white-balancing. White-balancing is performed before
%     demosaicing, and will affect both `I_bayer`, and `I_color`.
%
% align -- Bayer pattern description
%   A four-character character vector, specifying the Bayer tile pattern.
%   For example, 'gbrg'. `align` is not required if `ops.demosaic` is
%   `false`.
%
%   `align` has the same form as the `sensorAlignment` input argument of
%   `demosaic()`.
%
% wb -- White balancing multipliers
%   A 3-element vector, containing the multipliers to be applied to the
%   Red, Green, and Blue pixels, respectively. If the image file contains
%   the tag 'AsShotNeutral', the multipliers from the tag will take
%   precedence over those in `wb`. `wb` is not required if `ops.wb` is
%   `false`. White-balancing multipliers will be normalized by the
%   multiplier of the Green component (i.e. the second element of `wb`).
%   Note that white-balancing can result in clipping of the Red and Blue
%   channels, if these channels are given high multipliers relative to the
%   Green channel's multiplier.
%
% verbose -- Debugging flag
%   If true, graphical output will be generated for debugging purposes.
%
%   Defaults to false if not passed.
%
% ## Output Arguments
%
% I_bayer -- Colour filter array image
%   An image_height x image_width array, containing the colour filter array
%   data (Bayer pattern), subject to linearization and/or white-balancing,
%   as directed by `ops`.
%
% I_color -- Full-colour image
%   An image_height x image_width array x 3, containing the demosaicked
%   image, subject to the operations indicated in `ops`.
%
% ## Notes
% - The image is scaled to the range 0-1, following any linearization,
%   using on the black, and saturation values indicated in the image file's
%   tags. Presently, the code will ignore per-channel black, and saturation
%   values, and so will scale all channels using the same parameters.
%
% ## References
% - Processing RAW Images in MATLAB, by Rob Sumner (Department of
%   Electrical Engineering, UC Santa Cruz, 2014):
%   http://rcsumner.net/raw_guide/RAWguide.pdf
%
% See also demosaic

% Bernard Llanos
% Supervised by Dr. Y.H. Yang
% University of Alberta, Department of Computing Science
% File created August 9, 2017

if nargout > 1 && ~ops.demosaic
    error('Demosaicing was not requested, but `I_color` was requested as an output argument.')
end

if ops.convertColor && ~ops.demosaic
    error('Demosaicing was not requested, but is needed prior to colour space conversion, which was requested.')
end

if ops.wb
    narginchk(4, 5);
elseif ops.demosaic
    narginchk(3, 4);
else
    narginchk(2, 3);
end

wb = [];
align = [];
if ~isempty(varargin)
    if length(varargin) == 3
        align = varargin{1};
        wb = varargin{2};
        verbose = varargin{3};
    elseif length(varargin) == 2
        align = varargin{1};
        if ops.wb
            wb = varargin{2};
        else
            verbose = varargin{2};
        end
    elseif length(varargin) == 1
        if ops.demosaic
            align = varargin{1};
        else
            verbose = varargin{1};
        end
    end
else
    verbose = false;
end

%% Read image data

% warning off MATLAB:tifflib:TIFFReadDirectory:libraryWarning
if verbose
    [ ~, name, ext ] = fileparts(filename);
    display_name = [ name, ext ];
    display_name = replace(display_name, '_', '\_');
end
t = Tiff(filename,'r');
try
    offsets = getTag(t,'SubIFD'); % No such tag found on FLEA3 GigE images
    setSubDirectory(t,offsets(1));
catch ME
    if ~strcmp(ME.identifier,'MATLAB:imagesci:Tiff:tagRetrievalFailed')
        rethrow(ME)
    end
end
raw = read(t); % Create variable 'raw', the Bayer CFA data
close(t);
meta_info = imfinfo(filename);
% Crop to only valid pixels
if isfield(meta_info, 'SubIFDs')
    x_origin = meta_info.SubIFDs{1}.ActiveArea(2)+1; % +1 due to MATLAB indexing
    width = meta_info.SubIFDs{1}.DefaultCropSize(1);
    y_origin = meta_info.SubIFDs{1}.ActiveArea(1)+1;
    height = meta_info.SubIFDs{1}.DefaultCropSize(2);
    raw = double(raw(y_origin:y_origin+height-1,x_origin:x_origin+width-1));
end
raw = double(raw);

%% Linearization and intensity adjustment
% As I have disabled Gamma and LUT operations in FlyCapture2 when capturing
% images, this should not be necessary. The FLEA3 GigE Technical Reference,
% section 8.13, states that the sensor response is normally close to
% linear. In any case, `meta_info.SubIFDsf1g.LinearizationTable` does not
% exist for FLEA3 GigE images. This tag would be necessary for
% linearization. The tutorial says, "If this tag is empty (as for Canon
% cameras), you do not need to worry about this step."

linearized = false;
if isfield(meta_info, 'SubIFDs')
    if ops.linearize && isfield(meta_info.SubIFDs{1},'LinearizationTable')
        ltab = meta_info.SubIFDs{1}.LinearizationTable;
        raw = ltab(raw+1);
        linearized = true;
    end
    black = meta_info.SubIFDs{1}.BlackLevel(1);
    saturation = meta_info.SubIFDs{1}.WhiteLevel;
else
    black = meta_info.MinSampleValue;
    saturation = meta_info.MaxSampleValue;
end

lin_bayer = (raw-black)/(saturation-black);
lin_bayer = max(0,min(lin_bayer,1));

if verbose
    figure; imshow(lin_bayer);
    if linearized
        title(sprintf('Linearized Bayer image\n''%s''', display_name));
    else
        title(sprintf('Raw Bayer image (not linearized)\n''%s''', display_name));
    end
end

%% White Balancing

if ops.wb
    mask = bayerMask(size(lin_bayer,1),size(lin_bayer,2),align);
    if isfield(meta_info, 'AsShotNeutral')
        wb = (meta_info.AsShotNeutral).^-1;
    end
    wb = wb ./ wb(2);
    wbMask =...
        mask(:, :, 1) * wb(1) +...
        mask(:, :, 2) * wb(2) +...
        mask(:, :, 3) * wb(3);
    balanced_bayer = lin_bayer .* wbMask;
    
    if verbose
        figure; imshow(balanced_bayer);
        title(sprintf('White-balanced Bayer image\n''%s''', display_name));
    end
else
    balanced_bayer = lin_bayer;
end

%% Demosaicing
if ops.demosaic
    max_value = max(max(balanced_bayer));
    quantized = (balanced_bayer ./ max_value) * (2 ^ meta_info.BitDepth);
    switch meta_info.BitDepth
        case 8
            quantized = uint8(quantized);
        case 16
            quantized = uint16(quantized);
        case 32
            quantized = uint32(quantized);
        otherwise
            error('Unsupported bit depth.')
    end
    lin_rgb = double(demosaic(quantized, align)) * (max_value / (2 ^ meta_info.BitDepth));
    
    % Colour Space Conversion
    sRGB = false;
    if ops.convertColor && isfield(meta_info, 'ColorMatrix2')
        sRGB = true;
        
        % Not yet tested - This code may not work properly
        xyz2cam = meta_info.ColorMatrix2;
        xyz2cam = reshape(xyz2cam, 3, 3).'; % Convert from row-wise to column-wise order
        rgb2xyz = [
            0.4124564  0.3575761  0.1804375;
            0.2126729  0.7151522  0.0721750;
            0.0193339  0.1191920  0.9503041
            ];
        rgb2cam = xyz2cam * rgb2xyz; % Assuming previously defined matrices
        rgb2cam = rgb2cam ./ repmat(sum(rgb2cam,2),1,3); % Normalize rows to 1
        cam2rgb = inv(rgb2cam);
        
        r = cam2rgb(1,1)*lin_rgb(:,:,1)+cam2rgb(1,2)*lin_rgb(:,:,2)+cam2rgb(1,3)*lin_rgb(:,:,3);
        g = cam2rgb(2,1)*lin_rgb(:,:,1)+cam2rgb(2,2)*lin_rgb(:,:,2)+cam2rgb(2,3)*lin_rgb(:,:,3);
        b = cam2rgb(3,1)*lin_rgb(:,:,1)+cam2rgb(3,2)*lin_rgb(:,:,2)+cam2rgb(3,3)*lin_rgb(:,:,3);
        lin_srgb = cat(3,r,g,b);
        lin_srgb = max(0,min(lin_srgb,1)); % Always keep image clipped b/w 0-1
    else
        lin_srgb = lin_rgb;
    end
    
    % Display the result
    if verbose
        figure; imshow(lin_srgb);
        if sRGB
            title(sprintf('Colour image (after colour space conversion to sRGB)\n''%s''', display_name));
        else
            title(sprintf('Colour image (no colour space conversion to sRGB)\n''%s''', display_name));
        end
    end
end

%% Output
if ops.demosaic
    if nargout == 1
        I_bayer = lin_srgb;
    else
        I_bayer = balanced_bayer;
        I_color = lin_srgb;
    end
else
    I_bayer = balanced_bayer;
end