% returns the result of putting |mask| onto |image|
% for each pixel that is set in some channel of |mask|, saturates the pixel in
% the equivalent channel of |image|
function image = overlay_mask(image, mask, varargin)
    % parse options
    argc = size(varargin, 2);
    c = 1;
    while c <= argc
        arg = varargin{c};
        if strcmpi(arg, 'GrayScale')
            grayscale = 1;
        elseif strcmpi(arg, 'Saturation')
            if c + 1 > argc
                error('Saturation option should be followed by a double');
            end
            saturation = varargin{c + 1};
            c = c + 1;
        elseif strcmpi(arg, 'Lightness')
            if c + 1 > argc
                error('Lightness option should be followed by a double');
            end
            lightness = varargin{c + 1};
            c = c + 1;
        end
        c = c + 1;
    end

    % modify background image
    if exist('grayscale', 'var')
        gray_image = rgb2gray(image);
        image = cat(3, gray_image, gray_image, gray_image);
    end
    if exist('saturation', 'var')
        hsv_image = rgb2hsv(image);
        hsv_image(:,:,2) = hsv_image(:,:,2) * saturation;
        image = hsv2rgb(hsv_image);
    end
    if exist('lightness', 'var')
        hsv_image = rgb2hsv(image);
        hsv_image(:,:,3) = hsv_image(:,:,3) * lightness;
        image = hsv2rgb(hsv_image);
    end

    % lay mask onto background image
    num_channels = size(image, 3);
    channels = 1 : num_channels;
    for c = 1 : num_channels
        channel = image(:,:,c);
        mask_pixels = find(mask(:,:,c) == 1);
        channel(mask_pixels) = 255;
        image(:,:,c) = channel;
        for d = setdiff(channels, c)
            channel = image(:,:,d);
            channel(mask_pixels) = 0;
            image(:,:,d) = channel;
        end
    end
end
