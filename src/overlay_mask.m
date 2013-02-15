% creates a grayscale version of |image| and darkens it by |gray_level|
% returns the result of overlaying |mask| to the grayscale image
function image = overlay_mask(image, mask, gray_level)
    if nargin < 3
        gray_level = 0.3;
    end

    grayscale = rgb2gray(image) * gray_level;
    image = cat(3, grayscale, grayscale, grayscale);
    for c = 1 : 3
        channel = image(:,:,c);
        channel(mask(:,:,c) == 1) = 255;
        image(:,:,c) = channel;
    end
end
