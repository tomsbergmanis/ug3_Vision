% returns the result of putting |mask| onto |image|
% for each pixel that is set in some channel of |mask|, saturates the pixel in
% the equivalent channel of |image|
function image = overlay_mask(image, mask)
    [~, ~, num_channels] = size(image);
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
