% returns image with a cross centered on pixel (|x|, |y|) drawn in |channel|
% [0, 0] is the top left corner of the image
function image = overlay_cross(image, channel, x, y)
    [h w ~] = size(image);
    
    if channel == 1
        other1 = 2;
        other2 = 3;
    end
    if channel == 2
        other1 = 1;
        other2 = 3;
    end
    if channel == 3
       other1 = 1;
       other2 = 2;
    end
    
    if y + 1 < h 
        image(x, y + 1, channel) = 1;
        image(x, y + 1, other1) = 0;
        image(x, y + 1, other2) = 0;
    end
    
    if y - 1 > 0 
        image(x, y - 1, channel) = 1;
        image(x, y - 1, other1) = 0;
        image(x, y - 1, other2) = 0;
    end
    
    if x + 1 < w
        image(x + 1, y, channel) = 1;
        image(x + 1, y, other1) = 0;
        image(x + 1, y, other2) = 0;
    end   
    
    if x - 1 > 0
        image(x - 1, y, channel) = 1;
        image(x - 1, y, other1) = 0;
        image(x - 1, y, other2) = 0;
    end    
    
    image(x, y, channel) = 1;
    image(x, y, other1) = 0;
    image(x, y, other2) = 0;
    
end

