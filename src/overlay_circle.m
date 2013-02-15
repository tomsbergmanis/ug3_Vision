% returns image with a circle of |centre| = [x, y] and |radius|
% [0, 0] is the top left corner of the image
function image = overlay_circle(image, center, radius, thickness, color)
    if nargin < 4
        thickness = 1;
    end
    if nargin < 5
        color = [0 0 0];
    end

    for c = 0 : thickness
        image = overlay_circle_thin(image, center, radius - c, color);
    end
end


function image = overlay_circle_thin(image, center, radius, color)
    for c = 1 : length(color)
        channel = image(:,:,c);
        channel = overlay_circle_channel(channel, center, radius, color(c));
        image(:,:,c) = channel;
    end
end


function channel = overlay_circle_channel(channel, center, radius, color)
    Xc = center(1);
    Yc = center(2);

    for theta = 0 : 0.1 : 359
        x = round(Xc + radius * cos(theta));
        y = round(Yc + radius * sin(theta));
        channel(x, y) = color;
    end
end
