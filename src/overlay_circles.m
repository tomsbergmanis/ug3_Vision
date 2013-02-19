% draws the circles defined by the centres in |centers| and radiuses in |radii|
% onto |image| and returns the modified image
% [0, 0] is the top left corner of the image
% if |radii| is a number, draw all circles with that radius
function image = overlay_circles(image, centers, radii, color)
    if nargin < 4
        color = [255 255 255];
    end

    if length(radii) == 1
        radii = repmat(radii, length(centers), 1);
    end

    centers = round(centers);
    radii = round(radii);

    [~, ~, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        channel = overlay_circles_channel(channel, centers, radii, color(c));
        image(:,:,c) = channel;
    end
end


function channel = overlay_circles_channel(channel, centers, radii, color)
    [num_circles, ~] = size(centers);
    for c = 1 : num_circles
        Xc = centers(c,1);
        Yc = centers(c,2);
        radius = radii(c);

        for theta = 0 : 0.1 : 359
            x = round(Xc + radius * cos(theta));
            y = round(Yc + radius * sin(theta));
            channel(x, y) = color;
        end
    end
end
