% returns image with a lines drawn from the nth point in |points| to the n+1th
% [0, 0] is the top left corner of the image
function image = overlay_polygon(image, points, color)
    if nargin < 3
        color = [255 255 255];
    end
    
    points = round(points);

    [~, ~, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        channel = overlay_polygon_channel(channel, points, color(c));
        image(:,:,c) = channel;
    end
end


% Bresenham's line algorithm (simplified version)
% http://en.wikipedia.org/wiki/Bresenham's_line_algorithm#Simplification
function channel = overlay_polygon_channel(channel, points, color)
    for c = 1 : length(points) - 1
        start = points(c,:);
        stop = points(c + 1,:);
        x0 = start(1);
        y0 = start(2);
        x1 = stop(1);
        y1 = stop(2);

        dx = abs(x1 - x0);
        dy = abs(y1 - y0);

        if x0 < x1
            sx = 1;
        else
            sx = -1;
        end
        if y0 < y1
            sy = 1;
        else
            sy = -1;
        end

        err = dx - dy;

        while 1
            channel(x0, y0) = color;
            if x0 == x1 && y0 == y1
                break;
            end
            e2 = 2 * err;
            if e2 > -dy
                err = err - dy;
                x0 = x0 + sx;
            end
            if e2 < dx
                err = err + dx;
                y0 = y0 + sy;
            end
        end
    end
end
