% returns image with a line drawn from |start| = [x0, y0] to |stop| = [x1, y1]
% [0, 0] is the top left corner of the image
function image = overlay_line(image, start, stop, thickness, color)
    if nargin < 4
        thickness = 1;
    end
    if nargin < 5
        color = [0 0 0];
    end
    
    for c = 0 : thickness
        Xi = start(1) + ceil(c / 2) * (-1) ^ (c + 1);
        Xf = stop(1) + ceil(c / 2) * (-1) ^ (c + 1);
        image = overlay_line_thin(image, [Xi start(2)], [Xf stop(2)], color);
    end
end


function image = overlay_line_thin(image, start, stop, color)
    for c = 1 : length(color)
        channel = image(:,:,c);
        channel = overlay_line_channel(channel, start, stop, color(c));
        image(:,:,c) = channel;
    end
end


% Bresenham's line algorithm (simplified version)
% http://en.wikipedia.org/wiki/Bresenham's_line_algorithm#Simplification
function channel = overlay_line_channel(channel, start, stop, color)
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
