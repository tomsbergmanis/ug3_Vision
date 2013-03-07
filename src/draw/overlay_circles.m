% draws the circles defined by the centres in |centers| and radiuses in |radii|
% onto |image| and returns the modified image
% [0, 0] is the top left corner of the image
% if |radii| is a number, draw all circles with that radius
function image = overlay_circles(image, centers, radii, varargin)
    % parse options
    argc = size(varargin, 2);
    c = 1;
    while c <= argc
        arg = varargin{c};
        if strcmpi(arg, 'Color')
            if c + 1 > argc
                error('Color option should be followed by an integer tripplet');
            end
            colors = varargin{c + 1};
            c = c + 1;
        end
    end
    % option defaults
    if ~exist('colors', 'var')
        colors = [255 255 255];
    end

    num_circles = size(centers, 1);
    if size(radii, 1) == 1
        radii = repmat(radii, num_circles, 1);
    end
    if size(colors, 1) == 1
        colors = repmat(colors, num_circles, 1);
    end

    centers = round(centers);
    radii = round(radii);

    [~, ~, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        channel = overlay_circles_channel(channel, centers, radii, colors(:,c));
        image(:,:,c) = channel;
    end
end


function channel = overlay_circles_channel(channel, centers, radii, colors)
    [xmax, ymax] = size(channel);
    num_circles = size(centers, 1);
    for c = 1 : num_circles
        Xc = centers(c,1);
        Yc = centers(c,2);
        radius = radii(c);

        for theta = 0 : 0.1 : 359
            x = round(Xc + radius * cos(theta));
            y = round(Yc + radius * sin(theta));
            if x <= xmax && x > 0 && y <= ymax && y > 0
                channel(x, y) = colors(c,:);
            end
        end
    end
end
