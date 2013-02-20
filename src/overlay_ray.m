function image = overlay_ray(image, from, to, length, varargin)
    % parse options
    argc = size(varargin, 2);
    c = 1;
    while c <= argc
        arg = varargin{c};
        if strcmpi(arg, 'DrawEndpoint')
            draw_endpoint = 1;
        elseif strcmpi(arg, 'Color')
            if c + 1 > argc
                error('Color option should be followed by an integer tripplet');
            end
            color = varargin{c + 1};
            c = c + 1;
        end
        c = c + 1;
    end
    % option defaults
    if ~exist('color', 'var')
        color = [255, 255, 255];
    end
    if ~exist('draw_endpoint', 'var')
        draw_endpoint = 0;
    end

    x0 = from(1);
    y0 = from(2);
    x1 = to(1);
    y1 = to(2);
    dx = x0 - x1;
    dy = y0 - y1;
    lambda = min(sqrt(length ^ 2 / (dx ^ 2 + dy ^ 2)), ...
                -sqrt(length ^ 2 / (dx ^ 2 + dy ^ 2)));

    x = x0 + lambda * dx;
    y = y0 + lambda * dy;
    image = overlay_polygon(image, [x0, y0; x, y], color);
    if draw_endpoint
        image = overlay_circles(image, [x, y], length * 0.1, color);
    end
end
