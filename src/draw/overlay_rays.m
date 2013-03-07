function image = overlay_rays(image, from, to, length, varargin)
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
        c = c + 1;
    end
    % option defaults
    if ~exist('colors', 'var')
        colors = [255 255 255];
    end

    num_rays = size(from, 1);
    if size(colors, 1) == 1;
        colors = repmat(colors, num_rays, 1);
    end

    for c = 1 : num_rays
        if from(c,:) == to(c,:)
            continue;
        end
        x0 = from(c,1);
        y0 = from(c,2);
        x1 = to(c,1);
        y1 = to(c,2);
        dx = x0 - x1;
        dy = y0 - y1;
        lambda = min(sqrt(length ^ 2 / (dx ^ 2 + dy ^ 2)), ...
                    -sqrt(length ^ 2 / (dx ^ 2 + dy ^ 2)));

        x = x0 + lambda * dx;
        y = y0 + lambda * dy;

        image = overlay_polygon(image, [x0 y0; x y], colors(c,:));
    end
end
