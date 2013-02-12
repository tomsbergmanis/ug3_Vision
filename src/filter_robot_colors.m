function image_mask = filter_robot_colors(image, rgb_difference_threshold)
    if nargin < 2
        rgb_difference_threshold = 20;
    end

    [num_rows, num_cols, ~] = size(image);
    num_pixels = num_rows * num_cols;

    image = double(reshape(image, num_pixels, 3));
    image_mask = zeros(num_pixels, 3);

    for c = 1 : num_pixels
        r = image(c,1);
        g = image(c,2);
        b = image(c,3);
        h = rgb2hue(r, g, b);
        if is_red(h, r, g, b, rgb_difference_threshold)
            image_mask(c,1) = 1;
        elseif is_green(h, r, g, b, rgb_difference_threshold)
           image_mask(c,2) = 1;
        elseif is_blue(h, r, g, b, rgb_difference_threshold)
           image_mask(c,3) = 1;
        elseif is_cyan(h, r, g, b, rgb_difference_threshold)
            image_mask(c,3) = 1;
        end
    end

    image_mask = reshape(image_mask, num_rows, num_cols, 3);
end


function hue = rgb2hue(r, g, b)
    hue_rad = atan2(sqrt(3) * (g - b), 2 * r - g - b);
    hue_deg = mod(hue_rad * 180 / pi, 360);
    hue = round(hue_deg);
end


function rval = is_red(h, r, g, b, rgb_threshold)
    thresh = rgb_threshold * 2.0;
    hue_ok = h < 30 || h > 330;
    rgb_ok = abs(r - g) > thresh && abs(r - b) > thresh;
    rval = hue_ok && rgb_ok;
end


function rval = is_green(h, r, g, b, rgb_threshold)
    thresh = rgb_threshold * 0.8;
    hue_ok = h > 90 && h < 150;
    rgb_ok = abs(g - r) > thresh && abs(g - b) > thresh;
    rval = hue_ok && rgb_ok;
end


function rval = is_blue(h, r, g, b, rgb_threshold)
    thresh = rgb_threshold * 1.0;
    hue_ok = h > 210 && h < 270;
    rgb_ok = abs(b - r) > thresh && abs(b - g) > thresh;
    rval = hue_ok && rgb_ok;
end


function rval = is_cyan(h, r, g, b, rgb_threshold)
    thresh = rgb_threshold * 1.0;
    hue_ok = h > 150 && h < 210;
    rgb_ok = abs(g - r) > thresh && abs(b - r) > thresh;
    rval = hue_ok && rgb_ok;
end
