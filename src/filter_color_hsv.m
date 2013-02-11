function color_mask = filter_color_hsv(image, color)

[num_rows, num_cols, ~] = size(image);
color_mask = zeros(num_rows, num_cols);

if strcmpi(color, 'red')
    hue_min = 0;
    hue_max = 19;
elseif strcmpi(color, 'green')
    hue_min = 80;
    hue_max = 150;
elseif strcmpi(color, 'blue')
    hue_min = 171;
    hue_max = 264;
end

red = double(image(:,:,1));
green = double(image(:,:,2));
blue = double(image(:,:,3));

for e = 1 : numel(color_mask)
    r = red(e);
    g = green(e);
    b = blue(e);
    hue = rgb2hue(r, g, b);
    if hue >= hue_min && hue <= hue_max
        color_mask(e) = 1;
    end
end


function hue = rgb2hue(r, g, b)
    hue_rad = atan2(sqrt(3) * (g - b), 2 * r - g - b);
    hue_deg = mod(hue_rad * 180 / pi, 360);
    hue = round(hue_deg);
