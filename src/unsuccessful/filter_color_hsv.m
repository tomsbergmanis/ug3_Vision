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

    hsv = rgb2hsv(image);
    hue = hsv(:,:,1);

    for e = 1 : numel(color_mask)
        h = hue(e) * 360;
        if h >= hue_min && h <= hue_max
            color_mask(e) = 1;
        end
    end
end
