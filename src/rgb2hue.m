function hue = rgb2hue(r, g, b)
    hue_rad = atan2(sqrt(3) * (g - b), 2 * r - g - b);
    hue_deg = mod(hue_rad * 180 / pi, 360);
    hue = round(hue_deg);
end
