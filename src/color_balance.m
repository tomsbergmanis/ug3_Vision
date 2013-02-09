% simple color balancing algorihtm as described in
% http://www.ipol.im/pub/art/2011/llmps-scb/
function balanced_image = color_balance(image, s1, s2)

[num_rows, num_cols, ~] = size(image);
N = num_rows * num_cols;

red = double(image(:,:,1));
green = double(image(:,:,2));
blue = double(image(:,:,3));

% get cumulative histogram of color intensities
red_hist = cumsum(histc(red(:), 0:255));
green_hist = cumsum(histc(green(:), 0:255));
blue_hist = cumsum(histc(blue(:), 0:255));

% find vmin, vmax
red_vmin = find(red_hist <= N * s1 / 100.0, 1, 'first');
red_vmax = find(red_hist > N * (1 - s2 / 100.0), 1, 'last');
red_dv = red_vmax - red_vmin;
green_vmin = find(green_hist <= N * s1 / 100.0, 1, 'first');
green_vmax = find(green_hist > N * (1 - s2 / 100.0), 1, 'last');
green_dv = green_vmax - green_vmin;
blue_vmin = find(blue_hist <= N * s1 / 100.0, 1, 'first');
blue_vmax = find(blue_hist > N * (1 - s2 / 100.0), 1, 'last');
blue_dv = blue_vmax - blue_vmin;

% saturate pixels
red(red < red_vmin) = red_vmin;
green(green < green_vmin) = green_vmin;
blue(blue < blue_vmin) = blue_vmin;
red(red > red_vmax) = red_vmax;
green(green > green_vmax) = green_vmax;
blue(blue > blue_vmax) = blue_vmax;

% rescale pixels
red = arrayfun(@(x) (x - red_vmin) * 255 / red_dv, red);
green = arrayfun(@(x) (x - green_vmin) * 255 / green_dv, green);
blue = arrayfun(@(x) (x - blue_vmin) * 255 / blue_dv, blue);

% recombine image
balanced_image = cat(3, uint8(red), uint8(green), uint8(blue));
