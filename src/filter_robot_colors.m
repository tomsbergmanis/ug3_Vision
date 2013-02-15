function image_mask = filter_robot_colors(image)
    [num_rows, num_cols, ~] = size(image);
    num_pixels = num_rows * num_cols;
    image_mask = zeros(num_pixels, 3);

    rgb = double(reshape(image, num_pixels, 3));

    rgbN = double(reshape(normalise_rgb(image, 'approximate'), num_pixels, 3));
    rN_sdev = std(rgbN(:,1));
    gN_sdev = std(rgbN(:,2));
    bN_sdev = std(rgbN(:,3));
    rN_mean = mean(rgbN(:,1));
    gN_mean = mean(rgbN(:,2));
    bN_mean = mean(rgbN(:,3));

    hsv = reshape(rgb2hsv(image), num_pixels, 3);

    for c = 1 : num_pixels
        rN = rgbN(c,1);
        gN = rgbN(c,2);
        bN = rgbN(c,3);
        hue = hsv(c,1) * 360;
        % current pixel is red
        if      (hue >= 330 || hue <= 30) &&...
                (normal_prob(rN, rN_mean, rN_sdev) < 0.000001)
                    image_mask(c,1) = 1;
        % current pixel is green
        elseif  (hue >= 80 && hue < 150) &&...
                (normal_prob(gN, gN_mean, gN_sdev) < 0.0001)
                   image_mask(c,2) = 1;
        % current pixel is blue
        elseif  (hue >= 150 && hue <= 270) &&...
                (normal_prob(bN, bN_mean, bN_sdev) < 0.0000001)
                   image_mask(c,3) = 1;
        end
    end
   
    image_mask = reshape(image_mask, num_rows, num_cols, 3);

    image_mask = remove_noise(image_mask);
    image_mask = enforce_similar_channel_areas(image_mask);
end


function x = normal_prob(val, mu, sigma)
    x = 1.0 / (sigma * sqrt(2 * pi)) * exp(-(val - mu) ^ 2 / (2 * sigma ^ 2));
end


function image = remove_noise(image)
    [~, ~, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        channel = bwmorph(channel, 'majority', Inf);
        channel = bwmorph(channel, 'bridge', Inf);
        image(:,:,c) = channel;
    end
end


% we know that the robots are all about the same size - we can thus remove any
% channels in the mask that have a much smaller area than the other channels
% this catches some problems like the shadow of the blue robot in 
% data/1/00000095.jpg being detected as a red blob
function image = enforce_similar_channel_areas(image, area_proportion_threshold)
    if nargin < 2
        area_proportion_threshold = 0.125;
    end

    [~, ~, num_channels] = size(image);
    areas = zeros(num_channels, 1);
    for c = 1 : num_channels
        channel = image(:,:,c);
        region_props = regionprops(channel, 'Area');
        if length(region_props) > 0
            areas(c) = region_props.Area;
        end
    end
    avg_area = mean(areas(areas > 0));
    for c = 1 : num_channels
        if areas(c) < avg_area * area_proportion_threshold
            image(:,:,c) = image(:,:,c) * 0;
        end
    end
end
