function image_mask = filter_robot_colors(image)
    [num_rows, num_cols, ~] = size(image);
    num_pixels = num_rows * num_cols;
    image_mask = zeros(num_pixels, 3);

    rgb = double(reshape(image, num_pixels, 3));

    rgbN = double(reshape(normalise_rgb(image, 'approximate'), num_pixels, 3));
    rN_sdev = std(rgbN(:,1))
    gN_sdev = std(rgbN(:,2))
    bN_sdev = std(rgbN(:,3))
    rN_mean = mean(rgbN(:,1))
    gN_mean = mean(rgbN(:,2))
    bN_mean = mean(rgbN(:,3))

    hsv = reshape(rgb2hsv(image), num_pixels, 3);
    saturation_mean = mean(hsv(:,2));

    for c = 1 : num_pixels
        rN = rgbN(c,1);
        gN = rgbN(c,2);
        bN = rgbN(c,3);
        hue = hsv(c,1) * 360;
        % red
        if      (hue >= 330 || hue <= 30) &&...
                (normal_prob(rN, rN_mean, rN_sdev) < 0.0000001)
                    image_mask(c,1) = 1;
        % green
        elseif  (hue >= 90 && hue <= 150) &&...
                (normal_prob(gN, gN_mean, gN_sdev) < 0.001)
                   image_mask(c,2) = 1;
        % blue
        elseif  (hue >= 150 && hue <= 270) &&...
                (normal_prob(bN, bN_mean, bN_sdev) < 0.0000001)
                   image_mask(c,3) = 1;
        end
    end
   
    image_mask = reshape(image_mask, num_rows, num_cols, 3);
    %image_mask = clean(image_mask);
end


function x = normal_prob(val, mu, sigma)
    x = 1.0 / (sigma * sqrt(2 * pi)) * exp(-(val - mu) ^ 2 / (2 * sigma ^ 2));
end


function image = clean(image)
    [num_rows, num_cols, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        channel = bwmorph(channel, 'majority');
        image(:,:,c) = channel;
    end
end
