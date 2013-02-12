function normalised_image = normalise_rgb(image)
    red = double(image(:,:,1));
    green = double(image(:,:,2));
    blue = double(image(:,:,3));

    euclid_rgb = sqrt(red(:,:).^2 + green(:,:).^2 + blue(:,:).^2);
    norm_red = uint8(round(red(:,:) ./ euclid_rgb .* 255));
    norm_green = uint8(round(green(:,:) ./ euclid_rgb .* 255));
    norm_blue = uint8(round(blue(:,:) ./ euclid_rgb .* 255));

    normalised_image = cat(3, norm_red, norm_green, norm_blue);
end
