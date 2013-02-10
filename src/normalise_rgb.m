function normalised_image = normalise_rgb(image)

red = double(image(:,:,1));
green = double(image(:,:,2));
blue = double(image(:,:,3));

euclid_rgb = sqrt(red(:,:).^2 + green(:,:).^2 + blue(:,:).^2);
normalised_red = uint8(round(red(:,:)./euclid_rgb.*255));
normalised_green = uint8(round(green(:,:)./euclid_rgb.*255));
normalised_blue = uint8(round(blue(:,:)./euclid_rgb.*255));

normalised_image = cat(3, normalised_red, normalised_green, normalised_blue);
