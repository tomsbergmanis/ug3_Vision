function normalised_image = rgb_normalise(image)

red = double(image(:,:,1));
green = double(image(:,:,2));
blue = double(image(:,:,3));

euclid_rgb = sqrt(red(:,:).^2 + green(:,:).^2 + blue(:,:).^2);
normalised_red = red(:,:)./euclid_rgb;
normalised_green = green(:,:)./euclid_rgb;
normalised_blue = blue(:,:)./euclid_rgb;

normalised_image = cat(3, normalised_red, normalised_green, normalised_blue);
