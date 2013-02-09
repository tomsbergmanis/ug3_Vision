function norm_image = rgb_normalise(image)

R = double(image(:,:,1));
G = double(image(:,:,2));
B = double(image(:,:,3));

NormalizedRed = R(:,:)./sqrt(R(:,:).^2+G(:,:).^2+B(:,:).^2);
NormalizedGreen = G(:,:)./sqrt(R(:,:).^2+G(:,:).^2+B(:,:).^2);
NormalizedBlue = B(:,:)./sqrt(R(:,:).^2+G(:,:).^2+B(:,:).^2);

norm(:,:,1) = NormalizedRed(:,:);
norm(:,:,2) = NormalizedGreen(:,:);
norm(:,:,3) = NormalizedBlue(:,:);
norm_image = cat(3, NormalizedRed, NormalizedGreen, NormalizedBlue);
