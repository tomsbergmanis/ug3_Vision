% normalises the values of the red, green, and blue channels of |image| in order
% to eliminate illumination differences in the image
% formula used: {r, g, b} = {r, g, b} / sqrt(r^2 + g^2 + b^2)
% if 'approximate' is passed as an additional parameter, instead use
% {r, g, b} = {r, g, b} / (r + g + b) 
% this is approximately two times faster than the exact normalisation
function normalised_image = normalise_rgb(image, varargin)
    approximate = ~isempty(find(strcmpi(varargin, 'approximate')));

    red = double(image(:,:,1));
    green = double(image(:,:,2));
    blue = double(image(:,:,3));

    if approximate
        euclid_rgb = red(:,:) + green(:,:) + blue(:,:);
    else
        euclid_rgb = sqrt(red(:,:).^2 + green(:,:).^2 + blue(:,:).^2);
    end
    norm_red = round(red(:,:) ./ euclid_rgb .* 255);
    norm_green = round(green(:,:) ./ euclid_rgb .* 255);
    norm_blue = round(blue(:,:) ./ euclid_rgb .* 255);

    normalised_image = cat(3, norm_red, norm_green, norm_blue);
end
