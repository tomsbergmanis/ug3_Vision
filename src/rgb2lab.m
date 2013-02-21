function image = rgb2lab(image)
    image = applycform(image, makecform('srgb2lab'));
end
