function color_mask = threshold_color_rgb(image, color, threshold)

[num_rows, num_cols, ~] = size(image);
color_mask = zeros(num_rows, num_cols);

image = rgb_normalise(image);

red = image(:,:,1);
green = image(:,:,2);
blue = image(:,:,3);

if strcmpi(color, 'red')
    self = red;
    other1 = green;
    other2 = blue;
elseif strcmpi(color, 'green')
    self = green;
    other1 = red;
    other2 = blue;
elseif strcmpi(color, 'blue')
    self = blue;
    other1 = red;
    other2 = green;
end 

for e = 1 : numel(self)
    if self(e) > other1(e) + threshold && self(e) > other2(e) + threshold
        color_mask(e) = 1;
    end
end
