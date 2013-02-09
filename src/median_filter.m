function background = median_filter(folder)

reds = cell(100, 1);
greens = cell(100, 1);
blues = cell(100, 1);
for c = 1 : 100
    image = imread(sprintf('%s/%08d.jpg', folder, c));
    red = image(:,:,1);
    green = image(:,:,2);
    blue = image(:,:,3);
    reds{c} = red;
    greens{c} = green;
    blues{c} = blue;
end
dim = 2;
redMedian = median(cat(dim + 1, reds{:}), dim + 1);
greenMedian = median(cat(dim + 1, greens{:}), dim + 1);
blueMedian = median(cat(dim + 1, blues{:}), dim + 1);
background = cat(dim + 1, redMedian, greenMedian, blueMedian);
