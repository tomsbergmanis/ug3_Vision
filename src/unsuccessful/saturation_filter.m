function image_mask = saturation_filter(image, min_saturation, precision)
    if nargin < 2
        min_saturation = 0.95;
    end
    if nargin < 3
        precision = 0.01;
    end

    [num_rows, num_cols, ~] = size(image);
    range = 0:precision:1;

    hsv = rgb2hsv(image);
    saturation = hsv(:,:,2);

    saturation_hist = histc(saturation(:), range);
    saturation_percentage = cumsum(saturation_hist) ./ sum(saturation_hist);
    idx = find(saturation_percentage >= min_saturation);
    saturation_threshold = range(idx(1) + 1);

    image_mask = zeros(num_rows, num_cols);
    image_mask(saturation >= saturation_threshold) = 1;
