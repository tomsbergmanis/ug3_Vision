function color_mask = filter_color_hist(image, color)
    [num_rows num_cols ~] = size(image);
    color_mask = zeros(num_rows, num_cols);

    if strcmpi(color, 'red')
        image_channel = image(:,:,1);
    elseif strcmpi(color, 'green')
        image_channel = image(:,:,2);
    elseif strcmpi(color, 'blue')
        image_channel = image(:,:,3);
    end

    color_mask(image_channel >= right_threshold(image_channel)) = 1;
end


function thresh = right_threshold(image_channel)
    hist = smooth_histogram(histc(image_channel(:), 0 : 255), 50, 4);
    [len ~] = size(hist);

    peak = find(hist == max(hist), 1, 'first');

    xmaxr = -1;
    peak_right = -1;
    for ind = peak + 1 : len - 1
        value = hist(ind);
        if hist(ind - 1) < value && value >= hist(ind + 1)
            if value >= xmaxr
                xmaxr = value;
                peak_right = ind;
            end
        end
    end
    if peak_right == -1
        peak_right = len;
        xmaxr = 1;
    end

    xminr = max(hist) + 1;
    valley_right = -1;
    for ind = peak + 1 : peak_right - 1
        value = hist(ind);
        if hist(ind - 1) > value && value <= hist(ind + 1)
            if value <= xminr
                xminr = value;
                valley_right = ind;
            end
        end
    end
    if valley_right == -1
        valley_right = len - 1;
        xminr = 2;
    end

    thresh = valley_right;
end

function smooth_hist = smooth_histogram(hist, filterlen, sizeparam)
    [len ~] = size(hist);

    filter = gaussian_window(filterlen, sizeparam);
    filter = filter/sum(filter);

    convolution = conv(filter, hist);
    smooth_hist = convolution(1 + filterlen/2 : len + filterlen/2);
end

function x = gaussian_window(n, w)
    x = exp(-0.5 * (w / n * [-(n - 1) : 2  : n - 1]') .^ 2);
end
