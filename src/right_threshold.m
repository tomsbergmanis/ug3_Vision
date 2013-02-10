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
