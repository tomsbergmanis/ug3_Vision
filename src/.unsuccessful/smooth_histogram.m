function smooth_hist = smooth_histogram(hist, filterlen, sizeparam)
    if nargin < 2
        filterlen = 50;
    end
    if nargin < 3
        sizeparam = 4;
    end

    [len ~] = size(hist);

    filter = gaussian_window(filterlen, sizeparam);
    filter = filter/sum(filter);

    convolution = conv(filter, hist);
    smooth_hist = convolution(1 + filterlen/2 : len + filterlen/2);
end


function x = gaussian_window(n, w)
    x = exp(-0.5 * (w / n * [-(n - 1) : 2  : n - 1]') .^ 2);
end
