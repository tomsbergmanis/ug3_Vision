function smooth_hist = smooth_histogram(hist, filterlen, sizeparam)

[len ~] = size(hist);

filter = exp(-0.5*(sizeparam/filterlen*[-(filterlen-1) : 2 : filterlen-1]').^2);
filter = filter/sum(filter);

convolution = conv(filter, hist);
smooth_hist = convolution(1 + filterlen/2 : len + filterlen/2);
