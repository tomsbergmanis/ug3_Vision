function idx = threshold_histogram(hist)
    highest_peak_value = -Inf;
    highest_peak_index = 0;
    lowest_peak_value = Inf;
    lowest_peak_index = 0;
    for c = 2 : length(hist) - 1
        cur = hist(c);
        prv = hist(c - 1);
        nxt = hist(c + 1);
        if prv < cur && nxt < cur
            if cur > highest_peak_value
                highest_peak_value = cur;
                highest_peak_index = c;
            elseif cur < lowest_peak_value
                lowest_peak_value = cur;
                lowest_peak_index = c;
            end
        end
    end
    min_idx = min([lowest_peak_index, highest_peak_index]);
    max_idx = max([lowest_peak_index, highest_peak_index]);
    deepest_valley_value = Inf;
    deepest_valley_index = 0;
    for c = min_idx + 1 : max_idx - 1
        cur = hist(c);
        prv = hist(c - 1);
        nxt = hist(c + 1);
        if prv > cur && nxt > cur && cur < deepest_valley_value
            deepest_valley_value = cur;
            deepest_valley_index = c;
        end
    end
    idx = deepest_valley_index;
end
