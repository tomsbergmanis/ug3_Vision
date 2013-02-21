% extracts all the relevant information from an image
% 1) a mask where pixels in channel 1 are set to 1 iif they belong to the red
% robot (and similarly for channel 2/green robot and channel 3/blue robot)
% 2) an array with the centroids of the pixels in the mask
% 3) an array with the centroids of the convex hulls of the pixels in the mask
function [color_mask, varargout] = analyse_image(image)
    [num_rows, num_cols, num_channels] = size(image);

    color_mask_ugly = get_color_mask(image);
    color_mask = zeros(num_rows, num_cols, num_channels);

    centroids = zeros(num_channels, 2);
    convex_centroids = zeros(num_channels, 2);
    for c = 1 : num_channels
        channel = color_mask_ugly(:,:,c);
        if ~any(channel(:))
            continue;
        end
        props = regionprops(channel, 'Centroid', 'ConvexImage', 'BoundingBox');
        convex_props = regionprops(props.ConvexImage, 'Centroid');
        convex_centroid = convex_props.Centroid;
        convex_centroid = [convex_centroid(2) + props.BoundingBox(2), ...
                           convex_centroid(1) + props.BoundingBox(1)];
        convex_centroids(c,:) = convex_centroid;
        centroids(c,:) = [props.Centroid(2), props.Centroid(1)];
        convex_image = props.ConvexImage;
        [num_rows_convex, num_cols_convex] = size(convex_image);
        for row = 1 : num_rows_convex
            for col = 1 : num_cols_convex
                if convex_image(row, col) == 1
                    newrow = round(row + props.BoundingBox(2));
                    newcol = round(col + props.BoundingBox(1));
                    color_mask(newrow, newcol, c) = 1;
                end
            end
        end
        color_mask(:,:,c) = bwmorph(color_mask(:,:,c), 'remove');
    end
    color_mask = overlay_rays(color_mask, centroids, convex_centroids, 100, ...
                              'Color', [1 0 0; 0 1 0; 0 0 1]);
    varargout{1} = centroids;
    varargout{2} = convex_centroids;
    varargout{3} = color_mask_ugly;
end


function image_mask = get_color_mask(image)
    [num_rows, num_cols, ~] = size(image);
    num_pixels = num_rows * num_cols;
    image_mask = zeros(num_pixels, 3);

    rgb = double(reshape(image, num_pixels, 3));

    rgbN = double(reshape(normalise_rgb(image, 'approximate'), num_pixels, 3));
    rN_sdev = std(rgbN(:,1));
    gN_sdev = std(rgbN(:,2));
    bN_sdev = std(rgbN(:,3));
    rN_mean = mean(rgbN(:,1));
    gN_mean = mean(rgbN(:,2));
    bN_mean = mean(rgbN(:,3));

    hsv = reshape(rgb2hsv(image), num_pixels, 3);

    for c = 1 : num_pixels
        rN = rgbN(c,1);
        gN = rgbN(c,2);
        bN = rgbN(c,3);
        hue = hsv(c,1) * 360;
        value = hsv(c,3) * 100;
        % current pixel is red
        if      (hue >= 330 || hue <= 30) && ...
                (value >= 15) && ...
                (normal_prob(rN, rN_mean, rN_sdev) <  0.00001)
                
                    image_mask(c,1) = 1;
        % current pixel is green
        elseif  (hue >= 80 && hue < 150) && ...
                (value >= 20) && ...
                (normal_prob(gN, gN_mean, gN_sdev) <  0.005)
                   image_mask(c,2) = 1;
        % current pixel is blue
        elseif  (hue >= 150 && hue <= 270) && ...
                (value >= 26) && ...
                (normal_prob(bN, bN_mean, bN_sdev) < 0.0000075)
                   image_mask(c,3) = 1;
        end
    end

    image_mask = reshape(image_mask, num_rows, num_cols, 3);

    image_mask = remove_noise(image_mask);
    image_mask = remove_outliers(image_mask);
    image_mask = enforce_similar_channel_areas(image_mask);
end


function x = normal_prob(val, mu, sigma)
    x = 1.0 / (sigma * sqrt(2 * pi)) * exp(-(val - mu) ^ 2 / (2 * sigma ^ 2));
end


function image = remove_noise(image)
    [~, ~, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        channel = bwmorph(channel, 'majority', Inf);
        channel = bwmorph(channel, 'bridge', Inf);
        image(:,:,c) = channel;
    end
end


% finds the connected components in each channel of |image| and removes those
% that are far away from the centroid of the pixels in that channel
% here, 'far away' means more distant than |distance_proprtion_threshold| times
% the average distance of each connected component to the channel centroid
% this removes big areas of noise such as the big green blob inside of the black
% arrow of the red robot in data/1/00000006.jpg
function image = remove_outliers(image, distance_proportion_threshold)
    if nargin < 2
        distance_proportion_threshold = 2.0;
    end

    [~, ~, num_channels] = size(image);
    for c = 1 : num_channels
        channel = image(:,:,c);
        % skip empty channels
        if ~any(channel(:))
            continue;
        end
        channel_properties = regionprops(channel, 'Centroid');
        channel_centroid = channel_properties.Centroid;
        regions = bwconncomp(channel);
        regions_properties = regionprops(regions, 'Centroid', 'PixelIdxList');
        regions_centroids = {regions_properties.Centroid};
        distances = cellfun(@(x) norm(x - channel_centroid), regions_centroids);
        mean_distance = mean(distances);
        for d = 1 : length(distances)
            if distances(d) > mean_distance * distance_proportion_threshold
                idx = regions_properties(d).PixelIdxList;
                channel(idx) = 0;
            end
        end
        image(:,:,c) = channel;
    end
end


% we know that the robots are all about the same size - we can thus remove any
% channels in the mask that have a much smaller area than the other channels
% this catches some problems like the shadow of the blue robot in 
% data/1/00000095.jpg being detected as a red blob
function image = enforce_similar_channel_areas(image, area_proportion_threshold)
    if nargin < 2
        area_proportion_threshold = 0.5;
    end

    [~, ~, num_channels] = size(image);
    areas = zeros(num_channels, 1);
    for c = 1 : num_channels
        channel = image(:,:,c);
        if ~any(channel(:))
            continue;
        end
        region_props = regionprops(channel, 'Area');
        areas(c) = region_props.Area;
    end
    max_area = max(areas(areas > 0));
    for c = 1 : num_channels
        area = areas(c);
        if area == 0
            continue;
        end
        if  (area < max_area * area_proportion_threshold) || ...
            (area > max_area / area_proportion_threshold)
                image(:,:,c) = image(:,:,c) * 0;
        end
    end
end
