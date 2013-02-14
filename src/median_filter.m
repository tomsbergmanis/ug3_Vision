% generates a background image from a set of sample images
% subtracting the background from the sample images eases object detection
% does not work on these datasets because the blue/cyan robot does not move
% for most of the images i.e. will be considered part of the background
function background = median_filter(path, image_type, start_offset)
    if nargin < 2
        image_type = 'jpg';
    end
    if nargin < 3
        start_offset = 0;
    end

    dim = 2;

    files = dir(sprintf('%s/*.%s', path, image_type));
    filenames = {files.name};
    [~, num_files] = size(filenames);

    reds = cell(num_files - start_offset, 1);
    greens = cell(num_files - start_offset, 1);
    blues = cell(num_files - start_offset, 1);
    for c = 1 + start_offset : num_files
        image = imread(sprintf('%s/%s', path, filenames{c}));
        reds{c} = image(:,:,1);
        greens{c} = image(:,:,2);
        blues{c} = image(:,:,3);
    end
    red_median = median(cat(dim + 1, reds{:}), dim + 1);
    green_median = median(cat(dim + 1, greens{:}), dim + 1);
    blue_median = median(cat(dim + 1, blues{:}), dim + 1);

    background = cat(dim + 1, red_median, green_median, blue_median);
end
