function image = random_image(path, image_type)
    if nargin < 2
        image_type = 'jpg';
    end

    files = dir(sprintf('%s/*.%s', path, image_type));
    filenames = {files.name};
    [~, num_files] = size(filenames);

    random_index = randi([1 num_files]);
    filename = sprintf('%s/%s', path, filenames{random_index});
    image = imread(filename);
