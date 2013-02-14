%writes out images to test filter robot colors on whole dataset
function res=test_filter_robot_colors(path, topath, image_type, start_offset)
    if nargin < 2
        image_type = 'jpg';
    end
    if nargin < 3
        start_offset = 0;
    end

    files = dir(sprintf('%s/*.%s', path, image_type));
    filenames = {files.name};
    [~, num_files] = size(filenames);

   
    for c = 1 + start_offset : num_files
        image = imread(sprintf('%s/%s', path, filenames{c}));
        
        imwrite(filter_robot_colors(image),sprintf('%s/%s', topath, filenames{c}))
    end
   
end
