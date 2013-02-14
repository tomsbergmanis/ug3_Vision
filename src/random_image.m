% returns a random image of type |image_type| found at |path|
% if |path| is not specified, set it to a random folder in "ug3_Vision/data"
function image = random_image(path, image_type)
    if nargin < 1
        path = random_dir('ug3_Vision', 'data');
    end
    if nargin < 2
        image_type = '.jpg';
    end

    if ~strcmp(image_type(1), '.')
        image_type = strcat('.', image_type);
    end

    image = imread(random_file(path, image_type));
end


% looks for a directory |tl_dir| somewhere up from this file's location
% returns the absolute path to one of the directories in |tl_dir|/|branch_dir|
function path = random_dir(tl_dir, branch_dir)
    cur_path = mfilename('fullpath');
    tl_path = cur_path(1:strfind(cur_path, tl_dir) + length(tl_dir));
    branch_path = strcat(tl_path, branch_dir);
    branch_path_contents = dir(branch_path);
    branch_path_dirs = {};
    for c = 1 : length(branch_path_contents)
        elem = branch_path_contents(c);
        if ~elem.isdir || strcmp(elem.name, '.') || strcmp(elem.name, '..')
            continue;
        end
        branch_path_dirs{end + 1} = fullfile(branch_path, elem.name);
    end
    path = branch_path_dirs{randi([1 length(branch_path_dirs)])};
end


% returns a random file ending with |extension| found at |path|
function path = random_file(path, extension)
    if nargin < 2
        extension = '';
    end

    file_type = strcat('*', extension);
    files = dir(fullfile(path, file_type));
    filenames = {files.name};

    path = fullfile(path, filenames{randi([1 length(filenames)])});
end
