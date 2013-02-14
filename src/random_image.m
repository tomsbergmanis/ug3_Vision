% returns a random image from some directory D and print the path to that image
% the function understands the following options:
%   'Quiet'             don't print the path to the image
%   'ImageType', C      look for images of type C (default = "jpg")
% any remaining parameters are taken to be the path to D
% if D is not specified, default to a random sub-directory of "ug3_Vision/data"
function image = random_image(varargin)
    TL_DIR = 'ug3_Vision';
    BRANCH_DIR = 'data';
    % parse options
    argc = length(varargin);
    c = 1;
    while c <= argc
        arg = varargin{c};
        if strcmpi(arg, 'Quiet')
            quiet = 1;
        elseif strcmpi(arg, 'ImageType')
            if c + 1 > argc
                error('ImageType option should be followed by a string');
            end
            image_type = varargin{c + 1};
            c = c + 1;
        else
            path = arg;
        end
        c = c + 1;
    end
    % option defaults
    if ~exist('quiet', 'var')
        quiet = 0;
    end
    if ~exist('image_type', 'var')
        image_type = '.jpg';
    end
    if ~exist('path', 'var')
        path = random_dir(TL_DIR, BRANCH_DIR);
    end

    if ~strcmp(image_type(1), '.')
        image_type = strcat('.', image_type);
    end

    image_path = random_file(path, image_type);
    image = imread(image_path);

    if ~quiet
        idx = strfind(image_path, TL_DIR);
        if isempty(idx)
            idx = 1;
        end
        fprintf(1, 'image = %s\n', image_path(idx : end));
    end
end


% looks for a directory |tl_dir| somewhere up from this file's location
% returns the absolute path to one of the directories in |tl_dir|/|branch_dir|
function path = random_dir(tl_dir, branch_dir)
    cur_path = mfilename('fullpath');
    tl_path = cur_path(1 : strfind(cur_path, tl_dir) + length(tl_dir));
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
