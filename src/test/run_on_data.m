clear all;
clc;

if ~exist('TL_DIR', 'var')
    TL_DIR = 'ug3_Vision';
end
if ~exist('IN_DIR', 'var')
    IN_DIR = 'data';
end
if ~exist('FILTER_FN', 'var')
    FILTER_FN = 'analyse_image';
end
if ~exist('OUT_DIR', 'var')
    OUT_DIR = fullfile('res', strrep(FILTER_FN, '_', '-'));
end
filter_fn = str2func(FILTER_FN);
disp(sprintf('function: %s\ninput: %s\n', FILTER_FN, fullfile(TL_DIR, IN_DIR)));

curpath = mfilename('fullpath');
tlpath = curpath(1 : strfind(curpath, TL_DIR) + length(TL_DIR));
inpath = strcat(tlpath, IN_DIR);
outpath = strcat(tlpath, OUT_DIR);
inpath_contents = dir(inpath);
inpath_dirs = {};
for c = 1 : length(inpath_contents)
    elem = inpath_contents(c);
    if ~elem.isdir || strcmp(elem.name, '.') || strcmp(elem.name, '..')
        continue;
    end
    inpath_dirs{end + 1} = fullfile(inpath, elem.name);
end

num_dirs = length(inpath_dirs);
times = [];
for c = 1 : num_dirs
    in_dir = inpath_dirs{c};
    files = dir(strcat(in_dir, filesep, '*.jpg'));
    file_names = {files.name};
    [~, dir_name] = fileparts(in_dir);
    out_dir = fullfile(outpath, strcat(IN_DIR, '-', dir_name));
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    num_files = length(file_names);
    for d = 1 : num_files
        file_name = file_names{d};
        input = fullfile(in_dir, file_name);
        disp(sprintf('[dir %d/%d] [file %d/%d]', c, num_dirs, d, num_files));
        disp(sprintf('\tinput = %s', input(strfind(input, TL_DIR) : end)));
        image = imread(input);
        timer = tic;
        [directions, ~, ~, hull] = filter_fn(image);
        elapsed = toc(timer);
        times(end + 1) = elapsed;
        for dim = 1 : size(hull, 3)
            hull(:,:,dim) = bwmorph(hull(:,:,dim), 'remove');
        end
        output = fullfile(out_dir, file_name);
        image = overlay_mask(image, directions, 'Saturation', 0.75);
        image = overlay_mask(image, hull);
        imwrite(image, output, 'jpg');
        disp(sprintf('\toutput = %s', output(strfind(input, TL_DIR) : end)));
        disp(sprintf('\tprocessing time = %fs', elapsed));
    end
end

disp(sprintf('\naverage processing time per image: %fs', mean(times)));
