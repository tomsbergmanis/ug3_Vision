function res = main(path, varargin)
    % parse options
    image_type = 'jpg';
    start_offset = 0;
    time_step = 0.05;
    show_robots = 1;
    argc = size(varargin, 2);
    c = 1;
    while c <= argc
        arg = varargin{c};
        if strcmpi(arg, 'ImageType')
            if c + 1 > argc
                error('ImageType options should be followed by a string');
            end
            image_type = varargin{c + 1};
            c = c + 1;
        elseif strcmpi(arg, 'StartOffset')
            if c + 1 > argc
                error('StartOffset options should be followed by an integer');
            end
            start_offset = varargin{c + 1};
            c = c + 1;
        elseif strcmpi(arg, 'TimeStep')
            if c + 1 > argc
                error('Timstep options should be followed by a number');
            end
            time_step = varargin{c + 1};
            c = c + 1;
        elseif strcmpi(arg, 'NoShowRobots')
            show_robots = 0;
        end
        c = c + 1;
    end

    curpath = mfilename('fullpath');
    curpath = fileparts(curpath);
    addpath(fullfile(curpath, 'algo'));
    addpath(fullfile(curpath, 'draw'));
    addpath(fullfile(curpath, 'test'));

    files = dir(sprintf('%s/*.%s', path, image_type));
    filenames = {files.name};
    [~, num_files] = size(filenames);

    [m n ~] = size(imread(sprintf('%s/%s', path, filenames{1})));
    track_mask = zeros(m, n);

    rc = zeros(num_files - start_offset, 2);
    gc = zeros(num_files - start_offset, 2);
    bc = zeros(num_files - start_offset, 2);

    reds = cell(num_files - start_offset, 1);
    greens = cell(num_files - start_offset, 1);
    blues = cell(num_files - start_offset, 1);

    for i = 1 + start_offset : num_files
        image = imread(sprintf('%s/%s', path, filenames{i}));
        [direction_mask, centroids, ~, robot_mask] = analyse_image(image);
        for dim = 1 : size(robot_mask, 3)
            robot_mask(:,:,dim) = bwmorph(robot_mask(:,:,dim), 'remove');
        end
        r = uint16(centroids(1,:));
        if r(1) > 0 && r(2) > 0
            track_mask = overlay_cross(track_mask, 1, r(1), r(2));
            rc(i,:) = [r(1), r(2)];
        end

        g = uint16(centroids(2,:));
        if g(1) > 0 && g(2) > 0
            track_mask = overlay_cross(track_mask, 2, g(1), g(2));
            gc(i,:) = [g(1), g(2)];
        end

        b = uint16(centroids(3,:));
        if b(1) > 0 && b(2) > 0
            track_mask = overlay_cross(track_mask, 3, b(1), b(2));
            bc(i,:) = [b(1), b(2)];
        end

        reds{i} = image(:,:,1);
        greens{i} = image(:,:,2);
        blues{i} = image(:,:,3);
        analysed_image = overlay_mask(image, direction_mask);
        if show_robots == 1
            analysed_image = overlay_mask(analysed_image, robot_mask);
        end
        imshow(analysed_image);

        pause(time_step);
    end

    rc = rc(any(rc,2),:);
    gc = gc(any(gc,2),:);
    bc = bc(any(bc,2),:);

    red_median = median(cat(3, reds{:}), 3);
    green_median = median(cat(3, greens{:}), 3);
    blue_median = median(cat(3, blues{:}), 3);

    bg = cat(3, red_median, green_median, blue_median);
    
    bg = overlay_polygon(bg, rc, [255, 255, 255]);
    bg = overlay_polygon(bg, gc, [255, 255, 255]);
    bg = overlay_polygon(bg, bc, [255, 255, 255]);
    bg = overlay_mask(bg, track_mask);
    imshow(bg)
end
