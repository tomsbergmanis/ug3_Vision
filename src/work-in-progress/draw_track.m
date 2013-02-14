%i'll improve on this, it's not finished

function res = draw_track(path, image_type, start_offset)
    if nargin < 2
        image_type = 'jpg';
    end
    if nargin < 3
        start_offset = 0;
    end

   

    files = dir(sprintf('%s/*.%s', path, image_type));
    filenames = {files.name};
    [~, num_files] = size(filenames);

    bg=median_filter(path, image_type, 50);
   
    for i = 1 + start_offset : num_files

        image = imread(sprintf('%s/%s', path, filenames{i}));
        r=identify_robots(image, 1);
        g=identify_robots(image, 2);
        b=identify_robots(image, 3);
        if numel(r)>1
            bg(uint16(r(2)),uint16(r(1)),1)=255;
            bg(uint16(r(2)),uint16(r(1)),2)=0;
            bg(uint16(r(2)),uint16(r(1)),3)=0;
        end
       
        if numel(g)>1
            bg(uint16(g(2)),uint16(g(1)),2)=255;
            bg(uint16(g(2)),uint16(g(1)),1)=0;
            bg(uint16(g(2)),uint16(g(1)),3)=0;
        end
        
        if numel(b)>1
            bg(uint16(b(2)),uint16(b(1)),3)=255;
            bg(uint16(b(2)),uint16(b(1)),1)=0;
            bg(uint16(b(2)),uint16(b(1)),2)=0;
        end
        
    end
    imshow(bg)
end