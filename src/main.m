function res = main(path, image_type, start_offset, time_step)
    if nargin < 2
        image_type = 'jpg';
    end
    if nargin < 3
        start_offset = 0;
    end

    files = dir(sprintf('%s/*.%s', path, image_type));
    filenames = {files.name};
    [~, num_files] = size(filenames);
    
    [col row ~]=size(imread(sprintf('%s/%s', path, filenames{1})));
    track_mask=zeros(col,row);
    
    rc=zeros(num_files-start_offset,2);
    gc=zeros(num_files-start_offset,2);
    bc=zeros(num_files-start_offset,2);
    
    reds = cell(num_files - start_offset, 1);
    greens = cell(num_files - start_offset, 1);
    blues = cell(num_files - start_offset, 1);
       
    for i = 1 + start_offset : num_files
        image = imread(sprintf('%s/%s', path, filenames{i}));
       	[directions, centroids, ~, masks] = analyse_image(image);
                            
        r = uint16(centroids(1,:));
        if r(1)>0 && r(2)>0
            track_mask=draw_center(track_mask,1, r(1),r(2));
            rc(i,:)=[r(1),r(2)];
        end
        
        g = uint16(centroids(2,:));
        if g(1)>0 && g(2)>0
            track_mask=draw_center(track_mask,2, g(1),g(2));
            gc(i,:)=[g(1),g(2)];
        end

 	    b = uint16(centroids(3,:));
        if b(1)>0 && b(2)>0
            track_mask=draw_center(track_mask,3, b(1),b(2));
            bc(i,:)=[b(1),b(2)];
        end
        
        reds{i} = image(:,:,1);
        greens{i} = image(:,:,2);
        blues{i} = image(:,:,3);
        imshow(overlay_mask(image,directions));
        pause(time_step);
    end
    
    rc = rc(any(rc,2),:);
    gc = gc(any(gc,2),:);
    bc = bc(any(bc,2),:);
    
    red_median = median(cat(3, reds{:}), 3);
    green_median = median(cat(3, greens{:}), 3);
    blue_median = median(cat(3, blues{:}), 3);

    bg = cat(3, red_median, green_median, blue_median);
    
    bg=overlay_polygon(bg, rc, [255,255,255]);
    bg=overlay_polygon(bg, gc, [255,255,255]);
    bg=overlay_polygon(bg, bc, [255,255,255]);
    bg=overlay_mask(bg, track_mask);
    imshow(bg)
end


function image = draw_center(image,which, x, y)
    [h w ~] = size(image);
    
    if which==1
        other1=2;
        other2=3;
    end
    
    if which==2
        other1=1;
        other2=3;
    end
    
    if which==3
       other1=1;
       other2=2;
    end
    
    if y+1<h 
        image(x,y+1,which)=1;
        image(x,y+1,other1)=0;
        image(x,y+1,other2)=0;
    end
    
    if y-1>0 
        image(x,y-1,which)=1;
        image(x,y-1,other1)=0;
        image(x,y-1,other2)=0;
    end
    
    if x+1<w
        image(x+1,y,which)=1;
        image(x+1,y,other1)=0;
        image(x+1,y,other2)=0;
    end   
    
    if x-1>0
        image(x-1,y,which)=1;
        image(x-1,y,other1)=0;
        image(x-1,y,other2)=0;
    end    
    
    image(x,y,which)=1;
    image(x,y,other1)=0;
    image(x,y,other2)=0;
    
end
