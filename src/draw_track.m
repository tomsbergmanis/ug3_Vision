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

    
    bg=median_filter(path, image_type, 1, 10);
    gbg=rgb2gray(bg);
    bg=cat(3,gbg,gbg,gbg);
    
    rc=zeros(num_files-start_offset,2);
    gc=zeros(num_files-start_offset,2);
    bc=zeros(num_files-start_offset,2);
       
    for i = 1 + start_offset : num_files
        image = imread(sprintf('%s/%s', path, filenames{i}));
       	[~, centroids,~,~,tri] = analyse_image(image);
                            
        r = uint16(centroids(1,:));
        if r(1)>0 && r(2)>0
            bg=draw_center(bg,1, r(1),r(2));
            rc(i,:)=[r(1),r(2)];
        end
        
        g = uint16(centroids(2,:));
        if g(1)>0 && g(2)>0
            bg=draw_center(bg,2, g(1),g(2));
            gc(i,:)=[g(1),g(2)];
        end

 	    b = uint16(centroids(3,:));
        if b(1)>0 && b(2)>0
            bg=draw_center(bg,3, b(1),b(2));
            bc(i,:)=[b(1),b(2)];
        end
    end
    
    rc = rc(any(rc,2),:);
    gc = gc(any(gc,2),:);
    bc = bc(any(bc,2),:);
   
    bg=overlay_polygon(bg, rc, [255,255,255]);
    bg=overlay_polygon(bg, gc, [255,255,255]);
    bg=overlay_polygon(bg, bc, [255,255,255]);
    
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
        image(x,y+1,which)=255;
        image(x,y+1,other1)=0;
        image(x,y+1,other2)=0;
    end
    
    if y-1>0 
        image(x,y-1,which)=255;
        image(x,y-1,other1)=0;
        image(x,y-1,other2)=0;
    end
    
    if x+1<w
        image(x+1,y,which)=255;
        image(x+1,y,other1)=0;
        image(x+1,y,other2)=0;
    end   
    
    if x-1>0
        image(x-1,y,which)=255;
        image(x-1,y,other1)=0;
        image(x-1,y,other2)=0;
    end    
    
    image(x,y,which)=255;
    image(x,y,other1)=0;
    image(x,y,other2)=0;
    
end
