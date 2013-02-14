
function result = identify_robots(image, which)
    masks=filter_robot_colors(image); 
    centroid=regionprops(masks(:,:,which),'Centroid' );
    result=cat(1, centroid.Centroid);
    
    
end