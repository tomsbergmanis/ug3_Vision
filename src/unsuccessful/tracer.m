function trace = tracer(path, image_type, start_offset)
%scritp I was playing around with. 
if nargin < 2
    image_type = 'jpg';
end
if nargin < 3
    start_offset = 0;
end



files = dir(sprintf('%s/*.%s', path, image_type));
filenames = {files.name};
[~, num_files] = size(filenames);
frame = imread(sprintf('%s/%s', path, filenames{1}));
[x y ~]=size(frame);
tracet=cell(num_files-start_offset-1);

for c = 1 + start_offset : num_files-1
    frame = imread(sprintf('%s/%s', path, filenames{c}));
    nextFrame = imread(sprintf('%s/%s', path, filenames{c+1}));
    tracet{c}=im2bw(rgb2gray(frame-nextFrame),0.25);
    
    % tracet{c}=frame-nextFrame;
end
trace=tracet{1};
for c = 2 + start_offset : num_files-1
    trace=trace+tracet{c};
end
    

