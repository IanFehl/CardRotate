% Group members: Ian Fehl and Mary Parker
clearvars
close all

% Change this variable to the folder of images you want to test
path = 'C:\Users\ianfe\Documents\MATLAB\Image Processing\Project2-Card_Rotate\ImageSet2';

dinfo = dir(fullfile(path)); % directory path
dinfo([dinfo.isdir]) = []; % removes "." and ".." from the dinfo struct
numfiles = length(dinfo); % Determines the amount of files in the folder

for j = 1 : numfiles % look at every file in the entire folder
    file = fullfile(path,dinfo(j).name);
    I = imread(file); % Read in the image file
    
    rot_crop_img = rotate_crop(I); % returns rotated and cropped image
    
    figure
    subplot(1,2,1); % plot original image
    imshow(I)
    title('Input image')
    
    subplot(1,2,2); % plot rotated and cropped image
    imshow(rot_crop_img)
    title('Output image - aligned and cropped')
end

%--------------------------------------------------------------------------
%                                Functions
%--------------------------------------------------------------------------

% rotate_crop(I) takes in the input image and calls all other functions in 
% the program to return the rotated and cropped image
function rot_crop_img = rotate_crop(I)
    rot_img = rotate_image(I); % returns rotated image
    thresh_img = thresh(rot_img); % threshold the rotated image
    [ymin,ymax,xmin,xmax] = final_crop(thresh_img); % get coordinates for  cropping
    rot_crop_img = rot_img(ymin:ymax,xmin:xmax); % crop the rotated image and return it
end

% edge_det(I) uses sobel kernels to perform edge detection on the input image
function [image, dir] = edge_det(I)
    hor1 = [-1 -2 -1;
             0  0  0;
             1  2  1];
         
    hor2 = [ 1  2  1;
             0  0  0;
            -1 -2 -1];
    
    vert1 = [-1 0 1; 
             -2 0 2; 
             -1 0 1];
        
    vert2 = [1 0 -1;
             2 0 -2;
             1 0 -1];

    % convolve input image with all 4 filters and convert images from
    % double to uint8
    hor_im1 = conv2(I, hor1, 'same');
    hor_im1 = uint8(hor_im1);
    hor_im2 = conv2(I, hor2, 'same');
    hor_im2 = uint8(hor_im2);

    vert_im1 = conv2(I, vert1, 'same');
    vert_im1 = uint8(vert_im1);
    vert_im2 = conv2(I, vert2, 'same');
    vert_im2 = uint8(vert_im2);

    % call imfuse to combine the resulting images to get images of the
    % horizontal and vertical edges. imfuse returns a rgb image so it must
    % be grayscaled afterwards
    hor = rgb2gray(imfuse(hor_im1, hor_im2));
    vert = rgb2gray(imfuse(vert_im1, vert_im2));

    % image that contains all vertical and horizontal edges in the input
    image = rgb2gray(imfuse(hor,vert)); 
    
    % convert from uint8 to double to use for atan2d
    hor = double(hor); 
    vert = double(vert);
    
    dir = atan2d(vert,hor); % direction of the gradient
end

% crop(I) uses regionprops to find and return the coordinates of the
% bounding box that surrounds the card 
function [ymin,ymax,xmin,xmax] = crop(I)
    image_bi = imbinarize(I, 'global');
    
    box = regionprops(image_bi, 'BoundingBox');
    for k = 1 : length(box) % look at all the bounding boxes
         BB = box(k).BoundingBox;
         if BB(3) > 50 && BB(4) > 50 % only uses large enough bounding boxes
            xmin = ceil(BB(1)); % BB(1) = x coordinate of top left pixel in box
            ymin = ceil(BB(2)); % BB(2) = y coordinate of top left pixel in box
            xmax = xmin + BB(3) - 1; % BB(3) = the width of the bounding box
            ymax = ymin + BB(4) - 1; % BB(4) = the height of the bounding box
         end
    end
end

% rotate_image(I) calculates the rotation angle needed to straighten the
% card vertically and rotates the input image accordingly
function rot_img = rotate_image(I)
    [ymin,ymax,xmin,xmax] = crop(I); % find bounding box coordinates
    crop_img = I(ymin:ymax,xmin:xmax); % crop input image to those coordinates
    new_image_bi = imbinarize(crop_img, 'global');
    
    % first pixel in the first row of the cropped image, i.e. one corner of
    % the card
    first_row= find(new_image_bi(1,:),1,'first');
%     last_row= find(new_image_bi(size(new_image_bi,1),:),1,'last');

    % first pixel in the first column of the cropped image, i.e. another 
    % corner of the card
    first_col = find(new_image_bi(:,1),1,'first');
%     last_col= find(new_image_bi(:,size(new_image_bi,2)),1,'last');

    % euclidean distance between the two corners
    distance = sqrt( (first_row - 1)^2 + (first_col - 1)^2); 

    % calculate the angle of rotation needed based on the angle of either
    % the long side of the card or the short side
    if distance >= 450
        angle = 90 + atand( (1 - first_col) / (first_row - 1));
    else
        angle = atand( (1 - first_col) / (first_row - 1));
    end

    rot_img = imrotate(I, angle); % rotated image by the calculated angle
end

% final_crop(I) returns the coordinates needed to crop the card from the
% rotated and previously cropped card. this is similar to crop(I), but the
% main difference is final_crop(I) returns the coordinates of the biggest
% bounding box
function [ymin,ymax,xmin,xmax] = final_crop(I)    
    stats = regionprops(I, 'Centroid', 'Area', 'BoundingBox');
    box = stats.BoundingBox;
    
    xmin = ceil(box(1)); % box(1) = x coordinate of top left pixel in box
    ymin = ceil(box(2)); % box(2) = y coordinate of top left pixel in box
    xmax = xmin + box(3) - 1; % box(3) = the width of the bounding box
    ymax = ymin + box(4) - 1; % box(4) = the height of the bounding box
end

 % thresh(I) removes the grayish background around the card. all gray
 % values around the card were found to have an intensity value lower than
 % 130.
function thresh_img = thresh(I)
    temp_img = I;
    temp_img(I < 130) = 0;
    thresh_img = edge_det(temp_img);
end