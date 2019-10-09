clc;
clear;
close all;

im_background = im2double(imread('./background.jpeg'));
im_object = im2double(imread('./target.png'));

% get source region mask from the user
objmask = get_mask(im_object);

% align im_s and mask_s with im_background
[im_s, mask_s] = align_source(im_object, objmask, im_background);

clone = im_s.*(mask_s==1)+im_background.*(mask_s~=1);
imshow(clone);
imwrite(clone,'output1.png');

% blend
disp('start');
im_blend = poisson_blend(im_s, mask_s, im_background);
disp('end');

imwrite(im_blend,'output2.png');
figure(), hold off, imshow(im_blend);
