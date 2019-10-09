function imgout = poisson_blend(im_s, mask_s, im_t, split)
% im_t = im_background;
% -----Input
% im_s     source image (object)
% mask_s   mask for source image (1 meaning inside the selected region)
% im_t     target image (background)
% -----Output
% imgout   the blended image
if nargin < 4
   split = 0 ;
end
[imh, imw, nb] = size(im_s);

%TODO: consider different channel numbers
if nb == 1 && split ==0
    imgend = reshape( poisson_blend(im_s(:,:,1), mask_s, im_t(:,:,1), 1), [imh,imw] );
    imgout = imgend.*(mask_s==1)+ im_t.*(mask_s~=1);
    return
end
if nb == 3
    img_R = reshape( poisson_blend(im_s(:,:,1), mask_s, im_t(:,:,1), 1), [imh,imw] );
    img_G = reshape( poisson_blend(im_s(:,:,2), mask_s, im_t(:,:,2), 1), [imh,imw] );
    img_B = reshape( poisson_blend(im_s(:,:,3), mask_s, im_t(:,:,3), 1), [imh,imw] );
    imgend(:,:,1) = img_R;
    imgend(:,:,2) = img_G;
    imgend(:,:,3) = img_B;
    imgout = imgend.*(mask_s==1)+ im_t.*(mask_s~=1);
    return
end

%TODO: initialize counter, A (sparse matrix) and b.
%Note: A don't have to be k¡Ák,
%      you can add useless variables for convenience,
%      e.g., a total of imh*imw variables
i = zeros(1);  
j = zeros(1);
v = zeros(1); 

%TODO: fill the elements in A and b, for each pixel in the image
counter = 1;
row = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if mask_s( (n-1)*imh+m )  == 1
            i(row) = counter;
            j(row) = (n-1)*imh+m;
            v(row) = 4;
            counter = counter + 1;
            row = row + 1;
        end
    end
end
counter = 1;
row = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if n==1 || n==imw || m==1 || m== imh
            i(row) = counter;
            j(row) = (n-1)*imh+m;
            v(row) = 1;
        else
            i(row) = counter;
            j(row) = (n-1)*imh+m;
            v(row) = 4;
        end
        counter = counter + 1;
        row = row + 1;
    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if n==1 || n==imw || m==1 || m== imh 
        else
            i(row) = counter;
            j(row) = (n-1)*imh+m+1;
            v(row) = -1;
            row = row + 1;
        end
        counter = counter + 1;
    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if n==1 || n==imw || m==1 || m== imh 
        else
            i(row) = counter;
            j(row) = (n-1)*imh+m-1;
            v(row) = -1;
            row = row + 1;
        end
        counter = counter + 1;

    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if n==1 || n==imw || m==1 || m== imh 
        else
            i(row) = counter;
            j(row) = (n-2)*imh+m;
            v(row) = -1;
            row = row + 1;
        end
        counter = counter + 1;

    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if n==1 || n==imw || m==1 || m== imh 
        else
            i(row) = counter;
            j(row) = (n)*imh+m;
            v(row) = -1;
            row = row + 1;
        end
        counter = counter + 1;
    end
end

A = sparse(i,j,v);
b = A*(reshape(im_s(:,:,1),[imh*imw,1]));

% imshow( reshape( b,[imh,imw]) );

im_t_new = im_t(:,:,1);
for n = 1:(imw)
    for m = 1:(imh)
        if mask_s( (n-1)*imh+m )  == 1
            b( (n-1)*imh+m ) = b( (n-1)*imh+m ) + im_t_new( (n-2)*imh+m )*( mask_s( (n-2)*imh+m )==0 ) + im_t_new( (n)*imh+m )*( mask_s( (n)*imh+m )==0 ) + im_t_new( (n-1)*imh+m+1 )*( mask_s( (n-1)*imh+m+1 )==0 ) + im_t_new( (n-1)*imh+m-1 )*( mask_s( (n-1)*imh+m-1 )==0 );       
        end
    end
end

% figure;
% imshow( reshape( b,[imh,imw]) );

i = zeros(1);  
j = zeros(1);
v = zeros(1);

counter = 1;
row = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if n==1 || n==imw || m==1 || m== imh
            i(row) = counter;
            j(row) = (n-1)*imh+m;
            v(row) = 1;
        else
            i(row) = counter;
            j(row) = (n-1)*imh+m;
            v(row) = 4; % *( mask_s( (n-1)*imh+m )==1 )
        end
        counter = counter + 1;
        row = row + 1;
    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if mask_s( (n-1)*imh+m )  == 1
            i(row) = counter;
            j(row) = (n-1)*imh+m+1;
            v(row) = -1; % *( mask_s( (n-1)*imh+m+1 )==1 )
            row = row + 1;
        end
        counter = counter + 1;
    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if mask_s( (n-1)*imh+m )  == 1
            i(row) = counter;
            j(row) = (n-1)*imh+m-1;
            v(row) = -1; %*( mask_s( (n-1)*imh+m-1 )==1 )
            row = row + 1;
        end
        counter = counter + 1;
    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if mask_s( (n-1)*imh+m )  == 1
            i(row) = counter;
            j(row) = (n-2)*imh+m;
            v(row) = -1; % *( mask_s( (n-2)*imh+m )==1 )
            row = row + 1;
        end
        counter = counter + 1;
    end
end
counter = 1;
for n = 1:(imw)
    for m = 1:(imh)
        if mask_s( (n-1)*imh+m )  == 1
            i(row) = counter;
            j(row) = (n)*imh+m;
            v(row) = -1; % *( mask_s( (n)*imh+m )==1 )
            row = row + 1;
        end
        counter = counter + 1;
    end
end

A2 = sparse(i,j,v);

%TODO: add extra constraints (if any)
%-----
%-----


%TODO: solve the equation
%use "lscov" or "\", please google the matlab documents
solution = A2\b;
% error = sum(abs(A*solution-b));
% disp(error)
imgout = solution;
%TODO: copy those variable pixels to the appropriate positions
%      in the output image to obtain the blended image
%-----
%-----
