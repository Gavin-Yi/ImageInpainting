clc;
clear;
close all;

imgin = im2double(imread('./test.png'));

[imh, imw, nb] = size(imgin);
assert(nb==1);
% the image is grayscale

V = zeros(imh, imw);
V(1:imh*imw) = 1:imh*imw;
% V(y,x) = (y-1)*imw + x
% use V(y,x) to represent the variable index of pixel (x,y)
% Always keep in mind that in matlab indexing starts with 1, not 0

%TODO: initialize counter, A (sparse matrix) and b.
i = zeros(1);  
j = zeros(1);
v = zeros(1); 

%TODO: fill the elements in A and b, for each pixel in the image
% Inner
n = 1:imh*imw;
i(n) = 1:imw*imh;
j(n) = 1:imw*imh;
v(n) = 1 + 1*( (1<n)&(n<imh) | ((imw-1)*imh<n)&(n<imw*imh) | (mod(n,imh)==1&n~=1&n~=((imw-1)*imh+1)) | mod(n,imh)==0&n~=imh&n~=(imw*imh) )  ...
    + 3*( ~( (1<=n)&(n<=imh) | ((imw-1)*imh<=n)&(n<=imw*imh) | (mod(n,imh)==1&n~=1&n~=((imw-1)*imh+1)) | mod(n,imh)==0&n~=imh&n~=(imw*imh) ) ) ;
% top
n = (imh*imw+1):2*imh*imw;
i(n) = 1:imh*imw;
j(n) = (n - (imw*imh)-1 ).*( (1<mod(n - (imw*imh),imh)) )+1*(~( (1<mod(n - (imw*imh),imh)) ));
v(n) = -1 * ( (1<mod(n - (imw*imh),imh)) );
% bottom
n = (2*imh*imw+1):3*imh*imw;
i(n) = 1:imh*imw;
j(n) = (n - (2*imw*imh)+1 ).*( 1<(mod(n - (2*imw*imh),imh)) )+1*(~( 1<(mod(n - (2*imw*imh),imh)) ));
v(n) = -1 * ( 1<(mod(n - (2*imw*imh),imh)) );
% left
n = (3*imh*imw+1):4*imh*imw;
i(n) = 1:imh*imw;
j(n) = ( n - (3*imw*imh)-imh ).*( ((n - (3*imw*imh)) >imh)&((n - (3*imw*imh)) <(imw-1)*imh+1) )+1*(~( ((n - (3*imw*imh)) >imh)&((n - (3*imw*imh)) <(imw-1)*imh+1) ));
v(n) = -1 * ( ((n - (3*imw*imh)) >imh)&((n - (3*imw*imh)) <(imw-1)*imh+1) );
% right
n = (4*imh*imw+1):5*imh*imw;
i(n) = 1:imh*imw;
j(n) = ( n - (4*imw*imh)+imh ).*( ((n - (4*imw*imh)) < ((imw-1)*imh+1))&((n - (4*imw*imh)) > imh) )+1*(~( ((n - (4*imw*imh)) < ((imw-1)*imh+1))&((n - (4*imw*imh)) > imh) ));
v(n) = -1 * ( ((n - (4*imw*imh)) < ((imw-1)*imh+1))&((n - (4*imw*imh)) > imh) );

% A1 and b
A = sparse(i,j,v);
b = A*(reshape(imgin,[imw*imh,1]));

%TODO: add extra constraints
% right - bottom brightness
b(imh*imw) = b(imh*imw);
% right - top brightness
b(imh*(imw-1)+1) = b(imh*(imw-1)+1);
% left -bottom brightness
b(imh) = b(imh)+0.5;
% left - top brightness
b(1) = b(1)+0.5;

%TODO: solve the equation
%use "lscov" or "\", please google the matlab documents
solution = A\b;
error = sum(abs(A*solution-b));
disp(error)
imgout = reshape(solution,[imh,imw]);

imwrite(imgout,'output4.png');
figure(), hold off, imshow(imgout);
figure(), imshow(imgin);

