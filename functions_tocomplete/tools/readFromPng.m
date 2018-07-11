function [Mp,n,m] = readFromPng(filename)
% This function reads the image from the png file, and displays it

M = imread('../images/0-brd.png'); % Read picture
figure; imshow(M); % Display picture

M1 = M(:,:,1); M2 = M(:,:,2); M3 = M(:,:,3); % Split RGB
n = size(M,1); m = size(M,2); % Get the size

% Transform the matrices into vectors and put them as a 3-MRHS
Mp = [M1(:),M2(:),M3(:)];
Mp = double(Mp)/double(intmax('uint16')); % And pass form uint16 to float format

end
