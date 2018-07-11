function [M,n,m] = readFromMat(filename)
% This function reads the image from the matrix file, and displays it

Mp = load(filename);
n = Mp.n; m = Mp.m; M = Mp.Mp;  % Recover the matricial data

Mb1 = reshape(M(:,1),[n,m]); % Re-put the picture in the rgb-matrix format
Mb2 = reshape(M(:,2),[n,m]);
Mb3 = reshape(M(:,3),[n,m]);
Mb = zeros(n,m,3);
Mb(:,:,1) = Mb1; Mb(:,:,2) = Mb2; Mb(:,:,3) = Mb3;
figure; imshow(Mb); % And display it

end
