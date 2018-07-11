function A = operator(n,m,sigma,t)
% Create the Gaussian blurring opearator

A = zeros(n*m);

tic
Ax = zeros(n*m); Ay = zeros(n*m);
j = 1:m; l = 1:m;
for i=1:n
   ji = (j-1)*n + i;
   kmin = max(1,i-t); kmax = min(n,i+t);
   for k=kmin:kmax
      lk = (l-1)*n + k;
      Ax(ji,lk) = exp(-(k-i)^2/(2*sigma^2))*eye(m);
   end
end

i = 1:n; k = 1:n;
for j=1:m
   ji = (j-1)*n + i;
   lmin = max(1,j-t); lmax = min(m,j+t);
   for l=lmin:lmax
      lk = (l-1)*n + k;
      Ay(ji,lk) = exp(-(l-j)^2/(2*sigma^2))*eye(n);
   end
end

A = Ax*Ay;
for z=1:size(A,1)
   A(z,:) = A(z,:)/sum(A(z,:)); % Normalize
end
toc

end
