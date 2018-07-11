function [Dx,Dy] = derivative(n,m)
% Create a derivative norm operator

Dx = sparse(n*(m-1),n*m); Dy = sparse(m*(n-1),n*m);

j = 1:m;
for i=1:n-1
   Dx( (j-1)*(n-1) + i, (j-1)*n + i )   = 1*eye(m);
   Dx( (j-1)*(n-1) + i, (j-1)*n + i+1 ) = -1*eye(m);
end

i = 1:n;
for j=1:m-1
   Dy((j-1)*n + i,(j-1)*n + i) = 1*eye(n);
   Dy((j-1)*n + i,j*n + i)     = -1*eye(n);
end

end
