function [Itere,residual,regulari] = myorthodir(A,b,niter,L)
% Homemade OrhtoDir ('cause I want to store all the solutions and all the norms)

ndofs = size(b,1); ncase = size(b,2);

Itere = zeros( ndofs, ncase*(niter+1) );
p     = zeros( ndofs, ncase*(niter+1) );
q     = zeros( ndofs, ncase*(niter+1) );
Delta = zeros( ncase*(niter+1), ncase );
% Yep, I'm storing all those guys. My name is Leak. Memory Leak

residual = zeros(niter+1,ncase);
regulari = zeros(niter+1,ncase);

%Non-zero initialization would be there
atimesItere = A*Itere(:,[1:ncase]);

Res = b - atimesItere;
p(:,[1:ncase]) = Res;

residual(1,:) = sqrt( diag( Res'*Res ) ) ;
regulari(1,:) = sqrt( diag( Itere(:,[1:ncase])'*L*Itere(:,[1:ncase]) ) );

q(:,[1:ncase]) = A*p(:,[1:ncase]);

for iter = 1:niter
   multindex   = [ncase*(iter-1)+1:ncase*iter];
   multindexp1 = multindex + ncase;

%   if rank(Delta) ~= size(Delta,1)
%      warning('Your test-cases are of lower dimension as their number. It means that you should not do a multiple Rhs.');
%   end

   Delta(multindex,:)   = q(:,multindex)'*q(:,multindex);
   gammai               = q(:,multindex)'*Res;
   alphai               = pinv(Delta(multindex,:))*gammai;

   Itere(:,multindexp1) = Itere(:,multindex) + p(:,multindex)*alphai;
   Res                  = Res - q(:,multindex)*alphai;

   residual(iter+1,:) = sqrt( diag( Res'*Res )) ;
   regulari(iter+1,:) = sqrt( diag ( ...
                               Itere(:,multindexp1)'*L*Itere(:,multindexp1) ) );

   Ari = A*Res;

   %% Orthogonalization
   p(:,multindexp1) = Res;
   q(:,multindexp1) = Ari;

   for jter=1:iter
      multjndex = [ncase*(jter-1)+1:ncase*jter];
      phiij  = q(:,multjndex)'*Ari;
      betaij = pinv(Delta(multjndex,:))*phiij; % You could spare a little CPU
      % by storing pinv(Delta) but not so much as Delta is size 3 xP
      p(:,multindexp1) = p(:,multindexp1) - p(:,multjndex) * betaij;
      q(:,multindexp1) = q(:,multindexp1) - q(:,multjndex) * betaij;
  end

end

end
