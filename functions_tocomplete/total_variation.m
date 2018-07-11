close all;
clear all;
addpath(genpath('./tools'));

%[ Mp, n, m ] = readFromPng('../images/0-brd.png');
[ Mp, n, m ] = readFromMat('../images/0-brd.mat'); % Load the data

A = operator(n,m,5,max(n,m)); % Build the blur operator
sA = size(A,1); % His size

%% Let the fun begin
tic; % Initialize the chronometer
% There will be no L-curve for this one
% (unless someone's ready to wait 'till the end of times to get a result)

%% Sooo... Now the cowboyish part:
% People usually use weired quad optimization algorithms and I didn't want
% to dig into that, so here is the plan :
% 1) initialize with an under-regularized quadratic version
% 2) Transform a bit the absolute value into sqrt{x^2+\epsilon}
% 3) Use a steepest descent to minimize from the initialization, 
% BUUUUUUT:
% The algorithm is accelerated by the parameter kappa, that makes
% the first iterates to be faster, and relaxed by the parameter step.
% At each iteration, if the cost-function did not decrease, the iteration is
% canceled and kappa is decreased. If kappa is already 1, step is then decreased
% what's more, all that is done in multiple right hand side

epsi = 1e-12; % Parameter to smooth the absolute value
niter = 100000; % Nb of iterations
mu = % Regularization parameter
step = 1; % Initialize the realxation parameter
mu0 = 1e-14; % Regularization parameter for the first quadratic computation
kappa = 2^10; % Better let a power of 2
[Dx,Dy] = derivative(n,m); Isum = ones(size(Dx,1),3); Isu1 = ones(size(Dx,1),1);

% Underregularizad pre-computation
L = Dx'*Dx + Dy'*Dy;
a0 = (A+mu0*L*norm(A,'fro')/norm(L,'fro'))\(Mp);

Mr = a0; % Initialization
fc   = zeros(niter+1,3);
stag = zeros(niter+1,1);

AnAn = A'*A; AnMp = A'*Mp; MpMp = Mp'*Mp; % Pre-compute what can be

% And the cost-function
fc(1,:) = diag(.5*Mr'*AnAn*Mr - Mr'*AnMp + .5*MpMp + ...
               mu*(Isum'*sqrt((Dx*Mr).^2+epsi)+Isum'*sqrt((Dy*Mr).^2+epsi)));

dir = zeros(sA,3);% First direction of research

% Again some pre-computation
ix = Dx*Mr; iy = Dy*Mr; absx = sqrt(ix.^2+epsi); absy = sqrt(iy.^2+epsi);

for i=1:niter % Multiple Rhs (life without risk is not real life)
   Mro = Mr;
   Sx = ix./(absx); Sy = iy./(absy);
   res = AnAn*Mr - AnMp + kappa*mu*( Dx'*(Sx.*Isum) + Dy'*(Sy.*Isum) );
   Mr = Mr - step * res;
   ix = Dx*Mr; iy = Dy*Mr; absx = sqrt(ix.^2+epsi); absy = sqrt(iy.^2+epsi);

   fc(i+1,:) = diag(.5*Mr'*AnAn*Mr - Mr'*AnMp + .5*MpMp + mu*(Isum'*absx+Isum'*absy));

   if norm(fc(i+1,:)) < norm(fc(i,:))
      % Just continue
   else % Cancel the iteration
      fc(i+1,:) = fc(i,:);
      Mr = Mro;
      if kappa > 1 % Decrease parameters
         kappa = .5*kappa;
      else
         step = .5*step;
      end
   end
end

% Plot the evolution of the cost-function for the parameters
figure; hold on;
plot(log10(fc(:,1)),'Color','red');
plot(log10(fc(:,2)),'Color','green');
plot(log10(fc(:,3)),'Color','blue');

%% And reconstruct the picture with the original format
Mc = zeros(n,m,3);
Mc(:,:,1) = reshape(Mr(:,1),[n,m]);
Mc(:,:,2) = reshape(Mr(:,2),[n,m]);
Mc(:,:,3) = reshape(Mr(:,3),[n,m]);
figure; imshow(Mc); % Display picture
toc; % Display the computation time
