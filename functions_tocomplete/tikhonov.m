close all;
clear all;
addpath(genpath('./tools'));

%[ Mp, n, m ] = readFromPng('../images/0-brd.png');
[ Mp, n, m ] = readFromMat('../images/0-brd.mat'); % Load the data

A = operator(n,m,5,max(n,m)); % Build the blur operator
sA = size(A,1); % His size

%% Let the fun begin
tic; % Initialize the chronometer
dolcurve = 0;

%L  = eye(sA); % The regularizing operator
[Dx,Dy] = derivative(n,m); L = Dx'*Dx + Dy'*Dy;

if dolcurve == 0

   mu = % Regularization parameter
   Mr = % Actually invert the regularized system

else

   % We're doing a L-curve: first compute the solution for many values of the 
   % regularization parameter, and then plot the curve and choose the best one.

   mu = % List of all the tested values for mu0
   smu = max(size(mu)); % mu0 is supposed to be a vector

   res = zeros(smu,1); % Initialize the residual
   nor = zeros(smu,1); % Initialize the norm
   Mrr = zeros(sA,3*smu); % Yep, I'll store 1 rgb-picture per value of mu
   % But after all, we're already storing a full operator of size nb pixels^2...

   for i = 1:smu
      ind = [3*i-2,3*i-1,3*i]; % Multinidex
      Mrr(:,ind) = % Invert the system
      res(i) = % Compute the residual
      % By the way, if you want to use a norm for a multiple vector, use
      % the Frobenius norm (unless of course you're a maniac)
      nor(i) = % And the norm of the solution
   end

   zebest = findCorner (res, nor, 3); % And find the corner
   % By the way, there exist an entire library for inverse problems called
   % Regularization Tools, by Christian Hansen (who's some sort of semi-god of
   % inverse problems). We're not using it right here because this library has
   % depandancies and I didn't want to be fucked by a freaking 
   % depandancy problem.
   % If interested, go to http://www.imm.dtu.dk/~pcha/Regutools/

   % Plot the L-curve
   figure; hold on;
   loglog(res,nor,'+-','Color','blue');
   loglog(res(zebest),nor(zebest),'o','Color','red');
   legend('L-curve'); 

   Mr = % Select the chosen solution

end

%% And reconstruct the picture with the original format
Mc = zeros(n,m,3);
Mc(:,:,1) = reshape(Mr(:,1),[n,m]);
Mc(:,:,2) = reshape(Mr(:,2),[n,m]);
Mc(:,:,3) = reshape(Mr(:,3),[n,m]);
figure; imshow(Mc); % Display picture
toc; % Display the computation time
