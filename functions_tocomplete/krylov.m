close all;
clear all;
addpath(genpath('./tools'));

%[ Mp, n, m ] = readFromPng('../images/0-brd.png');
[ Mp, n, m ] = readFromMat('../images/0-brd.mat'); % Load the data

A = operator(n,m,5,max(n,m)); % Build the blur operator
sA = size(A,1); % His size

%% Let the fun begin
tic; % Initialize the chronometer
dolcurve = 0; % No Picard plot for today as it would require to build the 
% Ritz basis, and that is not straightforward for GMRES

% Remark: we could use a preconditionner, it would be the inverse of the 
% regularization matrix we would use in Tikhonov
% Anyway, nothing prevents us to use it for the computation of the L-curve

L  = eye(sA); % The regularizing operator
%[Dx,Dy] = derivative(n,m); L = Dx'*Dx + Dy'*Dy;

if dolcurve == 0

   nmax = ; % Max iterations

   Itere = myorthodir (A, Mp, nmax, L);
   Mr = % Select the last one (in Mrhs)

else

   % We're doing a L-curve: 1 point per iteration, it's not more expensive
   nmax = ; % Iterations
   [Itere,res,nor] = myorthodir (A, Mp, nmax, L);

   nor(1,:) = nor(2,:)/2; % Erase the 0 here (as it will be logarithmified)

   zebest1 = findCorner2 (res(:,1), nor(:,1), 3, 1, 7); % And find the corners
   zebest2 = findCorner2 (res(:,2), nor(:,2), 3, 1, 7);
   zebest3 = findCorner2 (res(:,3), nor(:,3), 3, 1, 7);
   % /!\ The L-curve is not browsed in the same direction as for Tikhonov /!\
   % Hence, the function is not totally the same

   % Plot the L-curve
   figure; hold on;
   loglog(res(:,1),nor(:,1),'+-','Color','red');
   loglog(res(zebest1,1),nor(zebest1,1),'o','Color','red');
   loglog(res(:,2),nor(:,2),'+-','Color','blue');
   loglog(res(zebest2,2),nor(zebest2,2),'o','Color','blue');
   loglog(res(:,3),nor(:,3),'+-','Color','green');
   loglog(res(zebest3,3),nor(zebest3,3),'o','Color','green');
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
