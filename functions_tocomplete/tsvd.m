close all;
clear all;
addpath(genpath('./tools'));

%[ Mp, n, m ] = readFromPng('../images/0-brd.png');
[ Mp, n, m ] = readFromMat('../images/0-brd.mat'); % Load the data

A = operator(n,m,5,max(n,m)); % Build the blur operator
sA = size(A,1); % His size

%% Let the fun begin
tic; % Initialize the chronometer
dopicard = 1;
nmax = % Max nb of singular values computed
ntot = % Nb of singular values used (inferior to nmax)

% Remark: we could use a preconditionner (GSVD), it would be the invert of the 
% regularization matrix we would use in Tikhonov

%% Actually estimate the nmax first eigenvalues 
%(don't worry too much about that part)
Za    = zeros(sA,sA); OAATO = [ Za, A ; A', Za ];
[UV, thetatheta] = eigs(OAATO, 2*nmax);
[theta,Ind3] = sort(diag(thetatheta),'descend');
UV = UV(:,Ind3);
% From the eigenvalues of that fat matrix, we will deduce the sing values of A

U     = sqrt(2)*UV(1:sA,1:2*nmax);
V     = sqrt(2)*UV(sA+1:end,1:2*nmax);

indneg = find(theta<0);
theta(indneg) = []; V(:,indneg) = []; U(:,indneg) = [];
% Remove float abherrations and mirror values
% Mirror O Mirror on the wall, tell me who is the most beautyful of them all
tm1 = diag(1./theta);
Sigma = diag(theta(1:ntot));

if dopicard == 0

   V = ; tm1 = ; U = ; % Truncate
   Mr = ; % And solve

else

   rhs = U'*Mp; % The right hand side

   rplp1 = rhs(:,1)./theta;
   rplp2 = rhs(:,2)./theta;
   rplp3 = rhs(:,3)./theta;

   me1 = mean(abs(rplp1))/1e5; arplo1 = max(me1,abs(rplp1)); % Remove 0 values
   me2 = mean(abs(rplp2))/1e5; arplo2 = max(me2,abs(rplp2)); % 'cause of the log
   me3 = mean(abs(rplp3))/1e5; arplo3 = max(me3,abs(rplp3));

   % Plot the Picard Plot
   figure; hold on;
   semilogy(arplo1,'Color','red');
   semilogy(arplo2,'Color','green');
   semilogy(arplo3,'Color','blue');
   legend('Red','Green','Blue (bet you didnt guess)');

   indm1 = findPicard2 (log10(arplo1), ceil(sA/7), 1, 3); % Find the best number
   indm2 = findPicard2 (log10(arplo2), ceil(sA/7), 1, 3); % For each term
   indm3 = findPicard2 (log10(arplo3), ceil(sA/7), 1, 3);
   % Rem : don't worry too much for the arguments in findPicard2
   % Don't worry also for the singular matrix warnings, it's nothing

   V1 =; tm11 =; U1 =; % Truncate
   V2 =; tm12 =; U2 =;
   V3 =; tm13 =; U3 =;
   Mr(:,1) = ; % And solve
   Mr(:,2) = ;
   Mr(:,3) = ;

end

%% And reconstruct the picture with the original format
Mc = zeros(n,m,3);
Mc(:,:,1) = reshape(Mr(:,1),[n,m]);
Mc(:,:,2) = reshape(Mr(:,2),[n,m]);
Mc(:,:,3) = reshape(Mr(:,3),[n,m]);
figure; imshow(Mc); % Display picture
toc; % Display the computation time
