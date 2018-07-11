function [index] = findCorner2 (x, y, varargin)
 % This function finds the corner of an L-curve
 % Yep, I did copy/paste the entire function just because I didn't want to code
 % a clean way to integrate the (-) at line 87. @dirtycoding

 if size(x,1) <= 2
    error(['Not enough points on L-curve : it has no corner', ...
      ' Rem : if you actually have lots of points, try to transpose the vectors'])
 end

 d = 2;
 if numel(varargin)>0
     d = cell2mat(varargin(1));
 end
 
 dolog = 1;
 if numel(varargin)>1
     dolog = cell2mat(varargin(2));
 end

 % put into loglog
 if dolog == 1
    x = log10(x); y = log10(y);
 end

 % Normalize
 x = (x-min(x))/(max(x)-min(x));
 y = (y-min(y))/(max(y)-min(y));

 sx = size(x,1);     % size of x (hopefully too the size of y)
 n  = floor(sx/d)+1; %max(sx-2,floor(sx/2)+1);   % degree of polynoms
 t  = (1:1:sx)';     % coarse mesh (unused)
 %tp = (1:.1:sx)';    % fine mesh (useless for now)

 %% Compute the curvilign abscissa
 % Differentiation matrix
 Dd = zeros(sx-1,sx);
 for i=1:sx-1
    Dd(i,i) = 1; Dd(i,i+1) = -1;
 end
 dx = Dd*x; dx = [0;dx];
 dy = Dd*y; dy = [0;dy];
 dl = sqrt(dx.^2 + dy.^2);
 t = cumsum(dl);
 %%
 if numel(varargin)>2 % This prescripted order erases the other one
     n = cell2mat(varargin(3));
 end

 % First, interpolate x and y
 px = polyfit(t,x,n)';
 py = polyfit(t,y,n)';

 % build the derivation matrix
 Md = zeros(n+1);
 for i=2:n+1
    Md(i,i-1) = n-i+2;
 end
 
 % Then derivate
 px1 = Md*px;
 px2 = Md*px1;
 py1 = Md*py;
 py2 = Md*py1;

 tt = zeros(n+1,sx);
 for j=1:n+1
    tt(j,:) = t.^(n+1-j);
 end
 xx = px'*tt; xx1 = px1'*tt; xx2 = px2'*tt;
 yy = py'*tt; yy1 = py1'*tt; yy2 = py2'*tt;
% figure
% hold on
% plot(xx,yy,'Color','red');
% plot(x,y,'Color','blue');

 % Find the point with the smallest curve radius
 R = 0;
 index = 0;
 Ga = zeros(sx,1);
 for i=2:sx-1  % First and last are forbitten (because of interpolation bound effects)
    
    denom = (xx1(i)^2 + yy1(i)^2)^(3/2);
    
    if denom ~= 0
       Ga(i) =  (yy2(i)*xx1(i) - xx2(i)*yy1(i)) / denom; % Yes, right there
    else
       Ga(i) = R+1.;  % eliminate this candidate
    end
    
    % Test if the candidate is smaller than the current curvature
    if R == 0 || Ga(i) < R % Ga < 0
       %if xx1(i+1)>0 % residual has to decrease
          R = Ga(i);
          index = i;
       %end
    end
 end

end
