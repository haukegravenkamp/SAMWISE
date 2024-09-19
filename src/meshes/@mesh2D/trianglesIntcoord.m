function [coordGl,coordLc] = trianglesIntcoord(vCoord,p,v)
%% compute coordinates of internal nodes
%
% Hauke Gravenkamp, gravenkamp.research@gmail.com, 2019
% m: edge order
% vCoord: vertex coordinates
% v: 1d nodal positions

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% Nodal positions according to Pozrikidis textbook, Section 5.8

% nModesDom = sum(1:(p-2));                                                   % number of nodes in the interior
% 
% v = (v+1)/2;                                                                % scale to 0->1
% 
% xi  = zeros(nModesDom,1);                                                   % local coordinates on triangle (between 0 and 1!)
% eta = zeros(nModesDom,1);
% 
% c = 1;                                                                      % counter
% for i = 2:p                                                                 % compute local coordinates
%     for j = 2:(p+1-i)
%         k=p+3-i-j;
%         xi(c)  = 1/3*(1 + 2*v(i) -   v(j) - v(k));
%         eta(c) = 1/3*(1 -   v(i) + 2*v(j) - v(k));
%         c = c+1;
%     end
% end

[xi,eta] = triIntCoordsLocal(p,v);

% vertex coordinates
x1 = vCoord(1,1);
x2 = vCoord(2,1);
x3 = vCoord(3,1);
y1 = vCoord(1,2);
y2 = vCoord(2,2);
y3 = vCoord(3,2);

% map to Cartesian coordinates
x = x1 + (x2-x1)*xi + (x3-x1)*eta;
y = y1 + (y2-y1)*xi + (y3-y1)*eta;

% coordinates of internal nodes
coordGl = [x,y];
coordLc = [xi,eta];

   



end
