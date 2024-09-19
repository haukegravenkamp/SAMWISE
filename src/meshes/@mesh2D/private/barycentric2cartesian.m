
function [x,y,J] = barycentric2cartesian(u,v,w,vert)
% conversion of barycentric coordinates u,v,w to Cartesian
% coordinates x,y in triangle with vertex coordinates vert

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
x=[u,v,w]*vert(:,1);                                                        % x= u*x1 + v*x2 + w*x3
y=[u,v,w]*vert(:,2);                                                        % y= u*y1 + v*y2 + w*y3

xdu=vert(1,1)-vert(3,1);                                                    % x= u*x1 + v*x2 + (1-u-v)*x3
xdv=vert(2,1)-vert(3,1);
%         xdw=vert(3,1);

ydu=vert(1,2)-vert(3,2);
ydv=vert(2,2)-vert(3,2);
%         ydw=vert(3,2);

J=[xdu, ydu; xdv ydv];                                          % Jacobian

end
