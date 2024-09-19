function [coordGlobal,J,Jdet]=coordinateTransform(~,~,coordNodes,eta,~,~,~,~)
% cylindrical coordinates are chosen as (x,theta,r)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% endpoints of 1D element 
r1 = coordNodes(1,3);
r2 = coordNodes(2,3);

L = (r2-r1);                                                                % element length
coordGlobal = zeros(numel(eta),3);
coordGlobal(:,3) = r1 + (eta+1)*L/2;                                        % vertical coordinate z

J = zeros(3,3,numel(eta));

xdx = 1;
zdeta = L/2;
ydtheta = -(r1+(1+eta)/2*L);

J(1,1,:) = xdx;
J(2,2,:) = ydtheta;
J(3,3,:) = zdeta;
Jdet = L/2*r1 + L^2/4*(eta+1);


end
