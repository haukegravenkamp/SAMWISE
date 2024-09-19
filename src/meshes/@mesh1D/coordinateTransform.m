function [coordGlobal,J,Jdet]=coordinateTransform(~,~,coordNodes,eta,~,~,~,~)
% coordinate transformation for plate structures
% we always use linear interpolation between the first and last node of
% each layer, assuming the interface is a straight line

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% endpoints of 1D element 
z1 = coordNodes(1,3);
z2 = coordNodes(2,3);

% global coordinates at Gauss points
coordGlobal = zeros(numel(eta),3);
coordGlobal(:,3) = (z2+z1)/2 + (z2-z1)/2*eta;

% Jacobian
J = zeros(3,3,numel(eta));
zdeta2 = (z2-z1)/2;
J(1,1,:) = 1;
J(2,2,:) = 1;
J(3,3,:) = zdeta2;
Jdet = ones(numel(eta),1)*zdeta2;


end
