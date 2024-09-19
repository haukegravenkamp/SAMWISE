function [coordGlobal,J,Jdet]=coordinateTransformQuad(coordNodes,eta12,N,Ndeta1,Ndeta2,useLinear)

% all nodes
yv = coordNodes(:,2);
zv = coordNodes(:,3);

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% if shape functions not provided or explicitly asked for linear
% interpolation
if (nargin < 4) || useLinear
    % local coordinates
    e1 = eta12(:,1);
    e2 = eta12(:,2);

    % linear functions for interpolation of geometry
    N      = 0.25*[(1-e1).*(1-e2), (1+e1).*(1-e2), (1+e1).*(1+e2), (1-e1).*(1+e2)];
    Ndeta1 = 0.25*[(  -1).*(1-e2), (1   ).*(1-e2), (1   ).*(1+e2), (  -1).*(1+e2)];
    Ndeta2 = 0.25*[(1-e1).*(  -1), (1+e1).*(  -1), (1+e1)        , (1-e1)        ];

    % use corner nodes only
    yv = yv(1:4);
    zv = zv(1:4);
end

% coordinate derivatives
ydeta1 = Ndeta1 * yv;
zdeta1 = Ndeta1 * zv;
ydeta2 = Ndeta2 * yv;
zdeta2 = Ndeta2 * zv;
J = zeros(3,3,numel(ydeta1));
J(1,1,:) = 1;
J(2,2,:) = ydeta1;
J(2,3,:) = zdeta1;
J(3,2,:) = ydeta2;
J(3,3,:) = zdeta2;
Jdet = zeros(size(J,3),1);
for i = 1:size(J,3)
    Jdet(i) = abs(det(J(:,:,i)));
end

coordGlobal = zeros(numel(ydeta1),3);
coordGlobal(:,2,:) = N*yv;
coordGlobal(:,3,:) = N*zv;

end
