function [N,Nd1,Nd2,Na,Nad1,Nad2] = shapeFuntionsTri(p,eta12,dof,intTab,nModes,nDom)
% compute shape functions and derivatives at all Gauss points
% We return not only the scalar shape functions N and derivatives Nd1, Nd2, 
% but also the matrices of all shape functions for all dofs Na, Nad1, Nad2
% This is convenient because we need the scalar shape functions for the
% geometry interpolation but the vector version for the FE matrices.
% If the PDE is scalar, both are the same though

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
nEta = size(eta12,1);                                                       % number of Gauss points

N = zeros(nEta,nModes);                                                     % allocate scalar shape functions
Nd1 = N;
Nd2 = N;

Na = zeros(dof,nModes*dof,nEta);                                            % allocate vectorial shape functions for all dofs
Nad1 = Na;
Nad2 = Na;

eta1=intTab(p(1)+1).GLL.xi;                                                 % nodes at GLL points
eta2=intTab(p(2)+1).GLL.xi;

[eta1,eta2] = meshgrid(eta1(2:end-1),eta2(2:end-1));                        % 2D grid of interior local coordinates
eta1=eta1(:);                                                               % convert to vector
eta2=eta2(:);
coordDom = [eta1,eta2];                                                     % local coordinates of interior nodes

[proriolOrders,iV] = prepareInnerModesTri(p(1),intTab);
% iT = blendingTri(p,coordDom,nModes,nDom,intTab,proriolOrders,iV);                               % blending operator for creating nodal shape functions
iT = [];

for i = 1:nEta                                                              % loop Gauss points

    [shp,shpDeta1,shpDeta2] = shapeFunctionsXnYtri(p,eta12(i,:),nModes,nDom,intTab,iT,proriolOrders,iV);            % compute XnY shape functions

    N(i,:)   = shp;                                                    % shape function matrix
    Nd1(i,:) = shpDeta1;                                                    % shape function derivatives
    Nd2(i,:) = shpDeta2;                                                    % shape function derivatives

    Na(:,:,i)   = kron(shp,eye(dof));                                  % shape function matrix for all dofs
    Nad1(:,:,i) = kron(shpDeta1,eye(dof));                                  % shape function derivatives
    Nad2(:,:,i) = kron(shpDeta2,eye(dof));                                  % shape function derivatives

end

end
