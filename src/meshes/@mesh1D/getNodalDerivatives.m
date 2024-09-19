function shpDz = getNodalDerivatives(p,dofsGl,nEleDofs,nDofs,dof,intTab,Jdet)
% compute derivatives at nodes
etaNodes  = intTab(p+1).GLL.xi;                                             % local coordinates at all nodes
shpDzC    = zeros(nEleDofs,nEleDofs);                                       % allocate dN/dz for current element, local dofs
shpDz = zeros(nDofs,nDofs);                                                 % allocate dN/dz for current element, global dofs

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
for iN = 1:numel(etaNodes)                                                  % loop nodes

    shp = shpLagrange(p+1,etaNodes(iN),intTab);                             % get shape functions and derivatives
    Jeti = Jdet(iN);                                                        % Jacobian determinant at current location

    Ndz = kron(shp(2,:),eye(dof))/Jeti;                                     % shape function derivatives wrt z

    shpDzC((iN-1)*dof+(1:dof),:) = Ndz;                                     % store in matrix for current element

end
shpDz(dofsGl,dofsGl) = shpDzC;                                              % store in cell for all elements
% This array is in global coordinates!

end
