function [N,Nd1,Nd2,Na,Nad1,Nad2] = getShapeFuntions(p,eta12,dof,intTab,nModes,nDom,nV)
% compute shape functions and derivatives at all Gauss points
% We return not only the scalar shape functions N and derivatives Nd1, Nd2, 
% but also the matrices of all shape functions for all dofs Na, Nad1, Nad2
% This is convenient because we need the scalar shape functions for the
% geometry interpolation but the vector version for the FE matrices.
% If the PDE is scalar, both are the same though

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if nV == 3
    [N,Nd1,Nd2,Na,Nad1,Nad2] = mesh2D.shapeFuntionsTri(p,eta12,dof,intTab,nModes,nDom);
elseif nV == 4
    [N,Nd1,Nd2,Na,Nad1,Nad2] = mesh2D.shapeFuntionsQuad(p,eta12,dof,intTab,nModes,nDom);
end





end
