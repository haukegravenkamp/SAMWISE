%% Mapping from local coordinates on triangle to square
% Eq. (5.7.1) in 
% Finite and spectral element methods using Matlab, Pozrikidis, C.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
function [xid,etad,J]=tri2square(xi,eta)

    xid=2*xi./(1-eta)-1;                                                    % map on standard square ("xi dash", "eta dash")
    etad=2*eta-1;

    xidDxi   = 2/(1-eta);
    xidDeta  = 2*xi./(1-eta).^2;
    etadDxi  = 0;
    etadDeta = 2;

    J=[xidDxi, etadDxi; xidDeta etadDeta];                                  % Jacobian


end
     
