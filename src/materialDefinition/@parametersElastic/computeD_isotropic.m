function D = computeD_isotropic(G, nu, dof,isHandle)

if nargin < 3                                                               % assume 3d by default
    dof = 3;
end

if isHandle                                                                 % parameters defined as functions
     if dof == 1   
        D = @(y,z) eye(2)*G(y,z);
    elseif dof == 2
        D = @(y,z) 2*G(y,z)/(1-2*nu(y,z))*...
            [1-nu(y,z)  nu(y,z)     0; ...
            nu(y,z)  1-nu(y,z)    0;...
            0   0 (1-2*nu(y,z))/2];
    elseif dof == 3
        D = @(y,z) 2*G(y,z)/(1-2*nu(y,z))*...
            [1-nu(y,z), nu(y,z) , nu(y,z) , 0 , 0 , 0; ...
            nu(y,z) ,1-nu(y,z), nu(y,z) , 0 , 0 , 0;...
            nu(y,z) , nu(y,z) ,1-nu(y,z), 0 , 0 , 0;...
            0 , 0 , 0 , (1-2*nu(y,z))/2, 0 , 0;...
            0 , 0 , 0 , 0  , (1-2*nu(y,z))/2, 0;...
            0 , 0 , 0 , 0  , 0 , (1-2*nu(y,z))/2];
    end

else

    if dof == 1
        D=eye(2)*G;
    elseif dof == 2
        D = 2*G/(1-2*nu)*...
            [1-nu  nu     0; ...
            nu  1-nu    0;...
            0   0 (1-2*nu)/2];
    elseif dof == 3
        D = 2*G/(1-2*nu)*...
            [1-nu, nu , nu , 0 , 0 , 0; ...
            nu ,1-nu, nu , 0 , 0 , 0;...
            nu , nu ,1-nu, 0 , 0 , 0;...
            0 , 0 , 0 , (1-2*nu)/2, 0 , 0;...
            0 , 0 , 0 , 0  , (1-2*nu)/2, 0;...
            0 , 0 , 0 , 0  , 0 , (1-2*nu)/2];
    end

end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
