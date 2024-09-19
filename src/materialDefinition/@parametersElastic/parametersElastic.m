classdef parametersElastic  < checkSetProperties
    % parametersElastic
    % constants of linear elastic material
    
    properties
        isotropy = 0                                                        % 0 isotropic; 1 transversely isotropic; 2 orthotropic; 3 general anisotropic
        rho                                                                 % mass density
        G                                                                   % shear modulus
        lambda                                                              % second Lamé parameter
        E                                                                   % Young's modulus
        nu                                                                  % Poisson's ratio
        cl                                                                  % longitudinal wave velocity
        cs                                                                  % shear wave velocity
        D                                                                   % elasticity matrix

        E1                                                                  % parameters of orthotropic material
        E2
        E3
        G23
        G31
        G12
        nu12
        nu13
        nu23
    end
    
    methods (Static)
        D = computeD_isotropic(G, nu, dof, isHandle)
        D2d = computeD2d(D,assumption2D)
        D = computeD_orthotropic(E11,E22,E33, v12,v13,v23, G23,G31,G12)

    end
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
