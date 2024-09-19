classdef bc < matlab.mixin.Heterogeneous & handle
    % boundary conditions 

    properties
        location                                                            % 'top', 'bottom'
    end

    properties (Hidden = true)
        materialNo                                                          % material number at the boundary condition                                                            
        coord                                                               % coordinates
        dofs                                                                % degree of freedom
        directions                                                          % direction(s) to apply b.c.
    end
    
    methods
        function obj = bc
        end
        
    end

    methods (Sealed)
        [obj, glb] = analyzeBoundaryConditions(obj,geo,msh,mat,opt,glb)
        glb = findFixedDofs(obj,msh,glb);
    end

end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
